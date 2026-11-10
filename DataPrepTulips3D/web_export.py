"""Writes a MESA_data_dict (the same dict convertMesaData() returns and
save_to_pickle()/save_to_texture() consume for the Blender addon) out as a
self-contained .tulips3dweb bundle - the format the TULIPS-3D-WEB browser
viewer (repos/TULIPS-3D-WEB) loads directly, no Blender, no EXR, no pickle
parsing in JS.

Why a bespoke tiny format at all, instead of shipping the pickle/EXR
textures directly: the pickle already holds every array needed in raw
numpy form, the EXR textures are a Blender-geometry-nodes-specific derived
artifact a browser has no use for, and pickle itself isn't safely
parseable in JS (numpy's pickling uses its own C reduce protocol). So:
re-serialize the handful of arrays the web viewer actually needs into
manifest.json (dtype/shape/byte-range per array) + one concatenated
data.bin, zipped together - directly analogous to how glTF separates JSON
metadata from a binary buffer.

Two bundle shapes, both sharing that one data.bin layout:

- export_for_web_dict() (format_version 1): ONE dataset, ONE
  time_scale_type (and, for a binary, ONE star) per bundle - the original
  shape, still produced by TULIPS-3D-WEB/tools/export_for_web.py's CLI for
  backfilling/one-off re-exports.
- export_for_web_combined() (format_version 2): every time_scale_type
  convertMesaData() was run for (Ben: "bundle the linear, log_to_end and
  model_number data ... into one ... bundle and give the user ... the
  choice"), and, for a binary, BOTH stars (Ben: "if it is a binary, just
  bundle the two stars into one ... bundle since the user can now switch
  between Profile and Binary view") - ONE file a user drops on the web
  app to get everything. Array names are namespaced per time_scale_type
  (and star, for a binary) so they all still fit in one flat
  manifest.arrays table - see TULIPS-3D-WEB's src/bundleLoader.js for how
  the browser side resolves one (time_scale_type[, star]) selection back
  down to the same shape export_for_web_dict() produces directly.
"""
import json
import os
import zipfile

import numpy as np

from . import blender_keys_prep as key

FORMAT_VERSION = 1
FORMAT_VERSION_COMBINED = 2


def _as_f32(arr):
    return np.ascontiguousarray(arr, dtype=np.float32)


def _collect_leaf(d, prefix, add_array):
    """Builds the per-(time_scale_type[, star]) manifest dict for `d`,
    writing its arrays via add_array(name, arr) with `prefix.` prepended
    to every array name (prefix="" for the unprefixed format_version-1
    case - see export_for_web_dict below)."""
    def pname(suffix):
        return f"{prefix}.{suffix}" if prefix else suffix

    manifest = {
        "binary_nr": int(d[key.binary_nr]) if d.get(key.is_binary) else None,
        "t_resolution": int(d[key.t_resolution]),
        "t_resolution_orig": int(d.get(key.t_resolution_orig, d[key.t_resolution])),
        "time_scale_type": d.get(key.time_scale_type, "unknown"),
        "r_resolution": int(d[key.r_resolution]),
        "r_label": d.get("r_label", "mass"),
        "prof_labels": list(d[key.prof_labels]),
        "chem_elem_labels": list(d.get(key.chem_elem_labels, [])),
    }

    # --- radial profile data (always present) ---
    add_array(pname(key.data_prof_t_r), d[key.data_prof_t_r])  # (n_labels, t, r)
    if key.data_chem_abun_id in d:
        add_array(pname(key.data_chem_abun_id), d[key.data_chem_abun_id])  # (t, r, n_theta_chem)
        manifest["n_theta_chem"] = int(d[key.data_chem_abun_id].shape[2])

    # Outer stellar radius in Rsun per frame - every DataPrepTulips3D
    # vintage stores this redundantly across all prof_labels rows; row 0
    # matches what the Blender addon itself reads (R_star_from_grid[0, :]).
    r_max_rsun = np.asarray(d[key.data_r_max_Rsun])
    add_array(pname("r_max_rsun"), r_max_rsun[0, :] if r_max_rsun.ndim == 2 else r_max_rsun)

    # The true (non-uniform) radius-in-Rsun at each (time, r-index) grid
    # point, where the grid itself is built on a uniform ENCLOSED-MASS
    # coordinate (r_label="mass") - this is what "Rsun mode" actually
    # samples; absent in older DataPrepTulips3D-vintage pickles (handled
    # client-side with a documented uniform-spacing fallback - see
    # TULIPS-3D-WEB's src/radialScale.js).
    manifest["has_radial_remap"] = "data_radial_remap_rsun" in d
    if manifest["has_radial_remap"]:
        add_array(pname("data_radial_remap_rsun"), d["data_radial_remap_rsun"])  # (t, r)

    # --- every data_t scalar-per-frame field, generically ---
    data_t_keys = []
    for k, arr in d[key.data_t].items():
        add_array(pname(f"data_t.{k}"), arr)
        data_t_keys.append(k)
    manifest["data_t_keys"] = data_t_keys

    # --- binary-only baked geometry (Roche lobes, potential, rings) ---
    manifest["has_rochelobe"] = key.data_rochelobe_shape_1 in d
    if manifest["has_rochelobe"]:
        add_array(pname(key.data_rochelobe_shape_1), d[key.data_rochelobe_shape_1])
        add_array(pname(key.data_rochelobe_shape_2), d[key.data_rochelobe_shape_2])
        manifest["nr_Th_RL"] = int(d[key.nr_theta_points_rochelobe])
        manifest["nr_Phi_RL"] = int(d[key.nr_phi_points_rochelobe])

    manifest["has_potential_grid"] = key.data_rochelobe_potential_grid in d
    if manifest["has_potential_grid"]:
        add_array(pname(key.data_rochelobe_potential_grid), d[key.data_rochelobe_potential_grid])
        manifest["nr_x_potential"] = int(d[key.nr_x_points_potential])
        manifest["nr_y_potential"] = int(d[key.nr_y_points_potential])
        manifest["potential_x_range"] = list(d[key.potential_x_range])
        manifest["potential_y_range"] = list(d[key.potential_y_range])

    manifest["has_equipotential_rings"] = key.data_rochelobe_l2_ring in d
    if manifest["has_equipotential_rings"]:
        add_array(pname(key.data_rochelobe_l2_ring), d[key.data_rochelobe_l2_ring])
        add_array(pname(key.data_rochelobe_l3_ring), d[key.data_rochelobe_l3_ring])
        manifest["nr_Th_equipot"] = int(d[key.nr_theta_points_equipotential])

    return manifest


def _write_bundle(manifest, arrays, buf, out_path):
    manifest["arrays"] = arrays
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    # ZIP_STORED (no compression), not ZIP_DEFLATED: this bundle is read
    # from a LOCAL file the user drops on the page, not downloaded over a
    # network, so on-disk size barely matters - but decompression speed
    # does. Node's zlib (what a quick `JSZip.loadAsync` timing test in
    # this project's own tooling uses) wraps native C zlib and inflated a
    # combined binary bundle's ~570MB data.bin in ~3.5s; the SAME file in
    # an actual browser - where JSZip has no native zlib to call into and
    # falls back to its own pure-JS inflate - took well over a minute and
    # was still not done (confirmed via a real Puppeteer test: every
    # debug hook past the loadBundle() call stayed unset for 80+
    # seconds). Switching to STORED makes reading an entry a plain
    # byte-range copy instead, independent of this gap - essential once a
    # combined bundle (Ben: "bundle the linear, log_to_end and
    # model_number data ... into one") can be hundreds of MB.
    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_STORED) as zf:
        zf.writestr("manifest.json", json.dumps(manifest))
        zf.writestr("data.bin", bytes(buf))
    size_mb = os.path.getsize(out_path) / 1e6
    print(f"Wrote {out_path} ({size_mb:.1f} MB, {len(arrays)} arrays)")


def export_for_web_dict(d, out_path):
    """Writes `d` (an in-memory MESA_data_dict, as returned by
    convertMesaData()/load_from_pickle()) to `out_path` as a single-
    time_scale_type, single-star .tulips3dweb bundle (format_version 1)."""
    arrays = {}  # name -> {dtype, shape, offset, length}
    buf = bytearray()

    def add_array(name, arr):
        a = _as_f32(arr)
        offset = len(buf)
        buf.extend(a.tobytes())
        arrays[name] = {"dtype": "float32", "shape": list(a.shape), "offset": offset, "length": a.nbytes}

    manifest = _collect_leaf(d, "", add_array)
    manifest["format_version"] = FORMAT_VERSION
    manifest["is_binary"] = bool(d.get(key.is_binary, False))
    manifest["binary_nr"] = int(d[key.binary_nr]) if d.get(key.is_binary) else None

    _write_bundle(manifest, arrays, buf, out_path)
    print(f"  is_binary={manifest['is_binary']}, t_resolution={manifest['t_resolution']}")
    return out_path


def export_for_web_combined(variants, out_path, is_binary=False):
    """Writes ONE .tulips3dweb bundle (format_version 2) covering every
    time_scale_type DataPrepTulips3D was run for, and - for a binary -
    both stars at once.

    `variants`: {time_scale_type: d} for a single star, or
    {time_scale_type: {1: d_star1, 2: d_star2}} for a binary (d = a
    MESA_data_dict, as returned by convertMesaData()/load_from_pickle()).
    """
    arrays = {}
    buf = bytearray()

    def add_array(name, arr):
        a = _as_f32(arr)
        offset = len(buf)
        buf.extend(a.tobytes())
        arrays[name] = {"dtype": "float32", "shape": list(a.shape), "offset": offset, "length": a.nbytes}

    time_scale_types = list(variants.keys())
    manifest = {
        "format_version": FORMAT_VERSION_COMBINED,
        "is_binary": bool(is_binary),
        "time_scale_types": time_scale_types,
        "variants": {},
    }

    for ts, d_or_pair in variants.items():
        if is_binary:
            d1, d2 = d_or_pair[1], d_or_pair[2]
            m1 = _collect_leaf(d1, f"{ts}.star1", add_array)
            m2 = _collect_leaf(d2, f"{ts}.star2", add_array)
            manifest["variants"][ts] = {"star1": m1, "star2": m2}
        else:
            manifest["variants"][ts] = _collect_leaf(d_or_pair, ts, add_array)

    _write_bundle(manifest, arrays, buf, out_path)
    print(f"  is_binary={is_binary}, time_scale_types={time_scale_types}")
    return out_path


def combine_from_prepped_dirs(variant_dirs, out_path, is_binary=False, pickle_filename="MESA_data_dict.pkl"):
    """Convenience wrapper around export_for_web_combined() that loads
    each variant from an already-prepped directory's pickle (as
    convertMesaData() writes it) rather than requiring the caller to keep
    every in-memory dict around.

    `variant_dirs`: {time_scale_type: save_to_dir} for a single star, or
    {time_scale_type: {1: save_to_dir_star1, 2: save_to_dir_star2}} for a
    binary.
    """
    from .mesa_data_interp import load_from_pickle  # local import: avoids a circular import with mesa_data_interp, which imports this module

    variants = {}
    for ts, d_or_dirs in variant_dirs.items():
        if is_binary:
            variants[ts] = {nr: load_from_pickle(os.path.join(dir_, pickle_filename)) for nr, dir_ in d_or_dirs.items()}
        else:
            variants[ts] = load_from_pickle(os.path.join(d_or_dirs, pickle_filename))

    return export_for_web_combined(variants, out_path, is_binary=is_binary)
