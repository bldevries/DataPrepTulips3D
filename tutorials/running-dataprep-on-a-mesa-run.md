# Running DataPrepTulips3D on a MESA run

This tutorial walks through turning a finished MESA run into the
pickle + EXR-texture data that the TULIPS-3D Blender addon reads, for
both a single star and a binary system.

## 1. What you need before you start

**A finished MESA run's `LOGS` directory**, containing:

- A history file (`history.data` by default, or a custom name - for a
  binary run MESA typically names these `history1.data`/`history2.data`,
  one per star).
- The saved profile files (`profile1.data`, `profile2.data`, ...) and the
  `profiles.index` file MESA writes alongside them. These are what MESA
  produces by default whenever profiles are being saved during the run
  (the normal case) - you don't need to change any run-time settings to
  get them.

You do **not** need eccentricity, rotation, or any other specific
`&controls`/history-column setting turned on - DataPrepTulips3D detects
what's actually present in your history file and skips anything that
isn't there (rotation arrows, Roche lobe geometry, mass-transfer data,
etc. are all optional and silently omitted if the underlying MESA columns
are missing - see section 5).

**A Python environment**. From this repo's own top-level directory (the
one containing `setup.py`):

```bash
pip install -e .
```

This pulls in everything DataPrepTulips3D actually needs
(`numpy`, `scipy`, `matplotlib`, `Pillow`, `mesaPlot`, `openexr_numpy` -
`mesaPlot` reads the MESA data itself, `Pillow`/`openexr_numpy` write the
output textures).

## 2. The one function you need: `convertMesaData`

Everything goes through a single entry point:

```python
from DataPrepTulips3D.mesa_data_interp import convertMesaData

d = convertMesaData(
    mesa_LOGS_directory="/path/to/your/run/LOGS",
    t_resolution=100,
    r_resolution=100,
    save_to_dir="/path/to/where/you/want/the/output",
    time_scale_type="log_to_end",
)
```

- `mesa_LOGS_directory`: your MESA run's `LOGS` folder.
- `t_resolution`: how many time steps to resample the run down to (the
  addon animates over these, so this is effectively "how many frames").
- `r_resolution`: how many radial steps each profile (temperature,
  density, composition, ...) gets resampled to.
- `save_to_dir`: where the output goes - gets created if it doesn't
  exist. This becomes the directory you later point the TULIPS-3D
  addon's sidebar at.
- `time_scale_type`: how the `t_resolution` time steps are spaced across
  the run - see section 4.

This single call reads the history + profile data, computes all the
derived quantities the addon visualizes (Teff coloring, chemical
composition, and - for a binary - Roche lobe shapes, the Roche potential
landscape, mass-transfer rate, L1/L2/L3 points), and writes out:

- `<save_to_dir>/MESA_data_dict.pkl` - the actual data, pickled.
- `<save_to_dir>/textures/` - the same data, ALSO baked into `.exr`
  images (one per quantity, one pixel column per time step) - this is
  the form the Blender addon actually samples at render time.

**Caching**: if `<save_to_dir>/MESA_data_dict.pkl` already exists,
`convertMesaData` does nothing and just loads and returns that existing
file - it will NOT re-process your MESA run. If you want to redo a run
(e.g. with a different `t_resolution`), delete the old `save_to_dir`
first, or use a different `save_to_dir`.

## 3. Single star

```python
import os
from DataPrepTulips3D.mesa_data_interp import convertMesaData

MESA_LOGS = "/path/to/your/run/LOGS"
OUT_DIR = "/path/to/prepped_data/my_star"

os.makedirs(OUT_DIR, exist_ok=True)

d = convertMesaData(
    mesa_LOGS_directory=MESA_LOGS,
    t_resolution=200,
    r_resolution=100,
    save_to_dir=OUT_DIR,
    time_scale_type="log_to_end",
)
print("Done:", list(d.keys()))
```

Point the TULIPS-3D addon's "Star options" sidebar directory field at
`OUT_DIR` and click "Create" - see the TULIPS-3D repo's own docs for the
addon side.

## 4. Binary system (two stars)

A binary run needs to be prepped **once per star**, into **two separate
output directories** - the addon then loads both and combines them into
one visualization. Pass `is_binary=True` and `binary_nr=1`/`2` so the
output is correctly tagged for which star it is:

```python
import os
from DataPrepTulips3D.mesa_data_interp import convertMesaData

BASE = "/path/to/your/binary/run"       # contains LOGS1/ and LOGS2/
OUT_BASE = "/path/to/prepped_data"

for nr, logs, hist in [(1, "LOGS1", "history1.data"), (2, "LOGS2", "history2.data")]:
    out_dir = os.path.join(OUT_BASE, f"my_binary_star{nr}")
    os.makedirs(out_dir, exist_ok=True)
    d = convertMesaData(
        mesa_LOGS_directory=os.path.join(BASE, logs),
        t_resolution=200,
        r_resolution=100,
        save_to_dir=out_dir,
        filename_history=hist,
        time_scale_type="log_to_end",
        is_binary=True,
        binary_nr=nr,
    )
    print(f"Star {nr} done:", list(d["data_t"].keys()))
```

`filename_history` is only needed if your history file isn't literally
called `history.data` (the usual case for a binary run's two stars, which
MESA names `history1.data`/`history2.data`).

Roche lobe shapes, the Roche potential landscape, L1/L2/L3 points, and
the mass-transfer rate are computed automatically whenever the history
file actually contains both stars' masses, the orbital separation, and
`r1` (i.e. whenever it genuinely is a binary history file) - there's
nothing extra to enable. In the Blender addon, point "Star 1" at
`my_binary_star1` and "Star 2" at `my_binary_star2`, then use "Create
binary".

## 5. Choosing `t_resolution` and `time_scale_type`

`t_resolution` trades off animation smoothness against processing time
and output size - data-prep re-reads MESA profile files per time step it
keeps, so cost scales roughly linearly with `t_resolution`. A few hundred
is a reasonable starting point; a full run's native resolution (every
row MESA saved, no resampling) can take many minutes and produce
gigabyte-sized output for a long run.

`time_scale_type` controls WHICH of the run's original time steps get
kept, not just how many:

- `"model_number"`: evenly spaced by MESA model number/row index. Simple,
  but if your run spends most of its models on a short, fast-evolving
  phase (common near a stellar-evolution end state), that phase will
  dominate the resampled output.
- `"linear"`: evenly spaced in real elapsed time (`star_age`). Tends to
  cluster samples toward whichever part of the run actually took the
  most real time (often early evolution).
- `"log_to_end"`: evenly spaced in `log10(final_age - age)` - biases
  samples toward the END of the run. This is the usual choice for
  visualizing a star's late-stage evolution (mass transfer, common
  envelope, a flash, etc.), since those short-but-dramatic phases would
  otherwise be represented by only a handful of frames under
  `"model_number"` or `"linear"`.

There's no wrong answer - try `"log_to_end"` first, and switch to
`"model_number"` or `"linear"` if a particular phase you care about isn't
represented well.

## 6. Other useful options

- `profiles`: by default, `["mass", "logT", "logRho", "he4", "en"]` gets
  loaded as radius-and-time-dependent profile data. Pass your own list
  of MESA profile-column names if you want different (or additional)
  quantities visualized.
- `r_grid_name` (default `"mass"`): which MESA profile column is used as
  the radial coordinate profiles get resampled against. Change this if
  your run's profile data is better indexed by a different column (e.g.
  `"radius"`).
- `n_theta_rochelobe`/`n_phi_rochelobe` (binary only, default 128/24) and
  `n_x_potential`/`n_y_potential` (default 64/64): mesh/grid resolution
  for the Roche lobe shape and the Roche potential landscape,
  respectively. Higher looks smoother in Blender but costs more to
  compute and store.

## 7. Next step

Once you have a prepped output directory (or two, for a binary), open
Blender with the TULIPS-3D addon installed, point its sidebar at the
directory/directories, and click "Create" (single star) or "Create
binary". See the TULIPS-3D repo for how to install and use the addon
itself.
