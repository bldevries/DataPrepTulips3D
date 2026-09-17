"""Particle-based mass-transfer stream: real test-particle trajectories in
the same dimensionless corotating Roche frame roche_lobe.py already uses
everywhere else in this project (donor/star 1 at the local origin, companion/
star 2 at (1,0,0), G(M1+M2)=1, separation a=1, rotation angular velocity
Omega=1 about the z axis through x_com=q/(1+q); q=M2/M1).

Ben: "I want the donor to emit particles towards the accretor... The
dynamics of the particles should be governed by the potential you
calculated." The equations of motion are the standard
circular-restricted-three-body-problem (CR3BP) equations for a massless
test particle:

    x'' = 2*y' - dPhi/dx   (Coriolis term, only when include_coriolis=True)
    y'' = -2*x' - dPhi/dy  (Coriolis term, only when include_coriolis=True)
    z'' = -dPhi/dz

where Phi = roche_lobe.roche_potential(x, y, z, q) (re-used directly, not
re-derived). dPhi/d{x,y,z} is the exact analytic gradient (derived and
cross-checked term-by-term against roche_lobe._dphi_dx_on_axis, not a
numerical finite-difference, since this needs to run cheaply inside a
tight RK4 loop over many steps x many particles): dPhi/dx =
x/((1+q)*r1**3) + q*(x-1)/((1+q)*r2**3) - (x-x_com), and similarly for y
(extra "-y" term from the centrifugal part) and z (no centrifugal
contribution, since the centrifugal term only involves x,y).

The Coriolis term is OFF by default (include_coriolis=False everywhere in
this module). It's not merely omitted for a shortcut: a particle actually
MOVING through the corotating frame (unlike the STATIC potential
landscape/equipotential surfaces built elsewhere in this project) does
physically feel it, and an earlier version of this module had it on by
default - real Roche-lobe-overflow streams curve and typically overshoot
the companion before looping back into a disk (the textbook "hook" shape)
because of this term, and turning it on reproduces exactly that (see the
prototype scripts/renders this module's git history references). But it
also made capture time HIGHLY sensitive to mass ratio q (fine at this
project's own q~0.625, but taking 5-10x longer, or needing a much wider -
and visually worse, rosette-destroying - capture radius, at other q even
within this SAME project's own dataset's range) - real, but non-trivial,
complexity. Ben, after seeing both versions: "let's take out the Coriolis
force, that should simplify things. The first thing I want is to show the
mass transfer between the stars. We can look at a more realistic particle
orbit later if need be." Without it, particles fall close to straight
from L1 toward the companion, and capture is fast and robust across
q=0.5-1.6 (verified directly) with no per-q tuning at all - a deliberately
simpler v1. include_coriolis stays a real, threaded-through parameter
(not deleted) specifically so the more physically realistic version can
be switched back on later without redoing this integrator.

Deliberately general, not hardcoded to "L1 -> companion": Ben also wants a
future stellar-wind emitter ("shoot particles in all directions... not
only towards the accretor"). generate_particle_pool() takes the particles'
initial positions/velocities as plain arrays, so a wind emitter is just a
different caller (isotropic launch points over a star's surface, radial
initial velocities) - none of the integration/capture/escape logic below
needs to change.

Time here is in the same dimensionless units as everywhere else in this
project (1/Omega; one full orbit = 2*pi). Quasi-static approximation
(deliberate, matching how every other per-frame quantity in this project -
the potential landscape, equipotential rings, Roche lobe shapes - is
already computed): a particle pool is generated using ONE MESA frame's
(q, mdot, ...) held fixed for the pool's whole (short, at most a few
orbits) flight time, since real stream/wind flight times are utterly
negligible next to the nuclear-evolution timescale between MESA frames
this project resamples to (t_resolution~1000 frames spanning a whole
stellar lifetime) - the same assumption already implicit in treating each
frame's potential landscape as its own static equilibrium picture rather
than integrating a genuinely time-dependent 3-body problem across frames.
"""
import numpy as np

from . import roche_lobe


def phi_gradient(x, y, z, q):
    '''Exact analytic gradient of roche_lobe.roche_potential (see module
    docstring for the derivation/cross-check). Vectorized/broadcastable
    over x, y, z, q exactly like roche_potential itself.'''
    r1 = np.sqrt(x**2 + y**2 + z**2)
    r2 = np.sqrt((x - 1.)**2 + y**2 + z**2)
    r1 = np.maximum(r1, 1e-8)  # avoid true division-by-zero at the exact stellar center
    r2 = np.maximum(r2, 1e-8)
    x_com = q / (1. + q)
    dphidx = x / ((1. + q) * r1**3) + q * (x - 1.) / ((1. + q) * r2**3) - (x - x_com)
    dphidy = y / ((1. + q) * r1**3) + q * y / ((1. + q) * r2**3) - y
    dphidz = z / ((1. + q) * r1**3) + q * z / ((1. + q) * r2**3)
    return dphidx, dphidy, dphidz


def _accel(state, q, include_coriolis=False):
    '''state: (..., 6) array of [x, y, z, vx, vy, vz]. Returns d(state)/dt,
    i.e. [vx, vy, vz, ax, ay, az] - the CR3BP equations of motion (gravity
    via phi_gradient, +Coriolis when include_coriolis, Omega=1). Vectorized
    over any leading shape.

    include_coriolis=False (Ben: "let's take out the Coriolis force, that
    should simplify things... the first thing I want is to show the mass
    transfer between the stars. We can look at a more realistic particle
    orbit later if need be") drops the two Coriolis terms below, leaving
    gravity + centrifugal only - i.e. particles just roll downhill on Phi,
    with inertia. This isn't only simpler to reason about: it's also the
    direct fix for the capture-radius/q-sensitivity tuning fight upstream
    of this change (see git history / mass_transfer_stream memory) - that
    entire problem came from the Coriolis-induced looping making capture a
    near-chaotic function of q; without it, particles fall roughly
    straight from L1 toward the companion and capture becomes robust
    across mass ratio with no per-q tuning at all (verified below).'''
    x, y, z, vx, vy, vz = (state[..., i] for i in range(6))
    dphidx, dphidy, dphidz = phi_gradient(x, y, z, q)
    if include_coriolis:
        ax = 2. * vy - dphidx
        ay = -2. * vx - dphidy
    else:
        ax = -dphidx
        ay = -dphidy
    az = -dphidz
    return np.stack([vx, vy, vz, ax, ay, az], axis=-1)


def integrate_particles(x0, y0, z0, vx0, vy0, vz0, q, dt=0.01, n_steps=2000,
                         capture_target="companion", capture_radius=0.02,
                         donor_radius=0.02, escape_radius=3.0, include_coriolis=False):
    '''Integrates a whole pool of test particles simultaneously via RK4,
    each until it's captured or escapes or n_steps runs out. All of
    x0/y0/z0/vx0/vy0/vz0 are 1D arrays of the same length N (one entry per
    particle) - a plain array-in, array-out API so the caller decides how
    particles are seeded (an L1 stream today; a whole-surface wind later -
    see module docstring).

    q may be a single scalar (one CR3BP frame - the quasi-static-per-frame
    assumption in the module docstring) OR an (N,) array giving each
    particle its OWN q - this is what lets bake_particle_pool below
    integrate an entire multi-MESA-frame pool (each particle tagged with
    its own birth frame's mass ratio) in a single vectorized call instead
    of one Python-level loop iteration per frame.

    capture_target: "companion" (the usual mass-transfer-stream case - the
    accretor sits at local (1,0,0)) or "donor" (for a future use where the
    "companion" role in this local frame is actually the mass-losing star,
    e.g. computing wind capture by a nearby star). capture_radius is
    dimensionless (units of the separation a=1) - Ben: "make this an
    option" - callers should pass capture_radius = (desired_radius_Rsun /
    separation_Rsun) for "hits a set radius", or the target's own physical
    radius converted the same way for "hits the star itself".

    Default 0.02 (with include_coriolis=False - see that flag's own
    docstring above). An EARLIER version of this function defaulted to
    include_coriolis=True, and capture time under that model turned out to
    be highly q-sensitive - fine (and produced a nice multi-loop rosette
    before capture) at this project's own q~0.625, but effectively never
    captured within any reasonable window at some other q even within
    this SAME project's own active dataset's range, forcing a much wider
    (and visually worse - it killed the rosette) capture_radius=0.04 to
    stay robust. Ben: "let's take out the Coriolis force... the first
    thing I want is to show the mass transfer between the stars" - without
    it, particles fall close to straight from L1 to the companion and
    capture is robust across q=0.5-1.6 (verified directly) even back at
    this original tight 0.02 radius, with no per-q tuning needed at all.

    donor_radius is always checked too (regardless of capture_target) -
    particles that fall back onto the star they were launched from are
    just destroyed, not tracked as a separate "capture" statistic, purely
    to avoid the 1/r**2 gravity term blowing up near a point mass when a
    particle grazes back past its own donor.

    escape_radius is measured from the system's true center of mass
    (x_com, 0, 0), not from the local origin (star 1) - Ben: "if the
    particle flies away from the system, it should be destroyed at ~3x
    the separation" - default 3.0 matches that.

    Returns a dict:
      trajectory : (N, n_steps+1, 3) array of (x,y,z), frozen (repeats the
                   last live position) after each particle's death - so
                   plain per-frame indexing by step number always gives a
                   sane position even for already-dead particles.
      alive_until : (N,) int array - the step index at which each particle
                    died (n_steps if it survived the whole integration
                    window without dying - "still flying" at the end,
                    which callers should probably treat as an escape too).
      death_reason : (N,) array of strings: "captured", "fell_back",
                     "escaped", or "survived".
    '''
    n = len(x0)
    q = np.broadcast_to(np.asarray(q, dtype=float), (n,)).copy()
    x_com = q / (1. + q)  # (N,) - each particle's own CM x-position, for the escape check below
    state = np.stack([
        np.asarray(x0, dtype=float), np.asarray(y0, dtype=float), np.asarray(z0, dtype=float),
        np.asarray(vx0, dtype=float), np.asarray(vy0, dtype=float), np.asarray(vz0, dtype=float),
    ], axis=-1)  # (N, 6)

    traj = np.zeros((n, n_steps + 1, 3))
    traj[:, 0, :] = state[:, :3]
    alive = np.ones(n, dtype=bool)
    alive_until = np.full(n, n_steps, dtype=int)
    death_reason = np.full(n, "survived", dtype=object)

    target_pos = np.array([1., 0., 0.]) if capture_target == "companion" else np.array([0., 0., 0.])
    donor_pos = np.array([0., 0., 0.]) if capture_target == "companion" else np.array([1., 0., 0.])

    for step in range(1, n_steps + 1):
        # classic RK4, vectorized over all still-alive particles at once
        k1 = _accel(state, q, include_coriolis)
        k2 = _accel(state + 0.5 * dt * k1, q, include_coriolis)
        k3 = _accel(state + 0.5 * dt * k2, q, include_coriolis)
        k4 = _accel(state + dt * k3, q, include_coriolis)
        new_state = state + (dt / 6.) * (k1 + 2. * k2 + 2. * k3 + k4)
        state = np.where(alive[:, None], new_state, state)

        pos = state[:, :3]
        d_target = np.linalg.norm(pos - target_pos, axis=-1)
        d_donor = np.linalg.norm(pos - donor_pos, axis=-1)
        cm_pos = np.stack([x_com, np.zeros_like(x_com), np.zeros_like(x_com)], axis=-1)  # (N, 3)
        d_cm = np.linalg.norm(pos - cm_pos, axis=-1)

        newly_captured = alive & (d_target < capture_radius)
        newly_fellback = alive & ~newly_captured & (d_donor < donor_radius)
        newly_escaped = alive & ~newly_captured & ~newly_fellback & (d_cm > escape_radius)

        death_reason[newly_captured] = "captured"
        death_reason[newly_fellback] = "fell_back"
        death_reason[newly_escaped] = "escaped"
        just_died = newly_captured | newly_fellback | newly_escaped
        alive_until[just_died] = step
        alive = alive & ~just_died

        traj[:, step, :] = np.where(alive[:, None], pos, traj[:, step - 1, :])

        if not alive.any():
            traj[:, step + 1:, :] = traj[:, step:step + 1, :]
            break

    return {"trajectory": traj, "alive_until": alive_until, "death_reason": death_reason}


DEFAULT_LAUNCH_SPEED = 0.02  # Ben, after seeing both tested speeds plotted: "start with a slow
# launch seed so that particles get captured (through a rosette) onto the accretor" - 0.08 sent
# most particles into a wide Coriolis-deflected loop that rains back onto the DONOR instead
# (a real CR3BP effect, not a bug, but not the intended visual); 0.02 reliably produces the
# clean rosette-onto-accretor look. See prototype_mass_transfer_particles.py's saved plot.


def launch_from_l1(q, n_particles, speed=0.02, spread=0.15, seed=None):
    '''Seeds a pool of particles at the L1 point with a small kick toward
    the companion (real Roche-lobe overflow leaves through L1 near rest in
    the corotating frame - a tiny push is needed to actually start the
    fall, matching e.g. a small thermal/gas-pressure velocity at the
    donor's photosphere) plus a small random transverse spread so the
    pool isn't a single degenerate line. speed/spread are dimensionless
    (units of the local circular orbital speed, since v=1 there in these
    units). Returns (x0, y0, z0, vx0, vy0, vz0), each length n_particles.'''
    rng = np.random.default_rng(seed)
    x_l1 = float(roche_lobe.find_L1(np.array([q]))[0])
    x0 = np.full(n_particles, x_l1)
    y0 = rng.normal(0., 1e-4, n_particles)  # negligible - just breaks exact degeneracy
    z0 = rng.normal(0., 1e-4, n_particles)
    vx0 = np.full(n_particles, speed)
    vy0 = rng.normal(0., speed * spread, n_particles)
    vz0 = rng.normal(0., speed * spread, n_particles)
    return x0, y0, z0, vx0, vy0, vz0


def bake_particle_pool(q_per_frame, n_particles_per_frame=30, dt=0.0005, n_local_steps=3000,
                        speed=DEFAULT_LAUNCH_SPEED, spread=0.15,
                        capture_radius=0.02, donor_radius=0.10, escape_radius=3.0, seed=0,
                        include_coriolis=False):
    '''Bakes a per-MESA-frame pool of candidate mass-transfer-stream
    trajectories, one launched-from-L1 pool per frame using THAT frame's
    own q (see module docstring's quasi-static-per-frame assumption) -
    all frames' pools integrated together in a single vectorized call
    (via integrate_particles' per-particle q support), not a Python loop
    over frames.

    q_per_frame : (n_frames,) array, one mass ratio (q1 = M2/M1) per MESA
    frame - e.g. mesa_data_interp.py's own q1 array.

    n_particles_per_frame is a fixed BAKED pool size, not how many
    particles are actually visible at any moment - that's meant to be a
    live, mdot-scaled fraction of this pool at Blender runtime (so the
    "proportionality constant" setting stays live-adjustable with no
    rebaking - see this module's own docstring / the design discussion
    this came out of). Bigger pools give smoother live thinning at the
    cost of more baked data - see bake cost measurements before picking a
    production value.

    Returns a dict:
      trajectory  : (n_frames, n_particles_per_frame, n_local_steps+1, 3)
      alive_until : (n_frames, n_particles_per_frame) int
      death_reason: (n_frames, n_particles_per_frame) object array of str
      local_dt    : dt (echoed back - callers need it to convert local
                    step index to dimensionless flight time)
    '''
    q_per_frame = np.asarray(q_per_frame, dtype=float)
    n_frames = len(q_per_frame)
    n = n_frames * n_particles_per_frame

    x0 = np.empty(n); y0 = np.empty(n); z0 = np.empty(n)
    vx0 = np.empty(n); vy0 = np.empty(n); vz0 = np.empty(n)
    q_flat = np.empty(n)
    for f in range(n_frames):
        sl = slice(f * n_particles_per_frame, (f + 1) * n_particles_per_frame)
        frame_seed = None if seed is None else (seed, f)
        xf, yf, zf, vxf, vyf, vzf = launch_from_l1(
            q_per_frame[f], n_particles_per_frame, speed=speed, spread=spread, seed=frame_seed)
        x0[sl], y0[sl], z0[sl] = xf, yf, zf
        vx0[sl], vy0[sl], vz0[sl] = vxf, vyf, vzf
        q_flat[sl] = q_per_frame[f]

    result = integrate_particles(x0, y0, z0, vx0, vy0, vz0, q_flat, dt=dt, n_steps=n_local_steps,
                                  capture_target="companion", capture_radius=capture_radius,
                                  donor_radius=donor_radius, escape_radius=escape_radius,
                                  include_coriolis=include_coriolis)

    shape4 = (n_frames, n_particles_per_frame, n_local_steps + 1, 3)
    shape2 = (n_frames, n_particles_per_frame)
    return {
        "trajectory": result["trajectory"].reshape(shape4),
        "alive_until": result["alive_until"].reshape(shape2),
        "death_reason": result["death_reason"].reshape(shape2),
        "local_dt": dt,
    }
