"""Roche lobe geometry: the classical Roche-lobe equipotential surface for a
binary system, computed from the mass ratio via the standard "cast a ray from
the star's own center in every direction, root-find where the corotating-
frame Roche potential first equals its critical value at the L1 saddle
point" method (the same approach behind published Roche-lobe renderers, e.g.
https://github.com/seanseungbeomlee/MESA-Visualizer, which uses
PyAstronomy's rochepot_dl/get_lagrange_1 for the same computation).

All lengths here are dimensionless, in units of the orbital separation a
(=1), with the star whose lobe is being computed (mass M1) fixed at the
local origin and its companion (mass M2) at local (1, 0, 0). Mass ratio
q = M2/M1 (companion over self). To get the OTHER star's own lobe shape (in
its own local frame, with ITS companion at ITS local +x), call the same
functions again with q' = 1/q = M1/M2 - the potential is written generically
in terms of "the star at the origin" and "the companion at distance 1", so
this role-swap is exact, not an approximation.

No scipy dependency (this needs to also run inside Blender's bundled Python,
which doesn't ship scipy) - both root-finds use a fixed-iteration bisection,
vectorized entirely in numpy so it stays fast even across many MESA frames x
many angular directions at once.
"""
import numpy as np


def roche_potential(x, y, z, q):
    '''Dimensionless Roche potential in the corotating frame (G(M1+M2)=1,
    separation a=1, star 1/mass M1 at the origin, star 2/mass M2 at
    (1,0,0), q=M2/M1). Rotation axis is z, through the center of mass at
    x_com = q/(1+q) on the line joining the stars. Vectorized/broadcastable
    over x, y, z, q.'''
    r1 = np.sqrt(x**2 + y**2 + z**2)
    r2 = np.sqrt((x - 1.)**2 + y**2 + z**2)
    x_com = q / (1. + q)
    return -1./((1.+q)*r1) - q/((1.+q)*r2) - 0.5*((x - x_com)**2 + y**2)


def _dphi_dx_on_axis(x, q):
    '''d(Phi)/dx along the line joining the stars (y=z=0) - used to find the
    L1 point, the saddle point of the potential restricted to this line.'''
    x_com = q / (1. + q)
    return 1./((1.+q)*x**2) - q/((1.+q)*(1.-x)**2) - (x - x_com)


def find_L1(q, n_iter=60):
    '''Solves d(Phi)/dx = 0 for x in (0,1): the L1 point's position along the
    line joining the stars, as a fraction of the separation, measured from
    star 1 (mass M1, at x=0) toward star 2 (mass M2, at x=1). q may be a
    scalar or a numpy array (e.g. one value per MESA frame) - fully
    vectorized fixed-iteration bisection.

    d(Phi)/dx -> +inf as x->0+ and -> -inf as x->1-, and is monotonically
    decreasing in between for any physical q > 0, so plain bisection always
    converges to the unique root.'''
    q = np.atleast_1d(np.asarray(q, dtype=float))
    lo = np.full_like(q, 1e-6)
    hi = np.full_like(q, 1. - 1e-6)
    for _ in range(n_iter):
        mid = 0.5*(lo+hi)
        f_mid = _dphi_dx_on_axis(mid, q)
        lo = np.where(f_mid > 0., mid, lo)
        hi = np.where(f_mid > 0., hi, mid)
    return 0.5*(lo+hi)


def roche_lobe_radius_grid(q, n_theta=24, n_phi=24, n_iter=50, n_samples=200):
    '''For the star at the local origin (companion at local (+1,0,0), mass
    ratio q=M2/M1), casts a ray in each direction of an (n_theta x n_phi)
    grid - theta = polar angle from the +x/companion axis, in [0, pi]; phi =
    azimuthal angle about that axis, in [0, 2pi) - and finds the radius r at
    which the Roche potential along that ray FIRST equals the critical
    potential at L1 (moving outward from the star): the classical Roche
    lobe surface. q may be an array (one value per MESA frame); returns
    (r_grid, x_l1) where r_grid has shape (len(q), n_theta, n_phi) [radii in
    units of the separation a] and x_l1 has shape (len(q),).

    The upper bound is pinned to just past x_L1 itself (not a fixed fraction
    of the separation) - the classical Roche lobe's largest extent in ANY
    direction is along the line of centers (toward L1), so this is a safe,
    tight bound in every direction.

    Special case: theta=0 (exactly toward the companion) is degenerate - by
    construction Phi(x_L1,0,0)=Phi_L1 exactly, so x_L1 is a *tangent* point
    of Phi(r) along that one ray (dPhi/dx=0 there, that's the literal
    definition of L1), staying below Phi_L1 on both sides rather than
    crossing it. That ring of the grid (theta=0, all phi - a single
    degenerate point in continuous terms) is set directly to x_L1 rather
    than root-found.

    IMPORTANT (found via a real, reproducible rendering artifact - a
    zigzagging Roche-lobe tube right at its tip, only visible once
    n_theta was raised from 24 to 128 for smoothness): for theta NEAR but
    not at the pole, the near-tangency at theta=0 doesn't vanish
    immediately - it echoes as a narrow, easy-to-miss "bump" where
    Phi(r) pokes briefly *above* Phi_L1 in a tiny range around r~x_L1,
    then drops back below it, before FINALLY crossing for good further
    out is NOT what happens either - empirically (q~2.15, theta~1.4 deg)
    the bump sits with bound (Phi<Phi_L1) regions on BOTH sides of it, and
    the classical/intended lobe boundary is the FIRST crossing (moving
    outward from the star) into that bump, not any later one. A single
    blind bisection over a WIDE bracket (the original approach here)
    implicitly assumes exactly one clean sign change between its two
    endpoints - for these near-pole rays, both endpoints can land on the
    SAME (negative) side of the bump, silently violating that assumption
    and converging to something meaningless (empirically, the upper bound
    itself) instead of raising an error. This went unnoticed at the
    previous coarse n_theta=24 (theta step ~7.8 deg skipped straight over
    the affected near-pole band entirely) - purely a resolution-dependent
    latent bug, not something the earlier Eggleton-formula cross-check
    (agreement <0.5% across q=0.1-5.0, at that coarser resolution) could
    have caught.

    Fixed with a coarse-sample-then-refine approach (same spirit as
    equipotential_boundary_xy's fix for an unrelated but analogous
    "blind bisection isn't safe here" problem): densely sample Phi(r)
    along each ray first, take the FIRST sample-to-sample sign change from
    negative to positive (the actual physical lobe boundary - the
    innermost point where the ray first leaves the bound region), then
    bisect-refine within just that narrow bracket for precision. Costs
    more (n_samples evaluations per ray up front) but is robust regardless
    of how close a ray passes to the near-pole degeneracy.'''
    q = np.atleast_1d(np.asarray(q, dtype=float))
    x_l1 = find_L1(q)  # (n_frames,)
    phi_l1 = roche_potential(x_l1, 0., 0., q)  # critical potential, (n_frames,)

    theta = np.linspace(0., np.pi, n_theta)
    phi = np.linspace(0., 2.*np.pi, n_phi, endpoint=False)
    dx = np.broadcast_to(np.cos(theta)[:, None], (n_theta, n_phi))
    dy = np.broadcast_to(np.sin(theta)[:, None] * np.cos(phi)[None, :], (n_theta, n_phi))
    dz = np.broadcast_to(np.sin(theta)[:, None] * np.sin(phi)[None, :], (n_theta, n_phi))

    n_frames = q.shape[0]
    r_grid = np.zeros((n_frames, n_theta, n_phi))

    for f in range(n_frames):
        q_f = q[f]
        phi_l1_f = phi_l1[f]
        hi_bound = x_l1[f] * 1.05
        r_samples = np.linspace(1e-5, hi_bound, n_samples)  # (n_samples,)

        X = r_samples[None, None, :] * dx[:, :, None]  # (n_theta, n_phi, n_samples)
        Y = r_samples[None, None, :] * dy[:, :, None]
        Z = r_samples[None, None, :] * dz[:, :, None]
        F = roche_potential(X, Y, Z, q_f) - phi_l1_f

        sign = np.sign(F)
        crosses = (sign[:, :, :-1] < 0.) & (sign[:, :, 1:] >= 0.)  # FIRST negative-to-positive transition
        idx = np.arange(n_samples - 1)
        idx_grid = np.where(crosses, idx[None, None, :], n_samples)
        first_idx = np.min(idx_grid, axis=-1)  # (n_theta, n_phi)
        no_cross = first_idx == n_samples
        first_idx_safe = np.where(no_cross, 0, first_idx)

        lo = r_samples[first_idx_safe]
        hi = r_samples[np.minimum(first_idx_safe + 1, n_samples - 1)]
        lo = np.where(no_cross, x_l1[f], lo)  # fallback (shouldn't trigger outside theta=0, which gets overridden below anyway)
        hi = np.where(no_cross, x_l1[f], hi)

        f_lo = roche_potential(lo*dx, lo*dy, lo*dz, q_f) - phi_l1_f
        for _ in range(n_iter):
            mid = 0.5*(lo+hi)
            f_mid = roche_potential(mid*dx, mid*dy, mid*dz, q_f) - phi_l1_f
            same_side_as_lo = np.sign(f_mid) == np.sign(f_lo)
            lo = np.where(same_side_as_lo, mid, lo)
            f_lo = np.where(same_side_as_lo, f_mid, f_lo)
            hi = np.where(same_side_as_lo, hi, mid)
        r_grid[f] = 0.5*(lo+hi)

    # theta=0 pole fix (see docstring) - the degenerate ray toward the
    # companion, radius is exactly x_L1 regardless of phi.
    r_grid[:, 0, :] = x_l1[:, None]

    return r_grid, x_l1


def l4_l5_positions():
    '''Exact, mass-ratio-INDEPENDENT positions of L4/L5 (the third vertex of
    an equilateral triangle with the two masses - a classical exact result,
    no root-finding needed, unlike L1/L2/L3): star 1 at the origin, star 2
    at (1,0,0). Returns ((x_L4, y_L4), (x_L5, y_L5)).'''
    return (0.5, np.sqrt(3.) / 2.), (0.5, -np.sqrt(3.) / 2.)


def roche_potential_softened(x, y, z, q, softening):
    '''Same as roche_potential, but with a fixed length added in quadrature
    to each point mass's distance (r1, r2) - the standard "gravitational
    softening" trick (as in N-body simulations) for turning an unplottable
    1/r singularity at each star into a smooth, finite peak. Negligible
    effect more than a softening length or so away from either star - L1
    through L5 (all at least several tenths of the separation from the
    nearest star, for any realistic mass ratio) are essentially unaffected;
    verified numerically (softening=0.1 changes L4's potential by <0.4%).
    Used only for the full potential-SURFACE plot (a landscape
    visualization) - the actual Roche lobe boundary (roche_lobe_radius_grid)
    needs the real, unsoftened potential and must not use this.'''
    r1 = np.sqrt(x**2 + y**2 + z**2 + softening**2)
    r2 = np.sqrt((x - 1.)**2 + y**2 + z**2 + softening**2)
    x_com = q / (1. + q)
    return -1./((1.+q)*r1) - q/((1.+q)*r2) - 0.5*((x - x_com)**2 + y**2)


def _dphi_dx_general(x, q):
    '''d(Phi)/dx along the line joining the stars, valid EVERYWHERE except
    exactly at x=0 or x=1 (the stars themselves) - unlike _dphi_dx_on_axis
    (which hardcodes the sign convention for 0<x<1 only), this handles any
    x by explicitly tracking sign(x) and sign(x-1) (r1=|x|, r2=|x-1| for
    points on this axis, and d(1/|t|)/dt = -sign(t)/t^2). Used to find L2
    (x>1, beyond star 2) and L3 (x<0, beyond star 1) - _dphi_dx_on_axis's
    formula happens to match this exactly for 0<x<1 (sign(x)=+1,
    sign(x-1)=-1), which is why find_L1 doesn't need to change.'''
    x_com = q / (1. + q)
    return np.sign(x)/((1.+q)*x**2) + q*np.sign(x-1.)/((1.+q)*(x-1.)**2) - (x - x_com)


def _bisect_sign_aware(f, lo, hi, n_iter=80):
    '''Generic vectorized bisection that works regardless of which of f(lo)/
    f(hi) is positive (unlike find_L1's bisection, which hardcodes "f>0 on
    the lo side") - needed for find_L2/find_L3, where that sign convention
    flips depending on which side of the binary is being searched.'''
    lo = np.asarray(lo, dtype=float).copy()
    hi = np.asarray(hi, dtype=float).copy()
    f_lo = f(lo)
    for _ in range(n_iter):
        mid = 0.5*(lo+hi)
        f_mid = f(mid)
        same_side_as_lo = (np.sign(f_mid) == np.sign(f_lo))
        lo = np.where(same_side_as_lo, mid, lo)
        f_lo = np.where(same_side_as_lo, f_mid, f_lo)
        hi = np.where(same_side_as_lo, hi, mid)
    return 0.5*(lo+hi)


def find_L2(q, n_iter=80):
    '''Solves d(Phi)/dx = 0 for x>1: the outer Lagrange point beyond star 2
    (mass M2, at x=1), along the line joining the stars. q may be scalar or
    array. Bracket (1+eps, 51) is generous enough for any q in a realistic
    stellar-binary range (verified q=0.1-10); dPhi/dx runs from +inf (at
    x=1+) to -inf (as x->inf), monotonically, so bisection converges to the
    unique root.'''
    q = np.atleast_1d(np.asarray(q, dtype=float))
    lo = np.full_like(q, 1. + 1e-6)
    hi = np.full_like(q, 51.)
    return _bisect_sign_aware(lambda x: _dphi_dx_general(x, q), lo, hi, n_iter)


def find_L3(q, n_iter=80):
    '''Solves d(Phi)/dx = 0 for x<0: the outer Lagrange point beyond star 1
    (mass M1, at the origin), along the line joining the stars. Mirror
    image of find_L2 (bracket (-51, -eps); dPhi/dx runs from +inf as
    x->-inf to -inf at x=0-, so bisection converges the same way).'''
    q = np.atleast_1d(np.asarray(q, dtype=float))
    lo = np.full_like(q, -51.)
    hi = np.full_like(q, -1e-6)
    return _bisect_sign_aware(lambda x: _dphi_dx_general(x, q), lo, hi, n_iter)


def _flood_fill_bounded(mask, seed):
    '''Flood-fills the connected component of `mask` (a 2D boolean array)
    containing `seed` (a (row,col) tuple) via iterative BFS/DFS (no scipy -
    same "must also run in Blender's bundled Python" constraint as the rest
    of this module, though in practice this only ever runs in
    DataPrepTulips3D's own environment). Returns (visited, touches_edge) -
    touches_edge is True if the component reaches the array boundary
    (meaning the region is NOT safely bounded within the grid - the caller
    should widen the domain and retry rather than trust the result).'''
    ny, nx = mask.shape
    visited = np.zeros_like(mask, dtype=bool)
    stack = [seed]
    visited[seed] = True
    touches_edge = False
    while stack:
        j, i = stack.pop()
        if j == 0 or j == ny-1 or i == 0 or i == nx-1:
            touches_edge = True
        for dj, di in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            nj, ni = j+dj, i+di
            if 0 <= nj < ny and 0 <= ni < nx and mask[nj, ni] and not visited[nj, ni]:
                visited[nj, ni] = True
                stack.append((nj, ni))
    return visited, touches_edge


def _moore_boundary_trace(mask):
    '''Traces the outer boundary of a single connected True-region in a 2D
    boolean grid as an ordered list of (row, col) pixel coordinates, via the
    standard Moore-neighbor tracing algorithm (walk the boundary in a fixed
    rotational order until returning to the start). Grid-resolution
    ("staircase") accuracy - not a smooth sub-pixel contour - which is fine
    for a wireframe visualization once resampled to a modest point count.
    Returns [] if `mask` has no True pixels.'''
    ny, nx = mask.shape
    start = None
    for j in range(ny):
        idx = np.nonzero(mask[j])[0]
        if len(idx):
            start = (j, int(idx[0]))
            break
    if start is None:
        return []

    # 8-connected offsets, clockwise starting "north"
    directions = ((-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1))

    def is_in(j, i):
        return 0 <= j < ny and 0 <= i < nx and mask[j, i]

    boundary = [start]
    current = start
    b_dir = 6  # "west" - nothing can be True to the west of the first True pixel in row-major scan
    max_steps = 6 * nx * ny  # safety cap, not expected to be hit for a well-formed single blob
    for _ in range(max_steps):
        found = False
        for k in range(1, 9):
            d = (b_dir + k) % 8
            dj, di = directions[d]
            nj, ni = current[0]+dj, current[1]+di
            if is_in(nj, ni):
                boundary.append((nj, ni))
                b_dir = (d + 5) % 8  # resume search from just past "where we came from"
                current = (nj, ni)
                found = True
                break
        if not found:
            break  # isolated single pixel - degenerate, avoid spinning
        if current == start and len(boundary) > 2:
            boundary.pop()
            break
    return boundary


def equipotential_boundary_xy(q, phi_target, x_range, y_range, seed_xy, n_grid=320, n_theta=48, smooth_window=9):
    '''For ONE frame (q, phi_target scalars): finds the closed boundary
    curve of the connected region {Phi(x,y,0;q) > phi_target} that contains
    seed_xy (a point known to lie safely inside that region for this
    level - e.g. an L4/L5 position, which is the domain's own potential
    maximum away from the two point-mass singularities, so it is "inside"
    for any level below Phi(L4)).

    IMPORTANT: `_moore_boundary_trace` walks pixel-to-pixel on the (n_grid x
    n_grid) mask, so its raw output is a "staircase" at grid-cell
    resolution (Ben, after the first version: "this still does not look
    good at all" - correctly; a screenshot showed clearly faceted/blocky
    rings). Resampling to more points by arc length (n_theta below) does
    NOT fix this - it just places more points along the same staircase,
    since it interpolates linearly between the already-blocky pixel
    boundary rather than smoothing it. Two real fixes are needed together:
    a much finer underlying grid (n_grid, up from an initial 140 - directly
    shrinks each stair-step) AND an explicit circular moving-average smooth
    on the final resampled points (removes the residual high-frequency
    zigzag a finite grid can't avoid, however fine). Neither alone was
    enough in testing at the previous n_grid=140/no-smoothing settings -
    the L1 Roche-lobe tube (make_roche_lobe_curve, a smooth analytic
    theta-parameterization, no grid/staircase involved at all) never had
    this problem, confirmed by Ben ("the L1 looks good though") - it's
    specific to this grid-tracing method.

    Why this (flood-fill + boundary-trace), not a simpler ray-cast from a
    fixed center: verified numerically that the merged equipotential around
    both stars is NOT star-shaped from either the center of mass or the
    geometric midpoint for realistic (non-unity) mass ratios - a plain
    "cast a ray and root-find the first/last crossing" approach silently
    misses large arcs of the boundary (confirmed against several q spanning
    0.1-10, using this project's actual MESA-derived mass-ratio range).
    Flood-fill from a guaranteed-inside seed + boundary trace has no such
    assumption and was verified correct (bounded, non-edge-touching, sane
    area) across the same q range.

    Meant for the L2/L3 "outer" equipotentials specifically: phi_target
    should be JUST BELOW (more negative than) Phi at the L2 or L3 point
    itself (a small backoff, e.g. Phi_L2 minus ~2-5% of |Phi_L2-Phi_L1| -
    the EXACT critical value is itself a degenerate/pinched contour, same
    reasoning as why roche_lobe_radius_grid needs the theta=0 pole
    workaround at L1) - not used for the L1 lobes themselves (those already
    have their own, more precise, per-star ray-cast method).

    Returns (n_theta, 2) array of (x,y) points forming a closed loop,
    resampled to n_theta points by ARC LENGTH (not by the raw, unevenly-
    spaced grid-boundary pixel count) for even spacing regardless of grid
    resolution. Raises ValueError if the region touches the grid edge
    (caller should widen x_range/y_range and retry) or is empty.'''
    xs = np.linspace(x_range[0], x_range[1], n_grid)
    ys = np.linspace(y_range[0], y_range[1], n_grid)
    X, Y = np.meshgrid(xs, ys)
    Phi = roche_potential(X, Y, 0., q)
    mask = Phi > phi_target

    j0 = int(np.argmin(np.abs(ys - seed_xy[1])))
    i0 = int(np.argmin(np.abs(xs - seed_xy[0])))
    if not mask[j0, i0]:
        raise ValueError(f"equipotential_boundary_xy: seed {seed_xy} is not inside "
                          f"the target region (Phi(seed)={Phi[j0,i0]:.4f} <= phi_target={phi_target:.4f})")

    visited, touches_edge = _flood_fill_bounded(mask, (j0, i0))
    if touches_edge:
        raise ValueError("equipotential_boundary_xy: region touches the grid edge - "
                          "widen x_range/y_range and retry")

    boundary_px = _moore_boundary_trace(visited)
    if len(boundary_px) < 3:
        raise ValueError("equipotential_boundary_xy: degenerate (near-empty) boundary")

    pts = np.array([(xs[i], ys[j]) for j, i in boundary_px])  # (n_px, 2), grid-resolution accuracy

    # Resample to n_theta points, evenly spaced by arc length around the closed loop.
    closed = np.vstack([pts, pts[0]])
    seg_lengths = np.sqrt(np.sum(np.diff(closed, axis=0)**2, axis=1))
    cum_len = np.concatenate([[0.], np.cumsum(seg_lengths)])
    total_len = cum_len[-1]
    target_len = np.linspace(0., total_len, n_theta, endpoint=False)
    x_resampled = np.interp(target_len, cum_len, closed[:, 0])
    y_resampled = np.interp(target_len, cum_len, closed[:, 1])
    resampled = np.stack([x_resampled, y_resampled], axis=1)

    # Circular moving-average smooth - removes the residual grid-staircase
    # zigzag that arc-length resampling alone doesn't (see the docstring
    # above). Wraps around the closed loop (np.roll, not a clipped/edge-
    # padded window) since there's no start/end to the ring. smooth_window
    # should stay well below n_theta (a handful of points either side, not
    # a large fraction of the loop) so this blurs out grid noise without
    # eating genuine shape detail.
    if smooth_window and smooth_window > 1:
        half = smooth_window // 2
        smoothed = np.zeros_like(resampled)
        for shift in range(-half, half + 1):
            smoothed += np.roll(resampled, shift, axis=0)
        resampled = smoothed / (2 * half + 1)

    return resampled


def roche_potential_surface_grid(q, x_range, y_range, n_x=48, n_y=48, softening=0.1):
    '''Evaluates the (softened) Roche potential Phi(x,y,0;q) over a regular
    Cartesian grid in the orbital plane - a "height field" of the whole
    potential landscape (as opposed to roche_lobe_radius_grid's single
    equipotential contour), for a literal 3D surface-plot visualization
    (e.g. matplotlib-style wireframe renders of the Roche potential). q may
    be an array (one value per MESA frame). Returns (phi_grid, xs, ys)
    where phi_grid has shape (len(q), n_y, n_x) and xs/ys (length n_x/n_y)
    are the fixed (frame-independent) grid coordinates, in units of the
    separation a.

    x_range/y_range should extend a bit past L4/L5 (l4_l5_positions()) to
    show the full landscape including the outer Lagrange points, not just
    the immediate vicinity of the two stars.'''
    q = np.atleast_1d(np.asarray(q, dtype=float))
    xs = np.linspace(x_range[0], x_range[1], n_x)
    ys = np.linspace(y_range[0], y_range[1], n_y)
    X, Y = np.meshgrid(xs, ys)  # (n_y, n_x)
    Phi = roche_potential_softened(X[None, :, :], Y[None, :, :], 0., q[:, None, None], softening)
    return Phi, xs, ys
