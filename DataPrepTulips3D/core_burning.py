"""Reference ignition curves for the classic MESA "TRho_Profile" burn-region
diagram - copied verbatim (rounded to 4 decimals) from MESA's own bundled
data/star_data/plot_info/{hydrogen,helium,carbon,oxygen}_burn.data files (MESA
24.08.1). Each curve gives the ignition boundary as logT(logRho): the star's own
core (log_center_Rho, log_center_T) trajectory over its lifetime crosses these
curves in order (H then He then C then O, by increasing ignition temperature)
- exactly the curves MESA's own pgstar TRho_Profile plot draws as fixed dashed
reference lines to let a user visually judge "what stage is this star at now".
This module automates that same visual comparison into a discrete label.
"""

HYDROGEN_BURN_LOGRHO = [4.7897, 4.7819, 4.7729, 4.7633, 4.7545, 4.7434, 4.7376, 4.7275, 4.713, 4.705, 4.6934, 4.6813, 4.6668, 4.6502, 4.637, 4.6235, 4.6097, 4.5955, 4.581, 4.5643, 4.5411, 4.525, 4.5066, 4.4915, 4.4706, 4.4535, 4.435, 4.4078, 4.3825, 4.364, 4.345, 4.3123, 4.2875, 4.2677, 4.2415, 4.2068, 4.186, 4.1437, 4.12, 4.0822, 4.0566, 4.03, 3.9891, 3.9466, 3.92, 3.873, 3.8484, 3.805, 3.76, 3.7282, 3.6896, 3.6539, 3.6037, 3.5531, 3.5177, 3.4763, 3.4245, 3.3723, 3.3317, 3.274, 3.22, 3.1728, 3.0916, 3.04, 2.9758, 2.9154, 2.8611, 2.815, 2.77, 2.725, 2.68, 2.635, 2.59, 2.5476, 2.5, 2.455, 2.41, 2.3674, 2.32, 2.275, 2.23, 2.1838, 2.14, 2.095, 2.05, 1.9938, 1.9477, 1.895, 1.8668, 1.8034, 1.7635, 1.7109, 1.6772, 1.6145, 1.58, 1.5278, 1.4902, 1.4286, 1.396, 1.3459, 1.3067, 1.2685, 1.2171, 1.1657, 1.127, 1.089, 1.0379, 0.9868, 0.9497, 0.91, 0.8586, 0.8067, 0.7678, 0.7277, 0.6738, 0.6186, 0.562, 0.5031, 0.43, 0.3489, 0.2566, 0.1787, 0.07, -0.0575, -0.1643, -0.3085, -0.4542, -0.601, -0.7488, -0.8981, -1.0692, -1.2319, -1.3811, -1.55, -1.7167, -1.8703, -2.5]
HYDROGEN_BURN_LOGT = [6.0, 6.01, 6.02, 6.03, 6.04, 6.05, 6.06, 6.07, 6.08, 6.09, 6.1, 6.11, 6.12, 6.13, 6.14, 6.15, 6.16, 6.17, 6.18, 6.19, 6.2, 6.21, 6.22, 6.23, 6.24, 6.25, 6.26, 6.27, 6.28, 6.29, 6.3, 6.31, 6.32, 6.33, 6.34, 6.35, 6.36, 6.37, 6.38, 6.39, 6.4, 6.41, 6.42, 6.43, 6.44, 6.45, 6.46, 6.47, 6.48, 6.49, 6.5, 6.51, 6.52, 6.53, 6.54, 6.55, 6.56, 6.57, 6.58, 6.59, 6.6, 6.61, 6.62, 6.63, 6.64, 6.65, 6.66, 6.67, 6.68, 6.69, 6.7, 6.71, 6.72, 6.73, 6.74, 6.75, 6.76, 6.77, 6.78, 6.79, 6.8, 6.81, 6.82, 6.83, 6.84, 6.85, 6.86, 6.87, 6.88, 6.89, 6.9, 6.91, 6.92, 6.93, 6.94, 6.95, 6.96, 6.97, 6.98, 6.99, 7.0, 7.01, 7.02, 7.03, 7.04, 7.05, 7.06, 7.07, 7.08, 7.09, 7.1, 7.11, 7.12, 7.13, 7.14, 7.15, 7.16, 7.17, 7.18, 7.19, 7.2, 7.21, 7.22, 7.23, 7.24, 7.25, 7.26, 7.27, 7.28, 7.29, 7.3, 7.31, 7.32, 7.33, 7.34, 7.35, 7.36]

HELIUM_BURN_LOGRHO = [8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4835, 8.4841, 8.4845, 8.485, 8.4855, 8.486, 8.4863, 8.4868, 8.4873, 8.4878, 8.4882, 8.4887, 8.4891, 8.4896, 8.4901, 8.4906, 8.4911, 8.4912, 8.4917, 8.4922, 8.4927, 8.4931, 8.4935, 8.494, 8.4944, 8.4949, 8.4953, 8.4958, 8.4958, 8.4962, 8.4967, 8.4971, 8.4975, 8.4979, 8.4983, 8.4987, 8.4991, 8.4995, 8.4999, 8.5, 8.5003, 8.5007, 8.5011, 8.5014, 8.5018, 8.5021, 8.5025, 8.5028, 8.5032, 8.5035, 8.5036, 8.5039, 8.5042, 8.5045, 8.5048, 8.5051, 8.5054, 8.5057, 8.5059, 8.5062, 8.5065, 8.5065, 8.5068, 8.507, 8.5073, 8.5075, 8.5077, 8.5079, 8.5081, 8.5083, 8.5085, 8.5087, 8.5088, 8.509, 8.5091, 8.5093, 8.5094, 8.5096, 8.5097, 8.5098, 8.5099, 8.51, 8.5101, 8.5101, 8.5101, 8.5101, 8.5099, 8.5083, 8.5043, 8.494, 8.4798, 8.4585, 8.44, 8.4173, 8.4008, 8.38, 8.3454, 8.3262, 8.2999, 8.2683, 8.2362, 8.2076, 8.1801, 8.14, 8.1128, 8.08, 8.0431, 8.0163, 7.9716, 7.9379, 7.8949, 7.8519, 7.8212, 7.769, 7.72, 7.6846, 7.6352, 7.6, 7.5528, 7.4972, 7.4597, 7.4077, 7.36, 7.3062, 7.2773, 7.225, 7.18, 7.12, 7.06, 6.9915, 6.8966, 6.8253, 6.7415, 6.6424, 6.5636, 6.4441, 6.3194, 6.1898, 6.04, 5.86, 5.6655, 5.4521, 5.2538, 5.051, 4.8449, 4.6491, 4.4614, 4.2638, 4.06, 3.8507, 3.663, 3.4757, 3.28, 3.094, 2.9112, 2.7337, 2.5505, 2.38, 2.2258, 2.0465, 1.9, 1.7389, 1.5874, 1.43, 1.2955, 1.1385, 1.0, 0.8678, 0.7404, 0.6, 0.4885, 0.3639, 0.2404, 0.118]
HELIUM_BURN_LOGT = [6.25, 6.26, 6.27, 6.28, 6.29, 6.3, 6.31, 6.32, 6.33, 6.34, 6.35, 6.36, 6.37, 6.38, 6.39, 6.4, 6.41, 6.42, 6.43, 6.44, 6.45, 6.46, 6.47, 6.48, 6.49, 6.5, 6.51, 6.52, 6.53, 6.54, 6.55, 6.56, 6.57, 6.58, 6.59, 6.6, 6.61, 6.62, 6.63, 6.64, 6.65, 6.66, 6.67, 6.68, 6.69, 6.7, 6.71, 6.72, 6.73, 6.74, 6.75, 6.76, 6.77, 6.78, 6.79, 6.8, 6.81, 6.82, 6.83, 6.84, 6.85, 6.86, 6.87, 6.88, 6.89, 6.9, 6.91, 6.92, 6.93, 6.94, 6.95, 6.96, 6.97, 6.98, 6.99, 7.0, 7.01, 7.02, 7.03, 7.04, 7.05, 7.06, 7.07, 7.08, 7.09, 7.1, 7.11, 7.12, 7.13, 7.14, 7.15, 7.16, 7.17, 7.18, 7.19, 7.2, 7.21, 7.22, 7.23, 7.24, 7.25, 7.26, 7.27, 7.28, 7.29, 7.3, 7.31, 7.32, 7.33, 7.34, 7.35, 7.36, 7.37, 7.38, 7.39, 7.4, 7.41, 7.42, 7.43, 7.44, 7.45, 7.46, 7.47, 7.48, 7.49, 7.5, 7.51, 7.52, 7.53, 7.54, 7.55, 7.56, 7.57, 7.58, 7.59, 7.6, 7.61, 7.62, 7.63, 7.64, 7.65, 7.66, 7.67, 7.68, 7.69, 7.7, 7.71, 7.72, 7.73, 7.74, 7.75, 7.76, 7.77, 7.78, 7.79, 7.8, 7.81, 7.82, 7.83, 7.84, 7.85, 7.86, 7.87, 7.88, 7.89, 7.9, 7.91, 7.92, 7.93, 7.94, 7.95, 7.96, 7.97, 7.98, 7.99, 8.0, 8.01, 8.02, 8.03, 8.04, 8.05, 8.06, 8.07, 8.08, 8.09, 8.1, 8.11, 8.12, 8.13, 8.14, 8.15, 8.16, 8.17, 8.18, 8.19, 8.2, 8.21, 8.22, 8.23]

CARBON_BURN_LOGRHO = [9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9729, 9.9728, 9.9728, 9.9728, 9.9727, 9.9726, 9.9726, 9.9724, 9.9723, 9.9721, 9.9718, 9.9715, 9.9712, 9.9704, 9.97, 9.9693, 9.9679, 9.9668, 9.9653, 9.9631, 9.9607, 9.957, 9.9548, 9.9497, 9.9437, 9.941, 9.9313, 9.924, 9.9137, 9.9062, 9.8922, 9.89, 9.8756, 9.8662, 9.8531, 9.8405, 9.8253, 9.8121, 9.8008, 9.7825, 9.7769, 9.7564, 9.7352, 9.7287, 9.7068, 9.6852, 9.675, 9.6518, 9.6235, 9.6212, 9.591, 9.5675, 9.5424, 9.5239, 9.4984, 9.4764, 9.4497, 9.42, 9.4062, 9.3671, 9.3362, 9.3113, 9.2741, 9.245, 9.207, 9.1673, 9.1236, 9.0838, 9.0421, 8.9917, 8.9423, 8.8915, 8.8312, 8.7768, 8.7064, 8.6355, 8.5621, 8.4738, 8.385, 8.2927, 8.1814, 8.0796, 7.973, 7.8481, 7.72, 7.56, 7.4018, 7.2187, 7.0412, 6.8309, 6.6113, 6.3612, 6.1075, 5.8398, 5.5178, 5.206, 4.9097, 4.6049, 4.2874, 3.9601, 3.6373, 3.3021, 2.9746, 2.6338, 2.3063, 1.9714, 1.6483, 1.3221]
CARBON_BURN_LOGT = [7.11, 7.12, 7.13, 7.14, 7.15, 7.16, 7.17, 7.18, 7.19, 7.2, 7.21, 7.22, 7.23, 7.24, 7.25, 7.26, 7.27, 7.28, 7.29, 7.3, 7.31, 7.32, 7.33, 7.34, 7.35, 7.36, 7.37, 7.38, 7.39, 7.4, 7.41, 7.42, 7.43, 7.44, 7.45, 7.46, 7.47, 7.48, 7.49, 7.5, 7.51, 7.52, 7.53, 7.54, 7.55, 7.56, 7.57, 7.58, 7.59, 7.6, 7.61, 7.62, 7.63, 7.64, 7.65, 7.66, 7.67, 7.68, 7.69, 7.7, 7.71, 7.72, 7.73, 7.74, 7.75, 7.76, 7.77, 7.78, 7.79, 7.8, 7.81, 7.82, 7.83, 7.84, 7.85, 7.86, 7.87, 7.88, 7.89, 7.9, 7.91, 7.92, 7.93, 7.94, 7.95, 7.96, 7.97, 7.98, 7.99, 8.0, 8.01, 8.02, 8.03, 8.04, 8.05, 8.06, 8.07, 8.08, 8.09, 8.1, 8.11, 8.12, 8.13, 8.14, 8.15, 8.16, 8.17, 8.18, 8.19, 8.2, 8.21, 8.22, 8.23, 8.24, 8.25, 8.26, 8.27, 8.28, 8.29, 8.3, 8.31, 8.32, 8.33, 8.34, 8.35, 8.36, 8.37, 8.38, 8.39, 8.4, 8.41, 8.42, 8.43, 8.44, 8.45, 8.46, 8.47, 8.48, 8.49, 8.5, 8.51, 8.52, 8.53, 8.54, 8.55, 8.56, 8.57, 8.58, 8.59, 8.6, 8.61, 8.62, 8.63, 8.64, 8.65, 8.66, 8.67, 8.68, 8.69, 8.7, 8.71, 8.72, 8.73, 8.74, 8.75, 8.76]

OXYGEN_BURN_LOGRHO = [10.5839, 10.5839, 10.57, 10.56, 10.5546, 10.542, 10.5409, 10.53, 10.5198, 10.51, 10.5127, 10.5, 10.493, 10.4777, 10.469, 10.46, 10.44, 10.4409, 10.4409, 10.4361, 10.4358, 10.44, 10.418, 10.4004, 10.391, 10.3923, 10.38, 10.3707, 10.3508, 10.337, 10.321, 10.32, 10.3006, 10.264, 10.2649, 10.2726, 10.2748, 10.26, 10.2634, 10.264, 10.2, 10.1519, 10.1417, 10.1413, 10.1417, 10.1519, 10.08, 10.0448, 10.02, 10.0029, 9.9, 9.78, 9.72, 9.66, 9.6, 9.54, 9.48, 9.42, 9.36, 9.3, 9.24, 9.18, 9.12, 9.0, 8.76, 8.52, 8.34, 8.16, 7.98, 7.8, 7.62, 7.5, 6.9, 6.6, 2.5]
OXYGEN_BURN_LOGT = [7.5, 8.21, 8.235, 8.2539, 8.2645, 8.288, 8.2899, 8.309, 8.326, 8.341, 8.3373, 8.3569, 8.3667, 8.3878, 8.3991, 8.41, 8.4331, 8.4321, 8.432, 8.4376, 8.438, 8.4331, 8.456, 8.472, 8.4811, 8.48, 8.4898, 8.498, 8.5129, 8.5237, 8.534, 8.5348, 8.5479, 8.57, 8.5695, 8.5653, 8.564, 8.5721, 8.5703, 8.57, 8.6034, 8.6252, 8.6298, 8.63, 8.6298, 8.6252, 8.6539, 8.666, 8.6746, 8.6797, 8.7098, 8.7392, 8.7518, 8.7637, 8.7748, 8.785, 8.7945, 8.8036, 8.8117, 8.8199, 8.8272, 8.834, 8.841, 8.8531, 8.8744, 8.8924, 8.9042, 8.9146, 8.9242, 8.9329, 8.9408, 8.946, 8.9681, 8.9776, 9.1]

import numpy as np


# Order matters: checked from hottest (most-evolved) to coolest, so the
# first curve a star's current (T, Rho) has climbed past (scanning from
# the top) is reported - a star that has already reached, say, O burning
# has necessarily also climbed past the H/He/C curves earlier in its life,
# so we want the MOST ADVANCED stage it currently satisfies, not the first.
_BURN_STAGES = [
    (4, "O", OXYGEN_BURN_LOGRHO, OXYGEN_BURN_LOGT),
    (3, "C", CARBON_BURN_LOGRHO, CARBON_BURN_LOGT),
    (2, "He", HELIUM_BURN_LOGRHO, HELIUM_BURN_LOGT),
    (1, "H", HYDROGEN_BURN_LOGRHO, HYDROGEN_BURN_LOGT),
]

BURN_STAGE_LABELS = {0: "None", 1: "H", 2: "He", 3: "C", 4: "O"}

# Pre-sorted (by logRho, ascending) once at import time for np.interp, which
# requires a non-decreasing x array - each curve as given by MESA is
# parameterized by a monotonic logT sweep, not logRho, so a couple of tiny
# near-vertical/kinked segments (confirmed directly: the helium curve holds
# logRho ~8.4835 across logT 6.25-6.5; the oxygen curve backtracks by
# ~0.001 dex near logT~8.43) become slightly ambiguous once resorted this
# way - immaterial for classifying real stellar cores (those corners sit
# far outside any temperature/density a real star's core actually passes
# through on its way up each curve), so not worth a more elaborate
# nearest-segment/polyline approach for a 4-way discrete classification.
def _prep_curve(logRho, logT):
    logRho = np.asarray(logRho, dtype=float)
    logT = np.asarray(logT, dtype=float)
    order = np.argsort(logRho)
    return logRho[order], logT[order]


_PREPPED = {label: _prep_curve(rho, T) for _, label, rho, T in _BURN_STAGES}


def classify_core_burning(log_center_T, log_center_Rho):
    '''Classifies a star's current core burning stage from its own
    (log_center_T, log_center_Rho) - Ben: "MESA tracks the core
    temperature and density and it plots this showing which type of
    burning is happening (H, He, C and O). Could you also extract this
    core burning information from MESA models?"

    Compares each frame's (log_center_T, log_center_Rho) against the same
    four fixed reference ignition curves MESA's own pgstar TRho_Profile
    plot draws (see this module's own docstring) - a real star's core
    trajectory generally moves up-and-right through this diagram over
    its life (contracting increases both T and Rho), crossing each curve
    once, in order. At a given log_center_Rho, "ignited" for a species
    means log_center_T has reached or exceeded that curve's own
    threshold logT there (interpolated linearly between the curve's
    tabulated points, clamped to the curve's own edge value beyond its
    tabulated logRho range - these curves already span an extremely wide
    density range, so real stellar cores essentially never fall outside
    it).

    Parameters
    ----------
    log_center_T : array-like
        log10(core temperature [K]) at each frame.
    log_center_Rho : array-like
        log10(core density [g/cm^3]) at each frame.

    Returns
    -------
    numpy.ndarray of FLOAT (integer-valued: 0.0-4.0), same length as the
    inputs - float, not int, to match every other data_t array's own dtype
    (int64 specifically breaks tulips3d_create_star_profile_pie's own
    generic "copy every top-level dict entry onto the object as a custom
    property" loop, confirmed directly: `TypeError: object of type
    'numpy.int64' has no len()`). 0 = below the H-burn threshold (e.g.
    still pre-main-sequence), 1 = H burning, 2 = He burning, 3 = C
    burning, 4 = O burning - the MOST ADVANCED stage each frame currently
    satisfies (a star doing O burning has, by definition, also climbed
    past the H/He/C thresholds earlier in its life).
    '''
    log_center_T = np.asarray(log_center_T, dtype=float)
    log_center_Rho = np.asarray(log_center_Rho, dtype=float)
    stage = np.zeros(log_center_T.shape, dtype=float)
    remaining = np.ones(log_center_T.shape, dtype=bool)
    for code, label, _, _ in _BURN_STAGES:
        curve_logRho, curve_logT = _PREPPED[label]
        threshold_logT = np.interp(log_center_Rho, curve_logRho, curve_logT)
        ignited = remaining & (log_center_T >= threshold_logT)
        stage[ignited] = code
        remaining &= ~ignited
    return stage
