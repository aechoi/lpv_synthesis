"""Main loop for control synthesis of an LPV system with IQC defined uncertainty.

Details taken from the paper "Robust Synthesis for Linear Parameter Varying Systems
Using Integral Quadratic Constraints" by S. Wang, H. Pfifer, and P. Seiler.
"""

import control as ct
from lpv_synthesis import iqc


def robust_synthesis(
    plant_dynamics: ct.StateSpace,
    convergence_tolerance: float = 1e-3,
    max_iterations=50,
) -> None:
    """Main loop for control synthesis

    Args:
        plant_dynamics: The state-space representation of the plant dynamics. Should be
            in some form that allows for extraction of i/o dimensions of interconnect
            components.
        convergence_tolerance: The tolerance for convergence of the RP level.
        max_iterations: The maximum number of iterations to perform before giving up.
    """

    last_rp_level = float("inf")
    iqc_filters = []
    rp_levels = []

    iqc_filter = iqc.get_initializing_filter(plant_dynamics)

    for _ in range(max_iterations):
        controller = synthesize_controller(iqc_filter)
        iqc_filter, rp_level = analyze(controller)
        iqc_filters.append(iqc_filter)
        rp_levels.append(rp_level)

        if abs(rp_level - last_rp_level) < convergence_tolerance:
            break
        last_rp_level = rp_level
    else:
        print(
            f"Warning: Synthesis did not converge after {max_iterations} iterations. Final RP level: {rp_level}. Final RP difference: {abs(rp_level - last_rp_level)}"
        )
    return controller, (iqc_filters, rp_levels)


if __name__ == "__main__":
    plant = ct.ss()
    controller, _ = robust_synthesis(plant)
