from mip import Model, BINARY, INTEGER, OptimizationStatus, xsum
import stim

def mip_circuit_distance(dem: stim.DetectorErrorModel, time_limit: int = 300, *, verbose: bool = False, solver_name: str = None):
    """
        Solve for the minimum number of DEM error mechanisms whose combined effect:
        - flips NO detectors (all detectors even parity)
        - flips AT LEAST ONE logical observable (some logical has odd parity)

    Formulation:
        x_j ∈ {0,1}  for each error mechanism j
        y_i ∈ ℤ≥0    detector parity slack vars
        s_r ∈ ℤ≥0    logical parity slack vars
        w_r ∈ {0,1}  logical parity bits (mod 2)

        For each detector i:
            sum_j D[i,j] * x_j = 2 * y_i

        For each logical r:
            sum_j L[r,j] * x_j = 2 * s_r + w_r

        At least one logical flipped:
            sum_r w_r ≥ 1

        Objective:
            minimize sum_j x_j

    Returns:
        {
            "status": mip.OptimizationStatus,
            "distance": int | None,         # minimal number of DEM errors, or None if infeasible
            "error_indices": list[int] | None,  # indices j of chosen DEM errors
        }
    """
    err_detectors, err_observables, num_detectors, num_observables = parse_dem_errors(dem)

    num_errors = len(err_detectors)
    if num_errors == 0:
        return {"status": None, "distance": None, "error_indices": None}

    # Build inverted maps: for each detector/logical, which errors touch it
    det_to_errors = [[] for _ in range(num_detectors)]
    for j, dets in enumerate(err_detectors):
        for d in dets:
            det_to_errors[d].append(j)

    obs_to_errors = [[] for _ in range(num_observables)]
    for j, obs in enumerate(err_observables):
        for r in obs:
            obs_to_errors[r].append(j)

    # MIP model
    m = Model("shortest_undetected_logical_error", solver_name=solver_name)

    if not verbose:
        m.verbose = 0

    # Decision variables: x_j ∈ {0,1}, select error mechanisms
    x = [m.add_var(var_type=BINARY, name=f"x_{j}") for j in range(num_errors)]

    # Upper bounds for slack vars (parity counts)
    max_pairs = num_errors // 2

    # Detector slack vars: y_i ∈ ℤ≥0, enforce even parity (no detector flips)
    y = [
        m.add_var(
            var_type=INTEGER,
            lb=0,
            ub=max_pairs,
            name=f"y_det_{i}",
        )
        for i in range(num_detectors)
    ]

    # Logical slack vars and parity bits
    s = [
        m.add_var(
            var_type=INTEGER,
            lb=0,
            ub=max_pairs,
            name=f"s_log_{r}",
        )
        for r in range(num_observables)
    ]
    w = [
        m.add_var(var_type=BINARY, name=f"w_log_{r}")
        for r in range(num_observables)
    ]

    # Detector parity constraints: sum_j D[i,j] x_j = 2 y_i
    for i in range(num_detectors):
        if det_to_errors[i]:
            m += (
                xsum(x[j] for j in det_to_errors[i]) - 2 * y[i] == 0,
                f"detector_{i}_parity",
            )
        else:
            # No errors touch this detector; it is always even, no constraint needed
            pass

    # Logical parity constraints: sum_j L[r,j] x_j = 2 s_r + w_r
    for r in range(num_observables):
        if obs_to_errors[r]:
            m += (
                xsum(x[j] for j in obs_to_errors[r]) - 2 * s[r] - w[r] == 0,
                f"logical_{r}_parity",
            )
        else:
            # This logical is never affected; enforce w_r = 0
            m += (w[r] == 0, f"logical_{r}_always_even")

    # At least one logical has odd parity
    if num_observables > 0:
        m += (xsum(w) >= 1, "some_logical_flipped")
    else:
        # No logicals defined => no undetected logical error makes sense
        return {"status": None, "distance": None, "error_indices": None}

    # Objective: minimize number of errors
    m.objective = xsum(x)

    if time_limit is not None:
        m.max_mip_gap = 0  # exact optimality
        m.max_seconds = time_limit

    m.optimize()

    status = m.status

    ok_statuses = {OptimizationStatus.OPTIMAL, OptimizationStatus.FEASIBLE}
    if status not in ok_statuses:
        # Infeasible or no solution found
        return {"status": status, "distance": None, "error_indices": None}

    # Extract solution
    chosen = [j for j in range(num_errors) if x[j].x >= 0.5]
    distance = len(chosen)

    return {
        "status": status,
        "distance": distance,
        "error_indices": chosen,
        "lower_bound": m.objective_bound,
    }


def parse_dem_errors(dem: stim.DetectorErrorModel):
    """
    Extract, for each error mechanism in the DEM:
        - which detectors it flips
        - which logical observables it flips

    Returns:
        err_detectors: list[list[int]]
        err_observables: list[list[int]]
        num_detectors: int
        num_observables: int
    """
    err_detectors = []
    err_observables = []

    # Running offset for DETECTOR ids (Stim uses RELATIVE_DET ids with shifts)
    running_detector_offset = 0

    for inst in dem:
        if inst.type == "shift_detectors":
            # This increases the base index for subsequent relative detector ids
            # Targets: [relative_shift, time_shift], we only care about the detector shift
            # det_shift = inst.targets[0].val  # adjust to .value if needed
            # running_detector_offset += det_shift
            raise NotImplementedError("DEM parsing for shift_detectors not implemented yet.")

        elif inst.type == "error":
            dets = []
            obs = []
            for t in inst.targets_copy():
                if t.is_relative_detector_id():
                    # absolute detector index = running offset + relative id
                    dets.append(running_detector_offset + t.val)
                elif t.is_logical_observable_id():
                    # observable (logical) index
                    obs.append(t.val)

            # Deduplicate and sort for consistency
            dets = sorted(set(dets))
            obs = sorted(set(obs))

            err_detectors.append(dets)
            err_observables.append(obs)

        else:
            # other DEM instructions (logical_include, detector, etc.) are irrelevant here
            continue

    num_detectors = dem.num_detectors
    num_observables = dem.num_observables

    return err_detectors, err_observables, num_detectors, num_observables
