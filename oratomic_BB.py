import argparse
import time

import stim
from mip import OptimizationStatus

from ilp_circuit_distance import mip_circuit_distance

TIME_LIMIT = 10000


def make_code() -> tuple[list[stim.PauliString], list[stim.PauliString], list[stim.PauliString]]:
    l = 31
    m = 4
    terms_a = [(0, 0), (6, 1), (27, 0)]
    terms_b = [(0, 2), (15, 3), (24, 0)]
    n = 2 * l * m

    def left(i: int, j: int) -> int:
        return (i % l) + l * (j % m)

    def right(i: int, j: int) -> int:
        return l * m + (i % l) + l * (j % m)

    stabilizers: list[stim.PauliString] = []
    for i in range(l):
        for j in range(m):
            x_stabilizer = stim.PauliString(n)
            for dx, dy in terms_a:
                x_stabilizer[left(i + dx, j + dy)] = "X"
            for dx, dy in terms_b:
                x_stabilizer[right(i + dx, j + dy)] = "X"
            stabilizers.append(x_stabilizer)

            z_stabilizer = stim.PauliString(n)
            for dx, dy in terms_b:
                z_stabilizer[left(i - dx, j - dy)] = "Z"
            for dx, dy in terms_a:
                z_stabilizer[right(i - dx, j - dy)] = "Z"
            stabilizers.append(z_stabilizer)

    completed_tableau = stim.Tableau.from_stabilizers(
        stabilizers,
        allow_redundant=True,
        allow_underconstrained=True,
    )
    obs_indices = [
        k
        for k in range(len(completed_tableau))
        if completed_tableau.z_output(k) not in stabilizers
    ]
    observable_xs = [completed_tableau.x_output(k) for k in obs_indices]
    observable_zs = [completed_tableau.z_output(k) for k in obs_indices]
    return stabilizers, observable_xs, observable_zs


def make_memory_circuit(
    stabilizers: list[stim.PauliString],
    observables: list[stim.PauliString],
    basis: str,
) -> stim.Circuit:
    if basis not in {"x", "z"}:
        raise ValueError(f"Unsupported basis: {basis!r}")

    num_qubits = len(stabilizers[0])
    circuit = stim.Circuit()
    error_gate = "Z_ERROR" if basis == "x" else "X_ERROR"
    stabilizer_basis = "X" if basis == "x" else "Z"

    if basis == "x":
        circuit.append("H", range(num_qubits))
    circuit.append(error_gate, range(num_qubits), 1e-3)
    for k, observable in enumerate(observables):
        circuit.append("MPP", stim.target_combined_paulis(observable))
        circuit.append("OBSERVABLE_INCLUDE", stim.target_rec(-1), k)
    for stabilizer in stabilizers:
        if stabilizer.pauli_indices(stabilizer_basis):
            circuit.append("MPP", stim.target_combined_paulis(stabilizer))
            circuit.append("DETECTOR", stim.target_rec(-1))
    return circuit


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--basis",
        choices=["x", "z"],
        default="z",
        help="Compute the X- or Z-type distance by choosing the memory basis.",
    )
    parser.add_argument(
        "--time-limit",
        type=int,
        default=TIME_LIMIT,
        help="Time limit in seconds for the ILP solve.",
    )
    args = parser.parse_args()

    stabilizers, observable_xs, observable_zs = make_code()
    print(f"n = {len(stabilizers[0])}")
    print(f"num stabilizers = {len(stabilizers)}")
    print(f"k = {len(observable_zs)}")
    print(f"basis = {args.basis}")

    observables = observable_xs if args.basis == "x" else observable_zs
    circuit = make_memory_circuit(stabilizers, observables, args.basis)

    for solver in ["GRB"]:
        print(f"\nSolver: {solver}")
        t0 = time.monotonic()
        try:
            result = mip_circuit_distance(
                circuit.detector_error_model(), 
                solver_name=solver, 
                time_limit=args.time_limit, 
                # verbose=True
            )
        except Exception as exc:
            print(f" | failed: {exc}")
            continue
        print(" | distance: ", result["distance"])
        print(" | optimal: ", result["status"] == OptimizationStatus.OPTIMAL)
        if result["status"] != OptimizationStatus.OPTIMAL:
            print(" | status: ", result["status"])
            print(" | lower bound: ", result["lower_bound"])
        print(" | time: ", time.monotonic() - t0)
