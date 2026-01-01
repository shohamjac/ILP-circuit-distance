import time
import stim
from ilp_circuit_distance import mip_circuit_distance
from mip import OptimizationStatus

# This is a verbatim implementation of the code posted by Craig Gideny here: 
# https://quantumcomputing.stackexchange.com/questions/37289/compute-the-exact-minimum-distance-of-a-qecc-with-integer-linear-programming-met


def make_code() -> tuple[list[stim.PauliString], list[stim.PauliString], list[stim.PauliString]]:
    offsets = {1, -1, 1j, -1j, 3 + 6j, -6 + 3j}
    w = 24
    h = 12

    def wrap(c: complex) -> complex:
        return c.real % w + (c.imag % h)*1j

    def index_of(c: complex) -> int:
        c = wrap(c)
        return int(c.real + c.imag * w) // 2

    stabilizers: list[stim.PauliString] = []
    for x in range(w):
        for y in range(h):
            if x % 2 != y % 2:
                continue  # This is a data qubit.
            m = x + 1j*y
            basis = 'XZ'[x % 2]
            sign = -1 if basis == 'Z' else +1
            stabilizer = stim.PauliString(w * h // 2)
            for offset in offsets:
                stabilizer[index_of(m + offset * sign)] = basis
            stabilizers.append(stabilizer)

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
    observable_xs: list[stim.PauliString] = [
        completed_tableau.x_output(k)
        for k in obs_indices
    ]
    observable_zs: list[stim.PauliString] = [
        completed_tableau.z_output(k)
        for k in obs_indices
    ]

    return stabilizers, observable_xs, observable_zs


def make_z_basis_circuit(stabilizers: list[stim.PauliString], obs_zs: list[stim.PauliString]):
    num_qubits = len(stabilizers[0])
    circuit = stim.Circuit()
    circuit.append("X_ERROR", range(num_qubits), 1e-3)
    for k, observable in enumerate(obs_zs):
        circuit.append("MPP", stim.target_combined_paulis(observable))
        circuit.append("OBSERVABLE_INCLUDE", stim.target_rec(-1), k)
    for stabilizer in stabilizers:
        if stabilizer.pauli_indices('Z'):
            circuit.append("MPP", stim.target_combined_paulis(stabilizer))
            circuit.append("DETECTOR", stim.target_rec(-1))
    return circuit


if __name__ == "__main__":
    stabilizers, observable_xs, observable_zs = make_code()
    circuit = make_z_basis_circuit(stabilizers, observable_zs)

    for solver in ["CBC", "GRB"]:
        print(f"\nSolver: {solver}")
        t0 = time.monotonic()
        result = mip_circuit_distance(circuit.detector_error_model(), solver_name=solver)
        print(" | distance: ", result["distance"])
        print(" | optimal: ", result["status"] == OptimizationStatus.OPTIMAL)
        if not result["status"] == OptimizationStatus.OPTIMAL:
            print(" | status: ", result["status"])
            print(" | lower bound: ", result["lower_bound"])
        print(" | time: ", time.monotonic() - t0)


