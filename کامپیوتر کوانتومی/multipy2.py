from qiskit import QuantumCircuit, transpile, Aer, execute
from qiskit.circuit.library import QFT

# Function to create a quantum circuit for multiplying two numbers
def multiply_quantum(a, b):
    # Number of qubits needed
    n = max(a.bit_length(), b.bit_length())
    
    # Create a quantum circuit with 2n qubits for the inputs and n qubits for the result
    qc = QuantumCircuit(2*n + n, n)
    
    # Encode the numbers a and b into the quantum circuit
    for i in range(n):
        if (a >> i) & 1:
            qc.x(i)
        if (b >> i) & 1:
            qc.x(n + i)
    
    # Apply QFT to the result qubits
    qc.append(QFT(n).to_gate(), range(2*n, 2*n + n))
    
    # Controlled addition
    for i in range(n):
        for j in range(n):
            qc.cp((2 * 3.14159) / (2**(i + j + 1)), n + j, 2*n + i)
    
    # Apply inverse QFT to the result qubits
    qc.append(QFT(n, inverse=True).to_gate(), range(2*n, 2*n + n))
    
    # Measure the result
    qc.measure(range(2*n, 2*n + n), range(n))
    
    return qc

# Example usage
a = 3  # First number
b = 2  # Second number
qc = multiply_quantum(a, b)

# Simulate the quantum circuit
simulator = Aer.get_backend('qasm_simulator')
compiled_circuit = transpile(qc, simulator)
result = execute(compiled_circuit, simulator).result()

# Get the result
counts = result.get_counts(qc)
# Convert the result to a single number
result_number = int(list(counts.keys())[0], 2)
print("Result:", result_number)
