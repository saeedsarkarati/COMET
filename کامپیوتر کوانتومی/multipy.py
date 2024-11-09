'''
# ~ from qiskit import QuantumCircuit, transpile, assemble, execute
# ~ from qiskit.visualization import plot_histogram
# ~ from qiskit_aer import AerSimulator
from qiskit import QuantumCircuit, transpile, assemble
# Import section
from qiskit_aer import AerSimulator
# Function to create a quantum circuit for multiplying two numbers
def multiply_quantum(a, b):
    # Create a quantum circuit with 4 qubits and 4 classical bits
    qc = QuantumCircuit(4, 4)
    
    # Encode the numbers a and b into the quantum circuit
    if a == 1:
        qc.x(0)
    if b == 1:
        qc.x(1)
    
    # Apply a CNOT gate to multiply the numbers
    qc.ccx(0, 1, 2)
    
    # Measure the result
    qc.measure([2], [2])
    
    return qc

# Example usage
a = 1  # First number
b = 1  # Second number
qc = multiply_quantum(a, b)

# Simulate the quantum circuit
simulator = AerSimulator()
compiled_circuit = transpile(qc, simulator)
qobj = assemble(compiled_circuit)
result = simulator.run(qc, backend=simulator).result()

# Get the result
counts = result.get_counts(qc)
print("Result:", counts)
plot_histogram(counts)
'''
from qiskit import QuantumCircuit, transpile
from qiskit.visualization import plot_histogram
from qiskit_aer import AerSimulator

# Function to create a quantum circuit for multiplying two numbers
def multiply_quantum(a, b):
    # Create a quantum circuit with 4 qubits and 4 classical bits
    qc = QuantumCircuit(4, 4)
    
    # Encode the numbers a and b into the quantum circuit
    if a == 1:
        qc.x(0)
    if b == 1:
        qc.x(1)
    
    # Apply a CNOT gate to multiply the numbers
    qc.ccx(0, 1, 2)
    print (qc)
    # Measure the result
    qc.measure([2], [2])
    
    return qc

# Example usage
a = 10  # First number
b = 40  # Second number
qc = multiply_quantum(a, b)
print (qc)
# Simulate the quantum circuit
simulator = AerSimulator()
compiled_circuit = transpile(qc, simulator)
result = simulator.run(compiled_circuit).result()

# Get the result
counts = result.get_counts(qc)
# Convert the result to a single number
result_number = counts.keys()
print("Result:", result_number)
plot_histogram(counts)




