import numpy as np
import json


def build_circuit(n, m, k):
    number_of_inputs = n*m + m*k
    number_of_outputs = n*k
    number_of_mid_wires = n*k*m     # One for each initial mult-gate
    number_of_mult_gates = number_of_mid_wires + number_of_outputs  # We get an extra *1 for each sum (one for each output)
    number_of_wires = number_of_inputs + number_of_mid_wires + number_of_outputs
    mult_gates = []

    for i in range(n):
        for j in range(k):
            for l in range(m):
                a_index = m * i + l
                b_index = n*m + k * l + j
                mult_gates.append(([a_index], [b_index]))

    for i in range(n):
        for j in range(k):
            start_index = number_of_inputs + i * m*k + j*m
            mult_wires = []
            for l in range(m):
                index = start_index + l
                mult_wires.append(index)
            mult_gates.append((mult_wires, []))  # Empty list means *1

    circuit = {
        "numberOfWires": number_of_wires,
        "numberOfInputWires": number_of_inputs,
        "numberOfMidWires": number_of_mid_wires,
        "numberOfOutputWires": number_of_outputs,
        "numberOfMultiplicationGates": number_of_mult_gates,
        "gates": mult_gates
    }
    
    return circuit

# Small test
def eval(circuit, input):
    wire_values = []
    for i in range(circuit["numberOfInputWires"]):
        wire_values.append(input[i])
    for i in range(circuit["numberOfMultiplicationGates"]):
        left_wires = circuit["gates"][i][0]
        right_wires = circuit["gates"][i][1]
        left = 0
        right = 0
        for index in left_wires:
            left += wire_values[index]
        for index in right_wires:
            right += wire_values[index]
        if len(right_wires) == 0:
            right = 1
        wire_values.append(left * right)

    return wire_values

# Small test
def matrix_eval(circuit, A, B):
    A_flat = A.flatten().tolist()
    B_flat = B.flatten().tolist()
    input = A_flat + B_flat
    print(input)
    return eval(circuit, input)

def output_json(circuit):
    file = open("../out/matrix.txt", "w")
    file.write(json.dumps(circuit))
    file.close()

def output_text(circuit):
    file = open("../out/matrix2.txt", "w")
    file.write(str(circuit["numberOfWires"]))
    file.write("\n")
    file.write(str(circuit["numberOfInputWires"]))
    file.write("\n")
    file.write(str(circuit["numberOfMidWires"]))
    file.write("\n")
    file.write(str(circuit["numberOfOutputWires"]))
    file.write("\n")
    file.write(str(circuit["numberOfMultiplicationGates"]))
    file.write("\n")
    for gate in circuit["gates"]:
        file.write(str(gate[0]))
        file.write(", ")
        file.write(str(gate[1]))
        file.write("\n")
    file.close()


if __name__ == "__main__":
    A = np.array([[1,2], [3,4]])
    B = np.array([[2,0], [0,2]])
    
    circuit = build_circuit(2,2,2)
    
    output_text(circuit)
