import numpy as np
import json


def build_circuit(n, m, k):
    number_of_inputs = n*m + m*k + 1
    number_of_outputs = n*k
    number_of_mid_wires = n*k*m     # One for each initial mult-gate
    number_of_wires = number_of_inputs + number_of_mid_wires + number_of_outputs
    mult_gates = []

    for i in range(n):
        for j in range(k):
            for l in range(m):
                a_index = m * i + l + 1
                b_index = n*m + k * l + j + 1
                mult_gates.append(([a_index], [b_index]))

    for i in range(n):
        for j in range(k):
            start_index = number_of_inputs + i * m*k + j*m
            mult_wires = []
            for l in range(m):
                index = start_index + l
                mult_wires.append(index)
            mult_gates.append((mult_wires, [0]))  # Empty list means *1

    circuit = {
        "numberOfWires": number_of_wires,
        "numberOfInputWires": number_of_inputs,
        "numberOfOutputWires": number_of_outputs,
        "gates": mult_gates
    }
    
    return circuit

# Small test
def eval(circuit, input):
    wire_values = []
    for i in range(circuit["numberOfInputWires"]):
        wire_values.append(input[i])
    for i in range(circuit["numberOfWires"] - circuit["numberOfInputWires"]):
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
    input = [1] + A_flat + B_flat
    print(input)
    return eval(circuit, input)

def output_json(circuit):
    file = open("../out/matrix.txt", "w")
    file.write(json.dumps(circuit))
    file.close()

def output_text(circuit, name):
    file = open("../out/%s.txt" % name, "w")
    file.write(str(circuit["numberOfWires"]))
    file.write(" ")
    file.write(str(circuit["numberOfInputWires"]))
    file.write(" ")
    file.write(str(circuit["numberOfOutputWires"]))
    file.write("\n")
    for gate in circuit["gates"]:
        file.write(str(len(gate[0])) + " ")
        file.write(str(len(gate[1])) + " ")
        file.write(" ".join(map(str, gate[0])))
        file.write(" ")
        file.write(" ".join(map(str, gate[1])))
        file.write("\n")
    file.close()


if __name__ == "__main__":
    A = np.array([[1,2], [3,4]])
    B = np.array([[2,0], [0,2]])
    
    n = 10
    m = 10
    k = 10
    circuit = build_circuit(n, m, k)
    
    output_text(circuit, "n=%d_m=%d_k=%d" %(n, m, k))

    # print(matrix_eval(circuit, A, B))
