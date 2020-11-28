import numpy as np


def build_circuit(n, m, k):
    circuit = []
    for i in range(n*m + m*k):
        circuit.append(["INPUT"])

    for i in range(n):
        for j in range(k):
            for l in range(m):
                a_index = m * i + l
                b_index = n*m + k * l + j
                circuit.append(["MULT", [a_index], [b_index]])

    mult_start_index = n*m + m*k
    for i in range(n):
        for j in range(k):
            start_index = mult_start_index + i * m*k + j*m
            mult_wires = []
            for l in range(m):
                index = start_index + l
                mult_wires.append(index)
            circuit.append(["OUTPUT", mult_wires])
    
    return circuit

def eval(circuit, input):
    wire_values = []
    for i in range(len(circuit)):
        if circuit[i][0] == "INPUT":
            wire_values.append(input[i])
        if circuit[i][0] == "OUTPUT":
            sum_wires = circuit[i][1]
            sum = 0
            for index in sum_wires:
                sum += wire_values[index]
            wire_values.append(sum)
        if circuit[i][0] == "MULT":
            left_wires = circuit[i][1]
            right_wires = circuit[i][2]
            left = 0
            right = 0
            for index in left_wires:
                left += wire_values[index]
            for index in right_wires:
                right += wire_values[index]
            wire_values.append(left * right)

    return wire_values

def matrix_eval(circuit, A, B):
    A_flat = A.flatten().tolist()
    B_flat = B.flatten().tolist()
    input = A_flat + B_flat
    print(input)
    return eval(circuit, input)


if __name__ == "__main__":
    A = np.array([[1,2], [3,4]])
    B = np.array([[2,0], [0,2]])
    
    circuit = build_circuit(2,2,2)
    for i in range(len(circuit)):
        wire = circuit[i]
        print(i, ":", wire)
    
    evaluated = matrix_eval(circuit, A, B)
    for i in range(len(evaluated)):
        value = evaluated[i]
        print(i, ":", value)

    print(np.dot(A,B))