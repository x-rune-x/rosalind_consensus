import numpy as np


class FastaFile:
    def __init__(self, fasta_id, seq):
        self.fasta_id = fasta_id
        self.seq = seq

    def get_id(self):
        return self.fasta_id

    def get_seq(self):
        return self.seq


def create_matrix(file_name):
    file = open(file_name, "r")

    current_id = ""
    current_seq = ""
    string_matrix = []

    for line in file:
        if line[0] == ">" and current_id == "":
            current_id = line[1:].rstrip()
        elif line[0] == ">" and current_id != "":
            # print(f"id is {current_id}, seq is {current_seq}")
            seq_list = []
            for char in current_seq:
                seq_list.append(char)
            string_matrix.append(seq_list)
            current_id = line[1:].rstrip()
            current_seq = ""
        else:
            current_seq += line.rstrip()
    else:
        # print(f"id is {current_id}, seq is {current_seq}")
        seq_list = []
        for char in current_seq:
            seq_list.append(char)
        string_matrix.append(seq_list)

    file.close()
    return string_matrix


def create_profile_matrix(matrix):
    profile_matrix = [["A:", "C:", "G:", "T:"]]
    list_length = len(matrix)
    string_length = len(matrix[0])

    for column in range(string_length):
        a_count = 0
        c_count = 0
        g_count = 0
        t_count = 0
        for index in range(list_length):
            base = matrix[index][column]
            if base == "A":
                a_count += 1
            elif base == "C":
                c_count += 1
            elif base == "G":
                g_count += 1
            elif base == "T":
                t_count += 1

        column_list = [a_count, c_count, g_count, t_count]
        profile_matrix.append(column_list)
        final_matrix = np.array(profile_matrix)
        final_matrix = final_matrix.transpose()

    return final_matrix


def get_consensus(profile_matrix):
    profile_matrix = profile_matrix.transpose()
    string_length = len(profile_matrix)
    consensus_seq = ""
    for column in range(1, string_length):
        mc_base = max(profile_matrix[column])

        if mc_base == profile_matrix[column][0]:
            consensus_seq += "A"
        elif mc_base == profile_matrix[column][1]:
            consensus_seq += "C"
        elif mc_base == profile_matrix[column][2]:
            consensus_seq += "G"
        elif mc_base == profile_matrix[column][3]:
            consensus_seq += "T"

    return consensus_seq


first_matrix = create_matrix("rosalind_cons.txt")
second_matrix = create_profile_matrix(first_matrix)
consensus_seq = get_consensus(second_matrix)

output_file = open("consensus_seq.txt", "w")
output_file.write(consensus_seq + "\n")
np.savetxt(output_file, second_matrix, fmt='%s')

output_file.close()
