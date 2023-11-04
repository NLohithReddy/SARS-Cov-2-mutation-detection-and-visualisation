# import all of the required modules
from dna import dna
import numpy as np
import matplotlib.pyplot as plt
from scov import numpy_image_dict
from helper import *

# Read the dna sequence file-1 previously downloaded from NCBI.
dict_seq_1 = read_dna_seq('USA.txt')
# Modify the sequence with dummy 'N' nucleotide.
dict_seq_1 = gene_mod(dict_seq_1)

# Read the dna sequence file-2 previously downloaded from NCBI.
dict_seq_2 = read_dna_seq('China.txt')
# Modify the sequence with dummy 'N' nucleotide.
dict_seq_2 = gene_mod(dict_seq_2)

# Create matplotlib subplots for each gene.
f, ax = plt.subplots(nrows=11, ncols=3, figsize=(25, 30))
gene_name = list(numpy_image_dict.keys())
row = 0
col = 0
mut_dict = {}
new_sequence = ""  # Initialize a new sequence

for i in gene_name:
    G = i[5:]
    # Loop through each gene in the Coronavirus nucleotide sequence.
    gene_us = dna(dict_seq_1['gene='+G][1])
    # Invoke the transcription method of the class dna
    gene_us.transcription()
    # Invoke the method that converts the gene sequence into a numpy array.
    numpfy_usa = gene_us.numpfy()
    # Reshape the numpy array with a predefined shape from the numpy_image_dict dictionary.
    numpfy_usa = numpfy_usa.reshape(numpy_image_dict['gene='+G][0])
    # sub-plot the numpy array with matplotlib pcolor method.
    ax[row][col].pcolor(numpfy_usa)
    ax[row][col].set_title(G+' Gene - USA')
    col += 1
    gene_china = dna(dict_seq_2['gene='+G][1])
    # Invoke the transcription method of the class dna
    gene_china.transcription()
    # Invoke the method that converts the gene sequence into a numpy array.
    numpfy_china = gene_china.numpfy()
    # Reshape the numpy array with a predefined shape from the numpy_image_dict dictionary.
    numpfy_china = numpfy_china.reshape(numpy_image_dict['gene='+G][0])
    # sub-plot the numpy array with matplotlib pcolor method.
    ax[row][col].pcolor(numpfy_china)
    ax[row][col].set_title(G+' Gene - CHINA')
    col += 1

    # To find the gene mutation subtract the numpy array from the base sequence with the newer sequence.
    # Here the Chinese sequence is the base sequence and the USA sequence is a newer sequence.
    mut = numpfy_china - numpfy_usa
    if mut.any():
        # Here we are looking for a non-zero value in the mutated numpy array
        # (result of subtracting the 2 numpy arrays).
        # Presence of a non-zero value means that there is a difference between the 2 numpy arrays and the gene has
        # mutations.
        mut_nec = np.nonzero(mut)
        x = mut_nec[0]
        y = mut_nec[1]
        l = 0
        mut_dict[G] = []
        for i in x:
            us_base = numpfy_usa[i][y[l]]
            ch_base = numpfy_china[i][y[l]]
            mut_base = mut[i][y[l]]
            info_list = [ch_base, us_base, mut_base, (i, y[l])]
            mut_dict[G].append(info_list)
            print("Mutated DNA Base {} in China and Base {} in USA at position {} For the Gene {}".format(
                ch_base, us_base, (i, y[l]), G))
            l += 1

    # Combine the sequences and generate a new sequence with mutations
    combined_sequence = dict_seq_1['gene='+G][1] + dict_seq_2['gene='+G][1]
    mutated_sequence = ""
    mutation_rate = 0.01  # Define your mutation rate here (adjust as needed)

    for base in combined_sequence:
        if np.random.rand() < mutation_rate:
            # Mutation occurred, randomly choose a new base
            new_base = np.random.choice(['A', 'C', 'G', 'T'])
            mutated_sequence += new_base
        else:
            # No mutation, use the original base
            mutated_sequence += base

    new_sequence += mutated_sequence  # Add the mutated gene to the new sequence

    # Giving a title to the matplotlib subplot
    ax[row][col].pcolor(mut)
    ax[row][col].set_title(G+' Gene - Mutation')
    row += 1
    col = 0

f.tight_layout()
print(new_sequence)

# Saving the matplotlib subplot as a jpg.
f.savefig('Sars_Cov-2_Gene_Mutation.jpg')


