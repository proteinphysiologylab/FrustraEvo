import numpy as np
from Bio import SeqIO
import sys

def shannon_entropy(sequence):
    unique_residues, counts = np.unique(list(sequence), return_counts=True)
    probabilities = counts / len(sequence)
    entropy = -np.sum(probabilities * np.log2(probabilities))
    return entropy

def calculate_entropy_for_msa(msa_file):
    # Leer el archivo de MSA usando BioPython
    records = list(SeqIO.parse(msa_file, "fasta"))
    num_sequences = len(records)
    alignment_length = len(records[0].seq)
    
    # Crear una matriz para almacenar las secuencias del MSA
    alignment_matrix = np.zeros((num_sequences, alignment_length), dtype='str')
    
    # Rellenar la matriz con las secuencias del MSA
    for i, record in enumerate(records):
        alignment_matrix[i] = list(record.seq)
    
    # Calcular la entropía de Shannon para cada posición en el MSA
    entropy_values = []
    for position in range(alignment_length):
        entropy_values.append(shannon_entropy(alignment_matrix[:, position]))
    
    return entropy_values

# Ruta del archivo de MSA en formato FASTA
msa_file_path ='./../FrustraEvo_'+sys.argv[1]+'/OutPutFiles/'+'MSA_'+sys.argv[1]+'.fasta'

# Calcular la entropía para el MSA
entropies = calculate_entropy_for_msa(msa_file_path)
out=open('./../FrustraEvo_'+sys.argv[1]+'/OutPutFiles/SeqIC_'+sys.argv[1]+'.tab','w')
out.write('Position\tEntropy\n')
# Imprimir la entropía de cada posición en el MSA
for position, entropy in enumerate(entropies):
    pos=position + 1
    out.write(str(pos)+'\t'+str(entropy)+'\n')
out.close()
