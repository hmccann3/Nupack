from nupack import *
import pandas as pd
import matplotlib.pyplot as plt

"""
Runs nupack on all windows of length window_size in a fasta and generates plots
Records nupack structures and scores for each window in outfile
"""
### Define physical model for Nupack
my_model = Model(material='rna', celsius=28)
fig, ax = plt.subplots()
### Fasta containing all robin segments from study
genome_fasta_1 = open('/oak/stanford/groups/horence/Holly/robin.fasta', 'r')
seq_array = []
window_size = 75
k = -1
for line in genome_fasta_1:
    strip_line = line.strip()
    if '>' in strip_line:
        seq_array.append('')
        k += 1
    elif len(strip_line) > 0:
        seq_array[k] += strip_line
genome_fasta_1.close()
low = 0
low_index = 0
low_struct = ''
outfile = open('nupack_robin_' + str(window_size) + '_windows', 'a')
for seq in seq_array:
    ### Convert DNA to mRNA
    rna_full = seq.replace('T', 'U')
    j = 0
    y = []
    x = list(range(1, len(rna_full) + 1 - window_size))
    while (j < len(rna_full) - window_size):
        rna_sequence = rna_full[j:j+window_size]
        j += 1
        a = Strand(rna_sequence, name='a')
        t1 = Tube(strands={a:5e-6}, complexes=SetSpec(max_size=1), name='Tube t1')
        tube_result = tube_analysis(tubes=[t1], compute=['mfe'], model=my_model)
        walker_result = tube_result['(a)']
        for i, s in enumerate(walker_result.mfe):
            if i == 0:
                ### Record Nupack results for each sequence in outfile
                outfile.write(rna_sequence + '\n')
                outfile.write('%s'%(s.structure) + '\n')
                outfile.write(str(s.energy) + '\n')
                y.append(s.energy)
    ax.plot(x, y)
outfile.close()
plt.savefig("Nupack_full_robin_" + str(window_size), bbox_inches = "tight")
