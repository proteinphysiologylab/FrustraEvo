#!/usr/bin/python
# Import packages
import os
import sys

file = sys.argv[1]
output_path = sys.argv[2]


f1 = open(output_path+ '/' +
          os.path.basename(file) + '.protcoverage70.filtered', 'a')
f1.write('qseqid\tsseqid\tstaxid\tevalue\tpident\tbitscore\tqcovs\tsallseqid\tqlen\tslen\tqstart\tqend\tsstart\tsend\tlength\tnident\tqseq\tsseq\n')

f1b = open(output_path+ '/' +
          os.path.basename(file) + '.filtered', 'a')
f1b.write('qseqid\tsseqid\tstaxid\tevalue\tpident\tbitscore\tqcovs\tsallseqid\tqlen\tslen\tqstart\tqend\tsstart\tsend\tlength\tnident\tqseq\tsseq\n')


f2 = open(output_path+ '/' +
          os.path.basename(file) + '.protcoverage70.filtered.seqnames', 'a')

f2b = open(output_path+ '/' +
          os.path.basename(file) + '.filtered.seqnames', 'a')
    # Append 'hello' at the end of file

for line in open(file, 'r'): 
    #print(line)
    ls = line.split('\t')

    if ls[2] not in ['32630']:
        if ls[16].count('-') < 10:
            if float(ls[3]) < 0.05:
                if 'X' not in ls[17]:
                    if (float(ls[14]) + float(ls[14]) * 0.3) >= float(ls[8]) and  (float(ls[14]) - float(ls[14]) * 0.3) <= float(ls[8]): 
                        sacc = ls[1]
                        f1.write(line)
                        f2.write(sacc + '\n')
                    
                    sacc = ls[1]
                    f1b.write(line)
                    f2b.write(sacc + '\n')

# Close the file
f1.close()
f1b.close()
f2.close()
f2b.close()