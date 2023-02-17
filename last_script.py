ratio_ok = 0.6

sample2name = {}
samplenamesfile = open('samplenames.txt')
for line in samplenamesfile:
    l = line.strip().split('\t')
    sample2name[l[0] + '_' + l[1]] = l[2]
samplenamesfile.close()

sample2seqns = {}
samples = []
tablefile = open('table.csv')
for line in tablefile:
    if 'Assembly' in line:
        l = line.strip().split(',')
        sample = l[1].split('.')[1]
        if not sample in samples:
            samples.append(sample)
            sample2seqns[sample] = [int(l[0])]
        else:
            sample2seqns[sample].append(int(l[0]))
tablefile.close()

sample2confidence = {}
for sample in samples:
    confidence = max(sample2seqns[sample]) / sum(sample2seqns[sample])
    sample2confidence[sample] = confidence

outfile = open('TOTAL.fasta', 'w')
outfile.write('')
outfile.close()
badfile = open('TO_CHECK.txt', 'w')
badfile.write('')
badfile.close()

seqfile = open('consensus.fasta')
outfile = open('TOTAL.fasta', 'a')
badfile = open('TO_CHECK.txt', 'a')
lines = []
for line in seqfile:
    lines.append(line)
    if len(lines) == 2:
        if not 'Assembly  ' in lines[0]:
            sample = lines[0].split('.')[1]
            if sample2confidence[sample] > ratio_ok:
                #good
                outfile.write('>' + sample2name[sample] + '\n')
                outfile.write(lines[1])
            else:
                #bad
                badfile.write('Check sample ' + sample + ' (')
                badfile.write(sample2name[sample] + '). Ratio: ')
                badfile.write(str(sample2confidence[sample]) + '\n')
        lines = []
seqfile.close()
outfile.close()
badfile.close()
