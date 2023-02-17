text = ['#!/bin/bash']

text.append('dos2unix 01-demultiplexed/*.fastq')
text.append('mkdir 01-demultiplexed')
text.append('mkdir 02-quality_control')
text.append('mkdir 03-fused')
text.append('mkdir 04-fasta')
text.append('mkdir 05-blast')
text.append('mkdir 06-definitive')

for p in range(5,9):
    for i in range(1,25):
        for j in range(1,25):
            prefix = 'pool_' + str(p) + '.' + str(i) + '_' + str(j)
            # qualitycheck
            text.append('fastqc 01-demultiplexed/PREFIX.R1.demux.fastq'.replace('PREFIX',prefix))
            text.append('fastqc 01-demultiplexed/PREFIX.R2.demux.fastq'.replace('PREFIX',prefix))
            text.append('mv 01-demultiplexed/*.html 02-quality_control/')
            text.append('mv 01-demultiplexed/*.zip 02-quality_control/')
            # merge
            text.append('./tests/bbmap/fuse.sh -in1=01-demultiplexed/PREFIX.R1.demux.fastq -in2=01-demultiplexed/PREFIX.R2.demux.fastq out=03-fused/PREFIX.fused.fastq -pad=1 -fusepairs=t'.replace('PREFIX',prefix))
            # blast
            text.append('./tests/bbmap/reformat.sh -in=03-fused/PREFIX.fused.fastq -out=04-fasta/PREFIX.fused.fasta'.replace('PREFIX',prefix))
            blastline = ''
            blastline += 'blastn -db ../coidb/TardiCOIdb.fasta -query 04-fasta/PREFIX.fused.fasta '
            blastline += '-out 05-blast/PREFIX.fused.blastn -task blastn '
            blastline += '-evalue 1e-5 -max_target_seqs 1 -outfmt 6 -dust yes'
            text.append(blastline.replace('PREFIX',prefix))

text.append('python3 blast_selector_fused.py')

script = open('script_pools_5_8_fused.sh','w', newline='\n')
script.write('\n'.join(text))
script.close()
