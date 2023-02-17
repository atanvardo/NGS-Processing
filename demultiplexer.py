# Script by Alejandro Lopez-Lopez
# 
# Aims of the script: We received sequences from an Illumina run
# in which the tags (and in some cases the primers) are incomplete.
# This prevents us to use other scripts in published pipelines, so
# we had to create this script. This script will read the sequences,
# select those that are identifiable by the tag, and them demultiplex
# them, saving them in different files corresponding to the different
# samples. We will also remove the tags and leave them ready to be
# processed.

primers = ['CCHGAYATRGCHTTYCCHCG', 'TCDGGRTGNCCRAARAAYCA']

# Definition of functions that will be used in the script:


def identifier(sequencelines, primerlist, tolerance):
    # This function checks in a sequence has any of the primers of the provided list in the possible positions (0 to 7).
    # If found, it will return a list composed by the tag and the lines forming the new trimmed sequence.
    # [ tag , [ seqname, seqtrimmed, +, qctrimmed ] ]
    # If the primer is not found, it returns False

    sequence = sequencelines[1].strip()
    sequence_qc = sequencelines[3].strip()

    char2nuc = {'A': ['A'],
                'T': ['T'],
                'C': ['C'],
                'G': ['G'],
                'R': ['A', 'G'],
                'Y': ['C', 'T'],
                'N': ['A', 'C', 'G', 'T'],
                'W': ['A', 'T'],
                'S': ['G', 'C'],
                'M': ['A', 'C'],
                'K': ['G', 'T'],
                'B': ['C', 'G', 'T'],
                'H': ['A', 'C', 'T'],
                'D': ['A', 'G', 'T'],
                'V': ['A', 'C', 'G']}

    for position in range(0, 8):
        for primer in primerlist:
            target = sequence[position:(position + len(primer))]
            difs = 0
            for nuc in range(len(target)):
                if target[nuc] not in char2nuc[primer[nuc]]:
                    difs += 1
            if difs <= tolerance:
                # The primer is at that position!
                tag = sequence[0: position]
                if len(tag) == 0:
                    # We create a safeguard to avoid problems later. A missing tag
                    # is trated as no primer found. Else, it will cause false
                    # identifications when we try to determine which sample it is.
                    return False
                if len(tag) > 6:
                    tag = tag[-6:]
                trimmed_sequence = sequence[position:]
                trimmed_qc = sequence_qc[position:]
                return [tag, [sequencelines[0],
                              trimmed_sequence + '\n',
                              '+\n',
                              trimmed_qc + '\n'],
                        primer]
    return False


def determine_sample(sample_tag_f, sample_tag_r):
    # This function takes the two tags of a sample, and returns a string composed by the numbers of the
    # F and R primers corresponding to that sample, separated by an underscore. If one of the tags is not
    # found in the dictionary, it returns the string "unknown".
    tag2name = {'AACCGA': '1', 'CCGGAA': '2', 'AGTGTT': '3', 'CCGCTG': '4', 'AACGCG': '5', 'GGCTAC': '6',
                'TTCTCG': '7', 'TCACTC': '8', 'GAACTA': '9', 'CACAGT': '10', 'CAATCG': '11', 'CCGTCC': '12',
                'GGGACA': '13', 'AGCTCA': '14', 'ACTGGG': '15', 'GATCGG': '16', 'CTAGGC': '17', 'TGAGGT': '18',
                'TCAACT': '19', 'TACACA': '20', 'GATGAC': '21', 'AGTAGA': '22', 'TCCTTT': '23', 'ATGAGG': '24'}

    tags = list(tag2name.keys())

    # If the tags are complete, they will have a length of 6 and the comparison is easy:
    if len(sample_tag_f) == 6 and len(sample_tag_r) == 6:
        if sample_tag_f in tags:
            f_label = tag2name[sample_tag_f]
        else:
            return 'unknown'
        if sample_tag_r in tags:
            r_label = tag2name[sample_tag_r]
        else:
            return 'unknown'
        return f_label + '_' + r_label
    # If the length of the tags is smaller (incomplete but theoretically identifiable) we will have to create
    # a new dictionary in which the tags are truncated to the length of our tag:
    elif len(sample_tag_f) < 4 or len(sample_tag_r) < 4:
        return 'unknown'
    else:
        new_dictionary = {}
        for tag in tag2name:
            new_dictionary[tag[len(sample_tag_f):]] = tag2name[tag]
        tags = list(new_dictionary.keys())
        if sample_tag_f in tags:
            f_label = tag2name[sample_tag_f]
        else:
            return 'unknown'
        new_dictionary = {}
        for tag in tag2name:
            new_dictionary[tag[len(sample_tag_r):]] = tag2name[tag]
        tags = list(new_dictionary.keys())
        if sample_tag_r in tags:
            r_label = tag2name[sample_tag_r]
        else:
            return 'unknown'
        return f_label + '_' + r_label


def classify(r1filename, r2filename, prefix):
    # This is the main function of the script. It will read the two fasta files, taking 4 lines (1 fastq sequence)
    # each time. Every group of four lines includes:
    #   Line 1: name of the sequence, preceded by "@"
    #   Line 2: the sequence
    #   Line 3: the symbol "+"
    #   Line 4: the quality score of the sequence
    # If both reads (R1 and R2) have identifiable tags, we will store them in the r1/2outfile files.
    # If we find a tag that is not in our list, we will store them in the r1/2 unknown files.
    # If we can't find even the primer, we will store them in the r1/2discarded files.
    # We also include "good" files to store all the good sequences to check that we get all.

    # We start initializing the files:
    r1file = open(r1filename, 'r')
    r2file = open(r2filename, 'r')
    r1discardedfile = open(prefix + '.discarded.R1.fastq', 'a')
    r1unknownfile = open(prefix + '.unknown.R1.fastq', 'a')
    r2discardedfile = open(prefix + '.discarded.R2.fastq', 'a')
    r2unknownfile = open(prefix + '.unknown.R2.fastq', 'a')
    r1goodfile = open(prefix + '.good.R1.fastq', 'a')
    r2goodfile = open(prefix + '.good.R2.fastq', 'a')

    # And now we start reading the lines of the r1file and r2file by groups of four.
    # We will store the groups of four lines in these lists:
    logfile = open('log.txt', 'a')
    sequence_lines_r1 = []
    sequence_lines_r2 = []
    for line_r1, line_r2 in zip(r1file, r2file):
        # We read one line of each file and add it to the list
        sequence_lines_r1.append(line_r1)
        sequence_lines_r2.append(line_r2)
        # If after this we have four items in each list, that means that we have the four lines
        # that define each sequence, and we can process them.
        if len(sequence_lines_r1) == 4:
            # First, we are going to determine if there is a primer in the R1 sequence.
            identity_r1 = identifier(sequence_lines_r1, primers, 3)
            # And in the R2 sequence
            identity_r2 = identifier(sequence_lines_r2, primers, 3)
            # If at least one of them is equal to False, we have a missing primer. Then, the sequences will
            # be written to the discarded files:
            if not identity_r1 or not identity_r2:
                r1discardedfile.write(''.join(sequence_lines_r1))
                r2discardedfile.write(''.join(sequence_lines_r2))
            else:
                sequence_r1 = ''.join(identity_r1[1])
                sequence_r2 = ''.join(identity_r2[1])
                if identity_r1[2] == 'CCHGAYATRGCHTTYCCHCG':
                    tag_f = identity_r1[0]
                    tag_r = identity_r2[0]
                elif identity_r1[2] == 'TCDGGRTGNCCRAARAAYCA':
                    tag_f = identity_r2[0]
                    tag_r = identity_r1[0]
                else:
                    tag_f = 'unk'
                    tag_r = 'unk'
                samplename = determine_sample(tag_f, tag_r)
                if samplename == 'unknown':
                    r1unknownfile.write(sequence_r1)
                    r2unknownfile.write(sequence_r2)
                else:
                    if int(samplename.split('_')[0]) > 8:
                        logfile.write('STRANGEFOUND\n')
                        logfile.write(identity_r1[0] + '\n')
                        logfile.write(''.join(sequence_lines_r1))
                        logfile.write(identity_r1[2] + '\n')
                        logfile.write(identity_r2[0] + '\n')
                        logfile.write(''.join(sequence_lines_r2))
                        logfile.write(identity_r2[0] + '\n\n\n')
                    r1outfilename = prefix + '.' + samplename + '.R1.demux.fastq'
                    r2outfilename = prefix + '.' + samplename + '.R2.demux.fastq'
                    r1outfile = open(r1outfilename, 'a')
                    r2outfile = open(r2outfilename, 'a')
                    r1outfile.write(sequence_r1)
                    r2outfile.write(sequence_r2)
                    r1outfile.close()
                    r2outfile.close()
                    r1goodfile.write(sequence_r1)
                    r2goodfile.write(sequence_r2)

            # And, after processing these lines, we reset the list to leave them clean, waiting to read
            # the next sequence:
            sequence_lines_r1 = []
            sequence_lines_r2 = []
    logfile.close()

    # Close all the files:
    r1file.close()
    r2file.close()
    r1discardedfile.close()
    r2discardedfile.close()
    r1unknownfile.close()
    r2unknownfile.close()


classify('tardi-pilot_S1_L001_R1_001.fastq', 'tardi-pilot_S1_L001_R2_001.fastq', 'tardi-pilot')
