import numpy as np
import pandas as pd 
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
import pysam
import argparse

def merge_and_trim_sequences(seq1, seq2, umi1, umi2):
    # Align the two sequences using local alignment
    alignments = pairwise2.align.localms(seq1, seq2, 2, -1, -3, -0.5, one_alignment_only=True)
    aligned_seq1, aligned_seq2, score, start, end = alignments[0]

    # Initialize the consensus sequence with the first aligned sequence
    consensus = list(aligned_seq1)

    # Merge the aligned sequences into the consensus sequence
    for i, (base1, base2) in enumerate(zip(aligned_seq1, aligned_seq2)):
        if base1 == '-':
            consensus[i] = base2
        elif base2 == '-':
            continue
        elif base1 == 'N' and base2 != 'N':
            consensus[i] = base2
        elif base2 == 'N' and base1 != 'N':
            consensus[i] = base1
        elif base1 != base2:
            # Choose the base from seq2 if a mismatch occurs, unless seq2 contains 'N'
            consensus[i] = base2 if base2 != 'N' else base1

    # Join consensus list into string and remove gaps
    consensus = ''.join(consensus).replace('-', '')

    # Find UMI1 and remove everything to the left
    umi1_index = consensus.find(umi1)
    #if umi1_index != -1:
    #    consensus = consensus[umi1_index + len(umi1):]

    # Find UMI2 and remove everything to the right (reverse complement)


    umi2_rc = reverse_complement(umi2)
    umi2_index = consensus.rfind(umi2_rc)
    if (umi2_index != -1) and (umi1_index != -1):
        consensus = consensus[umi1_index + len(umi1):]
        consensus = consensus[:umi2_index]


    return consensus

def reverse_complement(seq):
    complement = str.maketrans("ATGC", "TACG")
    return seq.translate(complement)[::-1]

def extract_barcode(consensus, vector, barcode_length=12):
    # Align the consensus sequence to the known vector sequence
    alignments = pairwise2.align.localms(consensus, vector, 2, -1, -3, -0.5, one_alignment_only=True)
    if alignments:
        aligned_consensus_seq, aligned_vector_seq, score, start, end = alignments[0]
        barcode_start = end
        barcode_end = barcode_start + barcode_length
        barcode = aligned_consensus_seq[barcode_start:barcode_end]
        return barcode, score
    return "", 0

pre_barcode_vector1 = 'GTACCAAGGTCGGGCAGGAA' #9nt after
pre_barcode_vector2 = 'GTTTAAGAGCTAAGCTGGAA' #12nt after


Reference = 'NNNNNNNTGACTAACCAATTCAGTCGACTGGATCCGGTACCAAGGTCGGGCAGGAANNNNNNNNNGAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGGCCCAGACTGAGCACGTGAGTTTAAGAGCTAAGCTGGAANNNNNNNNNNNNAAAAAGTGGCACCGAGTCGGTGCGCATTACAGTACGAAGATGCATAGCTTTTTTTAAGCTTGGGCCGCTCGAGGTACCTCTCTATGACTNNNNNNN'

def extract_Promoter(consensus, barcode1, guide='GGCCCAGACTGAGCACGTGA'):

    start_index = consensus.find(barcode1)+len(barcode1) #Guide Index

    #Now get alignment of the guide
    alignments = pairwise2.align.localms(consensus, guide, 2, -1, -3, -0.5, one_alignment_only=True)
    aligned_consensus_seq, aligned_vector_seq, score, start_guide, end_guide = alignments[0]
    if score < 32: #Indicates vector is not mapping anywhere as well as it should
        return 'NG'
    elif start_index >= start_guide:
        return 'NB'
    else:
    #Return the sequence between the Barcode and the guide
        return consensus[start_index:start_guide]

def process_bam_file(input_bam):
    samfile = pysam.AlignmentFile(input_bam, "rb")
    Barcodes_Promoter = []
    
    # Track the last processed read name to detect and skip additional occurrences
    last_processed_read_name = None
    pending_read = None

    # Read through the BAM file and process reads in pairs
    try:
        for read in samfile:
            if read.query_name == last_processed_read_name:
                # Skip this read if we have already processed a pair with this read name
                continue

            if pending_read is None:
                # If there's no pending read, store the current read
                pending_read = read
            else:
                # We have a pair to process
                read1, read2 = pending_read, read

                # Reset the pending read
                pending_read = None

                # Process the pair
                UMI = read1.query_name.split('_')[-1]
                umi1, umi2 = UMI[:7], UMI[7:14]
                consensus_seq = merge_and_trim_sequences(read1.seq, read2.seq, umi1, umi2)
                barcode1, score1 = extract_barcode(consensus_seq, pre_barcode_vector1, barcode_length=9)
                barcode2, score2 = extract_barcode(consensus_seq, pre_barcode_vector2, barcode_length=12)
                if (score1 or score2) < 32: #Indicates vector is not mapping anywhere as well as it should
                    #Excluded from the analysis
                    #Either no vector or no barcode. 
                    Barcodes_Promoter.append((UMI,len(consensus_seq), '', '', 'NV'))
                    continue
                elif (Reference.find(barcode1) != -1) or (Reference.find(barcode2) != -1):
                    Barcodes_Promoter.append((UMI,len(consensus_seq), '', '', 'NB'))
                    continue
                    
                Promoter = extract_Promoter(consensus_seq, barcode1, guide='GGCCCAGACTGAGCACGTGA')
                Barcodes_Promoter.append((UMI,len(consensus_seq), barcode1, barcode2, Promoter))

                # Update the last processed read name
                last_processed_read_name = read1.query_name

    except StopIteration:
        # End of file reached
        pass
    
    samfile.close()
    return Barcodes_Promoter

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BAM file and export to CSV.")
    parser.add_argument("input_bam", help="Input BAM file path")
    parser.add_argument("output_file", help="Output CSV file path")
    args = parser.parse_args()
    results = process_bam_file(args.input_bam)
    df = pd.DataFrame(results, columns = ['UMI','Length', 'Barcode1', 'Barcode2', 'U6'])
    df.to_csv(args.output_file,  index = False, header = False)
