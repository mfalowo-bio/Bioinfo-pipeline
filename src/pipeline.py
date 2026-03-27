"""A modular, reproducible bioinformatics pipeline for processing and analyzing sequencing data."""


# simple quality control to count reads, average reads and gc content

# code
from Bio import SeqIO

fastq_file = "data/sample.fastq"

total_reads = 0
total_bases = 0
gc_count = 0

for record in SeqIO.parse(fastq_file, "fastq"):
    total_reads += 1
    seq = record.seq
    total_bases += len(seq)
    gc_count += seq.count("G") + seq.count("C")

avg_read_length = total_bases / total_reads
gc_content = gc_count / total_bases * 100

print(f"Total reads: {total_reads}")
print(f"Average read length: {avg_read_length:.2f}")
print(f"GC content: {gc_content:.2f}%")



def filter_reads_by_length(reads, min_length=50):
    """
    Remove reads shorter than a minimum length threshold.
    """
    filtered_reads = []

    for read in reads:
        if len(read) >= min_length:
            filtered_reads.append(read)

    return filtered_reads

def run_length_filtering(reads, length_threshold):
    print(f"Filtering reads shorter than {length_threshold} bases")
    return filter_reads_by_length(reads, min_length=length_threshold)