import functions as func
import os
import argparse

# initialize flag arguments
parser = argparse.ArgumentParser(description="Seqtrim is a quality trimming and adaptor removal tool compatible with genomic sequences stored in fastq files")

parser.add_argument("--read1", "-r1", help="Fastq file containing reads from pair 1")
parser.add_argument("--read2", "-r2", help="Fastq file containing reads from pair 2")
parser.add_argument("--adaptor", "-a", help="Full sequence of the adaptors used for library prep")
parser.add_argument("--max", "-m", help="Max mismatches for adaptor matching (0 = perfect match)")
parser.add_argument("--window", "-w", help="Length of sliding window for quality trimming")
parser.add_argument("--quality", "-q", help="Threshold to trim for mean phred quality in window ")
parser.add_argument("--out1", "-o1", help="Output path and file name for first paired reads")
parser.add_argument("--out2", "-o2", help="Output path and file name for second paired reads")

# store flag arguments as variables
args = parser.parse_args()

file1 = args.read1
file2 = args.read2
adaptor = args.adaptor
max_mismatch = int(args.max)
window_length = int(args.window)
threshold = int(args.quality)
out1 = args.out1
out2 = args.out2

# open both untrimmed and output files
read1 = open(file1, "r")
read2 = open(file2, "r")

out1 = open(out1, "a")
out2 = open(out2, "a")

entry = []
# read in untrimmed entries then write trimmed entries to output
for line1, line2 in zip(read1, read2):
    if len(entry) == 8:
        adaptor_trimmed_entry = func.adaptorTrim(entry, adaptor, max_mismatch)
        entry[2] = adaptor_trimmed_entry[0]
        entry[3] = adaptor_trimmed_entry[1]
        entry[6] = adaptor_trimmed_entry[2]
        entry[7] = adaptor_trimmed_entry[3]
        qual_trimmed_entry = func.qualTrim(entry, window_length, threshold)
        id1 = entry[0].strip()
        id2 = entry[1].strip()
        seq1 = qual_trimmed_entry[0].strip()
        seq2 = qual_trimmed_entry[1].strip()
        link1 = entry[4].strip()
        link2 = entry[5].strip()
        qual1 = qual_trimmed_entry[2].strip()
        qual2 = qual_trimmed_entry[3].strip()
        out1.write(id1 + "\n")
        out1.write(seq1 + "\n")
        out1.write(link1 + "\n")
        out1.write(qual1 + "\n")
        out2.write(id2 + "\n")
        out2.write(seq2 + "\n")
        out2.write(link2 + "\n")
        out2.write(qual2 + "\n")
    else:
        entry.append(line1)
        entry.append(line2)

# close files
read1.close()
read2.close()
out1.close()
out2.close()