# Overview
Seqtrim is quality trimming and adaptor removal tool built in Python. Both the quality trimming and adaptor removal functions of this program use a sliding window algorithm. For quality trimming, if the mean base quality inside the window falls beneath the user set threshold, the read is trimmed from the window start to the end of the read. For adaptor removal, if the sequence in the window matches the adaptor sequence (offset by "-m max_mismatches") the portion of sequence in the window is trimmed.

# Data Format
SeqTrim takes uncompressed FASTQ. It's output is also in an uncompressed format.

# Usage
```
python {installation_path}/main.py -r1 read1.fastq -r2 read2.fastq -a adaptor_sequence -m max_mismatches -w window_length -q phred-quality_threshold -o1 read1_output-path -o2 read2_output-path
```
