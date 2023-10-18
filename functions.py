import statistics as stat
##########################################################################################################################################
def adaptorTrim(entry,adaptor, max_mismatch):
    # initialize sliding window and functionv ariable
    window_start = 0
    window_end = len(adaptor)
    seq1 = entry[2].strip()
    seq2 = entry[3].strip()
    qual1 = entry[6].strip()
    qual2 = entry[7].strip()
    cut1_start = 0
    cut1_end = 0
    cut2_start = 0
    cut2_end = 0
    # slide window across seq1 until there is an adaptor match
    while window_end <= len(seq1):
        scores = [1 if adaptor[index] == seq1[window_start:window_end][index] else 0 for index in range(len(adaptor))]
        if sum(scores) == (len(adaptor) - max_mismatch):
            cut1_start = cut1_start + window_start
            cut1_end = cut1_end + window_end
            break
        else:
            window_start += 1
            window_end += 1
    # slide window across seq2 until there is an adaptor match
    while window_end <= len(seq2):
        scores = [1 if adaptor[index] == seq2[window_start:window_end][index] else 0 for index in range(len(adaptor))]
        if sum(scores) == (len(adaptor) - max_mismatch):
            cut2_start = cut2_start + window_start
            cut2_end = cut2_end + window_end
            break
        else:
            window_start += 1
            window_end += 1
    # trim both the sequences and quality scores where the adaptor lies
    seq1_cut = seq1[: cut1_start] + seq1[cut1_end :]
    seq2_cut = seq2[: cut2_start] + seq2[cut2_end :]
    qual1_cut = qual1[: cut1_start] + qual1[cut1_end :]
    qual2_cut = qual2[: cut2_start] + qual2[cut2_end :]
    # ensure both read pairs are of the same length, if different return original sequences
    if len(seq1_cut) == len(seq2_cut):
        return [seq1_cut, seq2_cut, qual1_cut, qual2_cut]
    else:
        return [seq1, seq2, qual1, qual2]
##########################################################################################################################################
def qualTrim(entry, window_length, threshold):
    # initialize sliding window and function variables
    window_start = 0
    window_end = (window_length)
    seq1 = entry[2].strip()
    seq2 = entry[3].strip()
    qual1 = entry[6].strip()
    qual2 = entry[7].strip()
    cut1 = 0
    cut2 = 0
    final_cut = 0
    # convert ascii quality values into phred33
    phred1 = [(ord(value) - 33) for value in qual1]
    phred2 = [(ord(value) - 33) for value in qual2]
    # slide window along phred1 values until mean of window falls below threshold
    while window_end <= len(phred1):
        mean = stat.mean(phred1[window_start:window_end])
        if mean < threshold:
            cut1 = cut1 + (window_start - 1)
            window_start = 0
            window_end = window_length
            break
        else:
            window_start += 1
            window_end += 1
    # slide window along phred2 values until mean of window falls below threshold
    while window_end <= len(phred2):
        mean = stat.mean(phred2[window_start:window_end])
        if mean < threshold:
            cut2 = cut2 + (window_start - 1)
            break
        else:
            window_start += 1
            window_end += 1
    # calculate range between two cut points then cut both reads at the median of the range
    if cut1 > cut2:
        cut_range = [num for num in range(cut2,cut1 + 1)]
        final_cut = (final_cut + round(stat.median(cut_range),0)) - 1 
    else:
        cut_range = [num for num in range(cut1,cut2 + 1)]
        final_cut = (final_cut + round(stat.median(cut_range),0)) - 1
    # trim sequence and quality string
    seq1 = seq1[0:final_cut]
    seq2 = seq2[0:final_cut]
    qual1 = qual1[0:final_cut]
    qual2 = qual2[0:final_cut]
    return [seq1, seq2, qual1, qual2]
##########################################################################################################################################



    

