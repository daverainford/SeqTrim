def seqtrim(min_qual, min_length, adaptor, infile, outfile):
    #Import libraries
    from tqdm import tqdm as tq
    import matplotlib.pyplot as plt
    #Define files, lists, and counts
    infile = open(infile, 'r')
    outfile = open(outfile, 'a')
    record = []
    base_qual_bef = []
    qual_count1 = 0
    base_count = 0
    #Iterate over each line in fastq file
    print('Quality trimming in progress. This may take several minutes...')
    for line in tq(infile):
        line = line.strip()
        record.append(line)
        #Halt iteration when the record list contains a full read record
        if len(record) == 4:
            #Initialize count for bases trimmed and define sequence and quality list/string objects
            count = 0
            seq_string = str(record[1])
            seq_list = [base for base in record[1]]
            qual_list = [ord(score) - 33 for score in record[3]]
            #Push per base quality to base_qual list for downstream visualiation
            if len(base_qual_bef) > 0:
                qual_count1 = qual_count1 + 1
                for index in range(len(qual_list)):
                    base_qual_bef[index] = base_qual_bef[index] + qual_list[index]
            else:
                qual_count1 = qual_count1 + 1
                for score in qual_list:
                    base_qual_bef.append(score)
            #If sequence contains adaptor sequence, remove adaptor index range from sequence and quality lists
            if seq_string.find(adaptor) != -1:
                start = seq_string.find(adaptor)
                end = start + len(adaptor)
                del seq_list[start:end]
                del qual_list[start:end]
                base_count = base_count + len(adaptor)
            #Create list of low quality base indicies that are below the minimum quality
            remove_index = [qual_list.index(score) for score in qual_list if score < 20]
            remove_index = list(set(remove_index))
            #Remove all bases from remove_index index list
            for index in remove_index:
                index = index - count
                seq_list.pop(index)
                qual_list.pop(index)
                count = count + 1
                base_count = base_count + 1
            qual_list = [chr(score + 33) for score in qual_list]
            #If sequence length meets minimum read length criteria, write trimmed sequence to outfile
            if len(seq_list) >= min_length:
                identifier = record[0] + '\n'
                seq_final = ''.join(seq_list) + '\n'
                plus = record[2] + '\n'
                qual_final = ''.join(qual_list) + '\n'
                outfile.write(identifier + seq_final + plus + qual_final)
            #If seuqence length does not meet minimum read length criteria, write original sequence to outfile
            else:
                identifier = record[0] + '\n'
                seq_final = record[1] + '\n'
                plus = record[2] + '\n'
                qual_final = record[3] + '\n'
                outfile.write(identifier + seq_final + plus + qual_final)
            record = []
            if qual_count1 == 1000000:
                break
    #Create line graph of per base quality before quality trimming
    for index in range(len(base_qual_bef)):
        base_qual_bef[index] = base_qual_bef[index] / qual_count1
    print(f'Quality trimming complete. {base_count} bases quality trimmed.')
    base_bef = [index + 1 for index in range(len(base_qual_bef))]
    plt.plot(base_bef, base_qual_bef, color ='blue')
    plt.xlabel("Base Position in Read")
    plt.ylabel("Average Quality Score")
    plt.title("Per Base Quality Before Trimming")
    plt.show()

#Call seqtrim function
seqtrim(20, 50, 'gtac', 'infile.fastq', 'outfile.fastq')
            
            



 
            
    