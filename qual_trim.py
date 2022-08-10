#CALL LIBRARIES
import os as os
from tqdm import tqdm as tq

#DEFINE FUNCTION TO QUALITY TRIM FASTQ FILES
def seqtrim(file_path, adaptor, min_length, min_qual):
    #SET WORKING DIRECTORY AND LIST FILES IN DIRECTORY
    os.chdir(file_path)
    files = os.listdir(file_path)
    for file in files:
        #CREATE EMPTY LIST TO STORE LINES FROM IN FILE (FASTQ) AND CREATE OUT FILE (FASTA) TO APPEND LINES TO
        record_list = []
        trimmed = open('path to outfile folder' + file[0:-5] + 'fasta', 'a', newline = '')
        file = open(file, 'r')
        for line in tq(file):
            line = line.strip()
            #ONCE A FULL RECORD HAS BEEN APPENDED TO (4 LINES) DUPLICATE LIST AND START SEQUENCE TRIMMING
            if len(record_list) == 4:
                record = record_list
                #TRIM BASES WITH A LOWER QUALITY THAT THE SET MINIMUM QUALITY SCORE
                index = [record[3].index(i) for i in list(record[3]) if (ord(i) - 33) < min_qual]
                list(set(index))
                count = 0
                for i in index:
                    i = i - count
                    record[1] = record[1][:i] + record[1][i+1:]
                    count = count + 1
                #REMOVE N BASES AND ADAPTOR SEQUENCES
                record[1] = record[1].replace('N', '')
                record[1] = record[1].replace(adaptor, '')
                #IF THE LENGTH OF THE SEQUENCE IS GREATER THAN THE MINIMUM SEQUENCE LENGTH POST TRIMMING, WRITE TRIMMED RECORD TO OUT FILE
                if len(record[1]) >= min_length:
                    trimmed.write('>' + record[0][1:-1] + '\n' + record[1] + '\n')
                    record_list = []
                    record = []
                #IF THE LENGTH OF THE SEQUENCE IS LESS THAN THE MINIMUM SEQUENCE LENGTH POST TRIMMING, WRITE ORIGINAL RECORD TO OUT FILE
                else:
                    trimmed.write('>' + record_list[0][1:-1] + '\n' + record_list[1] + '\n')
                    record_list = []
                    record = []
            record_list.append(line)

seqtrim('path to fastq file folder', 'adaptor sequence', minimum read length(integer), minimum quality score(integer))
