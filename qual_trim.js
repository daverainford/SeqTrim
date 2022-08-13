// Call modules
const pythonBridge = require('python-bridge');
const lineReader = require('line-reader');
const fs = require('fs')

// Function to trim low quality bases and adaptor sequences
function SeqTrim(fastq_path, temp_path, fasta_path, adaptor, min_qual, min_length) {
	// Define seq list and chunk of sequence records to be written to file
	const record = []
	const chunk = []
	// Iterate over lines in fastq file
	lineReader.eachLine(fastq_path,(line,last)=>{
		record.push(line)
		if (record.length == 4) {
			// Once seq list full, split list into one string to speed up python analysis then push to data chunk with linebreak
			let data = record.join('split')
			chunk.push(data + '\n')
			record.length = 0
		}
		// If data chunk contains 1000000 seq lists, join the list and write to temp file
		if (chunk.length == 1000000) {
			fs.appendFile(temp_path, chunk.join(''), (err) => {
				if (err) {
					console.log(err);
				} 
			})
			chunk.length = 0
		}
	})

	// Define Python intepreter and python code
	let python = pythonBridge();
	python.ex
	`
	#Import library
	import time

	#Start performance timer
	startTime = time.time()
	time.sleep(12)
	
	#Define file paths and start count for indexing
	temp = open(${temp_path}, 'r')
	trimmed = open(${fasta_path}, 'a')
	count = 0
	
	#Iterate over lines in temp file
	for line in temp:
		#Split record back into seq list
		record = line.split('split')
		name = record[0]
		seq = record[1]
		qual = record[3]
		#Trim adaptor if found in sequence
		if ${adaptor} in seq:
			adaptor_index.append(seq.find(${adaptor}))
			adaptor_len = len(${adaptor})
			seq = seq[:adaptor_index] + seq[adaptor_index + adaptor_len:]
			qual = qual[:adaptor_index] + qual[adaptor_index + adaptor_len:]
		#Get indexes of low quality and N bases then trim them
		base_index = [seq.index(i) for i,j in zip(seq, qual) if i == 'N' or ord(j) - 33 < ${min_qual}]
		for index in base_index:
			index = index - count
			seq = seq[:index] + seq[index + 1:]
			qual = qual[:index] + qual[index + 1:]
			count+=1
		#If sequence length is greater than the minimum length, write trimmed sequence to fasta
		if len(seq) >= ${min_length}:
			trimmed.write(name.replace('@', '>') + '''\n''' + seq + '''\n''')
		#If sequence length is less than the minimum length, write original sequence to fasta
		else:
			trimmed.write(name.replace('@', '>') + '''\n''' + record[0] + '''\n''')
	
	#Stop performance timer and print time elapsed since start
	executionTime = (time.time() - startTime)
	print('Execution time in seconds: ' + str(executionTime))
	`
	python.end()
}

// Call function
SeqTrim('./fastq.fastq', './temp.fastq', 'trimmed.fasta', 'CTGTCTCTTATACACATCT', 20, 50);
