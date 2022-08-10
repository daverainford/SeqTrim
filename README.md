# SeqTrim
SeqTrim is a fast and streamlined function to trim low-quality bases and remove adaptor sequences. A ~10GB fastq file can take just minutes to trim if ran through python, and if using the python-bridge module in Node.js that time is cut in half. I have included both the pure python and the javascript files in the repo.

It takes four inputs: the path to the folder containing your fastq files, your adaptor sequence, the minimum read length you deem acceptable, and the minimum quality score you deem acceptable.

This function will eventually be absorbed into the sequence analysis pipeline I am creating as well as be integrated into my website. I also am currently coding some tests to validate the function and will push those once complete.

Thank you for checking out my project and feel free to shoot me an email with any questions!
