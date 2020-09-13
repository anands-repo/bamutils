# bamutils
Some tools for manipulating BAM files that I didn't find in standard tools

split_bams.py: Allows you to split a BAM file into disjoint parts

split_fastq.py: Allows you to split a FASTQ file into disjoint parts. Expects the fastq reads to have one-to-one correspondence in the case of a paired-end fastq file (e.g., if samtools with multiple threads are used to convert bam to fastq, such ordering is not necessarily met)
