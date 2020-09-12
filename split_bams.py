import pysam
import argparse
import random
from collections import defaultdict
import logging
import numpy as np
import intervaltree


class Categorical:
	"""
	A categorical distribution for unequal proportion
	of outputs being piped
	"""
	def __init__(self, proportions):
		cumulative = np.cumsum([0] + proportions)
		assert(cumulative[-1]) == 1
		self.intervals = intervaltree.IntervalTree()
		for j, i in enumerate(range(len(cumulative) - 1)):
			self.intervals.addi(cumulative[i], cumulative[i + 1], j)

	def sample(self):
		sampling_point = random.uniform(0.0, 1.0)
		sampled_item = list(self.intervals[sampling_point]).pop()
		return sampled_item[2]


if __name__ == "__main__":
	options = argparse.ArgumentParser(description="Split a given BAM file into multiple BAM files")

	options.add_argument(
		"--bam",
		action="store",
		dest="bam",
		help="BAM file source",
		required=True
	)

	options.add_argument(
		"--prefix",
		action="store",
		dest="prefix",
		help="Prefix of the output file",
		required=True
	)

	options.add_argument(
		"--split_proportions",
		help="Comma-separated values of proportions of reads to send to each output file",
		required=True,
	)

	options.add_argument(
		"--seed",
		action="store",
		type=int,
		dest="seed",
		help="Seed for random sampling",
		default=13
	)

	opts = options.parse_args()
	random.seed(opts.seed)
	logging.basicConfig(level=logging.INFO, format='%(asctime)-15s %(message)s')

	proportions = list(map(float, opts.split_proportions.split(",")))
	logging.info("Will split bamfile in the following proportions: %s" % str(proportions))

	# Open the input bam file
	samfile = pysam.Samfile(opts.bam, 'r')

	# Open output bamfiles
	names = [opts.prefix + "_" + str(i) + ".bam" for i, _ in enumerate(proportions)]
	opsams = {name: pysam.AlignmentFile(name, 'wb', template=samfile) for name in names}

	# Keep a record of all outstanding templates. This is to track mate-pairs
	# so they go into the same file. Mate-pairs have the same template name
	outstanding_templates = dict()

	# Create a probability distribution with the available splits
	distribution = Categorical(proportions)

	# Traverse input bam and pipe to respective output files
	for i, read in enumerate(samfile.fetch()):
		query_name = read.query_name

		# Only primary alignments are considered
		if read.is_secondary:
			continue

		if (query_name is None) or (not read.is_paired):
			opname = names[distribution.sample()]
		elif query_name in outstanding_templates:
			opname = outstanding_templates[query_name]
			del outstanding_templates[query_name]
		else:
			opname = names[distribution.sample()]
			outstanding_templates[query_name] = opname

		target_bam = opsams[opname]

		if target_bam is not None:
			target_bam.write(read)

		if (i+1) % 100000 == 0:
			logging.info("Processed %d reads ..." %(i+1))
			logging.info("Number of outstanding templates %d" % len(outstanding_templates))

	# Close files
	map(lambda x: x.close(), opsams.values())
	samfile.close()
