import pysam
import argparse
import random

options = argparse.ArgumentParser();

options.add_argument("--bam", action="store", dest="bam", help="BAM file", required=True);
options.add_argument("--prefix", action="store", dest="prefix", help="Prefix of the output file", required=True);
options.add_argument("--num_splits", action="store", type=int, dest="num_splits", help="Number of splits", required=True);
options.add_argument("--conseve_mate_pairs", action="store_true", dest="conserve_mate_pairs", help="Keep mate pairs in the same output bam file", default=False);
options.add_argument("--seed", action="store", type=int, dest="seed", help="Seed for random sampling", default=13);

opts = options.parse_args();

random.seed(opts.seed);

# Open the input bam file
samfile = pysam.Samfile(opts.bam, 'r');

# Open output bamfiles
names  = [opts.prefix + "_" + str(i) + ".bam" for i in range(opts.num_splits)];
opsams = {name:pysam.AlignmentFile(name, 'wb', template=samfile) for name in names};

outstanding_templates = {name:{} for name in names};

# Traverse input bam
for i, read in enumerate(samfile.fetch()):
	query_name = read.query_name;
	target_bam = None;

	if (query_name is None) or (opts.conserve_mate_pairs is False):
		target_bam = opsams[names[random.randint(0, len(names)-1)]];
	else:
		query_name_found = False;

		for name in outstanding_templates:
			outstanding_template_dict = outstanding_templates[name];

			if query_name in outstanding_template_dict:
				query_name_found = True;
				target_bam       = opsams[name];

				# Not more than paired-end - save space
				del outstanding_template_dict[query_name];

				break;

		if not query_name_found:
			target_bam_name = names[random.randint(0, len(names)-1)];
			target_bam      = opsams[target_bam_name];

			outstanding_templates[target_bam_name][query_name] = 1;

	if target_bam is not None:
		target_bam.write(read);

	if (i+1) % 100000 == 0:
		print "Processed %d reads ..." %(i+1);

# Close files
bclose = lambda x : x.close();

map(bclose, opsams.values());
bclose(samfile);
