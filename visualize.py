import pysam
import ast
import sys
from optparse import OptionParser

def readNormalSeq(referenceFile, referenceRegion):
    variantGenotypes = {};
    variantAlleles   = {};
    referenceSeq     = [];
    normalSeq        = [];

    ##vcfRegion = 'chr' + referenceRegion;
    #vcfRegion = 'G15512.prenormal.sorted';
    #vcf       = pysam.VCF();
    #vcf.connect(vcfFile);

    #1. Read the reference sequence
    reference = pysam.FastaFile(referenceFile);
    for base in reference.fetch(referenceRegion):
        referenceSeq.append(base);        

    ##2. For each variant position, record the reference base as allele 0, and others starting from allele 1
    #for variant in vcf.fetch(referenceRegion):
    #    variantAlleles[variant.pos] = [];
    #    variantAlleles[variant.pos].append(referenceSeq[variant.pos]);
    #    for alt in variant.alt:
    #        variantAlleles[variant.pos].append(variant.alt);

    #    #Record the variant genotype
    #    variantGenotypes[variant.pos] = variant[vcfRegion]['GT'][0]

    ##3. Put them together based on the genotype at every location
    #for pos, base in enumerate(referenceSeq):
    #    positionalEntry = [];
    #    if pos in variantAlleles:
    #        genotype0 = variantGenotypes[pos][0];
    #        genotype1 = variantGenotypes[pos][2];

    #        if isinstance(variantAlleles[pos][genotype0], list):
    #            allele0 = variantAlleles[pos][genotype0][0];
    #        else:
    #            allele0 = variantAlleles[pos][genotype0];

    #        if isinstance(variantAlleles[pos][genotype1], list):
    #            allele1 = variantAlleles[pos][genotype1][0];
    #        else:
    #            allele1 = variantAlleles[pos][genotype1];

    #        positionalEntry.append(allele0);

    #        if allele0 != allele1:
    #            positionalEntry.append(allele1);

    #    else:
    #        positionalEntry.append(referenceSeq[pos]);

    #    normalSeq.append(positionalEntry);

    return referenceSeq, referenceSeq;

options = OptionParser();

options.add_option('--pos', action="store", type="int", dest="variantPosition");
options.add_option("--ref", action="store", type="string", dest="referenceFile");
options.add_option("--bam", action="store", type="string", dest="alignmentFile");
options.add_option("--chr", action="store", type="string", dest="referenceRegion");
options.add_option("--one_start", action="store_true", dest="one_start", default=False);

(opts, args) = options.parse_args();

variantPosition = opts.variantPosition;
referenceFile   = opts.referenceFile;
alignmentFile	= opts.alignmentFile; #"dream3.chr19.bqsr.ind.normal.bam" if sys.argv[-1] == 'normal' else 'dream3.chr19.bqsr.ind.tumor.bam';
referenceRegion = opts.referenceRegion;

if (variantPosition is None) or (referenceFile is None) or (alignmentFile is None) or (referenceRegion is None):
	sys.exit("Usage: python " + sys.argv[0] + " --pos <position whose pileup should be displayed> --ref <reference fasta> --bam <alignment file> --chr <chromosome> [--one_start, if position is 1-based]");

variantPosition = variantPosition-1 if opts.one_start else variantPosition;

sequenceRange   = [];

germlineSeq, referenceSeq = readNormalSeq(referenceFile,referenceRegion);
bamfile         = pysam.Samfile(alignmentFile);

margin = 10;

for pileupcolumn in bamfile.pileup(referenceRegion, variantPosition-margin/2, variantPosition+margin/2):
    sequenceRange.append(pileupcolumn.pos);

startPosition = min(sequenceRange);
endPosition   = max(sequenceRange);
dataset       = [];
referenceMap  = range(startPosition,endPosition+1);

for i in range(startPosition,endPosition+1):
    referenceMap[i-startPosition] = {'base':referenceSeq[i], 'columns':[], 'quality':[], 'spaces':0};

readContainer = [];
pileupreadList = [];

#Collect all the reads spanning the region and sort them based on start position
for pileupcolumn in bamfile.pileup(referenceRegion, startPosition, endPosition):
    for pileupread in pileupcolumn.pileups:
        if startPosition <= pileupread.alignment.reference_start <= pileupread.alignment.reference_end <= endPosition:
            readWithStartPos = (pileupread.alignment.reference_start, pileupread);
            read             = pileupread.alignment.query_sequence;
            if read not in pileupreadList:
                readContainer.append(readWithStartPos);
                pileupreadList.append(read);

print "Number of reads in container: " + str(len(readContainer));

readContainer = sorted(readContainer);
readsLabeled  = [r for (i,r) in readContainer];

#Initialize the reference map for each column for all reads
for i in range(startPosition,endPosition+1):
    for r in readsLabeled:
        referenceMap[i-startPosition]['columns'].append(" ");
        referenceMap[i-startPosition]['quality'].append(-1);

mapping_quality = [];

num_forward_strand_reads = 0;
num_reverse_strand_reads = 0;

strand = [];

#Iterate through each read and fill up the reference map
for index, pileupread in enumerate(readsLabeled):
    cigar     = pileupread.alignment.cigartuples;
    read      = pileupread.alignment.query_alignment_sequence;
    quality   = pileupread.alignment.query_alignment_qualities;
    position  = pileupread.alignment.reference_start;
    readPos   = 0;

    mapping_quality.append(pileupread.alignment.mapping_quality);

    if pileupread.alignment.is_reverse:
        num_reverse_strand_reads += 1;
    else:
        num_forward_strand_reads += 1;

    strand.append(1 if pileupread.alignment.is_reverse is True else 0);

    for desc in cigar:
        #Match
        if (desc[0] == 0) or (desc[0] == 7) or (desc[0] == 8):
            for j, i in enumerate(range(position,position+desc[1])):
                referenceMap[i-startPosition]['columns'][index] = read[j+readPos];
                referenceMap[i-startPosition]['quality'][index] = quality[j+readPos];

            position += desc[1];
            readPos  += desc[1];
            next;

        #Insert - index moves along read
        if desc[0] == 1:
            referenceMap[i-startPosition]['columns'][index] = read[readPos:readPos+desc[1]];
            referenceMap[i-startPosition]['quality'][index] = quality[readPos:readPos+desc[1]];
            readPos  += desc[1];
            next;
    
        #Delete - do nothing except advance along the reference - deletions are marked by default
        if desc[0] == 2:
            position += desc[1];
            next;
             
        #Never mind hard-clipped or soft-clipped or padded portions

#The list of strings for each position
lines = ["" for r in readsLabeled];

#Finalize the number of positions given to each base
for pos in referenceMap:
    pos['spaces'] = max([len(i) for i in pos['columns']]);

    for i, vertical in enumerate(pos['columns']):
        lines[i] += vertical;
 
        if pos['spaces'] - len(vertical) > 0:
            for j in range(pos['spaces'] - len(vertical)):
                lines[i] += " ";

lines.insert(0,"");

#Prepend the reference sequence and the base positions
for i, pos in enumerate(referenceMap):
    lines[0] += pos['base'];
    
    if pos['spaces'] > 1:
        for j in range(pos['spaces'] - 1):
            lines[0] += " ";

for line in lines:
    print line;


#Now print the column specific data: index: variantPosition - startPosition
data = [(v,s) for (v,s) in zip(referenceMap[variantPosition-startPosition]['columns'], referenceMap[variantPosition-startPosition]['quality']) if s > 0];
print "\nColumn data for column (" + referenceMap[variantPosition-startPosition]['base'] + ")" + str(variantPosition-startPosition) + " is : " + str(data);


print "\n Mapping qualities are : " + str(mapping_quality);

print "\nNum forward strand reads : " + str(num_forward_strand_reads);
print "Num reverse strand reads: " + str(num_reverse_strand_reads);

print "Strand : " + str(strand);
