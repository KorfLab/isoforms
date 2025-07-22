import argparse
import gzip

parser = argparse.ArgumentParser(description='gathers intron, donor, ' 
	'and acceptor sequences for PWM training')
parser.add_argument('genome')
parser.add_argument('gff')
parser.add_argument('--c', action='store_true', 
	help='gather only canonical GT/AG splice sites') 
parser.add_argument('--flank', required=False, type=int, 
	default=5, help='gather n amount of nucleotides on either side '
	'of intron boundary [%(default)i]')

args = parser.parse_args()

# i used to chatgpt to write this grep command...
# gets pattern match to 3rd field in tab delimited file
'''
grep -P "^([^\t]*\t){2}intron\b" annotation.gff3
'''

def revcomp(seq):

	comps = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

	rev = ''
	for i in range(1, len(seq) + 1):
		rev += comps[seq[-i]]

	return rev
	
if args.gff.endswith('.gz'):
	gff_fp = gzip.open(args.gff, 'rt')
else:
	gff_fp = open(args.gff, 'rt')
	
introns = []
exons = []
for line in gff_fp:
	line = line.split('\t')
	# need to search for either WormBase or RNAseq_splice
	# this just gets either, bad
	# coordinate could appear twice from different transcripts
	# then you accidently double count
	# the counts on RNAseq_splice are associated with the PROMOTER
	# not the splice site
	# in region you are looking at, what is the one count?
	# normalize the counts to region you are looking at
	# small genes are already regionized
	# highest count gets one count, rest gets partial count
	# weight from 0-1
	# compare all different methods in smallgenes set instead of 
	# whole genome. much easier
	if len(line) < 9: continue
	if line[2] == 'intron':
		intron = (line[0], int(line[3]), int(line[4]), line[6])
		introns.append(intron)
		#print(line)

gff_fp.close()

if args.genome.endswith('.gz'):
	genome_fp = gzip.open(args.genome, 'r')
else:
	genome_fp = open(args.genome, 'r')
	
# grep whole annotation.gff3 to keep only things of interest
# makes file smaller
	
chroms = {}
#or line in genome_fp.readlines():
for line in genome_fp:
	try:
		line = line.decode('utf-8')
		line = line.rstrip()
	except:
		line = line.rstrip()
	if line.startswith('>'):
		chrom = line.split(' ')[0][1:]
	else:
		if chrom not in chroms:
			chroms[chrom] = line
		else:
			chroms[chrom] += line
	
genome_fp.close()

'''
## collect subsequences
introns = []
exons = []
with open(args.gff, 'r') as gp:
	for line in gp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		if len(line) < 9: continue
		if line[2] == 'intron':
			intron = (line[0], int(line[3]), int(line[4]), line[6])
			introns.append(intron)
		if line[2] == 'exon':
			exon = (line[0], int(line[3]), int(line[4]), line[6])
			exons.append(exon)

chroms = {}
with open(args.genome, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			chrom = line.split(' ')[0][1:]
		else:
			if chrom not in chroms:
				chroms[chrom] = line
			else:
				chroms[chrom] += line
'''
accs = []
dons = []
iseqs = []
for i in introns:
	beg = i[1]-1
	end = i[2]
	if i[3] == '-':
		miseq = chroms[i[0]][beg-5:end+5]
		iseq = revcomp(miseq)
		print(iseq[5:7])
	else:
		iseq = chroms[i[0]][beg-5:end+5]
	aseq = iseq[-6:]
	dseq = iseq[:5]
	if args.c:
		if aseq.endswith('AG'):
			accs.append(aseq)
		if dseq.startswith('GT'):
			dons.append(dseq)
	else:
		accs.append(aseq)
		dons.append(dseq)
	iseqs.append(iseq)

with open('introns.txt', 'w') as ifp:
	for iseq in iseqs:
		ifp.write(f'{iseq}\n')

'''
with open('donors.txt', 'w') as dfp:
	for dseq in dons:
		dfp.write(f'{dseq}\n')

with open('acceptors.txt', 'w') as afp:
	for aseq in accs:
		afp.write(f'{aseq}\n')
'''





