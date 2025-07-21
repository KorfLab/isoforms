import argparse
import glob
from apc_analysis import SpliceSites

parser = argparse.ArgumentParser(
	description='creates PWMs from smallgenes')
parser.add_argument('smallgenes', type=str, metavar='<directory>', 
	help='path to smallgenes directory')
parser.add_argument('--don_len', required=False, type=int, default=5, 
	metavar='<integer>', help='donor site length')
parser.add_argument('--acc_len', required=False, type=int, default=5, 
	metavar='<integer>', help='acceptor site length')
parser.add_argument('--don_left', required=False, type=int, default=0, 
	metavar='<integer>', 
	help='get bases on donor site left flank [%(default)i]')
parser.add_argument('--acc_right', required=False, type=int, default=0, 
	metavar='<integer>', 
	help='get bases on acceptor site right flank [%(default)i]')

args = parser.parse_args()	

# shorten arguments for indexing
DN = args.don_len
AN = args.acc_len
DL = args.don_left 
AR = args.acc_right 

# make sure parameters work
#print(f'DN: {DN} AN: {AN} DL: {DL} AR: {AR}')
	
if args.smallgenes.endswith('/'):
	pass 
else:
	args.smallgenes = args.smallgenes + '/'
	
annotated_splice_sites = []
rnaseq_splice_sites = {}
for file in glob.glob(f'{args.smallgenes}*.fa'):
	ff = file
	gf = ff[:-2] + 'gff3'
	seq = ''
	with open(ff, 'rt') as fp:
		for line in fp:
			line = line.rstrip()
			if line.startswith('>'): continue
			else:
				seq += line
	# there are - strand features mixed in for 13 genes
	# only get + strand features
	with open(gf, 'rt') as fp:
		splice_site_set = {}
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			
			if line[6] == '-': continue
			
			if line[1] == 'WormBase' and line[2] == 'intron':
				donor = seq[int(line[3])-DL-1:int(line[3])+DN-1]
				acceptor = seq[int(line[4])-AN:int(line[4])+AR]
				print(donor, acceptor)
				annotated_splice_sites.append((donor, acceptor))
				
			if line[1] == 'RNASeq_splice' and line[2] == 'intron':
				donor = seq[int(line[3])-3:int(line[3])+5]
				acceptor = seq[int(line[4])-8:int(line[4])]
				splice_site_set[line[5]] = (donor, acceptor) 
				
		# assign splice sites to score
		rnaseq_splice_sites[line[0]] = splice_site_set
		
#print(annotated_splice_sites)

#print(rnaseq_splice_sites)

for file in glob.glob(f'{args.smallgenes}*.fa'):
	ff = file
	gf = ff[:-2] + 'gff3'
	wb_pwm = SpliceSites(ff, gf, DN, DL, AN, AR, source='WormBase')
	rna_pwm = SpliceSites(ff, gf, DN, DL, AN, AR, source='RNASeq')
	for d, a in wb_pwm.splice_sites():
		print(d, a)
	for d, a, s in rna_pwm.splice_sites():
		print(d, a, s)
		

	
	


