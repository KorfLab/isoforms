import argparse
import glob
import apc_analysis as aa
from apc_analysis import SpliceSites
from grimoire.genome import Reader

parser = argparse.ArgumentParser(
	description='creates PWMs from smallgenes')
parser.add_argument('smallgenes', type=str, metavar='<directory>', 
	help='path to smallgenes directory')
parser.add_argument('--don_len', required=False, type=int, default=5, 
	metavar='<integer>', help='donor site length [%(default)i]')
parser.add_argument('--acc_len', required=False, type=int, default=6, 
	metavar='<integer>', help='acceptor site length [%(default)i]')
parser.add_argument('--don_left', required=False, type=int, default=0, 
	metavar='<integer>', 
	help='get bases on donor site left flank [%(default)i]')
parser.add_argument('--acc_right', required=False, type=int, default=0, 
	metavar='<integer>', 
	help='get bases on acceptor site right flank [%(default)i]')

args = parser.parse_args()	

'''
make sure to filter only one isoform from each gene
some genes have more than one isoform in smallgenes 
one-based counting is not good for PWMs
remove non-coding RNA genes from genomic PWM training sequences
'''

# shorten arguments for indexing
DN = args.don_len
AN = args.acc_len
DL = args.don_left 
AR = args.acc_right 
	
if args.smallgenes.endswith('/'):
	pass 
else:
	args.smallgenes = args.smallgenes + '/'
	

all_rel_dons = []	
all_rel_accs = []
for ff in glob.glob(f'{args.smallgenes}/*.fa'):
	gf = ff[:-2] + 'gff3'
	genome = Reader(gff=gf, fasta=ff)
	chrom = next(genome)
	gene = chrom.ftable.build_genes()[0]
	tx = gene.transcripts()[0]
	gseq = gene.seq_str()

	# started writing code for one-based counting
	'''
	exon_coors = []
	for exon in tx.exons:
		exon_coors.append((exon.beg, exon.end))
	
	for i, intron in enumerate(tx.introns):
		intron_coor = (intron.beg, intron.end)
		#print(exon_coors[i], intron_coor, intron.score)
	'''	
		
	total_score = 0
	don_acc_sites = []
	for f in chrom.ftable.features: 
		if f.source == 'RNASeq_splice':
			total_score += f.score
			don = gseq[f.beg-100-DL:f.beg-100+DN]
			acc = gseq[f.end-99-AN:f.end-99+AR]
			don_acc = (don, acc, f.score)
			don_acc_sites.append(don_acc)
			
	# get relative counts
	for sites in don_acc_sites:
		rdon = [sites[0], sites[2] / total_score]
		racc = [sites[1], sites[2] / total_score]
		all_rel_dons.append(rdon)
		all_rel_accs.append(racc)
	
rel_dpwm = aa.build_weighted_pwm(all_rel_dons, len(all_rel_dons[0][0]))
rel_apwm = aa.build_weighted_pwm(all_rel_accs, len(all_rel_accs[0][0]))
	
aa.print_pwm(rel_dpwm, 'smallgenes_rnaseq_relative_donpwm')
aa.print_pwm(rel_apwm, 'smallgenes_rnaseq_relative_accpwm')




# old method
# use grimoire instead to standardize inputs
'''
parent_txs = {}
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
	# smallgenes converts - strand gene features to +
	# should not be any - features in the g
	with open(gf, 'rt') as fp:
		splice_site_set = {}
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			
			if line[6] == '-': continue
			
			# for ce.3.58 (and probably others), the same intron has 3 
			# different parent transcripts...so different transcripts 
			# overlap the same intron
			if line[1] == 'WormBase' and line[2] == 'intron':
				parent_tx = line[8].split(':')[1]
				donor = seq[int(line[3])-DL-1:int(line[3])+DN-1]
				acceptor = seq[int(line[4])-AN:int(line[4])+AR]
				if parent_tx not in parent_txs:
					parent_txs[parent_tx] = []
					parent_txs[parent_tx].append((donor, acceptor))
				else:
					parent_txs[parent_tx].append((donor, acceptor))
			if line[1] == 'RNASeq_splice' and line[2] == 'intron':
				donor = seq[int(line[3])-3:int(line[3])+5]
				acceptor = seq[int(line[4])-8:int(line[4])]
				splice_site_set[line[5]] = (donor, acceptor) 
				
		# assign splice sites to score
		rnaseq_splice_sites[line[0]] = splice_site_set
		
# only keep sites from first isoform (parent transcript) seen
annotated_splice_sites = []
seen = []
for parent in parent_txs:
	base = parent.split('.')[0]
	if base in seen: continue
	seen.append(base)
	for sites in parent_txs[parent]:
		annotated_splice_sites.append(sites)

wb_dons = []
wb_accs = []
rna_dons = []
rna_accs = []
for file in glob.glob(f'{args.smallgenes}*.fa'):
	ff = file
	gf = ff[:-2] + 'gff3'
	wb_sites = SpliceSites(ff, gf, DN, DL, AN, AR, source='WormBase')
	rna_sites = SpliceSites(ff, gf, DN, DL, AN, AR, source='RNASeq')
	for d, a in wb_sites.splice_sites():
		wb_dons.append(d)
		wb_accs.append(a)
	for d, a, s in rna_sites.splice_sites():
		rna_dons.append([d, s])
		rna_accs.append([a, s])
		
wb_dons = []
for sites in annotated_splice_sites:
	print(sites[0])
	wb_dons.append(sites[0])		
		
wb_dpwm = aa.build_pwm(wb_dons, len(wb_dons[0]))		
wb_apwm = aa.build_pwm(wb_accs, len(wb_accs[0]))
rna_dpwm = aa.build_weighted_pwm(rna_dons, len(rna_dons[0][0]))
rna_apwm = aa.build_weighted_pwm(rna_accs, len(rna_accs[0][0]))
'''


'''
aa.print_pwm(wb_dpwm, 'smallgenes_wormbase_donor_pwm')
aa.print_pwm(wb_apwm, 'smallgenes_wormbase_acceptor_pwm')
aa.print_pwm(rna_dpwm, 'smallgenes_rnaseq_donor_pwm')
aa.print_pwm(rna_apwm, 'smallgenes_rnaseq_acceptor_pwm')
'''







