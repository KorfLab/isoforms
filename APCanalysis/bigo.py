import argparse
import glob
import isoform
from isoform import Locus
import statistics as stats
import random
import time

parser = argparse.ArgumentParser(
	description='what is the time complexity of the APC algorithm?')
parser.add_argument('data_dir', type=str, metavar='<directory>', 
	help='directory with input fasta files')

args = parser.parse_args()

if not args.data_dir.endswith('/'): 
	args.data_dir = args.data_dir + '/'

flank = 99
emin = 25
imin = 35

genes = {}
for fasta in glob.glob(f'{args.data_dir}**.fa'):
	name, seq = next(isoform.read_fasta(fasta))
	gid = name.split(' ')[0]
	dons, accs = isoform.gtag_sites(seq, flank, emin)
	dfreq = round((len(dons) / (len(seq)/2)), 2)
	afreq = round((len(accs) / (len(seq)/2)), 2)
	# need to add time for apc
	genes[gid] = {'seq_len': len(seq), 'dfreq': dfreq, 'afreq': afreq}

dfreqs = [item[1]['dfreq'] for item in genes.items()]
afreqs = [item[1]['afreq'] for item in genes.items()]

global_dfreq = stats.mean(dfreqs)
global_afreq = stats.mean(afreqs)	
	
kmer_counts = {}
for fasta in glob.glob(f'{args.data_dir}**.fa'):
	name, seq = next(isoform.read_fasta(fasta))
	for i in range(len(seq)):
		kmer = seq[i:i+2]
		if len(kmer) == 2:
			if kmer not in kmer_counts:
				kmer_counts[kmer] = 1
			else:
				kmer_counts[kmer] += 1 
				
kmers = [item[0] for item in kmer_counts.items()]
kcounts = [item[1] for item in kmer_counts.items()]

print(kmers)
print(kcounts)
				

	
'''
start = time.time()
locus = Locus(name, seq, model, countonly=True)
end = time.time()
elapsed = end - start
print(locus.isocount, elapsed)
print(len(seq), len(dons), len(accs))
'''

'''
def ranseq(length):
	
	seq = ''
	for i in range(length):
		seq += random.choice(['A', 'C', 'G', 'T'])

	return seq
	
	
print(ranseq(100))

max_len = 1000
inc = 100
for i in range(300, max_len + inc, inc):
	seq = ranseq(i)
	d, a = isoform.gtag_sites(seq, flank, minex)
	print(i, len(d), len(a))
	
# frequency of GT and AG sites in real data?
# build seqs in 2mers, with options of AC, AG, AT, CG, CT

alph = ['a', 'b', 'c']
for i in range(len(alph)):
	for j in range(len(alph)):
		print(alph[i], alph[j])
		
alph = ['A', 'C', 'G', 'T']
dinucs = []
for i in range(len(alph)):
	for j in range(len(alph)):
		dinucs.append(f'{alph[i]}{alph[j]}')

dinuc_counts = {x:0 for x in dinucs}

print(dinuc_counts)		
'''
	







