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
parser.add_argument('model', type=str, metavar='<file>',
	help='worm splice model file')
parser.add_argument('--max_len', required=False, type=str, default=1000, metavar='<int>', 
	help='max length of sequence to test [%(default)i]')
parser.add_argument('--inc', required=False, type=str, default=100, metavar='<int>', 
	help='increment sequence length by [%(default)i]')

args = parser.parse_args()

if not args.data_dir.endswith('/'): 
	args.data_dir = args.data_dir + '/'

flank = 99
emin = 25
imin = 35

big_sum = 0
all_dons = 0
all_accs = 0

genes = {}
k2_counts = {}
for fasta in glob.glob(f'{args.data_dir}**.fa'):
	name, seq = next(isoform.read_fasta(fasta))
	gid = name.split(' ')[0]
	dons, accs = isoform.gtag_sites(seq, flank, emin)
	
	for i in range(flank+emin, len(seq)-flank-emin):
		k2 = seq[i:i+2]
		if len(k2) == 2:
			if k2 not in k2_counts:
				k2_counts[k2] = 1
			else:
				k2_counts[k2] += 1
	
	dfreq = round((len(dons) / (len(seq)/2)), 2)
	afreq = round((len(accs) / (len(seq)/2)), 2)
	genes[gid] = {'seq_len': len(seq), 'dfreq': dfreq, 'afreq': afreq}
	print(gid, len(seq), len(dons), len(accs))
	big_sum += len(seq)
	all_dons += len(dons)
	all_accs += len(accs)
	
print(big_sum, all_dons, all_accs)	


#ce.1.340 1007 24 30
#ce.4.6 1024 26 44

# GTs per 100 bp
print((all_dons/big_sum) * 100)

dfreqs = [item[1]['dfreq'] for item in genes.items()]
afreqs = [item[1]['afreq'] for item in genes.items()]

global_dfreq = stats.mean(dfreqs)
global_afreq = stats.mean(afreqs)	
				
all_k2s = [item[0] for item in k2_counts.items()]
all_k2_counts = [item[1] for item in k2_counts.items()]

#all_k2_counts[9] = round(all_k2_counts[9] / 10)
#all_k2_counts[8] = round(all_k2_counts[8] / 10)

print((all_k2_counts[9]/sum(all_k2_counts)) * 100)

#print(all_k2_counts[9], all_k2_counts[8])

print('########')

# this section gets new don/acc freqs just for comparison
ran_dfreqs = []
ran_afreqs = []
for i in range(300, args.max_len, args.inc):
	ranseq = ''.join(random.choices(all_k2s, 
									weights = all_k2_counts, k = i))
	rdons, raccs = isoform.gtag_sites(ranseq, flank, emin)
	print(len(ranseq), len(rdons), len(raccs))
	rdfreq = round((len(rdons) / (len(ranseq)/2)), 2)
	rafreq = round((len(raccs) / (len(ranseq)/2)), 2)
	ran_dfreqs.append(rdfreq)
	ran_afreqs.append(rafreq)
	
# why doesn't simulated freqs match actual? 
# simulation used actual freqs?
'''
print(global_dfreq, global_afreq)
print(stats.mean(ran_dfreqs), stats.mean(ran_afreqs))
'''

'''
# now do apc
# takes too long
for i in range(300, args.max_len + args.inc, args.inc):
	name = f'{i}'
	ranseq = ''.join(random.choices(all_k2s, weights = all_k2_counts, 
					k = i))
	dons, accs = isoform.gtag_sites(ranseq, flank, emin)
	locus = Locus(name, ranseq, args.model, countonly=True)
	print(len(ranseq), locus.isocount, len(dons), len(accs))
'''

# average number of donor/acceptor sites by length of gene?




	







