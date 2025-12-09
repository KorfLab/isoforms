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
parser.add_argument('--max_len', required=False, type=int, 
	default=1000, metavar='<int>', 
	help='max length of sequence to test [%(default)i]')
parser.add_argument('--inc', required=False, type=int, default=100, 
	metavar='<int>', help='increment sequence length by [%(default)i]')
parser.add_argument('--reps', required=False, type=int, default=100,
	metavar='<int>', 
	help='number of times to repeat simulation [%(default)i]')

args = parser.parse_args()

if not args.data_dir.endswith('/'): 
	args.data_dir = args.data_dir + '/'

flank = 99
emin = 25
imin = 35

# find frequency of dinucleotides in real sequences
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
	
	#print(len(seq), len(dons), len(accs))
	genes[gid] = {'seq_len': len(seq), 'dcounts': len(dons), 
					'acounts': len(accs)}

k2s = [item[0] for item in k2_counts.items()]
k2_weights = [item[1] for item in k2_counts.items()]

# generate random sequences with GT/AG sites

with open('bigo.csv', 'w') as file:
	file.write(f'seq_len,iso_count,n_dons,n_accs,time\n')
	for i in range(300, args.max_len+args.inc, args.inc):
		for j in range(args.reps):
			ranseq = ''.join(random.choices(k2s, weights = k2_weights, 
							k = round(i/2)))
			rdons, raccs = isoform.gtag_sites(ranseq, flank, emin)
			#print(len(ranseq), len(rdons), len(raccs))
			# just use last instance of name from line 34
			print(f'working on {i}...')
			start = time.time()
			locus = Locus(name, ranseq, args.model, countonly=True)
			end = time.time()
			total_time = end - start
			file.write(f'{len(ranseq)},{locus.isocount},{len(rdons)},'
						f'{len(raccs)},{total_time}\n')
			file.flush()

'''
with open('bigo.out', 'w') as file:
	file.write(f'seq_len,iso_count,n_dons,n_accs\n')
	for i in range(300, args.max_len+args.inc, args.inc):
		ranseq = ''.join(random.choices(k2s, weights = k2_weights, 
						k = round(i/2)))
		rdons, raccs = isoform.gtag_sites(ranseq, flank, emin)
		#print(len(ranseq), len(rdons), len(raccs))
		# just use last instance of name from line 34
		print(f'working on {i}...')
		locus = Locus(name, ranseq, args.model, countonly=True)
		file.write(f'{len(ranseq)},{locus.isocount},{len(rdons)},'
					f'{len(raccs)}\n')
		file.flush()
'''


















'''
for i in range(300, args.max_len + args.inc, args.inc):
	name = f'{i}'
	ranseq = ''.join(random.choices(all_k2s, weights = all_k2_counts, 
					k = i))
	dons, accs = isoform.gtag_sites(ranseq, flank, emin)
	locus = Locus(name, ranseq, args.model, countonly=True)
	print(len(ranseq), locus.isocount, len(dons), len(accs))
'''
'''
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
'''	

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




	







