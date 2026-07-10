import argparse
import isoform

parser = argparse.ArgumentParser(description='look at the splice sites '
	'of poorly predicted genes')
parser.add_argument('apc_gff')
parser.add_argument('fa')
parser.add_argument('model')

args = parser.parse_args()

# 1.542

isos = {}
with open(args.apc_gff, 'rt') as fp:
	iso_count = 0
	for line in fp:
		line = line.rstrip()
		if line.startswith('#') or line == '': continue
		line = line.split('\t')
		if line[2] == 'gene': continue
		if line[2] == 'mRNA':
			iso_count += 1
		if iso_count not in isos:
			isos[iso_count] = {}
			isos[iso_count]['exon'] = []
			isos[iso_count]['intron'] = []
		if line[2] == 'exon':
			isos[iso_count]['exon'].append((line[3], line[4]))
		if line[2] == 'intron':
			isos[iso_count]['intron'].append((line[3], line[4]))
			
for i in isos.items():
	print(i)
	break

print('####')
	
seq = []	
with open(args.fa, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		if line.startswith('>'): continue
		for n in line:
			seq.append(n)
			
print(isos[1])

model = isoform.read_splicemodel(args.model)

print(model['don'])
print(model['acc'])


for intron in isos[1]['intron']:
	print(intron)
	beg = int(intron[0])
	dseq = ''.join(seq[beg-1:beg+5])
	dscore = isoform.score_pwm(model['don'], dseq, 0, memo=None)
	print(dseq, dscore)
			
print(seq[260:266])
		
	
print('####')
# need to score the splice sites again

# isoform.score_pwm(pwm, seq, pos, memo=None)
# isoform.read_splicemodel(file)



