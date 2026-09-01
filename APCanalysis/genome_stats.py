import argparse

parser = argparse.ArgumentParser(description='get genomic intron and '
	'exon length stats')
parser.add_argument('genome')
parser.add_argument('annotation')

args = parser.parse_args()

tx_ids = {}
with open(args.annotation, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		if line[1] == 'WormBase' and line[2] == 'gene':
			attributes = line[8].split(';')
			att = dict([x.split('=') for x in attributes])
			if att['biotype'] == 'protein_coding':
				wbg = line[8].split(';')[0].split(':')[1]
				# all wbg are unique?
				tx_ids[att['sequence_name']] = wbg
	
exls = {}
inls = {}				
with open(args.annotation, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		if line[1] == 'WormBase':
			type_id = line[8].split(';')[0]
			if not type_id.split(':')[0].split('=')[1] == 'Transcript': 
				continue
			tx_id = '.'.join(type_id.split(':')[1].split('.')[:2])
			if tx_id not in tx_ids: continue
				
			if line[2] == 'exon':
				if tx_ids[tx_id] not in exls:
					exls[tx_ids[tx_id]] = {}
				exl = int(line[4]) - int(line[3]) + 1
				if exl not in exls[tx_ids[tx_id]]:
					exls[tx_ids[tx_id]][exl] = 1
				else:
					exls[tx_ids[tx_id]][exl] += 1		
			
			if line[2] == 'intron':
				if tx_ids[tx_id] not in inls:
					inls[tx_ids[tx_id]] = {}
				inl = int(line[4]) - int(line[3]) + 1
				if inl not in inls[tx_ids[tx_id]]:
					inls[tx_ids[tx_id]][inl] = 1
				else:
					inls[tx_ids[tx_id]][inl] += 1
			
print(exls)
print(inls)

# WBGene00022221
# 5000bp intron in the 3' utr??

for item in inls.items():
	print(item)

#print(dict(sorted(exls.items())))

#print(f'exin,length,count')
#for ln in sorted(exls.items()):
#	print(f'exon,{ln[0]},{ln[1]}')
	
#for ln in sorted(inls.items()):
#	print(f'intron,{ln[0]},{ln[1]}')









