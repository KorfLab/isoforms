import sys
import glob

sg_dir = sys.argv[1]

seq_lens = {}
for fasta in glob.glob(f'{sg_dir}*.fa'):
	g_name = ''.join(fasta.split('.')[-2:-1])
	with open(fasta, 'rt') as fp:
		for line in fp:
			line = line.rstrip()
			if line.startswith('>'): continue
			if g_name not in seq_lens:
				seq_lens[g_name] = 0
				seq_lens[g_name] += len(line)
			else:
				seq_lens[g_name] += len(line)

print(seq_lens)
