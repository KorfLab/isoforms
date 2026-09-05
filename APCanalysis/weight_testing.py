import argparse
from copy import deepcopy

parser = argparse.ArgumentParser(description='remove certain weights from optiso')
parser.add_argument('optiso_res')

args = parser.parse_args()

no_emm_all = {}
no_imm_all = {}
no_exl_all = {}
no_inl_all = {}
header = ''
with open(args.optiso_res, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		if line.startswith('gid'):
			header = line
			continue
		line = line.split(',')
		
		# remove Markov model weights
		no_emm = deepcopy(line)
		no_imm = deepcopy(line)
		no_emm[4] = str(0)
		no_imm[5] = str(0)
		
		# remove length model weights
		no_exl = deepcopy(line)
		no_inl = deepcopy(line)
		no_exl[6] = str(0)
		no_inl[7] = str(0)
		
		no_emm_all[line[0]] = no_emm
		no_imm_all[line[0]] = no_imm
		no_exl_all[line[0]] = no_exl
		no_inl_all[line[0]] = no_inl

# write 4 different files
def write_csv(all_weights, header, fname):
	
	with open(fname, 'wt') as f:
		f.write(header+'\n')
		for item in all_weights.items():
			f.write(','.join(item[1])+'\n')

write_csv(no_emm_all, header, 'no_emm_optiso.csv')
write_csv(no_imm_all, header, 'no_imm_optiso.csv')
write_csv(no_exl_all, header, 'no_exl_optiso.csv')
write_csv(no_inl_all, header, 'no_inl_optiso.csv')









