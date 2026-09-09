import sys

optres = sys.argv[1]

ids = {}
with open(optres, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split(',')
		if line[0] == 'gid': continue
		idinfo = line[0].split('.')
		if idinfo[0] not in ids:
			ids[idinfo[0]] = [int(idinfo[1])]
		else:
			ids[idinfo[0]].append(int(idinfo[1]))
			
sorted_ids = {k: sorted(v) for k, v in ids.items()}

print(sorted_ids['1'])
