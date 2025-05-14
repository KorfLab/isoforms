import glob
import subprocess
import statistics
import sys
import time

if len(sys.argv) == 1: sys.exit('provide a limit')

limit = int(sys.argv[1])

params = '--apwm models/acc.pwm --dpwm models/don.pwm --emm models/exon.mm --imm models/intron.mm --elen models/exon.len --ilen models/intron.len'

t0 = time.time()
ratios = []
for n, ff in enumerate(glob.glob(f'apc/*.fa')):
	t0 = time.time()
	subprocess.run(f'isoformer {ff} {params} > /dev/null', shell=True)
	t1 = time.time()
	subprocess.run(f'geniso {ff} {params} > /dev/null', shell=True)
	t2 = time.time()
	e1 = t1 - t0
	e2 = t2 - t1
	print(ff, e1, e2, e2/e1)
	ratios.append(e2/e1)

	if n+1 == limit: break

print(statistics.mean(ratios), '+/-', statistics.stdev(ratios))
