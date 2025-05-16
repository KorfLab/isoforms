import files
import os
import random
import sys


genelist = sys.argv[1]
outdir = sys.argv[2]
samplesize = int(sys.argv[3])
samplecount = sys.argv[4]
scount = int(samplecount)

form = '0' + str(len(samplecount)) + 'd'

filelen = 0

os.makedirs(outdir, exist_ok=True)

with files.getfp(genelist) as glfp:
    for line in glfp:
        filelen += 1

for i in range(scount):
    ftracker = 0
    sampnums = random.sample(range(filelen), k=samplesize)
    sampnums.sort()
    nxline = sampnums.pop(0)
    print(nxline, sampnums)
    with files.getfp(genelist) as glfp, open(f'{outdir}/sample_{i + 1:{form}}.txt', 'w') as ofp:
        for line in glfp:
            ftracker += 1
            if ftracker == nxline:
                print(line.strip(), file=ofp)
                if len(sampnums) > 0:
                    nxline = sampnums.pop(0)
                else:
                    break

