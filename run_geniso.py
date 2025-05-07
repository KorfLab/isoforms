import argparse
import os
import subprocess
import multiprocessing as mp
from multiprocessing import Pool
import time

parser = argparse.ArgumentParser(description='wrapper for geniso3')
parser.add_argument('apc_dir', type=str, metavar='<directory>', 
    help='directory with APC .gff3 and .fa files')
parser.add_argument('model', type=str, metavar='<file>',
    help='splice model file')
parser.add_argument('--outdir', required=False, type=str, default='APCisos/',
    metavar='<outdir>', help='name of output directory [%(default)s]')
parser.add_argument('--limit', required=False, type=int, default=100,
	metavar='<int>', help='limit number of transcripts [%(default)i]')
parser.add_argument('--weights', type=str, metavar='<file>',
    help='file with model weights')
parser.add_argument('--program', required=False, type=str, default='geniso3', 
    metavar='<executable>', help='path to geniso3')

args = parser.parse_args()

fpaths = {}
for file in os.listdir(args.apc_dir):
    gID = file.split('.')[1]
    if gID not in fpaths:
        fpaths[gID] = [(), ()]
    if file.endswith('.fa'):
        fpaths[gID][0] = f'{args.apc_dir}{file}'
    if file.endswith('.gff3'):
        fpaths[gID][1] = f'{args.apc_dir}{file}'

weights = {}
with open(args.weights, 'r') as file:
    for line in file.readlines():
        line = line.rstrip()
        if line.startswith('%'): continue
        line = line.split(',')
        weights[line[0]] = [x for x in line[1:]]

def generate(
        prog, fasta, model, wacc, wdon, wexs, wins, wexl, winl, winf, limit
        ):

    cmd = (
        f'./{prog} {fasta} {model} --wacc {wacc} --wdon {wdon} '
        f'--wexs {wexs} --wins {wins} --wexl {wexl} --winl {winl} '
        f'--winf {winf} --limit {limit}'
    )
    cmd = cmd.split(' ')
    gid = cmd[1].split('.')[-2]
    output = subprocess.run(cmd, stdout=subprocess.PIPE, text=True).stdout
    print(f'working on {gid}...')
    return output.rstrip()

def worker(input):

    return generate(
        input[0], input[1], input[2], input[3], input[4],
        input[5], input[6], input[7], input[8], input[9], input[10]
    )
    
inputs = []
for gID in fpaths:
    wacc = weights[gID][1]
    wdon = weights[gID][2]
    wexs = weights[gID][3]
    wins = weights[gID][4]
    wexl = weights[gID][5]
    winl = weights[gID][6]
    winf = weights[gID][7]
    input = [
        args.program, fpaths[gID][0], args.model, wacc, wdon, wexs, wins,
        wexl, winl, winf, args.limit
    ]
    inputs.append(input)

#s = time.perf_counter()
with Pool(processes=mp.cpu_count()-1) as pool:
    result = pool.map(worker, inputs)
    #print(result)
#e = time.perf_counter()
#print('multi:', e-s)

start = time.perf_counter()
with Pool(processes=mp.cpu_count()-1) as pool:
    result = pool.map(worker, inputs)

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

for res in result:
    gID = res.split('\n')[0].split(' ')[2]
    f = open(f'{args.outdir}{gID}.APC.gff', 'w')
    for line in res.split('\n'):
        f.write(line+'\n')         
    f.close()

end = time.perf_counter()
print('time:', end-start)
