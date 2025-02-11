import argparse
import os
import subprocess
import multiprocessing as mp
from multiprocessing import Pool
import time

parser = argparse.ArgumentParser(description='wrapper for optiso2')
parser.add_argument('apc_dir', type=str, metavar='<directory>', 
    help='directory with APC .gff3 and .fa files')
parser.add_argument('model', type=str, metavar='<file>',
    help='splice model file')
parser.add_argument('--notes', required=False, type=str, metavar='<string>',
    help='any metadata you want to add')
parser.add_argument('--program', required=False, type=str, default='optiso2', 
    metavar='<executable>', help='path to optiso2')

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

cmds = []
for gID in fpaths:
    cmd = f'./{args.program} {fpaths[gID][0]} {fpaths[gID][1]} {args.model}'
    cmds.append(cmd)

s = time.perf_counter()
with open('results_optiso2_defaults.txt', 'w') as file:
    for c in cmds:
        c = c.split(' ')
        gid = c[1].split('.')[-2]
        print(f'working on {gid}...')
        result = subprocess.run(c, stdout=subprocess.PIPE, text=True)
        print(result)
e = time.perf_counter()
print('no multi:', e-s)

def optimize(prog, fasta, gff, model):

    cmd = f'./{prog} {fasta} {gff} {model}'
    cmd = cmd.split(' ')
    gid = cmd[1].split('.')[-2]
    output = subprocess.run(cmd, stdout=subprocess.PIPE, text=True).stdout.split()
    return [gid] + output

def worker(inputs):

    return optimize(inputs[0], inputs[1], inputs[2], inputs[3])

inputs = []
for gID in fpaths:
    input = [args.program, fpaths[gID][0], fpaths[gID][1], args.model]
    inputs.append(input)

if args.notes: print(f'# {args.notes}')

s = time.perf_counter()
with Pool(processes=mp.cpu_count()) as pool:
    result = pool.map(worker, inputs)
    print(result)
e = time.perf_counter()
print('multi:', e-s)

