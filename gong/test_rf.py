"""
validation test by using dev data to test rf lib

    result:
        the isoformtree do behave as a full split tree
        which that leaf nodes == splits + 1
        aka the random forest approach do success
"""

import random
import argparse

from lib import randomf
from lib.randomf import IsoformTree

parser = argparse.ArgumentParser(
    description='Tree test exe')
parser.add_argument('--trails', type=int, metavar='<number of trees>', default=100,
        help='input trail number [%(default)i]')
parser.add_argument('--ub', type=int, metavar = '<upper bound for number of dataset>', default=100,
        help='size of donor and acceptor [%(default)i]')
arg = parser.parse_args()

for _ in range(arg.trails):
    size = random.randint(0, arg.ub)
    dons, accs, pos2info = randomf.generate_dev_data(size)
    tree = IsoformTree(dons, accs, pos2info)
    splits = len(tree.rules)
    leaves = len(tree.output)

    try:
        assert( splits+1 == leaves)
    except:
        print(splits)
        print(tree.rules)
        print(leaves)
        print(tree.output)