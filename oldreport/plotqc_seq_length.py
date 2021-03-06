#! /usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

####################

parser = argparse.ArgumentParser(
    description='Plot sequence length distribution'
)

parser.add_argument(
    '--input', nargs='?', help='input text file (tab separated) with columns: [Length, Count] as columns, no header'
)
parser.add_argument(
    '--output', nargs='?', help='output chart in pdf format'
)

args = vars(parser.parse_args())

####################

data = pd.read_table(
    args['input'],
    header=None, comment='#'
)
data.columns = ['Length', 'Count']

data.loc[:,'LL'] = data['Length'].astype('str')\
                                 .str\
                                 .split('-')\
                                 .map(lambda x: x[0])\
                                 .astype(int)
if data.shape[0] != 1:
    width = (data['LL'][1:].values - data['LL'][:-1].values).mean()
else:
    width = data['LL'].mean() / 5
    data.loc[1] = [0, 0, 0]

####################

fig, axes = plt.subplots()
fig.suptitle('Sequence Length Distribution')

axes.bar(
    data['LL'], data['Count'], width, align='edge',
    color='#2166ac', edgecolor='white'
)

axes.set_xlabel('Sequence Length (bp)')
axes.set_ylabel('Count')

fig.savefig(args['output'], transparent=True)

################################################################################
