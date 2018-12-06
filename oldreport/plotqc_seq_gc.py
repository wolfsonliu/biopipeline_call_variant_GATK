#! /usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

####################

parser = argparse.ArgumentParser(
    description='Plot sequence mean GC content distribution'
)

parser.add_argument(
    '--input', nargs='?', help='input text file (tab separated) with columns: [GC Content, Count] , NO HEADER'
)
parser.add_argument(
    '--output', nargs='?', help='output chart in pdf format'
)

args = vars(parser.parse_args())

####################

data = pd.read_table(
    args['input'],
    header=None, comment='#'
).fillna(0)
data.columns = ['GC', 'Count']
meangc = (data['Count'] * data['GC']).sum() / (data['Count'].sum())


####################
ntcolor = {'A':'#a63603', 'C':'#006d2c', 'G':'#08519c', 'T':'#54278f'}


fig, axes = plt.subplots()
fig.suptitle('Sequence mean GC content distribution')

axes.plot(
    data['GC'], data['Count'], color='#0571b0'
)

axes.set_xticks(data['GC'][np.arange(data.shape[0]) % 10 == 0])
axes.set_xlabel('Mean GC Content (%)')
axes.set_ylabel('Count')
 

fig.savefig(args['output'], transparent=True)

################################################################################
