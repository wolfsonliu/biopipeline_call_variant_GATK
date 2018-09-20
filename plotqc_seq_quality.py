#! /usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

####################

parser = argparse.ArgumentParser(
    description='Plot sequence quality distribution'
)

parser.add_argument(
    '--input', nargs='?', help='input text file (tab separated) with columns: [Quality, Count] , NO HEADER'
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
data.columns = ['Quality', 'Count']

####################

width = (np.array(data['Quality'][1:]) - np.array(data['Quality'][:-1])).mean()

fig, axes = plt.subplots()
fig.suptitle('Sequence quality distribution')

axes.bar(
    data['LL'], data['Count'], width, align='edge',
    color='#2166ac', edgecolor='white'
)

axes.set_xlabel('Sequence Quality Score')
axes.set_ylabel('Count')


fig.savefig(args['output'], transparent=True)

################################################################################
