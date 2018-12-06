#! /usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

####################

parser = argparse.ArgumentParser(
    description='Plot pie chart of SNP types in VCF file'
)

parser.add_argument(
    '--input', nargs='?', help='input text file (tab separated) with columns: [Base, G, A, T, C] , NO HEADER'
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
data.columns = ['Base', 'G', 'A', 'T', 'C']
data.loc[:,'warning'] = np.logical_or(
    (data['A'] - data['T']).abs() > 10,
    (data['G'] - data['C']).abs() > 10
)

data.loc[:,'failure'] = np.logical_or(
    (data['A'] - data['T']).abs() > 20,
    (data['G'] - data['C']).abs() > 20
)


####################
ntcolor = {'A':'#a63603', 'C':'#006d2c', 'G':'#08519c', 'T':'#54278f'}

normalcolor = '#2166ac'
warningcolor = '#f4a582'
failurecolor = '#ca0020'


fig, axes = plt.subplots()
fig.suptitle('Content across all bases')


linents = dict()

for nt in ntcolor.keys():
    linents[nt] = axes.plot(
        np.arange(data.shape[0]), data[nt], color=ntcolor[nt]
    )


ybottom, ytop = axes.get_ylim()
# plot warning and failure bar
if data['warning'].sum() > 0:
    colors = np.array(
        [warningcolor, failurecolor]
    )[list(data.loc[data['warning'], 'failure'].astype(int))]
    axes.bar(
        np.arange(data.shape[0])[data['warning']],
        [ytop] * data['warning'].sum(),
        color=colors, edgecolor=None, alpha=.5
    )

# reserve place for x tick labels
box = axes.get_position()
axes.set_position(
    [box.x0, box.y0 + box.height * 0.15,
     box.width, box.height * 0.85]
)

axes.set_xticks(np.arange(data.shape[0]))
axes.set_xticklabels(
    data['Base'], {'rotation':'vertical', 'verticalalignment':'top'}
)
axes.set_xlabel('Position in reads')
axes.set_ylabel('Percent (%)')
axes.legend(
    (linents['A'][0], linents['T'][0], linents['G'][0], linents['C'][0]),
    ('A', 'T', 'G', 'C'),
    loc='upper center', ncol=4
)


fig.savefig(args['output'], transparent=True)

################################################################################
