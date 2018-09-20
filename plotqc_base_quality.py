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
    '--input', nargs='?', help='input text file (tab separated) with columns: [Base, Mean, Median, Lower_Quartile, Upper_Quartile, 10th_Percentile, 90th_Percentile] , NO HEADER'
)
parser.add_argument(
    '--output', nargs='?', help='output pie chart in pdf format'
)

args = vars(parser.parse_args())

####################

data = pd.read_table(
    args['input'],
    header=None, comment='#'
).fillna(0)
data.columns = ['Base', 'mean', 'med', 'q1', 'q3', 'whislo', 'whishi']
data.loc[:,'warning'] = np.logical_or(data['q1'] < 10, data['med'] < 25)
data.loc[:,'failure'] = np.logical_or(data['q1'] < 5, data['med'] < 20)

####################
normalcolor = '#2166ac'
warningcolor = '#f4a582'
failurecolor = '#ca0020'

plotdata = list(
    data[['med', 'q1', 'q3', 'whislo', 'whishi', 'mean']].apply(
        lambda x: dict(x), axis=1
    )
)

boxprops = dict(color=normalcolor)
whiskerprops = dict(color='#000000')
medianprops = dict(color='#000000')


fig, axes = plt.subplots()
fig.suptitle('Base Quality')

axes.hlines(
    [25, 20, 10, 5], xmin=0, xmax=len(plotdata)-1,
    colors=[warningcolor, failurecolor, warningcolor, failurecolor],
    linestyles=['dotted']*4
)

boxes = axes.bxp(
    plotdata, positions=np.arange(len(plotdata)),
    showfliers=False,
    medianprops=medianprops, capprops=medianprops,
    boxprops=boxprops, whiskerprops=whiskerprops
)

box = axes.get_position()

axes.set_position(
    [box.x0, box.y0 + box.height * 0.15,
     box.width, box.height * 0.85]
)


for i in range(len(plotdata)):
    if data['warning'][i]:
        boxes['boxes'][i].set_color(warningcolor)
    elif data['failure'][i]:
        boxes['boxes'][i].set_color(failurecolor)

axes.set_xticklabels(
    data['Base'], {'rotation':'vertical', 'verticalalignment':'top'}
)
axes.set_xlabel('Position')
axes.set_ylabel('Quality Score')

fig.savefig(args['output'], transparent=True)

################################################################################
