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
    '--input', nargs='?', help='input vcf file'
)
parser.add_argument(
    '--output', nargs='?', help='output pie chart in pdf format'
)

args = vars(parser.parse_args())

####################

ntcolor = {
    'A': {'A':'#a63603', 'C':'#e6550d', 'G':'#fd8d3c', 'T':'#fdae6b'},
    'C': {'C':'#006d2c', 'A':'#31a354', 'G':'#74c476', 'T':'#a1d99b'},
    'G': {'G':'#08519c', 'A':'#3182bd', 'C':'#6baed6', 'T':'#9ecae1'},
    'T': {'T':'#54278f', 'A':'#756bb1', 'C':'#9e9ac8', 'G':'#bcbddc'}
}

####################

data = pd.read_table(
    args['input'],
    header=None, comment='#'
)
data.columns = [
    'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '20'
]
data['reflen'] = data['REF'].map(len)
data['altlen'] = data['ALT'].map(len)
data['issnp'] = np.logical_and(data['reflen'] == 1, data['altlen'] == 1)
data['vtype'] = data[['REF', 'ALT']].apply(lambda x : x['REF'] + '->' + x['ALT'], axis=1)

snptype = pd.DataFrame(
    {'count':data.loc[data['issnp']].groupby(['vtype'])['issnp'].count()}
)
snptype['text'] = snptype.index + snptype['count'].map(lambda x:'\n({0})'.format(x))
snptype['ref'] = snptype.index.str.split('->').map(lambda x: x[0])
snptype['alt'] = snptype.index.str.split('->').map(lambda x: x[1])
snptype['outer_color'] = snptype.apply(lambda x: ntcolor[x['ref']][x['ref']], axis=1)
snptype['inner_color'] = snptype.apply(lambda x: ntcolor[x['ref']][x['alt']], axis=1)

####################

fig, axes = plt.subplots()
fig.suptitle('SNV Type Pie Chart')

axes.pie(
    x=snptype.groupby('ref')['count'].sum(),
    colors=snptype.groupby('ref')['outer_color'].unique().map(lambda x: x[0]),
    radius=1, wedgeprops=dict(width=0.05, edgecolor='w'),
    startangle=90, counterclock =False
)


axes.pie(
    x=snptype['count'], labels=snptype['text'], colors=snptype['inner_color'],
    autopct='%1.1f%%',
    radius=0.95, wedgeprops=dict(width=0.95, edgecolor='w'),
    startangle=90, counterclock =False
)

axes.axis('equal')

fig.savefig(args['output'], transparent=True)

################################################################################
