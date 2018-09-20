#! /usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

####################

parser = argparse.ArgumentParser(
    description='Plot bar chart of annotated vcf'
)

parser.add_argument(
    '--input', nargs='?', help='input annotation vcf file'
)
parser.add_argument(
    '--output', nargs='?', help='output bar chart in pdf format'
)

args = vars(parser.parse_args())

####################

# process vcf files
data = pd.read_table(
    args['input'],
    header=None,
    sep='\t',
    comment='#'
)

data.columns = [
    'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '20'
]
data.loc[:,'INFOdict'] = data['INFO'].str.split(';').map(
    lambda x: dict(a.split('=') for a in x if len(a.split('=')) == 2)
)
data.loc[:,'reflen'] = data['REF'].map(len)
data.loc[:,'altlen'] = data['ALT'].map(len)
data.loc[:,'issnp'] = np.logical_and(data['reflen'] == 1, data['altlen'] == 1)
data.loc[:,'vtype'] = data[['REF', 'ALT']].apply(lambda x : x['REF'] + '->' + x['ALT'], axis=1)

annoinfo = data[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'issnp']].copy()
annoinfo.loc[:,'refGene'] = data['INFOdict'].map(lambda x: x['Gene.refGene']).copy() != '.'
annoinfo.loc[:,'ExAC'] = data['INFOdict'].map(lambda x: x['ExAC_ALL']).copy() != '.'
annoinfo.loc[:,'ClinVar'] = data['INFOdict'].map(lambda x: x['CLNDBN']).copy() != '.'

annostat = annoinfo.groupby(
    ['issnp', 'refGene', 'ExAC', 'ClinVar']
)['CHROM'].count().unstack(level=0).fillna(0).reset_index()

# add 0 data
if not True in annostat.columns:
    annostat[True] = [0]*annostat.shape[0]
elif not False in annostat.columns:
    annostat[False] = [0]*annostat.shape[0]

for x in ['refGene', 'ExAC', 'ClinVar']:
    if not True in annostat[x].unique():
        annoadd = annostat.copy()
        annoadd.loc[:, x] = [True] * annoadd.shape[0]
        annoadd.loc[:, True] = [0]*annoadd.shape[0]
        annoadd.loc[:, False] = [0]*annoadd.shape[0]
        annostat = pd.concat([annostat, annoadd], axis=0, ignore_index=True)
    elif not False in annostat[x].unique():
        annoadd = annostat.copy()
        annoadd.loc[:, x] = [False] * annoadd.shape[0]
        annoadd.loc[:, True] = [0]*annoadd.shape[0]
        annoadd.loc[:, False] = [0]*annoadd.shape[0]
        annostat = pd.concat([annostat, annoadd], axis=0, ignore_index=True)

####################

# generate data used for plot
plotdata = {
    'SNVindb': (
        annostat.groupby(['refGene'])[True].sum()[True],
        annostat.groupby(['ExAC'])[True].sum()[True],
        annostat.groupby(['ClinVar'])[True].sum()[True]
    ),
    'Indelindb': (
        annostat.groupby(['refGene'])[False].sum()[True],
        annostat.groupby(['ExAC'])[False].sum()[True],
        annostat.groupby(['ClinVar'])[False].sum()[True]
    ),
    'SNVnotdb': (
        annostat.groupby(['refGene'])[True].sum()[False],
        annostat.groupby(['ExAC'])[True].sum()[False],
        annostat.groupby(['ClinVar'])[True].sum()[False]
    ),
    'Indelnotdb': (
        annostat.groupby(['refGene'])[False].sum()[False],
        annostat.groupby(['ExAC'])[False].sum()[False],
        annostat.groupby(['ClinVar'])[False].sum()[False]
    )
}


fig, axes = plt.subplots()
fig.suptitle('Database Annotation')

idx = np.arange(3)
width=0.35

bar1 = axes.bar(idx, plotdata['SNVindb'], width, color='#b2182b')
bar2 = axes.bar(
    idx, plotdata['Indelindb'], width, color='#2166ac',
    bottom=np.array(plotdata['SNVindb'])
)
bar3 = axes.bar(
    idx, plotdata['SNVnotdb'], width, color='#fddbc7',
    bottom=np.array(plotdata['SNVindb']) + np.array(plotdata['Indelindb'])
)
bar4 = axes.bar(
    idx, plotdata['Indelnotdb'], width,  color='#d1e5f0',
    bottom=np.array(plotdata['SNVindb']) +
    np.array(plotdata['Indelindb']) +
    np.array(plotdata['SNVnotdb'])
)

axes.set_xticks(idx, minor=False)
axes.set_xticklabels(('refGene', 'ExAC', 'ClinVar'))
axes.set_xlabel('Database')
axes.set_ylabel('Count')
# Shrink current axis's height by 10% on the bottom
box = axes.get_position()
axes.set_position(
    [box.x0, box.y0 + box.height * 0.15,
     box.width, box.height * 0.85]
)
axes.legend(
    (bar1, bar2, bar3, bar4),
    ('SNV in Database', 'INDEL in Database', 'SNV not in Database', 'INDEL not in Database'),
    loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2
)
fig.savefig(args['output'], transparent=True)

################################################################################
