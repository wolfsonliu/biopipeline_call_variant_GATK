#! /usr/bin/env python3
import numpy as np
import pandas as pd
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
####################

# plotqc_base_content
def plotqc_base_content(datapath, figpath):
    data = pd.read_table(datapath, header=None, comment='#').fillna(0)
    data.columns = ['Base', 'G', 'A', 'T', 'C']
    data.loc[:,'warning'] = np.logical_or(
        (data['A'] - data['T']).abs() > 10, (data['G'] - data['C']).abs() > 10
    )

    data.loc[:,'failure'] = np.logical_or(
        (data['A'] - data['T']).abs() > 20, (data['G'] - data['C']).abs() > 20
    )
    # set plot color
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
        [box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85]
    )
    axes.set_xticks(np.arange(data.shape[0]))
    axes.set_xticklabels(
        data['Base'], {'rotation':'vertical', 'verticalalignment':'top'}
    )
    axes.set_xlabel('Position in reads')
    axes.set_ylabel('Percent (%)')
    axes.legend(
        (linents['A'][0], linents['T'][0], linents['G'][0], linents['C'][0]),
        ('A', 'T', 'G', 'C'), loc='upper center', ncol=4
    )
    fig.savefig(figpath, transparent=True)
####################

# plotqc_base_quality
def plotqc_base_quality(datapath, figpath):
    data = pd.read_table(datapath, header=None, comment='#').fillna(0)
    data.columns = ['Base', 'mean', 'med', 'q1', 'q3', 'whislo', 'whishi']
    data.loc[:,'warning'] = np.logical_or(data['q1'] < 10, data['med'] < 25)
    data.loc[:,'failure'] = np.logical_or(data['q1'] < 5, data['med'] < 20)
    # set colors
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
        plotdata, positions=np.arange(len(plotdata)), showfliers=False,
        medianprops=medianprops, capprops=medianprops,
        boxprops=boxprops, whiskerprops=whiskerprops
    )
    # change box color
    box = axes.get_position()
    axes.set_position(
        [box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85]
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
    fig.savefig(figpath, transparent=True)
####################

# plotqc_seq_gc
def plotqc_seq_gc(datapath, figpath):
    data = pd.read_table(datapath, header=None, comment='#').fillna(0)
    data.columns = ['GC', 'Count']
    meangc = (data['Count'] * data['GC']).sum() / (data['Count'].sum())
    # set colors
    ntcolor = {'A':'#a63603', 'C':'#006d2c', 'G':'#08519c', 'T':'#54278f'}
    fig, axes = plt.subplots()
    fig.suptitle('Sequence mean GC content distribution')
    axes.plot(data['GC'], data['Count'], color='#0571b0')
    axes.set_xticks(data['GC'][np.arange(data.shape[0]) % 10 == 0])
    axes.set_xlabel('Mean GC Content (%)')
    axes.set_ylabel('Count')
    fig.savefig(figpath, transparent=True)
####################

# plotqc_seq_length
def plotqc_seq_length(datapath, figpath):
    data = pd.read_table(datapath, header=None, comment='#')
    data.columns = ['Length', 'Count']
    data.loc[:,'LL'] = data['Length'].astype(
        'str'
    ).str.split('-').map(lambda x: x[0]).astype(int)
    # adjust data with same length reads
    if data.shape[0] != 1:
        width = (data['LL'][1:].values - data['LL'][:-1].values).mean()
    else:
        width = data['LL'].mean() / 5
        data.loc[1] = [0, 0, 0]
    fig, axes = plt.subplots()
    fig.suptitle('Sequence Length Distribution')
    axes.bar(
        data['LL'], data['Count'], width, align='edge',
        color='#2166ac', edgecolor='white'
    )
    axes.set_xlabel('Sequence Length (bp)')
    axes.set_ylabel('Count')
    fig.savefig(figpath, transparent=True)
####################

# plotqc_seq_quality
def plotqc_seq_quality(datapath, figpath):
    data = pd.read_table(datapath, header=None, comment='#').fillna(0)
    data.columns = ['Quality', 'Count']
    if data.shape[0] != 1:
        width = (np.array(data['Quality'][1:]) - np.array(data['Quality'][:-1])).mean()
    else:
        width = data['Quality'].mean()
    fig, axes = plt.subplots()
    fig.suptitle('Sequence quality distribution')
    axes.bar(
        data['Quality'], data['Count'], width, align='edge',
        color='#2166ac', edgecolor='white'
    )
    axes.set_xlabel('Sequence Quality Score')
    axes.set_ylabel('Count')
    fig.savefig(figpath, transparent=True)
####################

# plotvcf_annovcf
def plotvcf_annovcf(datapath, figpath):
    data = pd.read_table(datapath, header=None, sep='\t', comment='#')
    data.columns = [
        'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '20'
    ]
    data.loc[:,'INFOdict'] = data['INFO'].str.split(';').map(
        lambda x: dict(a.split('=') for a in x if len(a.split('=')) == 2)
    )
    data.loc[:,'reflen'] = data['REF'].map(len)
    data.loc[:,'altlen'] = data['ALT'].map(len)
    data.loc[:,'issnp'] = np.logical_and(data['reflen'] == 1, data['altlen'] == 1)
    data.loc[:,'vtype'] = data[['REF', 'ALT']].apply(
        lambda x : x['REF'] + '->' + x['ALT'], axis=1
    )
    annoinfo = data[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'issnp']].copy()
    annoinfo.loc[:,'refGene'] = data['INFOdict'].map(
        lambda x: x['Gene.refGene']
    ).copy() != '.'
    annoinfo.loc[:,'ExAC'] = data['INFOdict'].map(lambda x: x['ExAC_ALL']).copy() != '.'
    annoinfo.loc[:,'ClinVar'] = data['INFOdict'].map(lambda x: x['CLNDN']).copy() != '.'
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
        ('SNV in Database', 'INDEL in Database',
         'SNV not in Database', 'INDEL not in Database'),
        loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2
    )
    fig.savefig(figpath, transparent=True)
####################

# plotvcf_piesnp
def plotvcf_piesnp(datapath, figpath):
    ntcolor = {
        'A': {'A':'#a63603', 'C':'#e6550d', 'G':'#fd8d3c', 'T':'#fdae6b'},
        'C': {'C':'#006d2c', 'A':'#31a354', 'G':'#74c476', 'T':'#a1d99b'},
        'G': {'G':'#08519c', 'A':'#3182bd', 'C':'#6baed6', 'T':'#9ecae1'},
        'T': {'T':'#54278f', 'A':'#756bb1', 'C':'#9e9ac8', 'G':'#bcbddc'}
    }
    data = pd.read_table(datapath, header=None, comment='#')
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
    fig.savefig(figpath, transparent=True)
####################
