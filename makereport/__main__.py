#! /usr/bin/env python3
####################

import os
import pandas as pd
import zipfile as zf
import argparse
from functools import reduce
####################

from repfuncs import write_list_to_file
from repfuncs import get_list_block
from repfuncs import get_list_startwith
from repfuncs import get_list_cut
from plotfuncs import plotqc_base_quality
from plotfuncs import plotqc_base_content
from plotfuncs import plotqc_seq_gc
from plotfuncs import plotqc_seq_length
from plotfuncs import plotqc_seq_quality
from plotfuncs import plotvcf_annovcf
from plotfuncs import plotvcf_piesnp
####################
parser = argparse.ArgumentParser(description='Generate Report TeX file.')
parser.add_argument(
    '--label', nargs='?', default='Exome'
)
parser.add_argument(
    '--report-title', nargs='?', default='Exome Variant Analysis Report', required=True
)
parser.add_argument(
    '--report-author', nargs='?', default='MS Health Care', required=True
)
parser.add_argument(
    '--fig-pipeline', nargs='?', default='', required=True
)
# parser.add_argument(
#     '--fq1', nargs='?', required=True,
#     help='fastq file path, input the fastq file path for single end sequencing' +
#     ' or the first fastq file path for paired end sequencing'
# )
# parser.add_argument(
#     '--fq2', nargs='?',
#     help='fastq file path, not required for single end sequencing' +
#     ' or the second fastq file path for paired end sequencing'
# )
parser.add_argument(
    '--qczip1', nargs='?', required=True, default='',
    help='FastQC zip file path, input the only file path for single end sequencing' +
    ' or the first fastq FastQC zip file path for paired end sequencing'
)
parser.add_argument(
    '--qczip2', nargs='?', default='',
    help='FastQC zip file path, not required for single end sequencing' +
    ' or the second fastq FastQC zip file path for paired end sequencing'
)
parser.add_argument(
    '--samstat', nargs='?', required=True,
    help='samtools stats output txt file path'
)
parser.add_argument(
    '--bcfstat', nargs='?', required=True,
    help='bcftools stats output txt file path'
)
parser.add_argument(
    '--vcf', nargs='?', required=True,
    help='vcf file annotated by ANNOVAR'
)
parser.add_argument(
    '--output-directory', nargs='?', help='output directory path'
)
parser.add_argument('--ref-genome', nargs='?', default='GRCh38')

argdict = vars(parser.parse_args())
argdict['label'] = argdict['label'].replace('_', '\\_')
argdict['report_title'] = argdict['report_title'].replace('_', '\\_')
argdict['report_author'] = argdict['report_author'].replace('_', '\\_')
argdict['ref_genome'] = argdict['ref_genome'].replace('_', '\\_')
argdict['paired'] = False

# paired end
# if argdict['fq1'] not None and argdict['fq2'] not None:
if argdict['qczip1'] != '' and argdict['qczip2'] != '':
    argdict['paired'] = True
# else:
#     raise ValueError('--qczip1 and --qczip2 should all be gaven for paired-end reads')

# make directory
outdirs = dict()
outdirs['base'] = argdict['output_directory']
outdirs['fig'] = os.path.join(outdirs['base'], 'fig')
if not os.path.exists(outdirs['fig']):
    try:
        os.makedirs(outdirs['fig'], mode=511)
    except:
        pass

fignames = dict()
fignames['qc'] = dict()
fignames['qc']['f'] = dict()
fignames['qc']['r'] = dict()
fignames['vcf'] = dict()

fignames['pipeline'] = argdict['fig_pipeline']
####################
# QC file process
# FastQC extract fastqc_data.txt
plotfuncs_dict = {
    'base_quality': plotqc_base_quality,
    'base_content': plotqc_base_content,
    'seq_length': plotqc_seq_length,
    'seq_quality': plotqc_seq_quality,
    'seq_gc': plotqc_seq_gc
}
qc_namedict = {
    'basic_stat': '>>Basic Statistics',
    'base_quality': '>>Per base sequence quality',
    'base_content': '>>Per base sequence content',
    'base_n': '>>Per base N content',
    'seq_quality': '>>Per sequence quality scores',
    'seq_gc': '>>Per sequence GC content',
    'seq_length': '>>Sequence Length Distribution',
    'seq_duplication': '>>Sequence Duplication Levels',
    'seq_overrepresented': '>>Overrepresented sequences',
    'seq_adapter': '>>Adapter Content'
}
qczip = dict()
qczip['f'] = dict()
qczip['r'] = dict()
qclabel = ['f']
if argdict['paired']:
    qclabel += ['r']
for x in qclabel:
    filepath = argdict['qczip1']
    if x == 'r':
        filepath = argdict['qczip2']
    qczip[x]['qczip'] = zf.ZipFile(filepath)
    qczip[x]['name'] = '.'.join(os.path.basename(filepath).split('.')[0:-1])
    qczip[x]['data'] = dict()
    qczip[x]['warn'] = dict()
    qczip[x]['datanames'] = dict()
    qczip[x]['fignames'] = dict()
    # fetch all the QC data
    with qczip[x]['qczip'].open(os.path.join(qczip[x]['name'], 'fastqc_data.txt'), 'r') as f:
        qczip[x]['data']['all'] = [line.decode().strip() for line in f]
    # write QC data to file
    for a, b in qc_namedict.items():
        qczip[x]['data'][a] = get_list_block(
            qczip[x]['data']['all'], b, '>>END_MODULE'
        )
        qczip[x]['datanames'][a] = os.path.join(
            outdirs['base'], qczip[x]['name'] + '_' + a + '.txt'
        )
        write_list_to_file(qczip[x]['data'][a], qczip[x]['datanames'][a])
        qczip[x]['warn'][a] = [
            line.split('\t')[1] for line in qczip[x]['data']['all']
            if line.find(qc_namedict[a]) == 0
        ][0]
    # plot the figures to be used in the report
    for a in ['base_quality', 'base_content', 'seq_length', 'seq_quality', 'seq_gc']:
        qczip[x]['fignames'][a] = os.path.join(
            outdirs['fig'], qczip[x]['name'] + '_' + a + '.pdf'
        )
        fignames['qc'][x][a] = qczip[x]['fignames'][a]
        plotfuncs_dict[a](qczip[x]['datanames'][a], qczip[x]['fignames'][a])


####################
# SAM stats file
samstat_dict = {
    'summary_number': ['SN', ''],
    'GC_first': ['GCF', 'gc\tcount'],
    'GC_last': ['GCL', 'gc\tcount'],
    'ACGT': ['GCC', 'cycle\tA\tC\tG\tT\tN\t0'],
    'ACGT_first': ['FBC', 'cycle\tA\tC\tG\tT\tN\t0'],
    'ACGT_last': ['LBC', 'cycle\tA\tC\tG\tT\tN\t0'],
    'read_length': ['RL', 'read_length\tcount'],
    'read_length_first': ['FRL', 'read_length\tcount'],
    'read_length_last': ['LRL', 'read_length\tcount'],
    'indel': ['ID', 'length\tinsertion\tdeletion'],
    'cycle_indel': ['IC', 'cycle\tinsertion_fwd\tinsertion_rev\tdeletion_fwd\tdelection_rev'],
    'coverage_distribution': ['COV', '']
}
samstat = dict()
samstat['data'] = dict()
samstat['datanames'] = dict()
with open(argdict['samstat'], 'r') as f:
    samstat['data']['all'] = [line.strip() for line in f]
for a, b in samstat_dict.items():
    tmp = get_list_startwith(samstat['data']['all'], b[0])
    if len(tmp) > 0:
        samstat['data'][a] = get_list_cut(tmp, '\t', '1-')
        samstat['datanames'][a] = os.path.join(
            outdirs['base'], 'samstat_' + a + '.txt'
        )
        if len(b[1]) > 0:
            outlist = [b[1]] + samstat['data'][a]
        else:
            outlist = samstat['data'][a]
        write_list_to_file(
            outlist, samstat['datanames'][a]
        )
samstat['summary'] = dict(
    (x.split('\t')[0].replace(':',''), x.split('\t')[1])
    for x in samstat['data']['summary_number']
)


####################
# VCF stats file
vcfstat_dict = {
    'summary_number': ['SN', 'id\tkey\tvalue'],
    'tstv': ['TSTV', 'id\tts\ttv\tts/tv\tts (1st ALT)\ttv (1st ALT)\tts/tv (1st ALT)'],
    'qual': ['QUAL', 'id\tQuality\tnumber of SNPs\tnumber of transitions (1st ALT)\tnumber of transversions (1st ALT)\tnumber of indels'],
    'indel_distribution': ['IDD', 'id\tlength (deletions negative)\tcount'],
    'substitution_types': ['ST', 'id\ttype\tcount'],
    'depth_distribution': ['DP', 'id\tbin\tnumber of genotypes\tfraction of genotypes (%)\tnumber of sites\tfraction of sites (%)']
}
vcfstat = dict()
vcfstat['data'] = dict()
vcfstat['datanames'] = dict()
with open(argdict['bcfstat'], 'r') as f:
    vcfstat['data']['all'] = [line.strip() for line in f]
for a, b in vcfstat_dict.items():
    tmp = get_list_startwith(vcfstat['data']['all'], b[0])
    if len(tmp) > 0:
        vcfstat['data'][a] = get_list_cut(tmp, '\t', '1-')
        vcfstat['datanames'][a] = os.path.join(
            outdirs['base'], 'bcfstat_' + a + '.txt'
        )
        if len(b[1]) > 0:
            outlist = [b[1]] + vcfstat['data'][a]
        else:
            outlist = vcfstat['data'][a]
        write_list_to_file(
            outlist, vcfstat['datanames'][a]
        )
vcfstat['summary'] = dict(
    (x.split('\t')[1].replace(':',''), x.split('\t')[2])
    for x in vcfstat['data']['summary_number']
)


####################
# VCF file process
fignames['vcf']['annotation'] = os.path.join(outdirs['fig'], 'vcf_annovcf.pdf')
plotvcf_annovcf(argdict['vcf'], fignames['vcf']['annotation'])
fignames['vcf']['snvtype'] = os.path.join(outdirs['fig'], 'vcf_snvtype.pdf')
plotvcf_piesnp(argdict['vcf'], fignames['vcf']['snvtype'])

# read vcf
vcf = pd.read_table(argdict['vcf'], header=None, comment='#')
vcf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '20']
vcf.loc[:,'INFOdict'] = vcf['INFO'].str.split(';').map(
    lambda x: dict(a.split('=') for a in x if len(a.split('=')) == 2)
)

clnvcf = vcf.loc[
    vcf['INFOdict'].map(lambda x: x['CLNDN'] != '.').values
].reset_index(drop=True)[['CHROM', 'POS', 'REF', 'ALT', 'INFOdict']].copy()
clnvcf['CLNDN'] = clnvcf['INFOdict'].map(lambda x: x['CLNDN'].replace('_', ' ').split('|'))
clnvcf['CLNSIG'] = clnvcf['INFOdict'].map(lambda x: x['CLNSIG'].replace('_', ' '))


####################
# make texvariable

texvariable = {
    'report_title': argdict['report_title'],
    'report_author': argdict['report_author'],
    'sum_total_seq': samstat['summary']['raw total sequences'],
    'sum_mapping_rate': round(
        float(samstat['summary']['reads mapped'])/
        float(samstat['summary']['raw total sequences']),
        2
    ),
    'sum_refgenome': argdict['ref_genome'],
    'sum_variants_number': vcfstat['summary']['number of records'],
    'sum_snv_number': vcfstat['summary']['number of SNPs'],
    'sum_indel_number': vcfstat['summary']['number of indels'],
    'seq_mode': ['Single-End', 'Paired-End'][argdict['paired']],
    'map_maptool': 'BWA MEM',
    'map_maptool_v': '0.7.17',
    'map_refgenome': 'GRCh38',
    'map_refgenome_v': 'UCSC hg38',
    'map_samtool': 'samtools',
    'map_samtool_v': '1.7',
    'map_mkdup': 'Picard',
    'map_mkdup_v': '2.18.11',
    'map_total_seq': samstat['summary']['raw total sequences'],
    'map_mapped_seq': samstat['summary']['reads mapped'],
    'map_unmapped_seq': samstat['summary']['reads unmapped'],
    'vc_total_variant': vcfstat['summary']['number of records'],
    'vc_snv_number': vcfstat['summary']['number of SNPs'],
    'vc_indel_number': vcfstat['summary']['number of indels'],
    'figpath': outdirs['fig'] + '/',
    'fig_pipeline': argdict['fig_pipeline'],
    'fig_vc_snvtype': os.path.basename(fignames['vcf']['snvtype']),
    'fig_vc_annovcf': os.path.basename(fignames['vcf']['annotation'])
}


for x in qclabel:
    if x == 'f':
        tag = ''
    else:
        tag = '_r'
    seq_length_data = pd.DataFrame(
        [
            (float(a.split('\t')[0]), float(a.split('\t')[1]))
            for a in samstat['data']['read_length']
        ]
    )
    seq_length_data.columns = ['length', 'count']
    seq_length_data['cumcount'] = seq_length_data['count'].cumsum()
    seq_length_fd = list()
    seq_length_data.apply(
        lambda x: seq_length_fd.extend([x['length']] * int(x['count'])), axis=1
    )
    seq_length_fd = pd.Series(seq_length_fd)
    texvariable['qc_seq_length_mean'+tag] = round(seq_length_fd.mean(), 2)
    texvariable['qc_seq_length_min'+tag] = round(seq_length_fd.min(), 2)
    texvariable['qc_seq_length_median'+tag] = round(seq_length_fd.median(), 2)
    texvariable['qc_seq_length_max'+tag] = round(seq_length_fd.max(), 2)
    texvariable['qc_seq_length_q1'+tag] = round(seq_length_fd.quantile(0.25), 2)
    texvariable['qc_seq_length_q3'+tag] = round(seq_length_fd.quantile(0.75), 2)
    seq_gc_data = pd.DataFrame(
        [
            (float(a.split('\t')[0]), float(a.split('\t')[1]))
            for a in qczip['f']['data']['seq_gc'][1:]
        ]
    )
    seq_gc_data.columns = ['gc', 'count']
    seq_gc_data['cumcount'] = seq_gc_data['count'].cumsum()
    seq_gc_fd = list()
    seq_gc_data.apply(lambda x: seq_gc_fd.extend([x['gc']] * int(x['count'])), axis=1)
    seq_gc_fd = pd.Series(seq_gc_fd)
    texvariable['qc_seq_gc_mean'+tag] = round(seq_gc_fd.mean(), 2)
    texvariable['qc_seq_gc_min'+tag] = round(seq_gc_fd.min(), 2)
    texvariable['qc_seq_gc_median'+tag] = round(seq_gc_fd.median(), 2)
    texvariable['qc_seq_gc_max'+tag] = round(seq_gc_fd.max(), 2)
    texvariable['qc_seq_gc_q1'+tag] = round(seq_gc_fd.quantile(0.25), 2)
    texvariable['qc_seq_gc_q3'+tag] = round(seq_gc_fd.quantile(0.75), 2)
    texvariable['fig_qc_base_quality_boxplot'+tag] = os.path.basename(
        fignames['qc'][x]['base_quality']
    )
    texvariable['fig_qc_base_content'+tag] = os.path.basename(
        fignames['qc'][x]['base_content']
    )
    texvariable['fig_qc_seq_length_distribution'+tag] = os.path.basename(
        fignames['qc'][x]['seq_length']
    )
    texvariable['fig_qc_seq_quality_distribution'+tag] = os.path.basename(
        fignames['qc'][x]['seq_quality']
    )
    texvariable['fig_qc_seq_gc'+tag] = os.path.basename(
        fignames['qc'][x]['seq_gc']
    )
    # 'qc2-seq-length-mean',
    # 'qc2-seq-length-min',
    # 'qc2-seq-length-median',
    # 'qc2-seq-length-max',
    # 'qc2-seq-length-q1',
    # 'qc2-seq-length-q3',
    # 'qc2-seq-gc-mean',
    # 'qc2-seq-gc-min',
    # 'qc2-seq-gc-median',
    # 'qc2-seq-gc-max',
    # 'qc2-seq-gc-q1',
    # 'qc2-seq-gc-q3',
    # 'fig-qc2-seq-length-distribution',
    # 'fig-qc2-base-quality-boxplot',
    # 'fig-qc2-seq-quality-distribution',
    # 'fig-qc2-base-content',
    # 'fig-qc2-seq-gc',


####################
# generate tex file
from texfuncs import tex_warning
from texfuncs import tex_clnvcf_table
from texfuncs import tex_tex

texblock = dict()
texblock['qc_warn_message'] = dict()

for x in qclabel:
    texblock['qc_warn_message'][x] = [
        ' '.join(
            ['The quality control result for', qczip[x]['name'], 'is:']
        ).replace('_', '\\_')
    ]+ tex_warning(
        [a for a,b in qczip[x]['warn'].items() if b == 'warn'],
        [a for a,b in qczip[x]['warn'].items() if b == 'fail']
    )

sigdescription = ['']
if clnvcf.shape[0] > 0:
    sigstat = clnvcf['CLNSIG'].value_counts()
    sigdescription = [
        '\\begin{{itemize}}',
        *[
            '\\item {0}: {1} {2}'.format(
            x, y, 'variants' if y > 1 else 'variant'
            ) for x,y in dict(sigstat).items()
        ],
        '\\end{{itemize}}'
    ]

reportvcf = clnvcf.copy()
if reportvcf.shape[0] > 20:
    reportvcf = reportvcf.loc[reportvcf['CLNSIG'].str.find('provided') == -1]
    if reportvcf.shape[0] > 20:
        reportvcf = reportvcf.loc[reportvcf['CLNSIG'].str.find('enign') == -1]

# calculate diseases
disvcf = reportvcf.loc[reportvcf['CLNSIG'].str.find('athogenic') > -1]
disvcf = reportvcf    # used for test
if disvcf.shape[0] > 0:
    disease = pd.Series(
        reduce(
            lambda x,y: x+y,
            [
                x for x in list(disvcf['CLNDN'])
            ]
        )
    ).value_counts()
    disease = disease.loc[disease.index != 'not provided']
    disdescription = [
        '\\begin{{itemize}}',
        *[ '\\item {0}: supported by {1} {2}'.format(
            x, y, 'variants' if y > 1 else 'variant'
        ) for x,y in dict(disease).items()],
        '\\end{{itemize}}'
    ]

if reportvcf.shape[0] > 0:
    texblock['clnvcf_table'] = [
        'The Clinical significance of variants summarized as follows:',
        '',
        *sigdescription,
        ''
    ]
    if disvcf.shape[0] > 0:
        texblock['clnvcf_table'] = texblock['clnvcf_table'] + [
            'According to the analysis of the sample provided, ',
            'the following diseases need consideration:',
            '',
            *disdescription,
            ''
        ]
    else:
        texblock['clnvcf_table'] = texblock['clnvcf_table'] + [
            'There is no likely pathogenic or pathogenic variants found in the sample.',
            ''
        ]
    texblock['clnvcf_table'] = texblock['clnvcf_table'] + [
        'The variants annotated with connection to diseases, ',
        'which need consideration, listed as follows (\\ref{{tab:vcfclinvar}}):'
    ] + tex_clnvcf_table(reportvcf)
else:
    if clnvcf.shape[0] > 0:
        texblock['clnvcf_table'] = [
            'The Clinical significance of variants listed as follows:',
            *sigdescription,
            '',
            'No known disease-associated variants found in the sample.'
        ]
    else:
        texblock['clnvcf_table'] = []

for x,y in tex_tex(argdict['paired']).items():
    texblock[x] = y

tex = texblock['dochead'] +\
      texblock['sec-summary'] +\
      texblock['sec-analysis_pipeline'] +\
      texblock['sec-analysis_result-head'] +\
      texblock['sec-analysis_result-qc-head'] +\
      [x for i in qclabel for x in texblock['qc_warn_message'][i]] +\
      texblock['sec-analysis_result-qc-body'] +\
      texblock['sec-analysis_result-mapping'] +\
      texblock['sec-analysis_result-cv'] +\
      texblock['sec-analysis_result-annotation'] + ['\n'] +\
      texblock['clnvcf_table'] +\
      texblock['docend']

outtex = '\n'.join(tex).format(**texvariable)

with open(os.path.join(outdirs['base'], 'report.tex'), 'w') as f:
    f.write(outtex)

################################################################################
