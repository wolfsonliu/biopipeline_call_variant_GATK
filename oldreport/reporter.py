#! /usr/bin/env python3

# tex = "\n".join([
#     "\\documentclass[10pt,oneside,a4paper]{{article}}",
#     "\\usepackage{{textcomp}}",
#     "\\usepackage{{graphicx}}",
#     "\\graphicspath{{{{{figpath}}}}}",
#     "\\usepackage[square,numbers]{{natbib}}",
#     "\\usepackage[colorlinks=false,bookmarks=true,pdfpagemode=FullScreen]{{hyperref}}",
#     "\\title{{{report_title}}}",
#     "\\author{{{report_author}}}",
#     "\\date{{\\today}}",
#     "\\begin{{document}}",
#     "\\maketitle",
#     "\\section{{Summary}}",
#     "\\label{{sec:summary}}",
#     "\\begin{{itemize}}",
#     "\\item Input FASTQ File with {sum_total_seq} reads in total.",
#     "\\item ${sum_mapping_rate}$ of reads mapped to reference genome ({sum_refgenome}).",
#     "\\item {sum_variants_number} variants called from data. \\\\",
#     "\\begin{{itemize}}",
#     "  \\item {sum_snp_number} SNVs (Single Nucleotide Variant)",
#     "  \\item {sum_indel_number} Indels (Insertion and Deletion)",
#     "  \\end{{itemize}}",
#     "\\end{{itemize}}",
#     "\\section{{Analysis Pipeline}}",
#     "\\label{{sec:pipeline}}",
#     "\\begin{{figure}}[h]",
#     "  \\centering",
#     "  \\includegraphics[width=0.75\\linewidth]{{{fig_pipeline}}}",
#     "  \\caption{{MS Variant Calling Pipeline}}",
#     "  \\label{{fig:pipeline}}",
#     "\\end{{figure}}",
#     "\\section{{Analysis Result}}",
#     "\\label{{sec:result}}",
#     "\\subsection{{Quality Control}}",
#     "\\label{{subsec:qc}}",
#     "FASTQ file statistics:",
#     "\\begin{{center}}",
#     "  \\begin{{tabular}}[h]{{l c}}",
#     "    \\hline",
#     "    Feature                & Statistics             \\\\",
#     "    \\hline",
#     "    Sequencing Mode        & {qc_seq_mode}          \\\\",
#     "    Reads number           & {qc_total_seq}         \\\\",
#     "    Pairs number           & {qc_total_pair}        \\\\",
#     "    Read Length Mean       & {qc_seq_length_mean}   \\\\",
#     "    Read Length Min        & {qc_seq_length_min}    \\\\",
#     "    Read Length Median     & {qc_seq_length_median} \\\\",
#     "    Read Length Max        & {qc_seq_length_max}    \\\\",
#     "    Read GC content Mean   & {qc_seq_gc_mean}       \\\\",
#     "    Read GC content Min    & {qc_seq_gc_min}        \\\\",
#     "    Read GC content Median & {qc_seq_gc_median}     \\\\",
#     "    Read GC content Max    & {qc_seq_gc_max}        \\\\",
#     "    \\hline",
#     "  \\end{{tabular}}",
#     "\\end{{center}}",
#     "\\subsubsection{{Sequence Length Distribution}}",
#     "\\label{{subsubsec:sequence_length}}",
#     "Read length distribution is a good measure of the sequence quality. In",
#     "most cases, the sequence length in the FASTQ file obtained from",
#     "sequencing service companies should be mostly similar after their",
#     "quality control. And the bad sequencing quality may lead to the",
#     "cutting of sequencing by the sequencing machine. So, if the sequences",
#     "mostly have same length or normally distributed, that means the",
#     "quality of sequences might be acceptable. However, if the sequence",
#     "length does not follow a normal distribution, the FASTQ file will need",
#     "more attention. And in the quality trim stage, the short sequences and",
#     "bad quality sequences should be discarded. The sequcne length",
#     "distribution is shown in Figure \\ref{{fig:seq_length_distribution}}.",
#     "\\begin{{center}}",
#     "  \\begin{{tabular}}[h]{{l c}}",
#     "    \\hline",
#     "    Read Length        & Nucleotides              \\\\",
#     "    \\hline",
#     "    Average Length     & {qc_seq_length_mean}        \\\\",
#     "    Min Length         & {qc_seq_length_min}    \\\\",
#     "    1st Quatile Length & {qc_seq_length_q1}     \\\\",
#     "    Median Length      & {qc_seq_length_median} \\\\",
#     "    3rd Quatile Length & {qc_seq_length_q3}     \\\\",
#     "    Max Length         & {qc_seq_length_max}    \\\\",
#     "    \\hline",
#     "  \\end{{tabular}}",
#     "\\end{{center}}",
#     "\\begin{{figure}}[h]",
#     "  \\centering",
#     "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_seq_length_distribution}}}}}",
#     "  \\caption{{Sequence Length Distribution}}",
#     "  \\label{{fig:seq_length_distribution}}",
#     "\\end{{figure}}",
#     "\\subsubsection{{Base Quality Boxplot}}",
#     "\\label{{subsubsec:base_quality_boxplot}}",
#     "Base quality boxplot is used to show the bases' quality",
#     "distribution. For the first several bases and the last several bases",
#     "of the reads on Illumina platform, the base quality is generaly of",
#     "large error rate, which should be trimmed before analysis. The base",
#     "quality is shown in Figure \\ref{{fig:base_quality_boxplot}}.",
#     "\\begin{{figure}}[h]",
#     "  \\centering",
#     "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_base_quality_boxplot}}}}}",
#     "  \\caption{{Base Quality Boxplot}}",
#     "  \\label{{fig:base_quality_boxplot}}",
#     "\\end{{figure}}",
#     "Phred score is used to evaluate the base quality, the larger the",
#     "better. when the quality score is higher than 30, that means the error",
#     "rate is lower than 0.001, which can be used as the empirical cutoff",
#     "for the base quality. For the good quality FASTQ file, most of the",
#     "bases have quality larger than 30. And for the bad quality FASTQ file,",
#     "the first quartiles may be below 30 for many sites. So, for the",
#     "following trimming stage, the bad FASTQ file should be trimmed and",
#     "filtered.",
#     "\\subsubsection{{Mean Sequence Quality Distribution}}",
#     "\\label{{subsubsec:mean_quality}}",
#     "The mean quality shows the quality of each sequence. The distribution",
#     "of mean sequence quality shows how the quality of the input FASTQ",
#     "file. For the good quality FASTQ file, most of the mean quality scores",
#     "are high and distributing in one small range. Meanwhile, for the bad",
#     "quality FASTQ files, except for the bad quality score, the range of",
#     "the quality score is also wider. And sometimes there will be multiple",
#     "peaks for the bad quality score distribution. Because there are",
#     "sequences having universally low quality values, which could be a",
#     "result of being poorly imaged in the sequencing stage. In the",
#     "following filter stage, the cutoff of bad quality sequence can be",
#     "chosen from the mean quality score distribution. The sequence quality",
#     "is shown in Figure \\ref{{fig:seq_quality_distribution}}.",
#     "\\begin{{figure}}[h]",
#     "  \\centering",
#     "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_seq_quality_distribution}}}}}",
#     "  \\caption{{Mean Sequence Quality Distributon}}",
#     "  \\label{{fig:seq_quality_distribution}}",
#     "\\end{{figure}}",
#     "\\subsubsection{{Base Content}}",
#     "\\label{{subsubsec:basecontent}}",
#     "The base GC content plot can show whether the AT or GC has different",
#     "Though GC content is species specific, the basewise G percentage of",
#     "per base should be consistent with the basewise C percentage, so as",
#     "with the A and T percentage basewise. Because in the randomly",
#     "generated library, the A and T bases are expected to be the same, so",
#     "as with the G and C bases. The discordance of base content for AT and",
#     "GC infers that the library or sequencing stage might have faults. The",
#     "base content is shown in Figure \\ref{{fig:base_content}}.",
#     "\\begin{{figure}}[h]",
#     "  \\centering",
#     "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_base_content}}}}}",
#     "  \\caption{{Base Content}}",
#     "  \\label{{fig:base_content}}",
#     "\\end{{figure}}",
#     "\\subsubsection{{Sequence GC Content Distribution}}",
#     "\\label{{subsubsec:sequencegc}}",
#     "The distribution of sequence GC content shows potential can provide",
#     "the sequence quality information as well. GC content is not a random",
#     "number in genomes of different species, which is driven by evolution",
#     "in fact and varies continuously.  For a random sequence library, the",
#     "GC content among difference bases should be similar. In ideal case,",
#     "the GC content of sequencing data should reflect the genome GC",
#     "content, but they should not show huge imbalance as well. However, the",
#     "using of random primer in sequencing preparation might introduce the",
#     "GC content bias at different bases.  The mean GC content of each",
#     "sequence can generate a density plot. Ideally, the distribution of",
#     "sequence GC content should be normally distributed. If the",
#     "distribution diverges from normal distribution a lot, that might be a",
#     "sign of low-quality or contamination of sequences from other",
#     "organisms. The sequence GC content distribution is shown in",
#     "Figure \\ref{{fig:seq_gc}}.",
#     "\\begin{{center}}",
#     "  \\begin{{tabular}}[h]{{l c}}",
#     "    \\hline",
#     "    Read GC content        & Nucleotides        \\\\",
#     "    \\hline",
#     "    Average GC content     & {qc_seq_gc_mean}   \\\\",
#     "    Min GC content         & {qc_seq_gc_min}    \\\\",
#     "    1st Quatile GC content & {qc_seq_gc_q1}     \\\\",
#     "    Median GC content      & {qc_seq_gc_median} \\\\",
#     "    3rd Quatile GC content & {qc_seq_gc_q3}     \\\\",
#     "    Max GC content         & {qc_seq_gc_max}    \\\\",
#     "    \\hline",
#     "  \\end{{tabular}}",
#     "\\end{{center}}",
#     "\\begin{{figure}}[h]",
#     "  \\centering",
#     "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_seq_gc}}}}}",
#     "  \\caption{{Sequence GC Content}}",
#     "  \\label{{fig:seq_gc}}",
#     "\\end{{figure}}",
#     "\\subsection{{Mapping}}",
#     "\\label{{subsec:mapping}}",
#     "Mapping stage finds",
#     "the reads position on the reference genome. The",
#     "result of mapping shows the coverage and depth of sequencing. The",
#     "following analysis depends on the accurate mapping result. The mapping",
#     "result is saved as SAM/BAM file. The possible PCR duplicated sequences",
#     "will be marked in the process. BWA is one of the opensource software",
#     "packages for mapping high-throughput sequencing data to reference",
#     "genome that generates accurate results in a fast speed.",
#     "\\subsubsection{{Mapping Process}}",
#     "\\label{{subsubsec:mappingprocess}}",
#     "\\begin{{center}}",
#     "  \\begin{{tabular}}[h]{{l c c}}",
#     "    \\hline",
#     "    Process              & Software        & Version             \\\\",
#     "    \\hline",
#     "    Mapping Software     & {map_maptool}   & {map_maptool_v}     \\\\",
#     "    Reference Genome     & {map_refgenome} & ({map_refgenome_v}) \\\\",
#     "    SAM/BAM Process      & {map_samtool}   & {map_samtool_v}     \\\\",
#     "    Mark PCR Duplicate   & {map_mkdup}     & {map_mkdup_v}       \\\\",
#     "    \\hline",
#     "  \\end{{tabular}}",
#     "\\end{{center}}",
#     "\\subsubsection{{Mapping Statistics}}",
#     "\\label{{subsubsec:mappingstat}}",
#     "\\begin{{center}}",
#     "  \\begin{{tabular}}[h]{{l r}}",
#     "    \\hline",
#     "    Item            & Count              \\\\",
#     "    \\hline",
#     "    Raw Total Reads & {map_total_seq}    \\\\",
#     "    Mapped Reads    & {map_mapped_seq}   \\\\",
#     "    Unmapped Reads  & {map_unmapped_seq} \\\\",
#     "    \\hline",
#     "  \\end{{tabular}}",
#     "\\end{{center}}",
#     "\\subsection{{Variant Calling}}",
#     "\\label{{subsec:variant}}",
#     "GATK 4 is used as the main variant calling tools in the analysis",
#     "pipeline. SNVs (Single Nucleotide Variant) and InDels (Insertion and",
#     "Deletion) are detected by the GATK software. In this pipeline, SNVs",
#     "are considered.",
#     "\\begin{{description}}",
#     "\\item[Variants Calling] GATK 4.0.8 % Setting: if different version",
#     "\\item[Variants Filter] SNP filters are applied. Good quality SNVs should satisfy the following conditions. \\\\",
#     "  \\begin{{itemize}}",
#     "  \\item Quality Depth (QD) $\\geq 2$",
#     "  \\item Mapping Quality (MQ) $\\geq 40$",
#     "  \\item MQRankSum $\\geq -12.5$",
#     "  \\item Fisher test of strand bias (FS) $\\leq 60$",
#     "  \\item ReadPosRankSum $\\geq -8$",
#     "  \\end{{itemize}}",
#     "\\end{{description}}",
#     "Result of variants calling is saved as VCF (variants calling format).",
#     "\\begin{{center}}",
#     "  \\begin{{tabular}}[h]{{l r}}",
#     "    \\hline",
#     "    Item           & Count              \\\\",
#     "    \\hline",
#     "    Total Variants & {vc_total_variant} \\\\",
#     "    SNV            & {vc_snv_number}    \\\\",
#     "    Indel          & {vc_indel_number}  \\\\",
#     "    \\hline",
#     "  \\end{{tabular}}",
#     "\\end{{center}}",
#     "The condition of SNV types is shown in Figure \\ref{{fig:snvtype}}.",
#     "\\begin{{figure}}[h]",
#     "  \\centering",
#     "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_vc_snvtype}}}}}",
#     "  \\caption{{SNV types}}",
#     "  \\label{{fig:snvtype}}",
#     "\\end{{figure}}",
#     "\\subsection{{Variant Annotation}}",
#     "\\label{{subsec:annotationg}}",
#     "Variant annotation links the variants with the existing variant",
#     "database information. The variants annotation uses ANNOVAR software.",
#     "\\begin{{description}}",
#     "\\item[RefSeq Gene] (refGene) RefSeq Gene is the NCBI database for gene",
#     "  information. Annotation with refGene shows the gene that variants",
#     "  belonging to.",
#     "\\item[ExAC] (exac) Exome Aggregation Consortium database, more than",
#     "  60K individuals SNV data.",
#     "\\item[dbSNP] (avsnp) NCBI Database of single nucleotide polymorphisms",
#     "  (SNPs) and other small variants.",
#     "\\item[dbNSFP] (dbnsfp) dbNSFP is a database of all non-synonymous",
#     "  single-nucleotide variants.",
#     "\\item[ClinVar] (clinvar\\_20160302) ClinVar is the NCBI database storeing",
#     "  variants associated with disease.",
#     "\\end{{description}}",
#     "\\subsubsection{{Annotation Statistics}}",
#     "\\label{{subsubsec:annotationstatistics}}",
#     "The annotation of each database is shown in Figure \\ref{{fig:snvannobar}}.",
#     "\\begin{{figure}}[h]",
#     "  \\centering",
#     "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_vc_annovcf}}}}}",
#     "  \\caption{{Variant Annotation}}",
#     "  \\label{{fig:snvannobar}}",
#     "\\end{{figure}}",
#     "\\end{{document}}"
# ])


doclist = dict()
doclist['dochead'] = [
    "\\documentclass[10pt,oneside,a4paper]{{article}}",
    "\\usepackage{{textcomp}}",
    "\\usepackage[table]{{xcolor}}",
    "\definecolor{{rowgray}}{{rgb}}{{0.9,0.9,0.9}}",
    "\\usepackage{{graphicx}}",
    "\\graphicspath{{{{{figpath}}}}}",
    "\\usepackage[square,numbers]{{natbib}}",
    "\\usepackage[colorlinks=false,bookmarks=true,pdfpagemode=FullScreen]{{hyperref}}",
    "\\title{{{report_title}}}",
    "\\author{{{report_author}}}",
    "\\date{{\\today}}",
    "\\begin{{document}}",
    "\\maketitle"
]
doclist['docend'] = ["\\end{{document}}"]
doclist['sec-summary'] = [
    "\\section{{Summary}}",
    "\\label{{sec:summary}}",
    "\\begin{{itemize}}",
    "\\item Input FASTQ File ({seq_mode}) with {sum_total_seq} reads (pairs) in total.",
    "\\item ${sum_mapping_rate}$ of reads mapped to reference genome ({sum_refgenome}).",
    "\\item {sum_variants_number} variants called from data. \\\\",
    "\\begin{{itemize}}",
    "  \\item {sum_snp_number} SNVs (Single Nucleotide Variant)",
    "  \\item {sum_indel_number} Indels (Insertion and Deletion)",
    "  \\end{{itemize}}",
    "\\end{{itemize}}"
]
doclist['sec-analysis_pipeline'] = [
    "\\section{{Analysis Pipeline}}",
    "\\label{{sec:pipeline}}",
    "\\begin{{figure}}[h]",
    "  \\centering",
    "  \\includegraphics[width=0.75\\linewidth]{{{fig_pipeline}}}",
    "  \\caption{{MS Variant Calling Pipeline}}",
    "  \\label{{fig:pipeline}}",
    "\\end{{figure}}"
]
doclist['sec-analysis_result-head'] = [
    "\\section{{Analysis Result}}",
    "\\label{{sec:result}}"
]
doclist['sec-analysis_result-subsec-qc-head-single-notrim'] = [
    "\\subsection{{Quality Control}}",
    "\\label{{subsec:qc}}",
    "FASTQ file statistics:",
    "\\begin{{center}}",
    "  \\rowcolors{{1}}{{}}{{rowgray}}",
    "  \\begin{{tabular}}[h]{{l c}}",
    "    \\hline",
    "    Feature                & Statistics             \\\\",
    "    \\hline",
    "    Read Length Mean       & {qc_seq_length_mean}   \\\\",
    "    Read Length Min        & {qc_seq_length_min}    \\\\",
    "    Read Length Median     & {qc_seq_length_median} \\\\",
    "    Read Length Max        & {qc_seq_length_max}    \\\\",
    "    Read GC content Mean   & {qc_seq_gc_mean}       \\\\",
    "    Read GC content Min    & {qc_seq_gc_min}        \\\\",
    "    Read GC content Median & {qc_seq_gc_median}     \\\\",
    "    Read GC content Max    & {qc_seq_gc_max}        \\\\",
    "    \\hline",
    "  \\end{{tabular}}",
    "\\end{{center}}"
]
doclist['sec-analysis_result-subsec-qc-body-single-notrim'] = [
    "\\subsubsection{{Sequence Length Distribution}}",
    "\\label{{subsubsec:sequence_length}}",
    "Read length distribution is a good measure of the sequence quality. In",
    "most cases, the sequence length in the FASTQ file obtained from",
    "sequencing service companies should be mostly similar after their",
    "quality control. And the bad sequencing quality may lead to the",
    "cutting of sequencing by the sequencing machine. So, if the sequences",
    "mostly have same length or normally distributed, that means the",
    "quality of sequences might be acceptable. However, if the sequence",
    "length does not follow a normal distribution, the FASTQ file will need",
    "more attention. And in the quality trim stage, the short sequences and",
    "bad quality sequences should be discarded. The sequcne length",
    "distribution is shown in Figure \\ref{{fig:seq_length_distribution}}.",
    "\\begin{{center}}",
    "  \\rowcolors{{1}}{{}}{{rowgray}}",
    "  \\begin{{tabular}}[h]{{l c}}",
    "    \\hline",
    "    Read Length        & Nucleotides              \\\\",
    "    \\hline",
    "    Average Length     & {qc_seq_length_mean}        \\\\",
    "    Min Length         & {qc_seq_length_min}    \\\\",
    "    1st Quatile Length & {qc_seq_length_q1}     \\\\",
    "    Median Length      & {qc_seq_length_median} \\\\",
    "    3rd Quatile Length & {qc_seq_length_q3}     \\\\",
    "    Max Length         & {qc_seq_length_max}    \\\\",
    "    \\hline",
    "  \\end{{tabular}}",
    "\\end{{center}}",
    "\\begin{{figure}}[h]",
    "  \\centering",
    "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_seq_length_distribution}}}}}",
    "  \\caption{{Sequence Length Distribution}}",
    "  \\label{{fig:seq_length_distribution}}",
    "\\end{{figure}}",
    "\\subsubsection{{Base Quality Boxplot}}",
    "\\label{{subsubsec:base_quality_boxplot}}",
    "Base quality boxplot is used to show the bases' quality",
    "distribution. For the first several bases and the last several bases",
    "of the reads on Illumina platform, the base quality is generaly of",
    "large error rate, which should be trimmed before analysis. The base",
    "quality is shown in Figure \\ref{{fig:base_quality_boxplot}}.",
    "\\begin{{figure}}[h]",
    "  \\centering",
    "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_base_quality_boxplot}}}}}",
    "  \\caption{{Base Quality Boxplot}}",
    "  \\label{{fig:base_quality_boxplot}}",
    "\\end{{figure}}",
    "Phred score is used to evaluate the base quality, the larger the",
    "better. when the quality score is higher than 30, that means the error",
    "rate is lower than 0.001, which can be used as the empirical cutoff",
    "for the base quality. For the good quality FASTQ file, most of the",
    "bases have quality larger than 30. And for the bad quality FASTQ file,",
    "the first quartiles may be below 30 for many sites. So, for the",
    "following trimming stage, the bad FASTQ file should be trimmed and",
    "filtered.",
    "\\subsubsection{{Mean Sequence Quality Distribution}}",
    "\\label{{subsubsec:mean_quality}}",
    "The mean quality shows the quality of each sequence. The distribution",
    "of mean sequence quality shows how the quality of the input FASTQ",
    "file. For the good quality FASTQ file, most of the mean quality scores",
    "are high and distributing in one small range. Meanwhile, for the bad",
    "quality FASTQ files, except for the bad quality score, the range of",
    "the quality score is also wider. And sometimes there will be multiple",
    "peaks for the bad quality score distribution. Because there are",
    "sequences having universally low quality values, which could be a",
    "result of being poorly imaged in the sequencing stage. In the",
    "following filter stage, the cutoff of bad quality sequence can be",
    "chosen from the mean quality score distribution. The sequence quality",
    "is shown in Figure \\ref{{fig:seq_quality_distribution}}.",
    "\\begin{{figure}}[h]",
    "  \\centering",
    "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_seq_quality_distribution}}}}}",
    "  \\caption{{Mean Sequence Quality Distributon}}",
    "  \\label{{fig:seq_quality_distribution}}",
    "\\end{{figure}}",
    "\\subsubsection{{Base Content}}",
    "\\label{{subsubsec:basecontent}}",
    "The base GC content plot can show whether the AT or GC has different",
    "Though GC content is species specific, the basewise G percentage of",
    "per base should be consistent with the basewise C percentage, so as",
    "with the A and T percentage basewise. Because in the randomly",
    "generated library, the A and T bases are expected to be the same, so",
    "as with the G and C bases. The discordance of base content for AT and",
    "GC infers that the library or sequencing stage might have faults. The",
    "base content is shown in Figure \\ref{{fig:base_content}}.",
    "\\begin{{figure}}[h]",
    "  \\centering",
    "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_base_content}}}}}",
    "  \\caption{{Base Content}}",
    "  \\label{{fig:base_content}}",
    "\\end{{figure}}",
    "\\subsubsection{{Sequence GC Content Distribution}}",
    "\\label{{subsubsec:sequencegc}}",
    "The distribution of sequence GC content shows potential can provide",
    "the sequence quality information as well. GC content is not a random",
    "number in genomes of different species, which is driven by evolution",
    "in fact and varies continuously.  For a random sequence library, the",
    "GC content among difference bases should be similar. In ideal case,",
    "the GC content of sequencing data should reflect the genome GC",
    "content, but they should not show huge imbalance as well. However, the",
    "using of random primer in sequencing preparation might introduce the",
    "GC content bias at different bases.  The mean GC content of each",
    "sequence can generate a density plot. Ideally, the distribution of",
    "sequence GC content should be normally distributed. If the",
    "distribution diverges from normal distribution a lot, that might be a",
    "sign of low-quality or contamination of sequences from other",
    "organisms. The sequence GC content distribution is shown in",
    "Figure \\ref{{fig:seq_gc}}.",
    "\\begin{{center}}",
    "  \\rowcolors{{1}}{{}}{{rowgray}}",
    "  \\begin{{tabular}}[h]{{l c}}",
    "    \\hline",
    "    Read GC content        & Percentage (\%)    \\\\",
    "    \\hline",
    "    Average GC content     & {qc_seq_gc_mean}   \\\\",
    "    Min GC content         & {qc_seq_gc_min}    \\\\",
    "    1st Quatile GC content & {qc_seq_gc_q1}     \\\\",
    "    Median GC content      & {qc_seq_gc_median} \\\\",
    "    3rd Quatile GC content & {qc_seq_gc_q3}     \\\\",
    "    Max GC content         & {qc_seq_gc_max}    \\\\",
    "    \\hline",
    "  \\end{{tabular}}",
    "\\end{{center}}",
    "\\begin{{figure}}[h]",
    "  \\centering",
    "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_seq_gc}}}}}",
    "  \\caption{{Sequence GC Content}}",
    "  \\label{{fig:seq_gc}}",
    "\\end{{figure}}"
]

doclist['sec-analysis_result-subsec-qc-head-single-trim'] = [
    "\\subsection{{Quality Control}}",
    "\\label{{subsec:qc}}",
    "FASTQ file statistics:",
    "\\begin{{center}}",
    "  \\rowcolors{{1}}{{}}{{rowgray}}",
    "  \\begin{{tabular}}[h]{{l c}}",
    "    \\hline",
    "    Feature                & Before Quality Control & After Quality Control  \\\\",
    "    \\hline",
    "    Read Length Mean       & {qc_seq_length_mean}   & {qc2_seq_length_mean}  \\\\",
    "    Read Length Min        & {qc_seq_length_min}    & {qc2_seq_length_min}   \\\\",
    "    Read Length Median     & {qc_seq_length_median} & {qc2_seq_length_median}\\\\",
    "    Read Length Max        & {qc_seq_length_max}    & {qc2_seq_length_max}   \\\\",
    "    Read GC content Mean   & {qc_seq_gc_mean}       & {qc2_seq_gc_mean}      \\\\",
    "    Read GC content Min    & {qc_seq_gc_min}        & {qc2_seq_gc_min}       \\\\",
    "    Read GC content Median & {qc_seq_gc_median}     & {qc2_seq_gc_median}    \\\\",
    "    Read GC content Max    & {qc_seq_gc_max}        & {qc2_seq_gc_max}       \\\\",
    "    \\hline",
    "  \\end{{tabular}}",
    "\\end{{center}}"
]

doclist['sec-analysis_result-subsec-qc-body-single-trim'] = [
    "\\subsubsection{{Sequence Length Distribution}}",
    "\\label{{subsubsec:sequence_length}}",
    "Read length distribution is a good measure of the sequence quality. In",
    "most cases, the sequence length in the FASTQ file obtained from",
    "sequencing service companies should be mostly similar after their",
    "quality control. And the bad sequencing quality may lead to the",
    "cutting of sequencing by the sequencing machine. So, if the sequences",
    "mostly have same length or normally distributed, that means the",
    "quality of sequences might be acceptable. However, if the sequence",
    "length does not follow a normal distribution, the FASTQ file will need",
    "more attention. And in the quality trim stage, the short sequences and",
    "bad quality sequences should be discarded. The sequcne length",
    "distribution is shown in Figure \\ref{{fig:seq_length_distribution}}.",
    "\\begin{{center}}",
    "  \\rowcolors{{1}}{{}}{{rowgray}}",
    "  \\begin{{tabular}}[h]{{l c}}",
    "    \\hline",
    "    Read Length        & Before Quality Control & After Quality Control   \\\\",
    "    \\hline",
    "    Average Length     & {qc_seq_length_mean}   & {qc2_seq_length_mean}   \\\\",
    "    Min Length         & {qc_seq_length_min}    & {qc2_seq_length_min}    \\\\",
    "    1st Quatile Length & {qc_seq_length_q1}     & {qc2_seq_length_q1}     \\\\",
    "    Median Length      & {qc_seq_length_median} & {qc2_seq_length_median} \\\\",
    "    3rd Quatile Length & {qc_seq_length_q3}     & {qc2_seq_length_q3}     \\\\",
    "    Max Length         & {qc_seq_length_max}    & {qc2_seq_length_max}    \\\\",
    "    \\hline",
    "  \\end{{tabular}}",
    "\\end{{center}}",
    "\\begin{{figure}}[h]",
    "  \\begin{{minipage}}[0.49\\linewidth]",
    "    \\centering",
    "    \\includegraphics[width=0.9\\linewidth]{{{{{fig_qc_seq_length_distribution}}}}}",
    "    \\subcaption{{Before Quality Control}}",
    "    \\label{{fig:seq_length_distribution_before}}",
    "  \\end{{minipage}}",
    "  \\begin{{minipage}}[0.49\\linewidth]",
    "    \\centering",
    "    \\includegraphics[width=0.9\\linewidth]{{{{{fig_qc2_seq_length_distribution}}}}}",
    "    \\subcaption{{After Quality Control}}",
    "    \\label{{fig:seq_length_distribution_after}}",
    "  \\end{{minipage}}",
    "  \\caption{{Sequence Length Distribution}}",
    "  \\label{{fig:seq_length_distribution_after}}",
    "\\end{{figure}}",
    "\\subsubsection{{Base Quality Boxplot}}",
    "\\label{{subsubsec:base_quality_boxplot}}",
    "Base quality boxplot is used to show the bases' quality",
    "distribution. For the first several bases and the last several bases",
    "of the reads on Illumina platform, the base quality is generaly of",
    "large error rate, which should be trimmed before analysis. The base",
    "quality is shown in Figure \\ref{{fig:base_quality_boxplot}}.",
    "\\begin{{figure}}[h]",
    "  \\begin{{minipage}}[0.49\\linewidth]",
    "    \\centering",
    "    \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_base_quality_boxplot}}}}}",
    "    \\subcaption{{Before Quality Control}}",
    "    \\label{{fig:base_quality_boxplot_before}}",
    "  \\end{{minipage}}",
    "  \\begin{{minipage}}[0.49\\linewidth]",
    "    \\centering",
    "    \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc2_base_quality_boxplot}}}}}",
    "    \\subcaption{{After Quality Control}}",
    "    \\label{{fig:base_quality_boxplot_after}}",
    "  \\end{{minipage}}",
    "  \\caption{{Base Quality Boxplot}}",
    "  \\label{{fig:base_quality_boxplot}}",
    "\\end{{figure}}",
    "Phred score is used to evaluate the base quality, the larger the",
    "better. when the quality score is higher than 30, that means the error",
    "rate is lower than 0.001, which can be used as the empirical cutoff",
    "for the base quality. For the good quality FASTQ file, most of the",
    "bases have quality larger than 30. And for the bad quality FASTQ file,",
    "the first quartiles may be below 30 for many sites. So, for the",
    "following trimming stage, the bad FASTQ file should be trimmed and",
    "filtered.",
    "\\subsubsection{{Mean Sequence Quality Distribution}}",
    "\\label{{subsubsec:mean_quality}}",
    "The mean quality shows the quality of each sequence. The distribution",
    "of mean sequence quality shows how the quality of the input FASTQ",
    "file. For the good quality FASTQ file, most of the mean quality scores",
    "are high and distributing in one small range. Meanwhile, for the bad",
    "quality FASTQ files, except for the bad quality score, the range of",
    "the quality score is also wider. And sometimes there will be multiple",
    "peaks for the bad quality score distribution. Because there are",
    "sequences having universally low quality values, which could be a",
    "result of being poorly imaged in the sequencing stage. In the",
    "following filter stage, the cutoff of bad quality sequence can be",
    "chosen from the mean quality score distribution. The sequence quality",
    "is shown in Figure \\ref{{fig:seq_quality_distribution}}.",
    "\\begin{{figure}}[h]",
    "  \\begin{{minipage}}[0.49\\linewidth]",
    "    \\centering",
    "    \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_seq_quality_distribution}}}}}",
    "    \\subcaption{{Before Quality Control}}",
    "    \\label{{fig:seq_quality_distribution_before}}",
    "  \\end{{minipage}}",
    "  \\begin{{minipage}}[0.49\\linewidth]",
    "    \\centering",
    "    \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc2_seq_quality_distribution}}}}}",
    "    \\subcaption{{After Quality Control}}",
    "    \\label{{fig:seq_quality_distribution_after}}",
    "  \\end{{minipage}}",
    "  \\caption{{Mean Sequence Quality Distributon}}",
    "  \\label{{fig:seq_quality_distribution}}",
    "\\end{{figure}}",
    "\\subsubsection{{Base Content}}",
    "\\label{{subsubsec:basecontent}}",
    "The base GC content plot can show whether the AT or GC has different",
    "Though GC content is species specific, the basewise G percentage of",
    "per base should be consistent with the basewise C percentage, so as",
    "with the A and T percentage basewise. Because in the randomly",
    "generated library, the A and T bases are expected to be the same, so",
    "as with the G and C bases. The discordance of base content for AT and",
    "GC infers that the library or sequencing stage might have faults. The",
    "base content is shown in Figure \\ref{{fig:base_content}}.",
    "\\begin{{figure}}[h]",
    "  \\begin{{minipage}}[0.49\\linewidth]",
    "    \\centering",
    "    \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_base_content}}}}}",
    "    \\subcaption{{Before Quality Control}}",
    "    \\label{{fig:base_content_before}}",
    "  \\end{{minipage}}",
    "  \\begin{{minipage}}[0.49\\linewidth]",
    "    \\centering",
    "    \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc2_base_content}}}}}",
    "    \\subcaption{{After Quality Control}}",
    "    \\label{{fig:base_content_after}}",
    "  \\end{{minipage}}",
    "  \\caption{{Base Content}}",
    "  \\label{{fig:base_content}}",
    "\\end{{figure}}",
    "\\subsubsection{{Sequence GC Content Distribution}}",
    "\\label{{subsubsec:sequencegc}}",
    "The distribution of sequence GC content shows potential can provide",
    "the sequence quality information as well. GC content is not a random",
    "number in genomes of different species, which is driven by evolution",
    "in fact and varies continuously.  For a random sequence library, the",
    "GC content among difference bases should be similar. In ideal case,",
    "the GC content of sequencing data should reflect the genome GC",
    "content, but they should not show huge imbalance as well. However, the",
    "using of random primer in sequencing preparation might introduce the",
    "GC content bias at different bases.  The mean GC content of each",
    "sequence can generate a density plot. Ideally, the distribution of",
    "sequence GC content should be normally distributed. If the",
    "distribution diverges from normal distribution a lot, that might be a",
    "sign of low-quality or contamination of sequences from other",
    "organisms. The sequence GC content distribution is shown in",
    "Figure \\ref{{fig:seq_gc}}.",
    "\\begin{{center}}",
    "  \\rowcolors{{1}}{{}}{{rowgray}}",
    "  \\begin{{tabular}}[h]{{l c}}",
    "    \\hline",
    "    Read GC content        & Before Quality Control & After Quality Control \\\\",
    "    \\hline",
    "    Average GC content     & {qc_seq_gc_mean}       & {qc2_seq_gc_mean}     \\\\",
    "    Min GC content         & {qc_seq_gc_min}        & {qc2_seq_gc_min}      \\\\",
    "    1st Quatile GC content & {qc_seq_gc_q1}         & {qc2_seq_gc_q1}       \\\\",
    "    Median GC content      & {qc_seq_gc_median}     & {qc2_seq_gc_median}   \\\\",
    "    3rd Quatile GC content & {qc_seq_gc_q3}         & {qc2_seq_gc_q3}       \\\\",
    "    Max GC content         & {qc_seq_gc_max}        & {qc2_seq_gc_max}      \\\\",
    "    \\hline",
    "  \\end{{tabular}}",
    "\\end{{center}}",
    "\\begin{{figure}}[h]",
    "  \\begin{{minipage}}[0.45\\linewidth]",
    "    \\centering",
    "    \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc_seq_gc}}}}}",
    "    \\subcaption{{Before Quality Control}}",
    "    \\label{{fig:seq_gc_before}}",
    "  \\end{{minipage}}",
    "  \\begin{{minipage}}[0.45\\linewidth]",
    "    \\centering",
    "    \\includegraphics[width=0.75\\linewidth]{{{{{fig_qc2_seq_gc}}}}}",
    "    \\subcaption{{After Quality Control}}",
    "    \\label{{fig:seq_gc_after}}",
    "  \\end{{minipage}}",
    "  \\caption{{Sequence GC Content}}",
    "  \\label{{fig:seq_gc}}",
    "\\end{{figure}}"
]

doclist['sec-analysis_result-subsec-mapping'] = [
    "\\subsection{{Mapping}}",
    "\\label{{subsec:mapping}}",
    "Mapping stage finds",
    "the reads position on the reference genome. The",
    "result of mapping shows the coverage and depth of sequencing. The",
    "following analysis depends on the accurate mapping result. The mapping",
    "result is saved as SAM/BAM file. The possible PCR duplicated sequences",
    "will be marked in the process. BWA is one of the opensource software",
    "packages for mapping high-throughput sequencing data to reference",
    "genome that generates accurate results in a fast speed.",
    "\\subsubsection{{Mapping Process}}",
    "\\label{{subsubsec:mappingprocess}}",
    "\\begin{{center}}",
    "  \\rowcolors{{1}}{{}}{{rowgray}}",
    "  \\begin{{tabular}}[h]{{l c c}}",
    "    \\hline",
    "    Process              & Software        & Version             \\\\",
    "    \\hline",
    "    Mapping Software     & {map_maptool}   & {map_maptool_v}     \\\\",
    "    Reference Genome     & {map_refgenome} & ({map_refgenome_v}) \\\\",
    "    SAM/BAM Process      & {map_samtool}   & {map_samtool_v}     \\\\",
    "    Mark PCR Duplicate   & {map_mkdup}     & {map_mkdup_v}       \\\\",
    "    \\hline",
    "  \\end{{tabular}}",
    "\\end{{center}}",
    "\\subsubsection{{Mapping Statistics}}",
    "\\label{{subsubsec:mappingstat}}",
    "\\begin{{center}}",
    "  \\begin{{tabular}}[h]{{l r}}",
    "    \\hline",
    "    Item            & Count              \\\\",
    "    \\hline",
    "    Raw Total Reads & {map_total_seq}    \\\\",
    "    Mapped Reads    & {map_mapped_seq}   \\\\",
    "    Unmapped Reads  & {map_unmapped_seq} \\\\",
    "    \\hline",
    "  \\end{{tabular}}",
    "\\end{{center}}"
]

doclist['sec-analysis_result-subsec-cv'] = [
    "\\subsection{{Variant Calling}}",
    "\\label{{subsec:variant}}",
    "GATK 4 is used as the main variant calling tools in the analysis",
    "pipeline. SNVs (Single Nucleotide Variant) and InDels (Insertion and",
    "Deletion) are detected by the GATK software. In this pipeline, SNVs",
    "are considered.",
    "\\begin{{description}}",
    "\\item[Variants Calling] GATK 4.0.8 % Setting: if different version",
    "\\item[Variants Filter] SNP filters are applied. Good quality SNVs should satisfy the following conditions. \\\\",
    "  \\begin{{itemize}}",
    "  \\item Quality Depth (QD) $\\geq 2$",
    "  \\item Mapping Quality (MQ) $\\geq 40$",
    "  \\item MQRankSum $\\geq -12.5$",
    "  \\item Fisher test of strand bias (FS) $\\leq 60$",
    "  \\item ReadPosRankSum $\\geq -8$",
    "  \\end{{itemize}}",
    "\\end{{description}}",
    "Result of variants calling is saved as VCF (variants calling format).",
    "\\begin{{center}}",
    "  \\rowcolors{{1}}{{}}{{rowgray}}",
    "  \\begin{{tabular}}[h]{{l r}}",
    "    \\hline",
    "    Item           & Count              \\\\",
    "    \\hline",
    "    Total Variants & {vc_total_variant} \\\\",
    "    SNV            & {vc_snv_number}    \\\\",
    "    Indel          & {vc_indel_number}  \\\\",
    "    \\hline",
    "  \\end{{tabular}}",
    "\\end{{center}}",
    "The condition of SNV types is shown in Figure \\ref{{fig:snvtype}}.",
    "\\begin{{figure}}[h]",
    "  \\centering",
    "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_vc_snvtype}}}}}",
    "  \\caption{{SNV types}}",
    "  \\label{{fig:snvtype}}",
    "\\end{{figure}}"
]

doclist['sec-analysis_result-subsec-annotation'] = [
    "\\subsection{{Variant Annotation}}",
    "\\label{{subsec:annotationg}}",
    "Variant annotation links the variants with the existing variant",
    "database information. The variants annotation uses ANNOVAR software.",
    "\\begin{{description}}",
    "\\item[RefSeq Gene] (refGene) RefSeq Gene is the NCBI database for gene",
    "  information. Annotation with refGene shows the gene that variants",
    "  belonging to.",
    "\\item[ExAC] (exac) Exome Aggregation Consortium database, more than",
    "  60K individuals SNV data.",
    "\\item[dbSNP] (avsnp) NCBI Database of single nucleotide polymorphisms",
    "  (SNPs) and other small variants.",
    "\\item[dbNSFP] (dbnsfp) dbNSFP is a database of all non-synonymous",
    "  single-nucleotide variants.",
    "\\item[ClinVar] (clinvar\\_20160302) ClinVar is the NCBI database storeing",
    "  variants associated with disease.",
    "\\end{{description}}",
    "\\subsubsection{{Annotation Statistics}}",
    "\\label{{subsubsec:annotationstatistics}}",
    "The annotation of each database is shown in Figure \\ref{{fig:snvannobar}}.",
    "\\begin{{figure}}[h]",
    "  \\centering",
    "  \\includegraphics[width=0.75\\linewidth]{{{{{fig_vc_annovcf}}}}}",
    "  \\caption{{Variant Annotation}}",
    "  \\label{{fig:snvannobar}}",
    "\\end{{figure}}"
]

doclist['sec-analysis_result-subsec-annotation-table'] = [
    "The variants annotated with connection to diseases have been list as follows (\\ref{{tab:vcfclinvar}}):"
]
####################

import os
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description='Generate Report TeX file.')
parser.add_argument('--output', nargs='?', type=str, default='result.tex')
parser.add_argument('--report-title', nargs='?', default='Exome Variant Analysis Report')
parser.add_argument('--report-author', nargs='?', default='MS Health Care')
parser.add_argument('--sum-total-seq', nargs='?', type=float, default=0)
parser.add_argument('--sum-mapping-rate', nargs='?', type=float, default=0)
parser.add_argument('--sum-refgenome', nargs='?', default='GRCh38')
parser.add_argument('--sum-variants-number', nargs='?', type=float, default=0)
parser.add_argument('--sum-snp-number', nargs='?', type=float, default=0)
parser.add_argument('--sum-indel-number', type=float, default=0)
parser.add_argument('--seq-mode', nargs='?', default='Single-End')
parser.add_argument('--total-seq', nargs='?', type=float, default=0)
parser.add_argument('--total-pair', nargs='?', type=float, default=0)
parser.add_argument('--qc-trimed', action='store_true', default=False)
parser.add_argument('--qc-warn', nargs='*', default=[])
parser.add_argument('--qc-fail', nargs='*', default=[])
parser.add_argument('--qc-warn2', nargs='*', default=[])
parser.add_argument('--qc-fail2', nargs='*', default=[])
parser.add_argument('--qc-seq-length-mean', nargs='?', type=float, default=0)
parser.add_argument('--qc-seq-length-min', nargs='?', type=float, default=0)
parser.add_argument('--qc-seq-length-median', nargs='?', type=float, default=0)
parser.add_argument('--qc-seq-length-max', nargs='?', type=float, default=0)
parser.add_argument('--qc-seq-length-q1', nargs='?', type=float, default=0)
parser.add_argument('--qc-seq-length-q3', nargs='?', type=float, default=0)
parser.add_argument('--qc-seq-gc-mean', nargs='?', type=float, default=0)
parser.add_argument('--qc-seq-gc-min', nargs='?', type=float, default=0)
parser.add_argument('--qc-seq-gc-median', nargs='?', type=float, default=0)
parser.add_argument('--qc-seq-gc-max', nargs='?', type=float, default=0)
parser.add_argument('--qc-seq-gc-q1', nargs='?', type=float, default=0)
parser.add_argument('--qc-seq-gc-q3', nargs='?', type=float, default=0)
parser.add_argument('--qc2-seq-length-mean', nargs='?', type=float, default=0)
parser.add_argument('--qc2-seq-length-min', nargs='?', type=float, default=0)
parser.add_argument('--qc2-seq-length-median', nargs='?', type=float, default=0)
parser.add_argument('--qc2-seq-length-max', nargs='?', type=float, default=0)
parser.add_argument('--qc2-seq-length-q1', nargs='?', type=float, default=0)
parser.add_argument('--qc2-seq-length-q3', nargs='?', type=float, default=0)
parser.add_argument('--qc2-seq-gc-mean', nargs='?', type=float, default=0)
parser.add_argument('--qc2-seq-gc-min', nargs='?', type=float, default=0)
parser.add_argument('--qc2-seq-gc-median', nargs='?', type=float, default=0)
parser.add_argument('--qc2-seq-gc-max', nargs='?', type=float, default=0)
parser.add_argument('--qc2-seq-gc-q1', nargs='?', type=float, default=0)
parser.add_argument('--qc2-seq-gc-q3', nargs='?', type=float, default=0)
parser.add_argument('--map-maptool', nargs='?', default='BWA MEM')
parser.add_argument('--map-maptool-v', nargs='?', default='0.7.17')
parser.add_argument('--map-refgenome', nargs='?', default='GRCh38')
parser.add_argument('--map-refgenome-v', nargs='?', default='UCSC hg38')
parser.add_argument('--map-samtool', nargs='?', default='samtools')
parser.add_argument('--map-samtool-v', nargs='?', default='1.7')
parser.add_argument('--map-mkdup', nargs='?', default='Picard')
parser.add_argument('--map-mkdup-v', nargs='?', default='2.18.11')
parser.add_argument('--map-total-seq', nargs='?', type=float, default=0)
parser.add_argument('--map-mapped-seq', nargs='?', type=float, default=0)
parser.add_argument('--map-unmapped-seq', nargs='?', type=float, default=0)
parser.add_argument('--vc-total-variant', nargs='?', type=float, default=0)
parser.add_argument('--vc-snv-number', nargs='?', type=float, default=0)
parser.add_argument('--vc-indel-number', nargs='?', type=float, default=0)
parser.add_argument('--vcf', nargs='?', type=str)
parser.add_argument('--figpath', nargs='?', default='fig/')
parser.add_argument('--fig-pipeline', nargs='?', default='')
parser.add_argument('--fig-qc-seq-length-distribution', nargs='?', default='')
parser.add_argument('--fig-qc2-seq-length-distribution',nargs='?', default='')
parser.add_argument('--fig-qc-base-quality-boxplot', nargs='?', default='')
parser.add_argument('--fig-qc2-base-quality-boxplot', nargs='?', default='')
parser.add_argument('--fig-qc-seq-quality-distribution', nargs='?', default='')
parser.add_argument('--fig-qc2-seq-quality-distribution', nargs='?', default='')
parser.add_argument('--fig-qc-base-content', nargs='?', default='')
parser.add_argument('--fig-qc2-base-content', nargs='?', default='')
parser.add_argument('--fig-qc-seq-gc', nargs='?', default='')
parser.add_argument('--fig-qc2-seq-gc', nargs='?', default='')
parser.add_argument('--fig-vc-snvtype', nargs='?', default='')
parser.add_argument('--fig-vc-annovcf', nargs='?', default='')


argdict = vars(parser.parse_args())
argdict['report_title'] = argdict['report_title'].replace('_', '\_')
argdict['report_author'] = argdict['report_author'].replace('_', '\_')

####################
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

clnvcf_table = "\n".join(
    [
        "\\begin{{table}}",
        "  \\centering",
        "  \\caption{{Disease associated variants}}",
        "  \\rowcolors{{1}}{{}}{{rowgray}}",
        "  \\begin{{tabular}}{{c | r | c | c | p{{0.15\\linewidth}} | p{{0.35\\linewidth}}}}",
        "    \\hline",
        "    \\multicolumn{{1}}{{c|}}{{Chrom}} & \\multicolumn{{1}}{{c|}}{{Position}} & \\multicolumn{{1}}{{c|}}{{Reference}} & \\multicolumn{{1}}{{c|}}{{Sample}} & \\multicolumn{{1}}{{c|}}{{Situation}} & \\multicolumn{{1}}{{c}}{{Disease}} \\\\",
        "    \\hline"
    ] + list(clnvcf.apply(
        lambda x: ' & '.join(
            str(i) for i in x[['CHROM', 'POS', 'REF', 'ALT', 'CLNSIG']]
        ) + ' & ' + ' | '.join(x['CLNDN']) + ' \\\\',
        axis=1
    )) + [
        "  \\label{{tab:vcfclinvar}}",
        "  \\end{{tabular}}",
        "\\end{{table}}"
    ]
)
 
####################

tex = doclist['dochead'] + \
      doclist['sec-summary'] + \
      doclist['sec-analysis_pipeline'] + \
      doclist['sec-analysis_result-head']

warnname = {
    'base_quality': 'Base Quality',
    'seq_quality': 'Sequence Mean Quality',
    'base_content': 'Base Content',
    'seq_gc': 'Sequence GC Content',
    'seq_length_distribution': 'Sequence Length Distribution'
}
qc_message_warn = {
    'base_quality': 'Any base is less than 10, or any median is less than 25. That maight be the drop of quality at the nd of the sequences. A quality trimming might be performed.',
    'seq_quality': 'If the mode sequence mean quality is below 27. A quality trimming might be performed.',
    'base_content': 'Difference between A and T, or G and C is greater than $10\%$. The sequences should be checked.',
    'seq_gc': 'More than $15\%$ reads deviated from normal distribution. There might be contamination.',
    'seq_length_distribution':'All sequences not the same length. Can be ignored.'
}
qc_message_fail = {
    'base_quality': 'Any base is less than 5, or any median is less than 20. That maight be the drop of quality at the end of the sequences. A quality trimming should be performed.',
    'seq_quality': 'If the mode sequence mean quality is below 20. A quality trimming should be performed.',
    'base_content': 'Difference between A and T, or G and C is greater than $20\%$. The sequences should be checked.',
    'seq_gc': 'More than $30\%$ reads deviated from normal distribution. There might be contamination.',
    'seq_length_distribution':'Any sequence with 0 length. Can be ignored.'
}

warning_messages = []
failure_messages = []


if len(argdict['qc_warn']) != 0:
    warning_messages = [
        "\\begin{{description}}",
    ] + [
        "\\item[\\textcolor{{orange}}{{WARNING}}] \\textbf{{" + warnname[x] + ":}} " + qc_message_warn[x] for x in argdict['qc_warn']
    ] + [
        "\\end{{description}}"
    ]

if len(argdict['qc_fail']) != 0:
    failure_messages = [
        "\\begin{{description}}",
    ] + [
        "\\item[\\textcolor{{red}}{{FAILURE}}] \\textbf{{" + warnname[x] + ":}} " + qc_message_fail[x] for x in argdict['qc_fail']
    ] + [
        "\\end{{description}}"
    ]


if argdict['seq_mode'] == 'Single-End':
    if argdict['qc_trimed']:
        tex = tex + doclist['sec-analysis_result-subsec-qc-head-single-trim'] + \
              warning_messages + failure_messages + \
              doclist['sec-analysis_result-subsec-qc-body-single-trim']
    else:
        tex = tex + doclist['sec-analysis_result-subsec-qc-head-single-notrim'] + \
              warning_messages + failure_messages + \
              doclist['sec-analysis_result-subsec-qc-body-single-notrim']
else:
    if len(argdict['qc_warn2']) != 0:
        warning_messages = warning_messages + [
            "And for the paired end input sequences:",
            "\\begin{{description}}",
        ] + [
            "\\item[" + x + "]" + qc_message_warn[x] for x in argdict['qc_warn2']
        ] + [
            "\\end{{description}}"
        ]
    if len(argdict['qc_fail']) != 0:
        failure_messages = failure_messages + [
            "And for the paired end input sequences:",
            "\\begin{{description}}",
        ] + [
            "\\item[" + x + "]" + qc_message_warn[x] for x in argdict['qc_fail2']
        ] + [
            "\\end{{description}}"
        ]
    if argdict['qc_trimed']:
        tex = tex + doclist['sec-analysis_result-subsec-qc-head-paired-trim'] + \
              warning_messages + failure_messages + \
              doclist['seq-analysis_result-subsec-qc-body-paired-notrim']
    else:
        tex = tex + doclist['sec-analysis_result-subsec-qc-head-paired-notrim'] + \
              warning_messages + failure_messages + \
              doclist['seq-analysis_result-subsec-qc-body-paired-notrim']

tex = tex + doclist['sec-analysis_result-subsec-mapping'] + \
      doclist['sec-analysis_result-subsec-cv'] + \
      doclist['sec-analysis_result-subsec-annotation'] + \
      doclist['sec-analysis_result-subsec-annotation-table'] + \
      clnvcf_table + \
      doclist['docend']

resulttex = '\n'.join(tex).format(**argdict)

with open(argdict['output'], 'w') as f:
    f.write(resulttex)

################################################################################
