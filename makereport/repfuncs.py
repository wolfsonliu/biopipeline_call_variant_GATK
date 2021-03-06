#! /usr/bin/env python3
# Functions in this file is used as assistant function for the makereport scripts.
# The functions solve small problems in the main script of makereport.

# a wrapper function to write data to file. 
def write_list_to_file(data, filepath):
    with open(filepath, 'w') as f:
        f.write('\n'.join(data))

# used for catching data from FastQC files.
def get_list_block(seqlist, start, end):
    result = list()
    inblock = False
    for x in seqlist:
        if x.find(start) >= 0:
            inblock = True
            continue
        elif x.find(end) >=0 and inblock:
            break
        if inblock:
            result.append(x)
    return result

# get lines starts with
def get_list_startwith(seqlist, start):
    return [line for line in seqlist if line.find(start) == 0]

# used to get data from samtools stat
def get_list_cut(seqlist, sep, cols):
    def parse_cols(colstring, maxcol):
        cslist = colstring.split(',')
        col_num = list()
        for i in cslist:
            if i.isnumeric():
                col_num.append(int(i))
            else:
                ilist = i.split('-')
                endcol = maxcol
                if ilist[1].isnumeric():
                    endcol = ilist[1]
                col_num.extend(list(range(int(ilist[0]), int(endcol))))
        return col_num
    max_col = len(seqlist[0].split(sep))
    columns = parse_cols(cols, max_col)
    return [sep.join([line.split(sep)[i] for i in columns]) for line in seqlist]
