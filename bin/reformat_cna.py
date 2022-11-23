import os
import pandas as pd

def reformat_cna(cna_file, sample):
    """
    :param cna_file: somatic CNA file from Sarek.
    :param sample:   string denoting sample name 'meta.id'
    """

    # Add to this dictionary when familiar with ASCAT, Control-FREEC
    index_keys = { 'cnvkit' : ['chromosome', 'start', 'end', 'probes', 'log2']}

    # Colnames expected by PCGR
    final_columns = ['Sample', 'Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean']

    # Stage file and guess tool
    df = pd.read_csv(cna_file, delim_whitespace=True)
    tool = guess_tool(df)

    # Log info
    print(f"Reading {cna_file} file.\n\nProcessing {len(df.index) -1} lines.\nCNA tool used:\n>{tool}\n")

    # Skip if test-dataset
    if tool != 'pcgr':

        # Subset columns corresponding to expected PCGR input.
        index = index_keys[tool]
        df = df[df.columns.intersection(index)]

        # Insert 'Sample' columns and rename cols
        df.insert(0, "Sample", sample)
        df.set_axis(final_columns, 1, inplace=True)

    # output TSV
    out = sample + "." + tool + '.tsv'
    df.to_csv(out, sep="\t", index=False)

def guess_tool(df):
    """
    Mainly for indexing, posterity by adding tool name to filename.
    'pcgr': only likely to occur using test data provided by PCGR
            i.e when using a test-dataset.
    """

    tool_keys = {
        "pcgr":   ['Sample', 'Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean'],
        "cnvkit": ['chromosome', 'start', 'end', 'gene', 'log2', 'depth', 'probes', 'weight', 'ci_lo', 'ci_hi']
        }

    # subset list preferred (returns string) to list comprehension which returns a set.
    tool = list(tool_keys.keys())[list(tool_keys.values()).index(list(df.axes[1]))]
    return tool


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
