#!/usr/bin/env python3

import os
import re
import pandas as pd

def pcgr_ready_vcf(sample):
    # pass meta.id to the script to remove .command* etc files from file search:
    sample_id = sample
    suffixes = (".cns", ".tbi")
    r = re.compile(f"{sample_id}*")
    sample_files = os.listdir("./")
    sample_files = list(filter(r.match, sample_files))
    sample_files = [ file for file in sample_files if not file.endswith(suffixes) ]

    # there will be a sample_keys.txt file in the work directory
    # load it as key_Df and remove from the VCF list before running vcf2tsvpy
    key = "".join([file for file in sample_files if file.endswith('.txt')])
    key_df = pd.read_table(key, header=None, sep="\t")
    key_df.index = pd.MultiIndex.from_arrays(key_df.values.T[(0,1,2,3),])
    sample_files.remove(key)

    # stage fields we are interested in capturing
    fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'TAF', 'TDP', 'NAF', 'NDP']

    # create dict of dataframes containing dataframes of the abpve fields.
    # strategy is to average TDP,NDP, TAF,NAF and append TAL from $meta.id_keys.txt

    vcf_dict = {}
    for idx, x in enumerate(sample_files):
        idx = idx + 1
        os.system(f'vcf2tsvpy --input_vcf {x} --out_tsv {idx}.tmp --skip_genotype_data --keep_rejected_calls')
        os.system(f'tail -n +2 {idx}.tmp > {idx}.tsv && rm {idx}.tmp')
        df = pd.read_table(f'{idx}.tsv', usecols=fields, low_memory=True, sep="\t")
        df.index = pd.MultiIndex.from_arrays(df.values.T[(0,1,3,4),])
        vcf_dict[idx] = df

    # average said columns
    # convert QUAL to floats so we can take the max. Do not want to discard quality scores if they exist.
    avg_df = pd.concat(list(vcf_dict.values())).assign(QUAL = lambda d: d['QUAL'].apply(lambda s: 0 if s == "." else float(s) )).groupby(level=[0,1,2,3]).agg({'TAF':'mean', 'TDP':'mean', 'NAF':'mean', 'NDP':"mean", 'QUAL':"max"})

    # round up TDP, NDP, vcf header specifies INT not FLOAT
    avg_df['TDP'] = avg_df['TDP'].apply(lambda x: str(round(x)))
    avg_df['NDP'] = avg_df['NDP'].apply(lambda x: str(round(x)))

    # removing duplicate lines should have the same affect on len() as averaging.
    # run a quick sanity check
    master_df = pd.concat(list(vcf_dict.values()))
    # Keep rs SNP ids if they are present
    master_df = master_df.sort_values('ID', ascending=False)
    master_df = master_df[~master_df.index.duplicated(keep='first')]

    assert len(avg_df) == len(master_df), 'averaged values dataframe does not match master dataframe with duplicate index vals removed'

    # add the rest to avg df. revert QUAL back to object
    avg_df = avg_df.assign( QUAL = lambda d: d['QUAL'].apply(lambda s: "." if s == 0.0 else float(s) ))
    avg_df['ID'] = master_df['ID']
    avg_df['FILTER'] = master_df['FILTER']

    # append TAL from keys
    key_df = key_df[~key_df.index.duplicated(keep='first')]
    avg_df['TAL'] = key_df[4]

    # convert index back to proper names
    avg_df = avg_df.reset_index()
    avg_df.rename(columns={'level_0':'#CHROM', 'level_1':'POS', 'level_2':'REF', 'level_3':'ALT'}, inplace=True)
    fields.append('TAL')
    fields[0] = '#CHROM'
    avg_df = avg_df[fields]


    # add other cols for merging.
    avg_df['GT'] = 'GT'
    avg_df['DPC'] = 'DPC'
    avg_df['DPT'] = 'DPT'
    avg_df['ADC'] = 'ADC'
    avg_df['ADT'] = 'ADT'
    avg_df['AL'] = 'AL'

    # genotype info
    avg_df['geno_GT'] = "0/1"
    avg_df['geno_DPT'] = avg_df['TDP']
    avg_df['geno_DPC'] = avg_df['NDP']
    avg_df['geno_ADT'] = ".,."
    avg_df['geno_ADC'] = ".,."
    avg_df['geno_AL'] = "1"

    avg_df['TAF'] = avg_df['TAF'].apply(lambda x: "{}{}".format('TAF=', x))
    avg_df['TDP'] = avg_df['TDP'].apply(lambda x: "{}{}".format('TDP=', x))
    avg_df['NAF'] = avg_df['NAF'].apply(lambda x: "{}{}".format('NAF=', x))
    avg_df['NDP'] = avg_df['NDP'].apply(lambda x: "{}{}".format('NDP=', x))
    avg_df['TAL'] = avg_df['TAL'].apply(lambda x: "{}{}".format('TAL=', x))

    avg_df['INFO'] = avg_df[['TAF', 'TDP', 'NAF', 'NDP', 'TDP', 'TAL']].apply(lambda x: ';'.join(x[x.notnull()]), axis = 1)
    avg_df = avg_df.drop(['TAF', 'TDP', 'NAF', 'NDP', 'TDP', 'TAL'], axis=1)

    avg_df['FORMAT'] = avg_df[['GT', 'DPC', 'DPT', 'ADC', 'ADT', 'AL']].apply(lambda x: ':'.join(x[x.notnull()]), axis = 1)
    avg_df = avg_df.drop(['GT', 'DPC', 'DPT', 'ADC', 'ADT', 'AL'], axis=1)

    avg_df[f'{sample_id}'] = avg_df[['geno_GT', 'geno_DPC', 'geno_DPT', 'geno_ADC', 'geno_ADT', 'geno_AL']].apply(lambda x: ':'.join(x[x.notnull()]), axis = 1)
    avg_df = avg_df.drop(['geno_GT', 'geno_DPC', 'geno_DPT', 'geno_ADC', 'geno_ADT', 'geno_AL'], axis=1)

    avg_df.to_csv('tmp.vcf', sep="\t", index=None, header=True)
    #os.system(f'bcftools view -h {sample_files[0]} | sed \'$d\' > header.txt')
    os.system(f'cat pcgr_header.txt tmp.vcf > {sample_id}.vcf')
    os.system(f'bgzip {sample_id}.vcf')
    os.system(f'tabix {sample_id}.vcf.gz')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
