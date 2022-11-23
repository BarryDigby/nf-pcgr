#!/usr/bin/env python

from pysam import VariantFile
import os
import time

# Inspired by: @gudeqing
# https://github.com/sigven/pcgr/issues/136#issuecomment-919273152
# Built upon by @BarryDigby to handle Strelka, Freebayes, Mutect2 VCF files.

def mutect2_vaf(record, sample_idx):
    VAF = record.samples[sample_idx]['AF'][0]
    VAF = VAF if VAF is not None else 0
    return VAF


def freebayes_vaf(record, sample_idx):
    AO = record.samples[sample_idx]['AO'][0]
    RO = record.samples[sample_idx]['RO']
    AO = AO if AO is not None else 0
    RO = RO if RO is not None else 0
    if (AO+RO) == 0:
        VAF = 0
    else:
        VAF = (AO/(AO+RO))
    return VAF

def strelka_snv_vaf(record, sample_idx):
    ref = str(record.ref + "U")
    alt = str(record.alts[0] + "U")
    tier1RefCounts = record.samples[sample_idx][ref][0]
    tier1AltCounts = record.samples[sample_idx][alt][0]
    tier1RefCounts = tier1RefCounts if tier1RefCounts is not None else 0
    tier1AltCounts = tier1AltCounts if tier1AltCounts is not None else 0
    if (tier1AltCounts+tier1RefCounts) == 0:
        VAF = 0
    else:
        VAF = (tier1AltCounts/(tier1AltCounts+tier1RefCounts))
    record.samples[sample_idx]['AD'] = [0, tier1AltCounts]
    return VAF

def strelka_indel_vaf(record, sample_idx):
    tier1RefCounts = record.samples[sample_idx]['TAR'][0]
    tier1AltCounts = record.samples[sample_idx]['TIR'][0]
    tier1RefCounts = tier1RefCounts if tier1RefCounts is not None else 0
    tier1AltCounts = tier1AltCounts if tier1AltCounts is not None else 0
    if (tier1AltCounts+tier1RefCounts) == 0:
        VAF = 0
    else:
        VAF = (tier1AltCounts/(tier1AltCounts+tier1RefCounts))
    record.samples[sample_idx]['AD'] = [0, tier1AltCounts]
    return VAF

vcf_formats = { "mutect2_vaf": ['AD', 'AF', 'DP', 'F1R2', 'F2R1', 'FAD', 'GQ', 'GT', 'PGT', 'PID', 'PL', 'PS', 'SB'],
                "freebayes_vaf": ['AD', 'AO', 'DP', 'GL', 'GQ', 'GT', 'MIN_DP', 'PL', 'QA', 'QR', 'RO'],
                "strelka_snv_vaf": ['AU', 'CU', 'DP', 'FDP', 'GU', 'SDP', 'SUBDP', 'TU'],
                "strelka_indel_vaf": ['BCN50', 'DP', 'DP2', 'DP50', 'FDP50', 'SUBDP50', 'TAR', 'TIR', 'TOR'] }


#####################################################################################################################
#####################################################################################################################

def reformat_vcf(vcf_file, out, reference):
    """
    BCFtools: split multi-allelic sites, normalize and remove sites where DP is missing or < 1. (divide by zero errors).
    Calculations need to be cross checked with someone with more experience in population genomics.
    BCFtools: reformatting sample names to NORMAL, TUMOR for downstream merging.
    """
    os.system(f'bcftools norm -f {reference} -m -both {vcf_file} | bcftools filter -e\'FORMAT/DP="." || FORMAT/DP<1\' -o tmp_.vcf')
    with VariantFile('tmp_.vcf') as fr:
        header = fr.header
        header.info.add('TDP', number=1, type='Integer', description='Tumor sample depth')
        header.info.add('NDP', number=1, type='Integer', description='Normal sample depth')
        header.info.add('TAF', number=1, type='Float', description='Tumor sample AF')
        header.info.add('NAF', number=1, type='Float', description='Normal sample AF')
        header.info.add('TAL', number='.', type='String', description='Algorithms that called the somatic mutation')
        samples = list(header.samples)
        formats = list(header.formats)
        fnc_str = list(vcf_formats.keys())[list(vcf_formats.values()).index(list(formats))]
        # append after algorithm-function mapping
        header.formats.add('AL', number='.', type='Integer', description='Codes for algorithms that produced the somatic call (1 = mutect2, 2 = freebayes, 3 = strelka)')
        if 'strelka' in fnc_str:
            header.formats.add('AD', number=2, type="Integer", description='AD flag for Strelka. Output as tuple so index rule for TAF does not need to be modified.')
        with VariantFile('tmp_1.vcf', 'w', header=header) as fw:
            tumor_is_first = 0
            tumor_is_second = 0
            algorithm = fnc_str.split('_', 1)[0]
            algorithm_code = 1 if algorithm == 'mutect2' else 2 if algorithm == 'freebayes' else 3
            for record in fr:
                VAF_sample0 = globals()[fnc_str](record, 0)
                VAF_sample1 = globals()[fnc_str](record, 1)
                if VAF_sample0 > VAF_sample1:
                    tumor_is_first += 1
                else:
                    tumor_is_second += 1
                tumor_idx = tumor_is_second >= tumor_is_first
                normal_idx = 1 - tumor_idx
                record.info['TDP'] = record.samples[tumor_idx]['DP']
                record.info['NDP'] = record.samples[normal_idx]['DP']
                record.info['TAF'] = round(record.samples[tumor_idx]['AD'][1]/record.samples[tumor_idx]['DP'], 3)
                record.info['NAF'] = round(record.samples[normal_idx]['AD'][1]/record.samples[normal_idx]['DP'], 3)
                record.info['TAL'] = algorithm
                record.samples[tumor_idx]['AL'] = algorithm_code
                record.samples[normal_idx]['AL'] = algorithm_code
                fw.write(record)
        ## write file for bcftools reheader.
        ## one per line, appearing in order of samples in VCF file
        normal = f'{samples[normal_idx]} NORMAL'
        tumor = f'{samples[tumor_idx]} TUMOR'
        ## order matters:
        if normal_idx == 0:
            with open("bcftools_reheader.txt", "w") as f:
                f.write(f'{normal}\n{tumor}')
        else:
            with open("bcftools_reheader.txt", "w") as f:
                f.write(f'{tumor}\n{normal}')

    print(f'we guess tumor sample is {samples[tumor_idx]} ')
    os.system(f'bcftools reheader -s bcftools_reheader.txt tmp_1.vcf > {out}')
    os.remove('tmp_.vcf')
    os.remove('tmp_1.vcf')
    os.system(f'bgzip {out}')
    os.system(f'tabix {out}.gz')

if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
