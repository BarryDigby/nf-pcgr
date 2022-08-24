from pysam import VariantFile
import os

# AUTHOR: @gudeqing
# https://github.com/sigven/pcgr/issues/136#issuecomment-919273152

def reformat_vcf(vcf_file, out, reference, tumor_sample=None):
    """
    [split and left-normalization] -> [re-calculate AF] -> [bgzip and tabix]
    :param vcf_file:
    :param out:
    :param reference:
    :param tumor_sample: tumor sample name or position, if the second sample is tumor, then it is 1 else 0
    :return:
    """
    os.system(f'bcftools norm -f {reference} -m -both {vcf_file} -o tmp_.vcf')
    with VariantFile('tmp_.vcf') as fr:
        header = fr.header
        header.info.add('TDP', number=1, type='Integer', description='Tumor sample depth')
        header.info.add('NDP', number=1, type='Integer', description='Normal sample depth')
        header.info.add('TAF', number=1, type='Float', description='Tumor sample AF')
        header.info.add('NAF', number=1, type='Float', description='Normal sample AF')
        samples = list(header.samples)

        if tumor_sample is not None:
            if tumor_sample not in [0, 1, '0', '1']:
                if tumor_sample in samples:
                    tumor_idx = samples.index(tumor_sample)
                    normal_idx = 1 - tumor_idx
                else:
                    raise Exception(f'{tumor_sample} is not in samples {samples} recorded in vcf')
            else:
                tumor_idx = int(tumor_sample)
                normal_idx = 1 - tumor_idx
        else:
            tumor_idx = guess_tumor_idx(vcf_file)
            normal_idx = 1 - tumor_idx

        with VariantFile(out, 'w', header=header) as fw:
            for record in fr:
                record.info['TDP'] = record.samples[tumor_idx]['DP']
                record.info['NDP'] = record.samples[normal_idx]['DP']
                # re-calculate AF since caller like sentieon may report AF that is not consistent with AD info
                record.info['TAF'] = round(record.samples[tumor_idx]['AD'][1]/record.samples[tumor_idx]['DP'], 3)
                record.info['NAF'] = round(record.samples[normal_idx]['AD'][1]/record.samples[normal_idx]['DP'], 3)
                fw.write(record)

    os.remove('tmp_.vcf')
    os.system(f'bgzip {out}')
    os.system(f'tabix {out}.gz')


def guess_tumor_idx(vcf_file):
    tumor_is_first = 0
    tumor_is_second = 0

    with VariantFile(vcf_file) as fr:
        samples = list(fr.header.samples)
        formats = list(fr.header.formats)
        if 'AF' not in formats:
            raise Exception('No AF in format info to detect tumor sample')
        for record in fr:
            if record.samples[0]['AF'][0] > record.samples[1]['AF'][0]:
                tumor_is_first += 1
            else:
                tumor_is_second += 1
    tumor_idx = tumor_is_second >= tumor_is_first
    print(f'we guess tumor sample is {samples[tumor_idx]} ')
    return tumor_idx


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
