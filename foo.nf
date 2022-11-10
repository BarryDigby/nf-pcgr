#!/usr/bin/env nextflow


mutect21 = Channel.of([
        'HCC2400T',
        [file("HCC2400T.mutect2.vcf.gz"),
        file("HCC2400T.mutect2.vcf.gz.tbi")],
        file("HCC2400T.cns")
    ])

strelka_snv  = Channel.of([
        'HCC1400T',
        [file("HCC1400T.strelka.somatic_snv.vcf.gz"),
        file("HCC1400T.strelka.somatic_snv.vcf.gz.tbi")],
        file("HCC1400T.cns")
    ])



foo = mutect21.mix(strelka_snv).flatten().collect().view()
