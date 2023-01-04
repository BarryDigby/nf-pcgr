#!/usr/bin/env nextflow


mutect2 = Channel.of( [ [ 'id:HCC2400T', 'sample:0', 'patient:200', 'lalala:foo' ], file("HCC2400T.mutect2.vcf.gz"), file("HCC2400T.mutect2.vcf.gz.tbi"), file("HCC2400T.cns") ])

mutect2.map{ meta, f1, f2, f3 -> return [ meta[0..2], f1 ] }.view()




