#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process foo {

    echo true

    input:
    tuple val(meta), path(files)

    output:
    stdout emit: output

    script:
    """
    echo $files
    """
}

workflow {

strelka_snv1  = Channel.of([
        [ id:'HCC2400T', tool:'strelka.somatic_snv' ],
        file("HCC2400T.strelka.somatic_snv.vcf.gz"),
        file("HCC2400T.strelka.somatic_snv.vcf.gz.tbi"),
        []
    ])

strelka_indel1 = Channel.of([
        [ id:'HCC2400T', tool:'strelka.somatic_indels' ],
        file("HCC2400T.strelka.somatic_indels.vcf.gz"),
        file("HCC2400T.strelka.somatic_indels.vcf.gz.tbi"),
        []
    ])

freebayes1 = Channel.of([
        [ id:'HCC2400T', tool:'freebayes' ],
        file("HCC2400T.freebayes.vcf.gz"),
        file("HCC2400T.freebayes.vcf.gz.tbi"),
        []
    ])

mutect21 = Channel.of([
        [ id:'HCC2400T', tool:'mutect2' ],
        file("HCC2400T.mutect2.vcf.gz"),
        file("HCC2400T.mutect2.vcf.gz.tbi"),
        []
    ])

strelka_snv  = Channel.of([
        [ id:'HCC1400T', tool:'strelka.somatic_snv' ],
        file("HCC1400T.strelka.somatic_snv.vcf.gz"),
        file("HCC1400T.strelka.somatic_snv.vcf.gz.tbi"),
        []
    ])

strelka_indel = Channel.of([
        [ id:'HCC1400T', tool:'strelka.somatic_indels' ],
        file("HCC1400T.strelka.somatic_indels.vcf.gz"),
        file("HCC1400T.strelka.somatic_indels.vcf.gz.tbi"),
        file("HCC1430T.cns")
    ])

freebayes = Channel.of([
        [ id:'HCC1400T', tool:'freebayes' ],
        file("HCC1400T.freebayes.vcf.gz"),
        file("HCC1400T.freebayes.vcf.gz.tbi"),
        file("HCC1420T.cns")
    ])

mutect2 = Channel.of([
        [ id:'HCC1400T', tool:'mutect2' ],
        file("HCC1400T.mutect2.vcf.gz"),
        file("HCC1400T.mutect2.vcf.gz.tbi"),
        file("HCC1400T.cns")
    ])

//id_caller = strelka_snv.mix(strelka_indel, mutect2, freebayes, strelka_indel1, strelka_snv1, mutect21, freebayes1).map{ it -> id = it[0].id; tool = it[0].tool; return [id, tool] }.groupTuple(by:0).map{ id, callers -> return [ id, callers.size() ] }

files = strelka_snv.mix(strelka_indel, mutect2, freebayes, strelka_indel1, strelka_snv1, mutect21, freebayes1)
//x = files.map{ it -> return it[1..-1] }.map{ it -> meta = [:]; meta.id = it[1].simpleName; return [ meta, it.toList() ] }.groupTuple(by:0).flatten().view()
// have info on how many vcf files we have per sample. if greater than 1, run the merge step..
x = files.map{ it -> return it[1..-1] }.flatten().unique().map{ it -> meta = it.simpleName; return [ meta, it ] }.groupTuple().view()
foo( x )
//foo.out.output

}
