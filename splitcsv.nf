#!/usr/bin/env nextflow

input = "results/constructed_samplesheet.csv"
foo = Channel.from(input)
bar = foo.map{ it ->
         sf = file(it)
         return sf }.view()

println(bar.getClass())

bar.splitCsv(header:true, sep:',').map{ row -> println(row) }
