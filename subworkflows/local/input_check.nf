//
// Check input and get vcf channels
// Sarek outputs VCF files in the form of:
// ${meta.id}.{tool}.vcf.gz..
// 'simpleName' should suffice here.
// Allow file or path as input, automatically check if VCF files are indexed
//

workflow INPUT_CHECK {
    take:
    samplesheet // samplesheet file or path to VCF files

    main:
    check_input(samplesheet)

    emit:
    vcf  // channel: [ val(meta), [ vcf.gz] , [ vcf.gz.tbi ] ]
}

def check_input(input){
    if(input.toString().endsWith('.csv')) {

        Channel.from(input)
            .splitCsv(header:false, sep:',')
            .map{ row ->

                fh = file(row[0])

                if(!file(fh).exists()){
                    log.error("ERROR: Please check input input -> VCF file does not exist: ${fh}")
                    System.exit(1)
                }else{
                    def meta = [:]
                    fh       = file(row[0])
                    meta.id  = fh.simpleName
                    tbi      = fh.toString() + '.tbi'
                    if(!file(tbi).exists()){
                        log.error("ERROR: Please make sure the VCF files are indexed with tabix. Offending file: ${tbi}")
                        System.exit(1)
                    }
                    return  [ meta, [ file(fh) ], [ file(tbi) ] ]
                }
            }
            .set { vcf }
    }else{

        Channel.fromPath( "${input}/*.vcf.gz" )
            .map{ it ->

                if(!file(it).exists()){
                    log.error("ERROR: VCF file does not exist at path: ${it}")
                    System.exit(1)
                }else{
                    def meta = [:]
                    meta.id  = file(it).simpleName
                    tbi      = it.toString() + '.tbi'
                    if(!file(tbi).exists()){
                        log.error("ERROR: Please make sure the VCF files are indexed with tabix. Offending file: ${tbi}")
                        System.exit(1)
                    }
                    return [ meta, [ file(it) ], [ file(tbi) ] ]
                }
            }
            .set { vcf }
    }
}
