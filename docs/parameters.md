# **nf-pcgr** parameters

## Global options

Mandatory parameters for running both CPSR/PCGR.

- `--genome` Genome assembly used to generate VCF files. Available: [grch37, grch38]
- `--database` Path to PCGR data bundle

## Input/output options

- `--input` Path to valid samplesheet CSV file. Please refer to documentation for valid samplesheet examples
- `--save_intermediates` Save tabixed, bgzipped and reformatted VCF files from a workflow run [default: false]
- `--outdir`  The output directory where the results will be saved. Must use absolute paths to storage on Cloud infrastructure [default:
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; results]

PCGR options
- `--assay` Type of DNA sequencing assay performed for input data (VCF) [default: WES]
- `--cpsr_report` CPSR report file (Gzipped JSON - file ending with 'cpsr.<genome_assembly>.json.gz' -  germline report of patient's blood/control sample.
- `--fasta` Reference fasta file used in Sarek analysis.
- `--tumor_site`                      [number]  Optional integer code to specify primary tumor type/site of query sample [default: 0]
  `--tumor_purity`                    [string]  Estimated tumor purity (between 0 and 1, (default: None)
  `--tumor_ploidy`                    [string]  Estimated tumor ploidy (default: None)
  `--cna_analysis`                    [boolean] Include somatic copy number alteration analysis.
  `--logr_gain`                       [number]  Log ratio-threshold (minimum) for segments containing copy number gains/amplifications (default: 0.8)
                                              [default: 0.8]
  `--logr_homdel`                     [number]  Log ratio-threshold (maximum) for segments containing homozygous deletions (default: -0.8) [default: -0.8]
  `--cna_overlap_pct`                 [number]  Mean percent overlap between copy number segment and gene transcripts for reporting of gains/losses in tumor
                                              suppressor genes/oncogenes, (default: 50) [default: 50]
  `--target_size_mb`                  [number]  For mutational burden analysis - approximate protein-coding target size in Mb of sequencing assay (default: 34
                                              WGS/WES) [default: 34]
  `--estimate_tmb`                    [boolean] Estimate tumor mutational burden from the total number of somatic mutations and target region size, default:
                                              False
  `--estimate_msi_status`             [boolean] Predict microsatellite instability status from patterns of somatic mutations/indels, default: False
  `--tmb_algorithm`                   [string]  Method for calculation of TMB, all coding variants (Chalmers et al., Genome Medicine, 2017), or non-synonymous
                                              variants only, default: all_coding [default: all_coding]
  `--estimate_signatures`             [boolean] Estimate relative contributions of reference mutational signatures in query sample and detect potential kataegis
                                              events, default: False
  `--min_mutations_signatures`        [number]  Minimum number of SNVs required for reconstruction of mutational signatures (SBS) by MutationalPatterns (default:
                                              200, minimum n = 100) [default: 200]
  `--all_reference_signatures`        [boolean] Use all reference mutational signatures (SBS, n = 67) in signature reconstruction rather than only those already
                                              attributed to the tumor type (default: False)
  `--include_artefact_signatures`     [boolean] Include sequencing artefacts in the collection of reference signatures (default: False)
  `--prevalence_reference_signatures` [number]  Minimum tumor-type prevalence (in percent) of reference signatures to be included in refitting procedure (default:
                                              5) [default: 5]
  `--include_trials`                  [boolean] Include relevant ongoing or future clinical trials, focusing on studies with molecularly targeted
                                              interventions
  `--tumor_dp_tag`                    [string]  Specify VCF INFO tag for sequencing depth (tumor, must be Type=Integer, default: TDP) [default: TDP]
  `--tumor_af_tag`                    [string]  Specify VCF INFO tag for variant allelic fraction (tumor,  must be Type=Float, default: TAF) [default:
                                              TAF]
  `--control_dp_tag`                  [string]  Specify VCF INFO tag for sequencing depth (control, must be Type=Integer, default: NDP) [default: NDP]
  `--control_af_tag`                  [string]  Specify VCF INFO tag for variant allelic fraction (control, must be Type=Float, default: NAF) [default:
                                              NAF]
  `--call_conf_tag`                   [string]  Specify VCF INFO tag for somatic variant call confidence (must be categorical, e.g. Type=String, default:
                                              _NA_) [default: TAL]
  `--tumor_dp_min`                    [number]  If VCF INFO tag for sequencing depth (tumor) is specified and found, set minimum required depth for inclusion in
                                              report (default: 0) [default: 0]
  `--tumor_af_min`                    [number]  If VCF INFO tag for variant allelic fraction (tumor) is specified and found, set minimum required AF for inclusion
                                              in report (default: 0) [default: 0]
  `--control_dp_min`                  [number]  If VCF INFO tag for sequencing depth (control) is specified and found, set minimum required depth for inclusion in
                                              report (default: 0) [default: 0]
  `--control_af_max`                  [number]  If VCF INFO tag for variant allelic fraction (control) is specified and found, set maximum tolerated AF for
                                              inclusion in report (default: 1) [default: 1]

PCGR (Tumor-only) options
  --tumor_only                      [boolean] Input VCF comes from tumor-only sequencing, calls will be filtered for variants of germline origin, (default:
                                              False) [default: false]
  --cell_line                       [boolean] Input VCF comes from tumor cell line sequencing (requires --tumor_only), calls will be filtered for variants of
                                              germline origin, (default: False)
  --pon_vcf                         [string]  VCF file with germline calls from Panel of Normals (PON) - i.e. blacklisted variants, (default: None)
  --exclude_pon                     [boolean] Exclude variants occurring in PoN (Panel of Normals, if provided as VCF (--pon_vcf), default: False) [default:
                                              false]
  --exclude_likely_hom_germline     [boolean] Exclude likely homozygous germline variants (100 pct allelic fraction for alternate allele in tumor, very unlikely
                                              somatic event, default: False)
  --exclude_likely_het_germline     [boolean] Exclude likely heterozygous germline variants (40-60 pct allelic fraction, AND presence in dbSNP + gnomAD, AND not
                                              existing as somatic event in COSMIC/TCGA, default: False)
  --exclude_dbsnp_nonsomatic        [boolean] Exclude variants found in dbSNP (only those that are NOT found in ClinVar(somatic origin)/DoCM/TCGA/COSMIC,
                                              defult: False)
  --exclude_nonexonic               [boolean] Exclude non-exonic variants, (default: False)
  --maf_onekg_eur                   [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes
                                              Project - European pop, default: 0.002) [default: 0.002]
  --maf_onekg_amr                   [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes
                                              Project - Ad mixed American pop, default: 0.002) [default: 0.002]
  --maf_onekg_afr                   [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes
                                              Project - African pop, default: 0.002) [default: 0.002]
  --maf_onekg_eas                   [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes
                                              Project -East Asian pop, default: 0.002) [default: 0.002]
  --maf_onekg_sas                   [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes
                                              Project - South Asian pop, default: 0.002) [default: 0.002]
  --maf_onekg_global                [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes
                                              Project - Global pop, default: 0.002) [default: 0.002]
  --maf_gnomad_nfe                  [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD -
                                              European (non-Finnish), default: 0.002) [default: 0.002]
  --maf_gnomad_asj                  [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD -
                                              Ashkenazi Jewish, default: 0.002) [default: 0.002]
  --maf_gnomad_fin                  [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD -
                                              European (Finnish), default: 0.002) [default: 0.002]
  --maf_gnomad_oth                  [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD -
                                              Other, default: 0.002) [default: 0.002]
  --maf_gnomad_amr                  [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD -
                                              Latino/Admixed American, default: 0.002) [default: 0.002]
  --maf_gnomad_afr                  [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD -
                                              African/African-American, default: 0.002) [default: 0.002]
  --maf_gnomad_eas                  [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD -
                                              East Asian, default: 0.002) [default: 0.002]
  --maf_gnomad_sas                  [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD -
                                              South Asian, default: 0.002) [default: 0.002]
  --maf_gnomad_global               [number]  Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD -
                                              global population, default: 0.002) [default: 0.002]

CSPR options
  --panel_id                        [number]  Select the virtual gene panel for cancer predisposition reports. [default: 0]
  --diagnostic_grade_only           [boolean] For panel_id's 1-42 (Genomics England PanelApp) - consider genes with a GREEN status only, default: False
  --ignore_noncoding                [boolean] Do not list non-coding variants in HTML report, default: False
  --pop_gnomad                      [string]  Population source in gnomAD used for variant frequency assessment (ACMG classification), default: nfe
                                              [default: nfe]
  --maf_upper_threshold             [number]  Upper MAF limit (gnomAD global population frequency) for variants to be included in the report, default: 0.9
                                              [default: 0.9]
  --classify_all                    [boolean] Provide CPSR variant classifications (TIER 1-5) also for variants with exising ClinVar classifications in output
                                              TSV, default: False
  --clinvar_ignore_noncancer        [boolean] Ignore (exclude from report) ClinVar-classified variants reported only for phenotypes/conditions NOT related to
                                              cancer, default: False

VCF filtering
  --filter_deepvariant              [string]  Apply filtering expression to DeepVariant VCF (germline) using bcftools filter. [default:
                                              "-i'FORMAT/DP>10'"]
  --filter_freebayes_germline       [string]  Apply filtering expression to FreeBayes (germline) using bcftools filter.
  --filter_freebayes_somatic        [string]  Apply filtering expression to FreeBayes (somatic/tumor-only) using bcftools filter.
  --filter_haplotypecaller          [string]  Apply filtering expression to GATK HaplotypeCaller (germline) using bcftools filter.
  --filter_mutect2                  [string]  Apply filtering expression to GATK Mutect2 (somatic/tumor-only) using bcftools filter.
  --filter_strelka_indels           [string]  Apply filtering expression to Illumina Strelka Indels (somatic) using bcftools filter.
  --filter_strelka_snvs             [string]  Apply filtering expression to Illumina Strelka snvs (somatic) using bcftools filter.
  --filter_strelka_variants         [string]  Apply filtering expression to Illumina Strelka variants (germline/tumor-only) using bcftools filter.

VEP options
  --vep_n_forks                     [number]  Number of forks (option '--fork' in VEP), default: 4 [default: 4]
  --vep_buffer_size                 [number]  Variant buffer size (variants read into memory simultaneously, option '--buffer_size' in VEP). Set to lower to
                                              reduce memory usage, default: 500 [default: 500]
  --vep_gencode_all                 [boolean] Consider all GENCODE transcripts with Variant Effect Predictor (VEP) (option '--gencode_basic' in VEP is used by
                                              default). [default: true]
  --vep_pick_order                  [string]  Comma-separated string of ordered transcript properties for primary variant pick                     ( option
                                              '--pick_order' in VEP), default: canonical,appris,biotype,ccds,rank,tsl,length,mane [default:
                                              canonical,appris,biotype,ccds,rank,tsl,length,mane]
  --vep_no_intergenic               [boolean] Skip intergenic variants during processing (option '--no_intergenic' in VEP), default: False

Generic options
  --multiqc_methods_description     [string]  Custom MultiQC yaml file containing HTML including a methods description.
