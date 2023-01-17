# **nf-pcgr** parameters

## Global options

Mandatory parameters for running both CPSR/PCGR.

- `--genome` Genome assembly used to generate VCF files. Available: [grch37, grch38]
- `--database` Path to PCGR data bundle
- `--fasta` Reference fasta file used in Sarek analysis.

## Input/output options

Define input paths, output directory for results and toggle saving intermediate analysis files.

- `--input` Path to valid samplesheet CSV file. Please refer to documentation for valid samplesheet examples
- `--save_intermediates` Save tabixed, bgzipped and reformatted VCF files from a workflow run [default: false]
- `--outdir`  The output directory where the results will be saved. Must use absolute paths to storage on Cloud infrastructure [default: results]

## PCGR options

Define `PCGR` parameters. To invoke copy number alteration analysis and tumor mutational burden analysis, set `--cna_analysis` and `--tmb_analysis` respectively, to `true`.

- `--assay` Type of DNA sequencing assay performed for input data (VCF) [default: WES]
- `--cpsr_report` CPSR report file (Gzipped JSON - file ending with 'cpsr.<genome_assembly>.json.gz' -  germline report of patient's blood/control sample.
- `--tumor_site` Optional integer code to specify primary tumor type/site of query sample [default: 0]

```console
0 = Any, 1 = Adrenal Gland, 2 = Ampulla of Vater, 3 = Biliary Tract, 4 = Bladder/Urinary Tract, 5 = Bone, 6 = Breast, 7 = Cervix, 8 = CNS/Brain, 9 = Colon/Rectum, 10 = Esophagus/Stomach, 11 = Eye, 12 = Head and Neck, 13 = Kidney, 14 = Liver, 15 = Lung, 16 = Lymphoid, 17 = Myeloid, 18 = Ovary/Fallopian Tube, 19 = Pancreas, 20 = Peripheral Nervous System, 21 = Peritoneum, 22 = Pleura, 23 = Prostate, 24 = Skin, 25 = Soft Tissue, 26 = Testis, 27 = Thymus, 28 = Thyroid, 29 = Uterus, 30 = Vulva/Vagina
```

- `--tumor_purity` Estimated tumor purity (between 0 and 1, [default: None]
- `--tumor_ploidy` Estimated tumor ploidy [default: None]
- `--cna_analysis` Include somatic copy number alteration analysis [default: false]
- `--logr_gain` Log ratio-threshold (minimum) for segments containing copy number gains/amplifications [default: 0.8]
- `--logr_homdel`  Log ratio-threshold (maximum) for segments containing homozygous deletions [default: -0.8]
- `--cna_overlap_pct` Mean percent overlap between copy number segment and gene transcripts for reporting of gains/losses in tumor suppressor genes/oncogenes [default: 50]
- `--target_size_mb` For mutational burden analysis - approximate protein-coding target size in Mb of sequencing assay [default: 34]
- `--estimate_tmb` Estimate tumor mutational burden from the total number of somatic mutations and target region size [default: false]
- `--estimate_msi_status` Predict microsatellite instability status from patterns of somatic mutations/indels, default: false
- `--tmb_algorithm` Method for calculation of TMB, all coding variants (Chalmers et al., Genome Medicine, 2017), or non-synonymous variants only. Available [all_coding, non_synom], [default: all_coding]
- `--estimate_signatures` Estimate relative contributions of reference mutational signatures in query sample and detect potential kataegis events [default: false]
- `--min_mutations_signatures` Minimum number of SNVs required for reconstruction of mutational signatures (SBS) by MutationalPatterns (default: 200, minimum n = 100) [default: 200]
- `--all_reference_signatures` Use all reference mutational signatures (SBS, n = 67) in signature reconstruction rather than only those already attributed to the tumor type [default: false]
- `--include_artefact_signatures` Include sequencing artefacts in the collection of reference signatures [default: false]
- `--prevalence_reference_signatures` Minimum tumor-type prevalence (in percent) of reference signatures to be included in refitting procedure [default: 5]
- `--include_trials` (Beta) Include relevant ongoing or future clinical trials, focusing on studies with molecularly targeted interventions [default: true]
- `--tumor_dp_min` If VCF INFO tag for sequencing depth (tumor) is specified and found, set minimum required depth for inclusion in report [default: 0]
- `--tumor_af_min` If VCF INFO tag for variant allelic fraction (tumor) is specified and found, set minimum required AF for inclusion in report [default: 0]
- `--control_dp_min` If VCF INFO tag for sequencing depth (control) is specified and found, set minimum required depth for inclusion in report [default: 0]
- `--control_af_max` If VCF INFO tag for variant allelic fraction (control) is specified and found, set maximum tolerated AF for inclusion in report [default: 1]

## PCGR (Tumor-only) options

- `--tumor_only` Input VCF comes from tumor-only sequencing, calls will be filtered for variants of germline origin [default: false]
- `--cell_line` Input VCF comes from tumor cell line sequencing (requires --tumor_only), calls will be filtered for variants of germline origin [default: false]
- `--pon_vcf` VCF file with germline calls from Panel of Normals (PON) - i.e. blacklisted variants [default: None]
- `--exclude_pon` Exclude variants occurring in PoN (Panel of Normals), if provided as VCF (--pon_vcf) [default: false]
- `--exclude_likely_hom_germline` Exclude likely homozygous germline variants (100 pct allelic fraction for alternate allele in tumor, very unlikely somatic event, [default: false]
- `--exclude_likely_het_germline` Exclude likely heterozygous germline variants (40-60 pct allelic fraction, AND presence in dbSNP + gnomAD, AND not existing as somatic event in COSMIC/TCGA) [default: false]
- `--exclude_dbsnp_nonsomatic` Exclude variants found in dbSNP (only those that are NOT found in ClinVar(somatic origin)/DoCM/TCGA/COSMIC [defult: false]
- `--exclude_nonexonic` Exclude non-exonic variants [default: false]
- `--maf_onekg_eur` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 GenomesProject - European pop) [default: 0.002]
- `--maf_onekg_amr` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes Project - Ad mixed American pop) [default: 0.002]
- `--maf_onekg_afr` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes Project - African pop) [default: 0.002]
- `--maf_onekg_eas` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes Project -East Asian pop) [default: 0.002]
- `--maf_onekg_sas` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes Project - South Asian pop) [default: 0.002]
- `--maf_onekg_global` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes Project - Global pop) [default: 0.002]
- `--maf_gnomad_nfe` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (gnomAD - European (non-Finnish)) [default: 0.002]
- `--maf_gnomad_asj` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (gnomAD - Ashkenazi Jewish) [default: 0.002]
- `--maf_gnomad_fin` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (gnomAD - European Finnish) [default: 0.002]
- `--maf_gnomad_oth` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (gnomAD - Other) [default: 0.002]
- `--maf_gnomad_amr` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (gnomAD - Latino/Admixed American) [default: 0.002]
- `--maf_gnomad_afr` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (gnomAD - African/African-American) [default: 0.002]
- `--maf_gnomad_eas` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (gnomAD - East Asian) [default: 0.002]
- `--maf_gnomad_sas` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (gnomAD - South Asian) [default: 0.002]
- `--maf_gnomad_global` Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (gnomAD - global population) [default: 0.002]

## CSPR options

-  `--panel_id` Select the virtual gene panel for cancer predisposition reports [default: 0]

```console
0 = CPSR exploratory cancer predisposition panel, 1 = Adult solid tumours cancer susceptibility, 2 = Adult solid tumours for rare disease, 3 = Bladder cancer pertinent cancer susceptibility, 4 = Brain cancer pertinent cancer susceptibility, 5 = Breast cancer pertinent cancer susceptibility, 6 = Childhood solid tumours cancer susceptibility, 7 = Colorectal cancer pertinent cancer susceptibility, 8 = Endometrial cancer pertinent cancer susceptibility, 9 = Familial Tumours Syndromes of the central & peripheral Nervous system, 10 = Familial breast cancer, 11 = Familial melanoma, 12 = Familial prostate cancer, 13 = Familial rhabdomyosarcoma, 14 = GI tract tumours, 15 = Genodermatoses with malignancies, 16 = Haematological malignancies cancer susceptibility, 17 = Haematological malignancies for rare disease, 18 = Head and neck cancer pertinent cancer susceptibility, 19 = Inherited MMR deficiency (Lynch syndrome) - Genomics England PanelApp, 20 = Inherited non-medullary thyroid cancer, 21 = Inherited ovarian cancer (without breast cancer), 22 = Inherited pancreatic cancer, 23 = Inherited polyposis, 24 = Inherited predisposition to acute myeloid leukaemia (AML) - Genomics England PanelApp, 25 = Inherited predisposition to GIST, 26 = Inherited renal cancer, 27 = Inherited phaeochromocytoma and paraganglioma, 28 = Melanoma pertinent cancer susceptibility, 29 = Multiple endocrine tumours, 30 = Multiple monogenic benign skin tumours, 31 = Neuroendocrine cancer pertinent cancer susceptibility, 32 = Neurofibromatosis Type 1, 33 = Ovarian cancer pertinent cancer susceptibility, 34 = Parathyroid Cancer, 35 = Prostate cancer pertinent cancer susceptibility, 36 = Renal cancer pertinent cancer susceptibility, 37 = Rhabdoid tumour predisposition, 38 = Sarcoma cancer susceptibility, 39 = Sarcoma susceptibility, 40 = Thyroid cancer pertinent cancer susceptibility, 41 = Tumour predisposition - childhood onset, 42 = Upper gastrointestinal cancer pertinent cancer susceptibility
```
-  `--diagnostic_grade_only` For panel_id's 1-42 - consider genes with a GREEN status only [default: false]
-  `--ignore_noncoding` Do not list non-coding variants in HTML report [default: false]
-  `--pop_gnomad` Population source in gnomAD used for variant frequency assessment (ACMG classification). Available: [afr,amr,eas,sas,asj,nfe,fin,global] [default: nfe]
-  `--maf_upper_threshold` Upper MAF limit (gnomAD global population frequency) for variants to be included in the report [default: 0.9]
-  `--classify_all` Provide CPSR variant classifications (TIER 1-5) also for variants with exising ClinVar classifications in output TSV [default: false]
-  `--clinvar_ignore_noncancer` Ignore (exclude from report) ClinVar-classified variants reported only for phenotypes/conditions NOT related to cancer [default: false]

## VCF filtering

`PCGR` and `CPSR` only report variants with `PASS` in the filter column. In the event you wish to apply your own filtering logic to variants, please construct a [valid bcftools filtering expression](https://samtools.github.io/bcftools/howtos/filtering.html). Include `-m x` in the expression to resets filters of sites which pass to "PASS".

- `--filter_deepvariant` Apply filtering expression to DeepVariant VCF (germline) using bcftools filter. [default: "-i'FORMAT/DP>10'"]
- `--filter_freebayes_germline` Apply filtering expression to FreeBayes (germline) using bcftools filter. [default: "-i'FORMAT/DP>10'"]
- `--filter_freebayes_somatic` Apply filtering expression to FreeBayes (somatic/tumor-only) using bcftools filter. [default: "-i'FORMAT/DP>10'"]
- `--filter_haplotypecaller` Apply filtering expression to GATK HaplotypeCaller (germline) using bcftools filter [default: "-i'FORMAT/DP>10'"]
- `--filter_mutect2` Apply filtering expression to GATK Mutect2 (somatic/tumor-only) using bcftools filt [default: "-i'FORMAT/DP>10'"]
- `--filter_strelka_indels` Apply filtering expression to Illumina Strelka Indels (somatic) using bcftools filt [default: "-i'FORMAT/DP>10'"]
- `--filter_strelka_snvs` Apply filtering expression to Illumina Strelka snvs (somatic) using bcftools filter [default: "-i'FORMAT/DP>10'"]
- `--filter_strelka_variants` Apply filtering expression to Illumina Strelka variants (germline/tumor-only) using [default: "-i'FORMAT/DP>10'"]

## VEP options

- `--vep_n_forks` Number of forks (option '--fork' in VEP) [default: 4]
- `--vep_buffer_size` Variant buffer size (variants read into memory simultaneously, option '--buffer_size' in VEP). Set to lower to reduce memory usage [default: 500]
- `--vep_gencode_all` Consider all GENCODE transcripts with Variant Effect Predictor (VEP) (option '--gencode_basic' in VEP is used by default). [default: true]
- `--vep_pick_order` Comma-separated string of ordered transcript properties for primary variant pick ( option'--pick_order' in VEP) [default: canonical,appris,biotype,ccds,rank,tsl,length,mane]
- `--vep_no_intergenic` Skip intergenic variants during processing (option '--no_intergenic' in VEP), [default: false]
