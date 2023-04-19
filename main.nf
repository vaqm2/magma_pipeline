#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def help_message() {
    log.info """
    Nextflow analysis pipeline to perform gene level analysis on sumstats file
    
    AUTHOR: Vivek Appadurai | vivek.appadurai@regionh.dk

    USAGE: nextflow run main.nf

    OPTIONS:

    --window [35,10] <Size of the annotation window around a gene>
    --genelist [gene.loc] <Gene locations file on the same build as snp-loc file>
    --out <output_prefix>
    --bfile [g1000_eur] PLINK file of genotypes to use as linkage disequilibrium reference
    --sumstats <sumstats file>
    --ncol <Column that specifies sample size in reference sumstats file>
    --genesets <Gene set annotations file>
    --help <prints this message log>
    """
}

if(params.help) {
    help_message()
    exit 0
}

log.info """
=================================================================
M A G M A - G E N E _ A N A L Y S I S _ P I P E L I N E
=================================================================
Reference GWAS:    $params.sumstats
SNP List:          $params.snplist
Gene List:         $params.genelist
LD Reference:      $params.bfile
N column:          $params.ncol
Gene Sets:         $params.genesets 
Annotation Window: $params.window
Output Prefix:     $params.out
=================================================================
"""

process create_snp_list {
    label 'low_mem'

    input: tuple path(sumstats),
    val(output_prefix)

    output: path(${output_prefix}.snps.list)
    
    script: """
    echo -e "SNP CHR BP" > ${output_prefix}.snps.list
    awk '{print \$1" "\$2" "\$3}' $sumstats >> ${output_prefix}.snps.list
    """
}

process annotate_snps {
    label 'low_mem'

    input: tuple path(snplist),
    path(genelist),
    val(annotation_window),
    val(output_prefix),
    path(magma_exe)

    output: path(${output_prefix}.genes.annot)

    script:
    """
    ./magma_exe --annotate $annotation_window \
        --snp-loc $snplist \
        --gene-loc $genelist \
        --output $output_prefix
    """
}

process gene_analysis {
    label 'mod_mem'
    publishDir launchDir

    input: tuple path(annotation_output),
    val(ld_reference_bfile),
    path(ld_reference_bem),
    path(ld_reference_bim),
    path(ld_reference_fam),
    path(sumstats),
    val(n_col_name),
    val(output_prefix),
    path(magma_exe)

    output: tuple path(${out_prefix}.genes.out), 
    path(${out_prefix}.genes.raw)

    script:
    """
    ./magma_exe --bfile $ld_reference_bfile \
        --gene-annot $annotation_output \
        --pval $sumstats ncol=$n_col_name \
        --out $output_prefix
    """
}

process gene_set_analysis {
    label 'mod_mem'
    publishDir launchDir

    input: tuple path(gene_analysis_output),
    path(gene_analysis_output_raw),
    path(gene_sets_file),
    val(output_prefix),
    path(magma_exe)

    output: path(${output_prefix}.sets.genes.out)

    script:
    """
    ./magma_exe --gene-results $gene_analysis_raw \
        --set-annot $gene_sets_file \
        --out $output_prefix
    """
}

workflow {
    Channel.fromPath(params.sumstats) \
    | combine(Channel.of(params.out)) \
    | create_snp_list \
    | combine(Channel.fromPath(params.genelist)) \
    | combine(Channel.of(params.window)) \
    | combine(Channel.of(params.out)) \
    | combine(Channel.fromPath(params.magma)) \
    | annotate_snps \
    | combine(Channel.of(params.bfile)) \
    | combine(Channel.fromPath('${params.bfile}.*'))
    | combine(Channel.fromPath(params.sumstats)) \
    | combine(Channel.of(params.ncol)) \
    | combine(Channel.of(params.out)) \
    | combine(Channel.fromPath(params.magma)) \
    | gene_analysis \
    | combine(Channel.fromPath(params.genesets)) \
    | combine(Channel.of(params.out)) \
    | combine(Channel.fromPath(params.magma)) \
    | gene_set_analysis
}