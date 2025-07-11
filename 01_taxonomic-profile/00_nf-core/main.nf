
println "projectDir: $projectDir"

log.info """\
    METAGENOMICS - N F   P I P E L I N E
    ===================================
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)

/*
 * Assembly of the high quality reads
 */
process ASSEMBLY {
    container "quay.io/biocontainers/megahit:1.2.9--h5b5514e_3"
    publishDir "${params.outdir}/assembly", mode:'copy'
    tag "Megahit on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}.contigs.fa"), emit: contigs_id
    path("${sample_id}/${sample_id}.contigs.fa"), emit: contigs
    path("${sample_id}/${sample_id}.log"), emit: contigs_logs

    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} \
    --out-prefix ${sample_id} \
    -o ${sample_id}
    """
}

/*
 * Calling genes on the contigs with Prodigal
 */
 process CALLGENES {
    container "quay.io/biocontainers/prodigal:2.6.3--h516909a_2"
    publishDir "${params.outdir}/prodigalGenes", mode:'copy'
    tag "Calling genes on $sample_id contigs"
    
    input:
    tuple val(sample_id), path(contigs_id)
    
    output:
    tuple val(sample_id), path("${sample_id}.gff"), emit: gene_annotations
    tuple val(sample_id), path("${sample_id}.fna"), emit: nucleotide_fasta
    tuple val(sample_id), path("${sample_id}.faa"), emit: amino_acid_fasta
    tuple val(sample_id), path("${sample_id}_all.txt"), emit: all_gene_annotations

    script:
    """
    prodigal -f gff -d ${sample_id}.fna \
    -o ${sample_id}.gff \
    -a ${sample_id}.faa \
    -s ${sample_id}_all.txt \
    -i ${contigs_id}
    """
}

/*
 * Definning the workflow
 */
workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .toSortedList( { a, b -> a[0] <=> b[0] } )
        .flatMap()
        .set { read_pairs_ch }
        read_pairs_ch.view()

    ASSEMBLY(read_pairs_ch)
    ASSEMBLY.out.contigs_id.view()
    ASSEMBLY.out.contigs.collect().view()
    ASSEMBLY.out.contigs_logs.collect().view()
    
    CALLGENES(ASSEMBLY.out.contigs_id)
    CALLGENES.out.gene_annotations.view()
    CALLGENES.out.nucleotide_fasta.view()
    CALLGENES.out.amino_acid_fasta.view()
    CALLGENES.out.all_gene_annotations.view()
}