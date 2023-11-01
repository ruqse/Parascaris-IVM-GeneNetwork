#!/usr/bin/env nextflow

nextflow.enable.dsl=2




/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data2/*_{R1,R2}*.fastq.gz"
params.sortmerna_ref = "$projectDir/data/smr_v4.3_fast_db.fasta"
params.adapters = "$projectDir/data/TruSeq3-PE-2.fa"
params.decoys = "$projectDir/data/decoys.txt"
params.gentrome = "$projectDir/data/gentrome.fa.gz"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"
params.sortmernaoutdir = "sortemerna"
params.trimmooutdir = "trimmomatic"
params.salmonoutdir = "salmon"
params.fastqcdir = "fastqc"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome}
    reads        : ${params.reads}
    genome       : ${params.genome}
    trim_adapters: ${params.adapters}
    outdir       : ${params.outdir}
    """
    .stripIndent()


/*
 * Remove remenant ribosomal RNA
 * given smr_v4.3_fast_db.fasta 
 */

process SORTMERNA {
    tag "SortmeRNA on $sample_id"

    publishDir(
        path: "${params.outdir}/${params.sortmernaoutdir}", mode:'copy'
    )
    input:
    path(sortmerna_ref)      
    tuple val(sample_id), path(read_pairs_ch)
   

    output:
    tuple val(sample_id), path("*.fastq.gz")  , emit: sortmerna_reads
    tuple val(sample_id), path("*.log")                   , emit: log

    script:
    """

    sortmerna --ref $sortmerna_ref --reads ${read_pairs_ch[0]} --reads ${read_pairs_ch[1]} --workdir . --aligned rRNA_reads --fastx --other non_rRNA_reads --paired_in --threads $task.cpus --out2
    mv non_rRNA_reads_fwd.fq.gz ${sample_id}_1.fastq.gz
    mv non_rRNA_reads_rev.fq.gz ${sample_id}_2.fastq.gz
    mv rRNA_reads.log ${sample_id}.sortmerna.log
    """

}


/*
 * Trim poor quality read ends including adapter sequences with TRIMMOMATIC
 * given adapters sequences templates in $adapters
 */

process TRIMMOMATIC {
    tag "Trimmomatic on $sample_id"

    publishDir(
        path: "${params.outdir}/${params.trimmooutdir}", mode:'copy'
    )
    input:
    path(adapters)      
    tuple val(sample_id), path(sortmerna_reads)
   

    output:
    tuple val(sample_id), path("*_paired.trim.fastq.gz")  , emit: trimmed_reads
    tuple val(sample_id), path("*_unpaired.trim.fastq.gz"), emit: unpaired_reads
    tuple val(sample_id), path("*.log")                   , emit: log

    script:
    """
    trimmomatic PE -threads $task.cpus -summary ${sample_id}.log \
    ${sortmerna_reads[0]} ${sortmerna_reads[1]} \
    ${sample_id}_1_paired.trim.fastq.gz ${sample_id}_1_unpaired.trim.fastq.gz \
    ${sample_id}_2_paired.trim.fastq.gz ${sample_id}_2_unpaired.trim.fastq.gz \
    ILLUMINACLIP:$adapters:2:30:10:2:True SLIDINGWINDOW:5:20 MINLEN:36
    """

}


/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process SALMON_INDEX {

    input:
    path gentrome
    path decoys
    output:
    path 'salmon_index'

    script:
    """
    salmon index -t $gentrome -d $decoys -p $task.cpus --keepDuplicates -i salmon_index
    
    """
}

process SALMON_QUANTIFICATION {
   tag "Salmon on $sample_id"

    publishDir(
        path:"${params.outdir}/${params.salmonoutdir}", mode:'copy'
    )

    input:
    path salmon_index
    tuple val(sample_id), path(trimmed_reads)

    output:
    path "$sample_id"

    script:
    """
    salmon quant --threads 20 --libType=A --gcBias -i $salmon_index -1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]} -o $sample_id
    """
}


process FASTQC {
    tag "FASTQC on $sample_id"

    publishDir(
        path: "${params.outdir}/${params.fastqcdir}", mode:'copy'
    )
    input:
    tuple val(sample_id), path(reads)
    tuple val(sample_id), path(sortmerna_reads)
    tuple val(sample_id), path(trimmed_reads)
    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads} ${sortmerna_reads} ${trimmed_reads}
    """
}

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}


workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    sortmerna_ch = SORTMERNA(params.sortmerna_ref, read_pairs_ch)  
    trim_ch = TRIMMOMATIC(params.adapters, SORTMERNA.out.sortmerna_reads)
    index_ch = SALMON_INDEX(params.gentrome, params.decoys)
    quant_ch = SALMON_QUANTIFICATION(index_ch, TRIMMOMATIC.out.trimmed_reads)
    fastqc_ch = FASTQC(read_pairs_ch, SORTMERNA.out.sortmerna_reads, TRIMMOMATIC.out.trimmed_reads)
    MULTIQC(quant_ch.mix(fastqc_ch).collect())

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
