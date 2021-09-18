#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.read_length = 250
params.accessions = "data/course_run_ids.txt"
params.min_contig_len = 1000
params.min_length = 100
params.min_quality = 20
params.trim_front = 5

process download {
    publishDir "${baseDir}/data/sra"
    cpus 4

    input:
    val(id)

    output:
    tuple val(id), path("*.fastq.gz")

    """
    fasterq-dump -t ${task.cpus} ${id}
    gzip *.fastq
    """
}

process preprocess {
    cpus 4
    publishDir "${baseDir}/data/preprocessed"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_filtered_R*.fastq.gz"), path("${id}_fastp.json"), path("${id}.html")

    """
    fastp -i ${reads[0]} -I ${reads[1]} \
        -o ${id}_filtered_R1.fastq.gz -O ${id}_filtered_R2.fastq.gz\
        --json ${id}_fastp.json --html ${id}.html \
        --trim_front1 ${params.trim_front} -l ${params.min_length} \
        -3 -M ${params.min_quality} -r -w ${task.cpus}
    """
}

process multiqc {
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    val(preprocess)

    output:
    path("multiqc_report.html")

    """
    multiqc ${baseDir}/data/preprocessed
    """
}


process megahit {
    cpus 4
    publishDir "${baseDir}/data/assembled", mode: "copy", overwrite: true

    input:
    tuple val(id), path(reads), path(json), path(report)

    output:
    tuple val(id), path("contigs/${id}.asm1.fna")

    """
    megahit -1 ${reads[0]} -2 ${reads[1]} -o contigs -t ${task.cpus} -m 0.5 \
            --min-contig-len ${params.contig_length} --out-prefix ${id}
    mv contigs/${id}.contigs.fa contigs/${id}.asm1.fna
    """
}

workflow {
    ids = Channel.fromPath("${baseDir}/${params.accessions}").splitText()

    ids | download | preprocess | megahit
    multiqc(preprocess.out.collect())
}
