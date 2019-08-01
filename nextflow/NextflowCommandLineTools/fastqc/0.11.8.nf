process fastqc {
    tag "${sample}_fastqc"
    publishDir "$params.outdir/$sample/fastQC", mode: 'copy'
    cpus 1
    penv 'threaded'
    memory '1 GB'
    time '1h'

    input:
    set val(sample), file(fastq: "*")

    output:
    file '*_fastqc.*'

    script:
    """
    $params.fastqc --noextract -t ${task.cpus} $fastq
    """

}
