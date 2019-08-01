process mips_trim_dedup {
    tag "${sample}_mips_trim_dedup"
    publishDir "$params.outdir/$sample/mapping", mode: 'copy'
    cpus 1
    penv 'threaded'
    memory '1 GB'
    time '1h'

    input:
    set val(sample), file(r1_fastqs: "*"), file(r2_fastqs: "*") from samples_R1_R2_fastq

    output:
    set val(sample), file('*_LMergedTrimmedDedup_R1_*.fastq.gz'), file('*_LMergedTrimmedDedup_R2_*.fastq.gz') into mips_trim_dedup_out

    script:
    def r1_args = r1_fastqs.collect{ "$it" }.join(" ")
    def r2_args = r2_fastqs.collect{ "$it" }.join(" ")

    """
    python $params.mips_trim_dedup -d $params.mip_design_file -r1 $r1_args -r2 $r2_args -l $params.mip_uuid_length -ur $params.mip_uuid_read > output.log 2> output.err
    """
}
