#!/usr/bin/env nextflow

params.samplesheet
params.outdir

Channel.fromPath( file(params.samplesheet) )
    .splitCsv(header: true, sep: '\t')
    .map{row ->
        def sample = row['Sample']
        def reads1 = row['R1'].tokenize( ',' ).collect { file(it) }
        def reads2 = row['R2'].tokenize( ',' ).collect { file(it) }
        return [ sample, reads1, reads2 ]
    }
    .tap{samples_R1_R2_fastq} // set of all fastq R1 R2 per sample
    .map { sample_ID, reads1, reads2 ->
        return [ sample_ID, [reads1, reads2 ].flatten() ]
    }
    .tap{samples_all_fastq} // set of all fastq per sample

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

process fastqc {
    tag "${sample}_fastqc"
    publishDir "$params.outdir/$sample/fastQC", mode: 'copy'
    cpus 1
    penv 'threaded'
    memory '1 GB'
    time '1h'

    input:
    set val(sample), file(fastq: "*") from samples_all_fastq

    output:
    file '*_fastqc.*' into fastqc_out

    script:
    """
    $params.fastqc --noextract -t ${task.cpus} $fastq
    """

}

process bwa_mem {
    tag "${sample}_bwa_mem"
    publishDir "$params.outdir/$sample/mapping", mode: 'copy'
    cpus 1
    penv 'threaded'
    memory '10 GB'
    time '1h'

    input:
    set val(sample), file(r1_fastq), file(r2_fastq) from mips_trim_dedup_out

    output:
    set val(sample), file("${sample}.bam"), file("${sample}.bai") into bwa_mem_out

    script:
    def barcode = r1_fastq.getName().split('_')[1]
    def bwa_readgroup = "\"@RG\\tID:${sample}_${barcode}\\tSM:${sample}\\tPL:ILLUMINA\\tLB:${sample}\\tPU:${barcode}\""

    """
    $params.bwa mem -t ${task.cpus} -c 100 -M -R $bwa_readgroup $params.genome $r1_fastq $r2_fastq | $params.samtools view -b | $params.samtools sort > ${sample}.bam
    $params.samtools index ${sample}.bam ${sample}.bai
    """
}

process gatk_UnifiedGenotyper {
    tag "${sample}_gatk_UG"
    publishDir "$params.outdir/fingerprint", mode: 'copy'
    cpus 2
    penv 'threaded'
    memory '5 GB'
    time '1h'

    input:
    set val(sample), file(input_bam), file(input_bai) from bwa_mem_out

    output:
    set val(sample), file("${sample}.vcf") into gatk_ug_out

    script:
    """
    module load Java/1.8.0_60
    java -jar $params.gatk -T UnifiedGenotyper --reference_sequence $params.genome --intervals $params.fingerprint_target --dbsnp $params.fingerprint_target --input_file $input_bam --out ${sample}.vcf
    """

}
