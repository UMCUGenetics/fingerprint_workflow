process UnifiedGenotyper {
    tag "${sample}_gatk_UG"
    publishDir "$params.outdir/fingerprint", mode: 'copy'
    cpus 2
    penv 'threaded'
    memory '5 GB'
    time '1h'

    input:
    set val(sample), file(input_bam), file(input_bai)

    output:
    set val(sample), file("${sample}.vcf")

    script:
    """
    module load Java/1.8.0_60
    java -jar $params.gatk -T UnifiedGenotyper --reference_sequence $params.genome --intervals $params.fingerprint_target --dbsnp $params.fingerprint_target --input_file $input_bam --out ${sample}.vcf
    """

}
