#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include mips_trim_dedup from 'NextflowCommandLineTools/mips/1.0.1.nf' params(
    params.mips_trim_dedup: params.mips_trim_dedup,
    params.design_file: params.mips_design_file,
    params.uuid_length: params.mips_uuid_length,
    params.uuid_read: params.mips_uuid_read,
    )
include fastqc from 'NextflowCommandLineTools/fastqc/0.11.8.nf' params(params)
include mem as bwa_mem from 'NextflowCommandLineTools/bwa/0.7.5a.nf' params(params)
include UnifiedGenotyper as gatk_UnifiedGenotyper from 'NextflowCommandLineTools/gatk/3.8-1-0-gf15c1c3ef.nf' params(
    process_outdir: 'fingerprint',
    gatk: params.gatk,
    genome: params.genome,
    intervals: params.fingerprint_target,
    dbsnp: params.fingerprint_target
    output_mode = 'EMIT_ALL_SITES'
    )

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

mips_trim_dedup(samples_R1_R2_fastq)
fastqc(samples_all_fastq)
bwa_mem(mips_trim_dedup.out)
gatk_UnifiedGenotyper(bwa_mem.out)
