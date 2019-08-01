#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include mips_trim_dedup from 'NextflowCommandLinetools/mips/1.0.1'
include fastqc from 'NextflowCommandLinetools/fastqc/0.11.8'
include mem as bwa_mem from 'NextflowCommandLinetools/bwa/0.7.5a'
include UnifiedGenotyper as gatk_UnifiedGenotyper from 'NextflowCommandLinetools/gatk/3.8-1-0-gf15c1c3ef'

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
