version 1.0

import "../wdl/structs.wdl"

task mem {
    input {
        String bwa_path
        String library_id
        String genome
        File r1_fastq
        File r2_fastq
        String output_name
    }

    command {
        ${bwa_path} mem -t 1 -c 100 -M -R "@RG\tID:${output_name}_${library_id}\tSM:${output_name}\tPL:ILLUMINA\tLB:${output_name}\tPU:${library_id}" ${genome} ${r1_fastq} ${r2_fastq} | samtools view -b | samtools sort > ${output_name}.bam
        samtools index ${output_name}.bam ${output_name}.bai
    }

    output {
        BamFile bam_file = object {
            file: "${output_name}.bam",
            index: "${output_name}.bai"
            }
    }

}
