version 1.0

import "../wdl/structs.wdl"

task UnifiedGenotyper {
    input {
        String gatk_path
        String genome
        File target_intervals
        File dbsnp
        BamFile input_bam_file
        String out_name
    }

    command {
        java -jar ${gatk_path} -T UnifiedGenotyper --reference_sequence ${genome} --intervals ${target_intervals} --dbsnp ${dbsnp} --input_file ${input_bam_file.file} --out ${out_name}.vcf
    }

    output {
        File vcf = "${out_name}.vcf"
    }
}
