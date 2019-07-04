version 1.0

task mips_trim_dedup {
    input {
        String mips_repo_path
        String design_file
        Array[File] r1_fastqs
        Array[File] r2_fastqs
        Int uuid_length
        String uuid_read
    }

    command {
        python ${mips_repo_path}/mips_trim_dedup.py -d ${design_file} -r1 ${sep=' ' r1_fastqs} -r2 ${sep=' ' r2_fastqs} -l ${uuid_length} -ur ${uuid_read} > output.log 2> output.err
    }

    output {
        Array[File] r1_fastq = glob("*_R1_*.fastq.gz")
        Array[File] r2_fastq = glob("*_R2_*.fastq.gz")
        File log = "output.log"
        File error = "output.err"
    }
}
