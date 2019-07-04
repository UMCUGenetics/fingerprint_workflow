version 1.0

task fastqc {
    input {
        String fastqc_path
        Array[File] r1_fastqs
        Array[File] r2_fastqs
    }

    command {
        ${fastqc_path} --noextract -t 1 -o ./ ${sep=' ' r1_fastqs} ${sep=' ' r2_fastqs}
    }

    output {
        Array[File] fastqc_reports = glob("*.zip")
    }
}
