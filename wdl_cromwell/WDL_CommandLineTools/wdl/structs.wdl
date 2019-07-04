version 1.0

struct Sample {
    String name
    String library_id
    Array[File] r1_fastqs
    Array[File] r2_fastqs
}

struct BamFile {
    File file
    File index
}
