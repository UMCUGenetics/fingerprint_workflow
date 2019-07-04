version 1.0

import "WDL_CommandLineTools/mips/1.0.1.wdl" as mips
import "WDL_CommandLineTools/bwa/0.7.5a.wdl" as bwa
import "WDL_CommandLineTools/gatk/3.8.1.wdl" as gatk
import "WDL_CommandLineTools/fastqc/0.11.8.wdl" as fastqc
import "WDL_CommandLineTools/wdl/structs.wdl"

workflow MipFingerprint {
    input {
        String genome
        Array[Sample] samples
    }

    scatter (sample in samples) {
        call fastqc.fastqc {
            input:
                r1_fastqs = sample.r1_fastqs,
                r2_fastqs = sample.r2_fastqs,
        }
        
        call mips.mips_trim_dedup {
            input:
                r1_fastqs = sample.r1_fastqs,
                r2_fastqs = sample.r2_fastqs,
        }

        call bwa.mem {
            input:
                r1_fastq = mips_trim_dedup.r1_fastq[0],
                r2_fastq = mips_trim_dedup.r2_fastq[0],
                genome = genome,
                library_id = sample.library_id,
                output_name = sample.name
        }

        call gatk.UnifiedGenotyper{
            input:
                genome = genome,
                input_bam_file = mem.bam_file,
                out_name = sample.name
        }
    }
}
