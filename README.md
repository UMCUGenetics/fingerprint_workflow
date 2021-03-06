# Fingerprint workflow
Comparison of workflow languages using the DX mip fingerprint pipeline.

- https://www.biocentric.nl/tag/framework/
- https://www.biocentric.nl/biocentric/bioinformatics-frameworks-part-4-nextflow/
- https://www.biocentric.nl/biocentric/bioinformatics-frameworks-part-5-cromwell-and-wdl/

#### Pipeline components:
- FastQC
- MIPS trim dedup
- BWA mem
- GATK UnifiedGenotyper

## Nextflow
- Website: https://www.nextflow.io/
- Docs: https://www.nextflow.io/docs/latest/index.html
- Useful resources:
  - https://github.com/nextflow-io/patterns
  - https://github.com/oliverSI/GATK4_Best_Practice
  - https://github.com/SciLifeLab/Sarek

#### Pros
- SGE / SLURM and cloud support.
- Continuous checkpoints and result caching on local file system.
- Singularity support.
- Easy to store final pipeline output in a structured way.
- Github support (running workflows from a github repo)

#### Cons
- Modules support in development, but coming soon.
    - Current strategy is to duplicate code.
    - https://www.nextflow.io/blog/2019/one-more-step-towards-modules.html

#### Install
```bash
mkdir /install/dir/location
cd /install/dir/location
module load Java/1.8.0_60
curl -s https://get.nextflow.io | bash
```
#### Running fingerprint pipeline
```bash
module load Java/1.8.0_60
../tools/nextflow run MipFingerprint.nf --samplesheet samples.tsv --outdir <output_dir_name>
```

## Cromwell / WDL
- Website: https://software.broadinstitute.org/wdl/
- Docs Cromwell: https://cromwell.readthedocs.io/en/stable/
- Docs WDL: https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md

#### Pros
- SGE / SLURM and cloud support.
- Module support.

#### Cons
- Result caching needs mysql database / server.
- Singularity support needs a lot of configuration but is possible.

#### Install
```bash
mkdir /install/dir/location
cd /install/dir/location
wget https://github.com/broadinstitute/cromwell/releases/download/38/cromwell-38.jar
```

#### Running fingerprint pipeline
```bash
java -Dconfig.file=cromwell.conf -jar ../tools/cromwell-38.jar run -i mip_fingerprint.inputs.json -o mip_fingerprint.options.json MipFingerprint.wdl
```
