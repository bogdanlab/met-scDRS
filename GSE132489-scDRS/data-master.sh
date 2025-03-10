# load PREAMBLE:
module load java/jdk-11.0.14
cd /u/project/geschwind/lixinzhe/nextflow-process/

# execute the nextflow script: 
nextflow run -resume \
    /u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE132489-scDRS/data-processing.nf