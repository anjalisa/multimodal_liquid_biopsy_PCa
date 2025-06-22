process FASTX_TRIMMER {
    tag "$meta.id"
    label 'process_medium'

    // FASTX-Toolkit: https://github.com/agordon/fastxtoolkit/releases/download/0.0.14/fastxtoolkit-0.0.14.tar.bz2

    conda (params.enable_conda ? "bioconda::fastx_toolkit=0.0.14" : null)
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
                "file://${projectDir}/singularity_container/fastx0.0.14.sif" : 'https://depot.galaxyproject.org/singularity/fastx_toolkit:0.0.14--h87f3376_10' }"
    
    /*
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastx_toolkit:0.0.14--h87f3376_10' :
        'quay.io/biocontainers/fastx_toolkit' }"
    */

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_{R1,R2}.fastq.gz'), emit: reads
   // path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """

      zcat ${reads[0]} | fastx_trimmer $args > ${prefix}_R1.fastq.gz

      zcat ${reads[1]} | fastx_trimmer $args > ${prefix}_R2.fastq.gz

    """

}

/*
# Define parameters
FASTX_TRIMMER=$1 #Path to fastx_trimer; '/home/jankef/.conda/pkgs/fastx_toolkit-0.0.14-0/bin/fastx_trimmer'
READ_LENGTH=$2 #Desired maximum read length; e.g., 100
INPUT=$3 #Path to input .fastq.gz file; e.g., '/omics/groups/OE0309/internal/janke/Literature_data/WGS/fastq/SRR17478151_1.fastq.gz'
OUTPUT=$4 #Path to ouput .fastq.gz file; e.g., '/omics/groups/OE0309/internal/janke/test/SRR17478151_1.fastq.gz'

# Trim .fastq.gz files to <READ_LENGTH>
zcat ${INPUT} | ${FASTX_TRIMMER} -l ${READ_LENGTH} -Q 33 -z > ${OUTPUT}
*/
