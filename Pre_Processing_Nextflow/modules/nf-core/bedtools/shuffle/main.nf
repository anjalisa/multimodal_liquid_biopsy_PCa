process BEDTOOLS_SHUFFLE {
    tag "$sampleID"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "file://${projectDir}/singularity_container/bedtools2-30-0.sif" :
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' }"
    // Alternative: 'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' 

    input:
    tuple val(sampleID), path(bed)

    output:
    tuple val(sampleID), path("*.bed")     , emit: bed
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sampleID}"
    """
    bedtools \\
        shuffle \\
        -i $bed \\
        -g hg19.genome \\
        $args \\
        > ${prefix}_shuffled_summits.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
