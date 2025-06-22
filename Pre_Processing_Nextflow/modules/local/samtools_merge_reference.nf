process SAMTOOLS_MERGE_REFERENCE {
    tag "$prefix"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
                "file://${projectDir}/singularity_container/samtools1.15.1.sif" : 'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' }"
    
    /*
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    */

    input:
    path(inputs)

    output:
    path '*.bam'                               , emit: merged_bam
    path '*.bam.bai'                           , emit: bai
    path '*.tsv'                              , emit: tsv
    path 'versions.yml'                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: ''
    //def bamlist= inputs.join(' ')
    """
    samtools \\
        merge \\
        -@ ${task.cpus-1} \\
        $args \\
        -o ${prefix}.bam \\
        ${inputs}
    
    samtools index ${prefix}.bam -@ ${task.cpus-1}

    coverage=\$(samtools view -c -@ 10 -q 10 -F 4 -F 2048 -F 256 -F 1024 ${prefix}.bam)
    echo -e "${prefix}.bam \t \${coverage}" >> readCount_${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}