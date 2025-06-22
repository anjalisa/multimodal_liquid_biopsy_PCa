process MULTIQC {
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
                "file://${projectDir}/singularity_container/multiqc1.12.sif" : 'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' }"
    
    /*
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
                'file://omics/odcf/analysis/OE0562_projects/early_detection_prostate/singularity_container/multiqc-1.9.sif' : 'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' }"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"
    */

    input: 
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    multiqc -f $args .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
