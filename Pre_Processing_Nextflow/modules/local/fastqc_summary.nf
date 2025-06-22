process FASTQC_SUMMARY {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "file://${projectDir}/singularity_container/rtools-amd64_2.0.sif" : 'docker://ariediger/rtools-amd64:2.0' }"


    input:
    path(fastqc_zip)

    output:
    path("fastqc_summary.tsv")                      , emit: fastqc_summary  
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in /bin/R/QC_summary.R
    // def prefix = task.ext.prefix ?: "${meta.id}"
    
    """ 
    echo ${fastqc_zip} | tr " " "\\n" > fastqc_files.txt
    
    Rscript ${projectDir}/bin/FASTQC_summary.R \\
        fastqc_files.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
    
    END_VERSIONS
    """
}