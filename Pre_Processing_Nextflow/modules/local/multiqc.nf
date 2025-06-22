process MULTIQC {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::multiqc=1.11" : null)
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
                "file://${projectDir}/singularity_container/multiqc1.12.sif" : 'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' }"
    
    /*
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
                'file:///omics/odcf/analysis/OE0562_projects/early_detection_prostate/singularity_container/multiqc-1.9.sif' : 'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' }"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"
    */
    
    input:
    path(multiqc_config)
    path(multiqc_custom_config)
    path('QC/final_bam/fastqc/*')
    path('QC/picard_multiple_metrics/metrics/*')
    path('trimming_cutadapt/metrics/*') 
    path('mapping_bowtie2/metrics/*')
    path('QC/mapping/samtools_stats/*')
    path('QC/mapping/samtools_stats/*')
    path('QC/mapping/samtools_stats/*')    
    path('QC/final_bam/samtools_stats/*')
    path('QC/final_bam/samtools_stats/*')
    path('QC/final_bam/samtools_stats/*')
    path('dedup_picard/picard_metrics/*')      
     
    //path(fastqc_final)
    //path(picard_multiplemetrics)
    //path(cutadapt)
    //path(bowtie)
    //path(flagstat)
    //path(idxstat)
    //path(stats)
    //path(flagstat_final)
    //path(idxstat_final)
    //path(stats_final)
    //path(markduplicates)
              
    output:
    path("*multiqc_report.html"), emit: report
    path("*_data")              , emit: data
    path("*_plots")             , optional:true, emit: plots
    path("versions.yml")        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc \\
        -f \\
        $args \\
        $custom_config \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}