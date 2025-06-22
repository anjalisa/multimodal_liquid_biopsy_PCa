process QC_SUMMARY {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "file://${projectDir}/singularity_container/rtools-amd64_2.0.sif" : 'docker://ariediger/rtools-amd64:2.0' }"


    input:
    path(saturation)
    path(cpg_coverage)
    path(enrichment)
    path(readcount)
    path(fastqc)
    val skip_methylation_qc
    // path(outdir)

    output:
    path("QC_summary.tsv")                      , emit: summary_tsv  
    path("Metric_summary.tsv")                  , emit: metrics_tsv
    path("Coverage_overview.svg")               , emit: coverage_svg
    path("Coverage_overview.png")               , emit: coverage_png
    path("Saturation_overview.png")             , emit: saturation_png
    path("Saturation_overview.svg")             , emit: saturation_svg
    path("CpG_enrichment.png")                  , emit: enrichment_png
    path("CpG_enrichment.svg")                  , emit: enrichment_svg
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in /bin/R/QC_summary.R
    // def prefix = task.ext.prefix ?: "${meta.id}"
    
    def  methyl_qc = skip_methylation_qc ? "SKIPPED" : "NOT_SKIPPED"
    """ 
    Rscript ${projectDir}/bin/QC_summary.R \\
        $saturation \\
        $cpg_coverage \\
        $enrichment \\
        $readcount \\
        $fastqc \\
        $methyl_qc
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        reshape2: \$(Rscript -e "library(reshape2); cat(as.character(packageVersion('reshape2')))")
        matrixStats: \$(Rscript -e "library(matrixStats); cat(as.character(packageVersion('matrixStats')))")

    END_VERSIONS
    """
}