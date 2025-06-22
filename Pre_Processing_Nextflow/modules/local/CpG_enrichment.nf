process CPG_ENRICHMENT {
    label 'process_intermediate'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "file://${projectDir}/singularity_container/rtools-amd64_2.0.sif" : 'docker://ariediger/rtools-amd64:2.0' }"

    input:
    path(finalbams)
    val single_end

    output:
    path("*.tsv")                     , emit: cpg_enrichment_summary
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in /bin/R/CpG_enrichment.R
    // def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_single = single_end ? "SINGLE" : "PAIRED"
    
    """
    echo ${finalbams} | tr " " "\\n" > finalbam_files.txt

    Rscript ${projectDir}/bin/CpG_enrichment.R \\
        finalbam_files.txt \\
        ${paired_single} 
   
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        MEDIPS: \$(Rscript -e "library(MEDIPS); cat(as.character(packageVersion('MEDIPS')))")
        BSgenome.Hsapiens.UCSC.hg19: \$(Rscript -e "library(BSgenome.Hsapiens.UCSC.hg19); cat(as.character(packageVersion('BSgenome.Hsapiens.UCSC.hg19')))")
        data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")
    
    END_VERSIONS
    """
}