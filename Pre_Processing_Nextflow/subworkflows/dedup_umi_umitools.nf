//
// UMI-tools dedup, index BAM file and run samtools stats, flagstat and idxstats
// nf-core workflow, but ADAPTED

include { UMITOOLS_DEDUP                } from '../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX                } from '../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS            } from './bam_stats_samtools'
include { FASTQC_BAM as FASTQC_FINAL    } from '../modules/local/fastqcbam'

workflow DEDUP_UMI_UMITOOLS {
    take:
    bam_bai                     // channel: [ val(meta), [ bam ], [ bai/csi ] ]
    umitools_dedup_stats        // boolean: true/false
    skip_samtools_stats_dedup   // boolean: true/false
    skip_fastqc_final           // boolean: true/false
    
    main:

    ch_versions = Channel.empty()

    //
    // UMI-tools dedup
    //
    UMITOOLS_DEDUP ( bam_bai, umitools_dedup_stats )
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions.first())

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX ( UMITOOLS_DEDUP.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    UMITOOLS_DEDUP.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    if (!params.skip_samtools_stats_dedup) {
        BAM_STATS_SAMTOOLS ( ch_bam_bai )
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions.first())
    }

    //  
    // FASTQC of final BAM files
    
    fastqc_final_html = Channel.empty()
    fastqc_final_zip  = Channel.empty()
    
    if (!skip_fastqc_final) {
            FASTQC_FINAL ( UMITOOLS_DEDUP.out.bam )
                fastqc_final_html  = FASTQC_FINAL.out.html
                fastqc_final_zip  = FASTQC_FINAL.out.zip
                ch_versions = ch_versions.mix(FASTQC_FINAL.out.versions)
    }

    emit:
    bam      = UMITOOLS_DEDUP.out.bam          // channel: [ val(meta), [ bam ] ]

    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    fastqc_final_html                          // channel: [ val(meta), [ html ] ]
    fastqc_final_zip                           // channel: [ val(meta), [ zip ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}