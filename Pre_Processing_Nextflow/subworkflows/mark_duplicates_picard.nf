//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES         } from '../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX                } from '../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS            } from './bam_stats_samtools'
include { FASTQC_BAM as FASTQC_FINAL    } from '../modules/local/fastqcbam'

workflow MARK_DUPLICATES_PICARD {
    take:
    bam                          // channel: [ val(meta), [ bam ] ]
    skip_samtools_stats_dedup    // boolean: true/false
    skip_fastqc_final            // boolean: true/false

    main:

    ch_versions = Channel.empty()

    //
    // Picard MarkDuplicates
    //
    PICARD_MARKDUPLICATES ( bam )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX ( PICARD_MARKDUPLICATES.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    PICARD_MARKDUPLICATES.out.bam
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
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)
    }

    //  
    // FASTQC of final BAM files
    
    fastqc_final_html = Channel.empty()
    fastqc_final_zip  = Channel.empty()
    
    if (!skip_fastqc_final) {
            FASTQC_FINAL ( PICARD_MARKDUPLICATES.out.bam )
                fastqc_final_html  = FASTQC_FINAL.out.html
                fastqc_final_zip  = FASTQC_FINAL.out.zip
                ch_versions = ch_versions.mix(FASTQC_FINAL.out.versions.first())
    }

    emit:
    bam      = PICARD_MARKDUPLICATES.out.bam     // channel: [ val(meta), [ bam ] ]
    metrics  = PICARD_MARKDUPLICATES.out.metrics // channel: [ val(meta), [ metrics ] ]

    bai      = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), [ csi ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), [ idxstats ] ]

    fastqc_final_html                            // channel: [ val(meta), [ html ] ]
    fastqc_final_zip                             // channel: [ val(meta), [ zip ] ]
    versions = ch_versions                       // channel: [ versions.yml ]
}