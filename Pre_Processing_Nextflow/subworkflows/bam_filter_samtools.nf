//
// Filter and index BAM file
//

include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FILTER     } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX                            } from '../modules/nf-core/samtools/index/main'

workflow BAM_FILTER_SAMTOOLS {
    take:
    ch_bam_bai          // channel: [ val(meta), [ bam ], [bai] ]
    fasta           // path(fasta)

    main:

    ch_versions = Channel.empty()

    
    SAMTOOLS_VIEW_FILTER ( ch_bam_bai , fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FILTER.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_VIEW_FILTER.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    SAMTOOLS_VIEW_FILTER.out.bam
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

    emit:
    bam      = SAMTOOLS_VIEW_FILTER.out.bam    // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}