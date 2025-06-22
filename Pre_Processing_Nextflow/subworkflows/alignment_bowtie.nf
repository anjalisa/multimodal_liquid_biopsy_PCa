//
// Alignment with Bowtie2, coordinate-sort and index BAM file
//

include { BOWTIE2_ALIGN                 } from '../modules/nf-core/bowtie2/align/main'
//include { BAM_SORT_SAMTOOLS           } from './bam_sort_samtools'
include { FASTQC_BAM as FASTQC_MAP      } from '../modules/local/fastqcbam'
include { SAMTOOLS_INDEX                } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT                 } from '../modules/nf-core/samtools/sort/main'
include { BAM_STATS_SAMTOOLS            } from '../subworkflows/bam_stats_samtools'

workflow ALIGNEMENT_BOWTIE {
    take: 
    reads                  // channel: [ val(meta), [ reads ] ]
    index                  // 
    set_readgroups         // boolean: true/false
    save_unaligned         // boolean: true/false
    //sort_bam               // boolean: true/false (true=sort ; false = view)
    skip_samtools_stats_bowtie2    // boolean: true/false
    skip_fastqc_map        // boolean: true/false

    main:

    ch_versions = Channel.empty()
     
    //
    // Map reads with BOWTIE2
    //
    
    BOWTIE2_ALIGN (reads, index, set_readgroups, save_unaligned)
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //
    // Sort and index BAM file (optional: and run samtools stats, flagstat and idxstats)
    //
    
    SAMTOOLS_SORT ( BOWTIE2_ALIGN.out.sam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    SAMTOOLS_SORT.out.bam
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
    
    if (!params.skip_samtools_stats_bowtie2) {
        BAM_STATS_SAMTOOLS ( ch_bam_bai )
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)
    }
    
    //
    // FASTQC of aligned BAM files
    
    fastqc_map_html = Channel.empty()
    fastqc_map_zip  = Channel.empty()
    
    if (!skip_fastqc_map) {
            FASTQC_MAP ( SAMTOOLS_SORT.out.bam )
                fastqc_map_html  = FASTQC_MAP.out.html
                fastqc_map_zip  = FASTQC_MAP.out.zip
                ch_versions = ch_versions.mix(FASTQC_MAP.out.versions.first())
    }
    
    emit:
    sam     = BOWTIE2_ALIGN.out.sam           // channel: [ val(meta), sam   ]
    log_summary = BOWTIE2_ALIGN.out.log       // channel: [ val(meta), log   ]
    fastq    = BOWTIE2_ALIGN.out.fastq         // channel: [ val(meta), fastq ]

    bam      = SAMTOOLS_SORT.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai      // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi      // channel: [ val(meta), [ csi ] ]
    bam_bai  = ch_bam_bai                  // channel: [ val(meta), [ bam ] , [ bai ] ]
    
    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    fastqc_map_html                                // channel: [ val(meta), [ html ] ]
    fastqc_map_zip                                 // channel: [ val(meta), [ zip ] ]

    versions = ch_versions                         // channel: [ versions.yml ]
}