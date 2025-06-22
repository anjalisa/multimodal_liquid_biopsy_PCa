//
// (Read QC), UMI extraction and trimming
//

// include { FASTQC           } from '../../modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT                      } from '../modules/nf-core/umitools/extract/main'
include { CUTADAPT                              } from '../modules/nf-core/cutadapt/main'
include { FASTQC as FASTQC_RAW                  } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM                 } from '../modules/nf-core/fastqc/main'

workflow UMITOOLS_CUTADAPT {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    skip_fastqc_raw     // boolean: true/false
    skip_fastqc_trim    // boolean: true/false
    skip_trimming       // boolean: true/false
    // umi_discard_read // integer: 0, 1 or 2

    main:

    ch_versions = Channel.empty()
    fastqc_raw_html = Channel.empty()
    fastqc_raw_zip  = Channel.empty()
        
    if (!skip_fastqc_raw) {
        FASTQC_RAW ( reads )
        fastqc_raw_html  = FASTQC_RAW.out.html
        fastqc_raw_zip  = FASTQC_RAW.out.zip
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
    }

    umi_reads = reads
    umi_log   = Channel.empty()
        
    UMITOOLS_EXTRACT ( reads ).reads.set{umi_reads}
   
    umi_log     = UMITOOLS_EXTRACT.out.log
    ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

                /*  // Discard R1 / R2 if required
                if (umi_discard_read in [1,2]) {
                    UMITOOLS_EXTRACT
                        .out
                        .reads
                        .map { meta, reads ->
                            if (!meta.single_end) {
                                meta['single_end'] = true
                                reads = reads[umi_discard_read % 2]
                            }
                            return [ meta, reads ]
                        }
                        .set { umi_reads } 
            } */
    
    umi_reads.view()

    trim_reads       = umi_reads       
    //trim_json        = Channel.empty()
    trim_log         = Channel.empty()
           
    fastqc_trim_html = Channel.empty()
    fastqc_trim_zip  = Channel.empty()
            
    if (!skip_trimming) {
        CUTADAPT ( umi_reads ).reads.set{trim_reads}
        //trim_json     = CUTADAPT.out.json
        trim_log      = CUTADAPT.out.log
        
        ch_versions   = ch_versions.mix(CUTADAPT.out.versions.first())
        
        if (!skip_fastqc_trim) {
            FASTQC_TRIM ( trim_reads )
            fastqc_trim_html  = FASTQC_TRIM.out.html
            fastqc_trim_zip  = FASTQC_TRIM.out.zip
            ch_versions = ch_versions.mix(FASTQC_TRIM.out.versions)
        }
    }

    emit:
    reads = trim_reads          // channel: [ val(meta), [ reads ] ]

    umi_log            // channel: [ val(meta), [ log ] ]

    trim_log           // channel: [ val(meta), [ txt ] ]
    //trim_json          // channel: [ val(meta), [ json ] ]
    
    fastqc_raw_html        // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip         // channel: [ val(meta), [ zip ] ]

    fastqc_trim_html        // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip         // channel: [ val(meta), [ zip ] ]

    versions = ch_versions  // channel: [ versions.yml ]

}