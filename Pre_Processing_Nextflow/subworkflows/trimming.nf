//
// (Read QC), UMI extraction and trimming
//

// include { FASTQC           } from '../../modules/nf-core/fastqc/main'
include { CUTADAPT                              } from '../modules/nf-core/cutadapt/main'
include { FASTQC as FASTQC_RAW                  } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM                 } from '../modules/nf-core/fastqc/main'

workflow TRIMMING_FASTQC {
    
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    skip_fastqc_raw     // boolean: true/false
    skip_fastqc_trim    // boolean: true/false

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
    
    //trim_json        = Channel.empty()
    trim_log         = Channel.empty()
    trim_reads       = Channel.empty()

    CUTADAPT ( reads )
    trim_reads    = CUTADAPT.out.reads
    //trim_json     = CUTADAPT.out.json
    trim_log      = CUTADAPT.out.log
            
    ch_versions   = ch_versions.mix(CUTADAPT.out.versions.first())

    fastqc_trim_html = Channel.empty()
    fastqc_trim_zip  = Channel.empty()        
    
    if (!skip_fastqc_trim) {
        FASTQC_TRIM ( trim_reads )
        fastqc_trim_html  = FASTQC_TRIM.out.html
        fastqc_trim_zip  = FASTQC_TRIM.out.zip
        ch_versions = ch_versions.mix(FASTQC_TRIM.out.versions.first())
    }     

    emit:
    reads = trim_reads          // channel: [ val(meta), [ reads ] ]

    trim_log           // channel: [ val(meta), [ txt ] ]
    //trim_json          // channel: [ val(meta), [ json ] ]
    
    fastqc_raw_html        // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip         // channel: [ val(meta), [ zip ] ]

    fastqc_trim_html        // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip         // channel: [ val(meta), [ zip ] ]

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

}