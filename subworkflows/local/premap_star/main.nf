//
// Reads pre-mapping with STAR
//

include { STAR_ALIGN              } from '../../../modules/nf-core/star/align'
include { PROCESS_PREMAPPED_READS } from '../../../modules/local/process_premapped_reads'

workflow PREMAP_STAR {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    index               // channel: [ val(meta), [ index ] ]
    gtf                 // channel: [ val(meta), [ gtf ] ]
    star_ignore_sjdbgtf // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
    seq_platform        // string : sequencing platform
    seq_center          // string : sequencing center
    premap_bed          // channel: path(premap_bed)

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with STAR
    //
    STAR_ALIGN ( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
    ch_orig_bam       = STAR_ALIGN.out.bam
    ch_log_final      = STAR_ALIGN.out.log_final
    ch_log_out        = STAR_ALIGN.out.log_out
    ch_log_progress   = STAR_ALIGN.out.log_progress
    ch_fastq          = STAR_ALIGN.out.fastq
    ch_versions       = ch_versions.mix(STAR_ALIGN.out.versions.first())

    //
    // Exclude reads overlapping with a BED file and turn the rest into a FASTQ file
    //
    ch_bam_and_unmapped = ch_orig_bam.join(ch_fastq, by: [0])
    PROCESS_PREMAPPED_READS ( ch_bam_and_unmapped, premap_bed )

    emit:
    excluded_bam    = PROCESS_PREMAPPED_READS.out.excluded_bam              // channel: [ val(meta), excluded_bam   ]
    premapped_fastq = PROCESS_PREMAPPED_READS.out.premapped_correct_reads   // channel: [ val(meta), included_fastq ]

    star_log_final      = ch_log_final     // channel: [ val(meta), log_final    ]
    star_log_out        = ch_log_out       // channel: [ val(meta), log_out      ]
    star_log_progress   = ch_log_progress  // channel: [ val(meta), log_progress ]

    versions            = ch_versions      // channel: [ versions.yml ]
}
