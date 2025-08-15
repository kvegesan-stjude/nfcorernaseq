process PROCESS_PREMAPPED_READS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(reads, stageAs: "input*/*")
    path(exclude_bed)

    output:
    tuple val(meta), path('*.pre-mapped.fq.gz')                         , emit: premapped_correct_reads
    tuple val(meta), path('*_pre-mapped.singletons.fq.gz')              , emit: premapped_signleton_reads, optional: true
    tuple val(meta), path('*_pre-mapped.broken.fq.gz')                  , emit: premapped_broken_reads   , optional: true
    tuple val(meta), path('*.excluded.bam'), path('*.excluded.bam.bai') , emit: excluded_bam
    path  "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def reads1 = [], reads2 = []
    meta.single_end ? [reads].flatten().each{reads1 << it} : reads.eachWithIndex{ v, ix -> ( ix & 1 ? reads2 : reads1) << v }

    // filtering parameters
    def fastq_out_primary = "", concatenate_cmd = ""
    if (meta.single_end) {
        fastq_out_primary = "-0 ${prefix}.pre-mapped.fq.gz"
        concatenate_cmd = "cat ${reads1.join(" ")} >> ${prefix}.pre-mapped.fq.gz"
    } else {
        fastq_out_primary = "-s ${prefix}_pre-mapped.singletons.fq.gz -0 ${prefix}_pre-mapped.broken.fq.gz -1 ${prefix}_1.pre-mapped.fq.gz -2 ${prefix}_2.pre-mapped.fq.gz"
        concatenate_cmd = "cat ${reads1.join(" ")} >> ${prefix}_1.pre-mapped.fq.gz && cat ${reads2.join(" ")} >> ${prefix}_2.pre-mapped.fq.gz"
    }
    def destination_dir = '/research_jude/rgs01_jude/groups/thomagrp/home/kvegesan/Projects/Z-DoTT/.temp/'

    """
    cp ${bam} ${destination_dir}
    # Sort and index the incoming BAM file
    samtools sort --threads ${task.cpus} -O bam -o ${prefix}.sorted.tmp.bam ${bam}
    samtools index --threads ${task.cpus} ${prefix}.sorted.tmp.bam

    # Separate reads into uniquely mappable to the exclude_bed and those that are not
    samtools view \\
        --threads ${task.cpus} \\
        --bam \\
        --exclude-flag 260 \\
        --targets-file ${exclude_bed} \\
        --output ${prefix}.excluded.bam \\
        ${prefix}.sorted.tmp.bam
    samtools index --threads ${task.cpus} ${prefix}.excluded.bam

    # Extract the read names of the excluded reads
    samtools view --threads ${task.cpus} ${prefix}.excluded.bam \\
        | cut -f1 \\
        | bgzip --threads ${task.cpus} --output ${prefix}.excluded-reads.txt.gz

    # Drop all excluded reads from the original bam and turn them into fastq
    samtools view --threads ${task.cpus} -h --qname-file ^<(zcat ${prefix}.excluded-reads.txt.gz) ${prefix}.sorted.tmp.bam | \\
    samtools collate --threads ${task.cpus} -O -r 100000 - collate-tmp | \\
    samtools fastq --threads ${task.cpus} -n ${fastq_out_primary} -

    rm ${prefix}.sorted.tmp.bam ${prefix}.sorted.tmp.bam.bai ${prefix}.excluded-reads.txt.gz

    # Concatenate extra reads to the end of the fastq files
    ${concatenate_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
