nextflow.enable.dsl=2

params.input_dir = "data/*.bam"
params.outdir = "results"
// Define where to store the persistent Bakta database
params.bakta_db_path = "${baseDir}/bakta_db"

workflow {
    input_ch = Channel.fromPath(params.input_dir)

    // Setup: Download Bakta DB if it doesn't exist
    db_ch = BAKTA_DOWNLOAD()

    // Main Pipeline
    filtered_reads = ADAPTER_FILT(input_ch)
    hifiasm_out    = HIFIASM(filtered_reads)
    oriented_fa    = ORIENTATION(hifiasm_out)
    
    // Annotation: Uses the oriented FASTA and the DB from the download step
    BAKTA_ANNOTATE(oriented_fa, db_ch.collect())
}

process BAKTA_DOWNLOAD {
    container 'staphb/bakta:1.12.0'
    storeDir params.bakta_db_path // Persistent storage across runs

    output:
    path "db-full", type: 'dir'

    script:
    """
    bakta_db download --output . --type full
    """
}


// ... Keep ADAPTER_FILT, HIFIASM, and ORIENTATION processes from previous script ...


process ADAPTER_FILT {
    container 'biocontainers/hifiadapterfilt:2.0.0_cv2'
    publishDir "${params.outdir}/${bam.baseName}/filtered_reads", mode: 'copy'
    
    input:
    path bam

    output:
    // hifiadapterfilt outputs a .filt.fastq.gz file by default
    path "${bam.baseName}.filt.fastq.gz"

    script:
    """
    # hifiadapterfilt.sh runs on all BAM/FASTQ in the work directory if no -p is provided.
    # We use -p to specify the specific file and -t for threads.
    hifiadapterfilt.sh -p ${bam.baseName} -t ${task.cpus}
    """
}

process HIFIASM {
    container 'quay.io/biocontainers/hifiasm:0.25.0--h5ca1c30_0'
    publishDir "${params.outdir}/${reads.baseName}/assembly", mode: 'copy'
    
    input:
    path reads

    output:
    tuple val(reads.baseName), path("${reads.baseName}.p_ctg.gfa"), path("${reads.baseName}.p_ctg.fa")

    script:
    """
    hifiasm -o ${reads.baseName} -t ${task.cpus} -l0 ${reads}
    awk '/^S/{print ">"\$2;print \$3}' ${reads.baseName}.bp.p_ctg.gfa > ${reads.baseName}.p_ctg.fa
    mv ${reads.baseName}.bp.p_ctg.gfa ${reads.baseName}.p_ctg.gfa
    """
}

process ORIENTATION {
    container 'quay.io/biocontainers/mappy:2.26--py310h50d9931_0'
    executor 'local'
    publishDir "${params.outdir}/${sample_id}/oriented", mode: 'copy'

    input:
    tuple val(sample_id), path(gfa), path(fasta)

    output:
    path "oriented_${sample_id}.fasta"

    script:
    """
    python ${baseDir}/orientation.py -i ${fasta} -output oriented_${sample_id}.fasta
    """
}


process BAKTA_ANNOTATE {
    container 'staphb/bakta:1.12.0'
    publishDir "${params.outdir}/${sample_id}/annotation", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)
    path db // Staged from the BAKTA_DOWNLOAD process

    output:
    path "${sample_id}_bakta/*"

    script:
    """
    bakta --db ${db} \
          --genus Salmonella \
          --species enterica \
          --complete \
          --gram - \
          --outdir "${sample_id}_bakta" \
          --prefix "${sample_id}" \
          --threads ${task.cpus} \
          ${fasta}
    """
}