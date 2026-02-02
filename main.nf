nextflow.enable.dsl=2

params.input_dir = "data/*.fastq.gz"
params.outdir = "results"

workflow {
    input_ch = Channel.fromPath(params.input_dir)

    // Assembly Step (now outputs both GFA and FASTA)
    hifiasm_out = HIFIASM(input_ch)

    // Orientation Step (uses the FASTA output)
    ORIENTATION(hifiasm_out)
}

process HIFIASM {
    container 'quay.io/biocontainers/hifiasm:0.25.0--h5ca1c30_0'
    // Save both the graph and the sequence to the results folder
    publishDir "${params.outdir}/${reads.baseName}/assembly", mode: 'copy'
    
    input:
    path reads

    output:
    // Pass both files to the next process as a tuple
    tuple val(reads.baseName), path("${reads.baseName}.p_ctg.gfa"), path("${reads.baseName}.p_ctg.fa")

    script:
    """
    # -l0: Disable duplication purging for homozygous/haploid genomes
    hifiasm -o ${reads.baseName} -t ${task.cpus} -l0 ${reads}

    # Convert GFA to FASTA using a standard awk command
    awk '/^S/{print ">"\$2;print \$3}' ${reads.baseName}.bp.p_ctg.gfa > ${reads.baseName}.p_ctg.fa
    
    # Rename primary GFA for cleaner output
    mv ${reads.baseName}.bp.p_ctg.gfa ${reads.baseName}.p_ctg.gfa
    """
}

process ORIENTATION {
    // This container includes both python and mappy (minimap2 python bindings)
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

