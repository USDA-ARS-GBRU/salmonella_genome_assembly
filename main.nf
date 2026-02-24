nextflow.enable.dsl=2

params.outdir = "/project/gbru_fy21_salm_compgen/annette/salmonella_genome_assembly/work"
params.reads = "/project/gbru_fy21_salm_compgen/annette/salmonella_genome_assembly/work/bam_data/*.fastq.gz"

// Define where to store the persistent Bakta database
params.bakta_db_path = "${projectDir}/bakta_db"
bakta_db_ch = Channel.value( file(params.bakta_db_path) )

workflow {
    // Reads were filtered manually using hifiadapfilter.
    filtered_reads = Channel.fromPath(params.reads)
    // Assemble 
    hifiasm_out    = HIFIASM(filtered_reads)
    // If the assembly is empty, shunt it to a different process and log it.
    valid_assemblies = hifiasm_out.filter { sample_id, gfa, fasta -> fasta.size() > 0 }
    empty_assemblies = hifiasm_out.filter { sample_id, gfa, fasta -> fasta.size() == 0 }

    LOG_EMPTY_ASSEMBLY(empty_assemblies)
    
    // Add the orientation script to each assembly tuple using map
    orientation_input = valid_assemblies.map { sample_id, gfa, fasta ->
        tuple(sample_id, gfa, fasta, file("${projectDir}/bin/orientation.py"))
    }

    // Non-empty assemblies continue in the process
    oriented_fa    = ORIENTATION(orientation_input)
    trna_out       = TRNASCAN(oriented_fa)

    // Annotation: Uses the oriented FASTA and the DB from the download step
    bakta_input = oriented_fa.combine(bakta_db_ch)
    BAKTA_ANNOTATE(bakta_input)
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

process LOG_EMPTY_ASSEMBLY {
    container 'quay.io/biocontainers/bakta:1.12.0--pyhdfd78af_0'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(gfa), path(fasta)

    script:
    """
    SIZE=\$(wc -c < "${fasta}")    
    echo "WARNING: Sample '${sample_id}' skipped — empty assembly (size: \$SIZE bytes)" >&2
    """
}

process ORIENTATION {
    // Container is defined in nextflow.config
    publishDir "${params.outdir}/${sample_id}/oriented", mode: 'copy'

    input:
    tuple val(sample_id), path(gfa), path(fasta), path(orientation_script)

    output:
    tuple val(sample_id), path("oriented_${sample_id}.fasta")

    script:
    """
    python3 ${orientation_script} -i ${fasta} -output oriented_${sample_id}.fasta
    """
}

process TRNASCAN {
    container 'quay.io/biocontainers/bakta:1.12.0--pyhdfd78af_0'
    publishDir "${params.outdir}/${sample_id}/trna", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}.trna.gff")

    script:
    """
    export TMPDIR=`pwd`
    tRNAscan-SE \
        -B \
        --detail \
        --thread ${task.cpus} \
        -o ${sample_id}.trna.gff \
        ${fasta}
    """
}


process BAKTA_ANNOTATE {
    container 'quay.io/biocontainers/bakta:1.12.0--pyhdfd78af_0'
    publishDir "${params.outdir}/${sample_id}/annotation", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta), path(bakta_db)

    output:
    path "${sample_id}_bakta"

    script:
    """
    bakta --db ${bakta_db} \
          --genus Salmonella \
          --species enterica \
          --complete \
          --gram - \
          --skip-trna \
          --output "${sample_id}_bakta" \
          --prefix "${sample_id}" \
          --threads ${task.cpus} \
          ${fasta}
    """
}
