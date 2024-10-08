process SNIFFLES {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::sniffles=2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sniffles:2.4--pyhdfd78af_0' :
        'quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)


    output:
    tuple val(meta), path("*_sniffles.vcf"), emit: sv_calls
    tuple val(meta), path("*_sniffles.snf"), emit: snf_calls
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sniffles \
        --input  $input \
        --vcf ${meta.id}_sniffles.vcf \
        --snf ${meta.id}_sniffles.snf \
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: \$(sniffles --help 2>&1 | grep Version |sed 's/^.*Version: //')
    END_VERSIONS
    """
}

