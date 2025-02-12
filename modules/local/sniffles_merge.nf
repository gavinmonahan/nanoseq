process SNIFFLES_MERGE {
    tag "sniffles_merge"
    label 'process_high'

    conda "bioconda::sniffles=2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sniffles:2.4--pyhdfd78af_0' :
        'quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(snf)


    output:
    path("sniffles_merged.vcf")  , emit: merged_svs
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def input_list = snf.collect{ "--input $it"}.join(' ') 
    
    """
    sniffles \
        $input_list \
        --vcf sniffles_merged.vcf \
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: \$(sniffles --help 2>&1 | grep Version |sed 's/^.*Version: //')
    END_VERSIONS
    """
}

