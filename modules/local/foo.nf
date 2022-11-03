process FOO {
    echo true
    tag "$meta.id"

    input:
    tuple val(meta), path(files)

    output:
    stdout emit: output

    script:
    """
    echo $files
    """
}
