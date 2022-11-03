process FOO {
    echo true
    tag "$meta"

    input:
    tuple val(meta), path(files)

    output:
    stdout emit: output

    script:
    """
    echo $files
    """
}
