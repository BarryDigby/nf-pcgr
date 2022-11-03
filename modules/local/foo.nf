process FOO {
    echo true
    tag "$meta"

    input:
    tuple val(meta), path '*'

    output:
    stdout emit: output

    script:
    """
    echo $files
    """
}
