'''
Authors: Aleksey Ermolaev, Danko Katerina

Script to filter fastq file.

Optional arguments:

    --min_length INT
    Minimal length for reads to pass the filter. Must be > 0.

    --keep_filtered
    Flag. Write not passed the filter reads in file.

    --gc_bounds INT INT
    Takes 1 or 2 values separated by sapce - lower and upper bounds of
    GC-content for read to pass the filter. If single value is get its
    counts as upper bound. Lower bound must lower of equal to upper.
    Both bounds must be >= 0.

    --output_base_name STR
    Prefix for output file(s).


Positional arguments:

    fastq file
    fastq file to process.


Example:
    python fastq_filtrator.py --min_length 50 --keep_filtered
    --gc_bounds 55 70 --output_base_name output_reads reads.fastq
'''
