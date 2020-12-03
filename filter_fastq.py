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

from collections import namedtuple
import sys
from pathlib import Path
import itertools

args = sys.argv[1:]

try:
    if not Path(args[-1]).exists():
        raise FileExistsError(f'File {args[-1]} does not exist')
    if not Path(args[-1]).is_file():
        raise IsADirectoryError(f'{args[-1]} is not a file')
except (FileExistsError, IsADirectoryError) as fastq_file_error:
    raise fastq_file_error
else:
    fastq_file = args[-1]
    del args[-1]

if '--min_length' in args:
    min_length_index = args.index('--min_length')
    try:
        if not args[min_length_index + 1].isdigit():
            raise ValueError('--min_length must get a numeric value')
        if int(args[min_length_index + 1]) <= 0:
            raise ValueError('Value for --min_length must be > 0')
    except (NameError, ValueError) as min_length_error:
        raise min_length_error
    else:
        min_length = int(args[min_length_index + 1])
        del args[min_length_index:min_length_index + 2]
else:
    min_length = None

if '--output_base_name' in args:
    output_base_name_arg = args.index('--output_base_name')
    try:
        if '--' in args[output_base_name_arg + 1]:
            raise ValueError("Argument --output_base_name must have a value")
    except ValueError as output_base_name_error:
        raise output_base_name_error
    else:
        output_base_name = args[output_base_name_arg + 1]
        del args[output_base_name_arg:output_base_name_arg + 2]
else:
    output_base_name = Path(fastq_file).stem

if '--keep_filtered' in args:
    arg_index = args.index('--keep_filtered')
    if arg_index == len(args) - 1:
        KEEP_FILTERED = True
        del args[arg_index]
    else:
        try:
            if '--' not in args[arg_index + 1]:
                raise ValueError(
                    'Argument --keep-filtered does not require value')
        except ValueError as keep_filtered_error:
            raise keep_filtered_error
        else:
            KEEP_FILTERED = True
            del args[arg_index]
else:
    KEEP_FILTERED = False

if '--gc_bounds' in args:
    arg_index = args.index('--gc_bounds')
    try:
        lower_bound = args[arg_index + 1]
        upper_bound = args[arg_index + 2]
        if (not args[arg_index + 1].isdigit() or
                not args[arg_index + 2].isdigit()):
            raise TypeError('Both bounds values must be numeric')
    except TypeError as gc_bound_error:
        raise gc_bound_error
    except IndexError:
        try:
            upper_bound = args[arg_index + 1]
            if not upper_bound.isdigit():
                raise TypeError('Bound value must be numeric') from None
        except TypeError as gc_bound_error:
            raise gc_bound_error
        else:
            bounds = (int(args[arg_index + 1]), 100)
    else:
        lower_bound = int(args[arg_index + 1])
        upper_bound = int(args[arg_index + 2])
        try:
            if lower_bound < 0 or upper_bound < 0:
                raise ValueError('Both bounds values must be >= 0')
            if lower_bound > upper_bound:
                raise ValueError('Upper bound must be >= lower bound')
        except ValueError as gc_bound_error:
            raise gc_bound_error
        else:
            bounds = (lower_bound, upper_bound)
else:
    bounds = None

args_nt = namedtuple(
    'args',
    'min_length keep_filtered gc_bounds output_base_name fastq_file')
args = args_nt(min_length=min_length,
               keep_filtered=KEEP_FILTERED,
               gc_bounds=bounds,
               output_base_name=output_base_name,
               fastq_file=fastq_file)

