[![Conda Version](https://img.shields.io/conda/vn/conda-forge/pygff.svg)](https://anaconda.org/conda-forge/pygff) 
[![PyPI version](https://badge.fury.io/py/pygff.svg)](https://badge.fury.io/py/pygff)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/betteridiot/pygff/blob/master/LICENSE) 
</br>

# Installation:

There are multiple options available for installation

```bash
# Recommended: from conda-forge
conda install -c conda-forge pygff

# From PyPI
pip install pygff

# From GitHub
git clone https://github.com/betteridiot/pygff.git
cd pygff
python3 ./setup.py install

```

Optionally, you may want to test the program. A very small testing suite has been
provided. **Note**: *`pytest` is required* for testing.

If you choose to test the program on your current environment/build:

```bash
# If installed from conda-forge or the PyPI
pytest --pyargs pygff

# Or, if installed from source
python ./setup.py test

```

---
# 'As-is' Warranty:
This program was written expressly for the purposes of education at the University
of Michigan's Department of Computational Medicine & Bioinformatics. It was only
intended to work for a specific subset of GFF/GTF files: **GFF3 files**. No promises are
made as to its wider functionality.

---
# Background:
The General Feature Format (GFF) was developed as a way to succinctly represent
genomic features (e.g. exons, introns, genes, etc). They are a 9-column, tab-delimited
plain text (or gzip compressed) file. The 9 columns are described as such:

| Column | Content | Description |
| :----- | :-----: | :---------- |
| 1 | seqid | The ID of the landmark used to establish the coordinate system for the current feature |
| 2 | source | The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this feature |
| 3 | type | The type of the feature (previously called the "method") |
| 4 | start | The start coordinate of the feature, given in positive 1-based integer coordinates, relative to the landmark given in column one |
| 5 | end | The end coordinate of the feature, given in positive 1-based integer coordinates, relative to the landmark given in column one |
| 6 | score | The score of the feature, a floating point number |
| 7 | strand | The strand of the feature. '+' for positive strand (relative to the landmark), '-' for minus strand, and '.' for features that are not stranded |
| 8 | phase | For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame |
| 9 | attributes | A list of feature attributes in the format tag=value |

**Note**: this information was taken, in part, from [The Sequence Ontology](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) group's GitHub.

---
# Description:
`pygff` was written to provide a helpful interface when handling GFF3 files. It 
allows the user to easily process the parsed GFF3 content, entry-by-entry, in a 
lazily generated way.

However, an additional functionality was written that produces a GFF3 index 
on the fly. This index allows for *pseudo-random access*. Due to how it is implemented, 
it is not recommended that this program is run from the command line directly since
the index is ephemeral by nature.

`pygff` exposes two classes for dealing with GFF3 data:

### `pygff.GffFile`:
Main class of the pygff package

Handles the opening, iterating, and closing GFF3 files. Can handle both
zipped and unzipped GFF3 files.

When iterated, it lazily returns a `pygff.GffEntry` object. These objects
can be compared against each other and programmatically accessed for all traits.

**Args**:
1. `filename` (`str`): /path/to/file.gff[.gz]
2. `periods` (`int`): For indexing purposes, the period to determine thresholding (default: 3)

**Raises**:
* `TypeError` if GFF file is not version 3

And exposes the `pygff.GffFile.fetch(seqid, start, stop)` method:

Generator that fetches all GFF entries within a given region. 

Also can only pull specific *types* of GFF entries (if supplied)

**Args**:
1. `seqid` (`str`): name of the chromosome of scaffold
2. `start` (`int`): start position of the feature (1-indexed)
3. `end` (`int`): end position of the feature (1-indexed)
4. `type` (`str`): GFF feature type (default: None)

**Yields**:
* (`pygff.GffEntry`): A given GFF entry from the region of interest

### `pygff.GffEntry`:
An object that represents a single GFF entry. 

This object also has the ability to perform total ordered comparison
operations (<, <=, ==, !=, >=, >) based on seqid first, then start
position, and finally the end position.

**Attributes**:
1. `seqid` (`str`): name of the chromosome of scaffold</br>
2. `source` (`str`): name of the program that generated the feature</br>
3. `type` (`str`): type of feature</br>
4. `start` (`int`): start position of the feature (1-indexed)</br>
5. `end` (`int`): end position of the feature (1-indexed)</br>
6. `score` (`float`): a quality score of the feature</br>
7. `strand` (`str`): either '+' (forward), '-'(reverse), or '.'</br>
8. `phase` (`int`): 0,1, or 2 that indicates that the first base of the is the first base of the codon</br>
9. `attributes` (`dict`): a dictionary of all tag/value pairs</br>

---
# Quickstart

### Importing

```python

import pygff

```

### Sequential Iteration

```python

with pygff.GffFile('/path/to/file.gff[.gz]') as gff:
    for entry in gff:
        do_something(entry)

```

### *Pseudo*-Random Access

```python
gff = pygff.GffFile('/path/to/file.gff[.gz]')
for entry in gff.fetch('chr1', 123040, 128040):
    do_something(entry)

```

### Output

```python

with open('outfile.gff', 'wb') as outfile:
    with pygff.GffFile('/path/to/file.gff[.gz]') as gff:
        for entry in gff:
            # Some filtering
            print(entry, file = outfile)

```

---
# Contributing & Code of Conduct:
This project is built on Open Science, Open Source, and Open Minds. To encourage
an environment of inclusivity and positivity, please see our [Code of Conduct](https://github.com/betteridiot/pygff/blob/master/CODE_OF_CONDUCT.md).

If you are interested in contributing to the project, please see the [CONTRIBUTING](https://github.com/betteridiot/pygff/blob/master/CONTRIBUTING.md) guidelines
