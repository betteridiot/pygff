#!/usr/bin/env python3
"""A simple GFF3 file parser

This program will lazily parse a GFF3 file and yield one entry at a time. It can
take both plain-text GFF3 files or GZIP-compressed files.

Usage:
    >>> import pygff as gfp
    >>> with gfp.GffFile('filename.gff.gz') as g3:
    ...     for line in g3:
    ...         print(line)

"""
import os
import sys
import gzip
import struct
import numpy as np
import pandas as pd
from bisect import bisect_right
from urllib.parse import unquote
from collections import OrderedDict
from functools import total_ordering
from collections.abc import MutableSet

_gzip_header = struct.Struct('<3B')
_true_header = (31, 139, 8)


def _is_version_3(filename, gzipped=False):
    if gzipped:
        gff = gzip.open(filename, 'rb')
    else:
        gff = open(filename, 'r')
    version = int(gff.readline().strip().split()[1])
    return True if version == 3 else False


def _is_zipped(filename):
    with open(filename, 'rb') as check:
        header = _gzip_header.unpack_from(check.peek(struct.calcsize('<3B')))
        return True if header == _true_header else False


def _parse_attrs(attribute_column):
    attr_dict = {}
    for attr in attribute_column.split(';'):
        tag, values = attr.split('=')
        for val in values.split(','):
            attr_dict.setdefault(tag, []).append(unquote(val))
    return attr_dict


class OrderedSet(MutableSet):

    def __init__(self, iterable=None):
        self.end = end = [] 
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:        
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)


@total_ordering
class GffEntry:
    """An object that represents a single GFF entry. 

    This object also has the ability to perform total ordered comparison
    operations (<, <=, ==, !=, >=, >) based on seqid first, then start
    position, and finally the end position.

    Attributes:
        seqid (str): name of the chromosome of scaffold
        source (str): name of the program that generated the feature
        type (str): type of feature
        start (int): start position of the feature (1-indexed)
        end (int): end position of the feature (1-indexed)
        score (float): a quality score of the feature
        strand (str): either '+' (forward), '-'(reverse), or '.'
        phase (int): 0,1, or 2 that indicates that the first base of the
                is the first base of the codon
        attributes (dict): a dictionary of all tag/value pairs

    """
    __slots__ = ['_features', '_attrs', '_attributes']

    __all__ = ['seqid', 'source', 'type', 'start', 'end', 'score',
            'strand', 'phase', 'attributes', 'has_tag', 'get_tag']

    def __init__(self, line):
        """Initialize the object

        Args:
            line(str): the string representation of GFF entry
        """
        *self._features, self._attrs = line.strip().split('\t')
        self._attributes = _parse_attrs(self._attrs)

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return "\t".join(self._features + [self._attrs])

    def __eq__(self, other):
        if self.seqid == other.seqid:
            if self.type == other.type:
                if self.start == other.start:
                    if self.end == other.end:
                        return True
        return False

    def __lt__(self, other):
        if self.seqid <= other.seqid:
            if self.start < other.start:
                return True
        return False
    
    @property
    def seqid(self):
        return self._features[0]

    @property
    def source(self):
        return self._features[1]

    @property
    def type(self):
        return self._features[2]

    @property
    def start(self):
        return int(self._features[3])

    @property
    def end(self):
        return int(self._features[4])

    @property
    def score(self):
        score = self._features[5]
        if score == '.':
            return score
        else:
            return float(self._features[5])

    @property
    def strand(self):
        return self._features[6]

    @property
    def phase(self):
        phase = self._features[7]
        if phase == '.':
            return phase
        else:
            return int(phase)

    @property
    def attributes(self):
        return self._attributes

    def has_tag(self, tag):
        """Checks to see if the tag is present in the attributes column

        Args:
            tag (str): a tag of interest

        Returns:
            (bool): whether or not the tag is present
        """
        return True if tag in self._attributes else False

    def get_tag(self, tag):
        """Retrieves the tag of interest

        Args:
            tag (str): the tag of interest

        Returns:
            (str): the value associated with the tag

        Raises:
            KeyError: if the tag isn't present in the attributes column
        """
        return self._attributes[tag]


def _get_average_diffs(seqid_starts, periods = 3):
    thresholds = {}
    for seqid in seqid_starts:
        unqiq_starts = list(seqid_starts[seqid])
        thresholds[seqid] = pd.Series(unqiq_starts).diff(periods = periods).dropna().mean().astype(int)
    return thresholds


def _get_thresholds(handle, periods = 3, gzipped = False):
    handle.seek(0)
    diffs = {}

    for line in handle:
        if gzipped:
            line = line.decode('ascii')
        if not line or line.startswith('#'):
            continue
        seqid, source, type, start, *others = line.strip().split('\t')
        diffs.setdefault(seqid, OrderedSet()).add(int(start))

    handle.seek(0)
    return _get_average_diffs(diffs, periods = periods)


def _gen_index(fileobject, periods = 3, gzipped = False):
    handle = fileobject
    thresholds = _get_thresholds(handle, periods, gzipped)
    index = OrderedDict()
    curr_pos = 0
    curr_idx = 0
    curr_threshold = None
    prev_start = None
    for line in handle:
        if gzipped:
            line = line.decode('ascii')
        if line.startswith('#') or not line:
            curr_pos += len(line)
            continue
        f_line = line.strip().split('\t')
        seqid = f_line[0]
        start = int(f_line[3])
        if seqid not in index:
            prev_start = None
            curr_idx = 0
            curr_threshold = thresholds[seqid]
        if prev_start is None or (start - prev_start) > curr_threshold:
            index.setdefault(seqid, {'start': OrderedDict(), 'pos': OrderedDict()})
            index[seqid]['start'].setdefault(curr_idx, start)
            index[seqid]['pos'].setdefault(curr_idx, curr_pos)
            curr_idx += 1
            prev_start = start
        curr_pos += len(line)
    return index


def _find_le(container, start):
    'Find rightmost value less than or equal to start'
    # Get the index/key position
    i = bisect_right(list(container['start'].values()), start)
    if i:
        # give the offset
        return container['pos'][i-1]
    raise ValueError("{} position not within the GFF start position range".format(start))


class GffFile:
    """Main class of the pygff package

    Handles the opening, iterating, and closing GFF3 files. Can handle both
    zipped and unzipped GFF3 files.

    When iterated, it lazily returns a `pygff.GffEntry` object. These objects
    can be compared against each other and programmatically accessed for all traits.

    Args:
        filename (str): /path/to/file.gff[.gz]
        periods (int): For indexing purposes, the period to determine thresholding (default: 3)

    Raises:
        TypeError if GFF file is not version 3
    """

    def __init__(self, filename, period = 3):
        """Initialize the fileobj

        Note: additional overhead applies due to a generating a small index of the GFF3 file

        Args:
            filename (str): /path/to/file.gff[.gz]
            periods (int): For indexing purposes, the period to determine thresholding (default: 3)

        Raises:
            TypeError if GFF file is not version 3
        """
        if _is_zipped(filename):
            self._handle = gzip.open(filename, 'rb')
            gzipped = True
        else:
            self._handle = open(filename, 'r')
            gzipped = False
        if not _is_version_3(filename, gzipped):
            raise TypeError('GFF file version is {}, but must be 3'.format(version))
        self._index = _gen_index(self._handle, gzipped)

    def fetch(self, seqid, start, end, type = None):
        """Generator that fetches all GFF entries within a given region. 

        Also can only pull specific *types* of GFF entries (if supplied)

        Args:
            seqid (str): name of the chromosome of scaffold
            start (int): start position of the feature (1-indexed)
            end (int): end position of the feature (1-indexed)
            type (str): GFF feature type (default: None)

        Yields:
            (`pygff.GffEntry`): A given GFF entry from the region of interest
        """
        start_offset = _find_le(self._index[seqid], start)
        self._handle.seek(start_offset)

        while True:
            entry = next(self)
            if entry.seqid == seqid:
                # Not there yet
                if entry.end < start:
                    continue
                # Went too far right
                elif entry.start >= end:
                    break
                # Goldilocks
                else:
                    if type is not None:
                        if entry.type != type:
                            continue
                    yield entry
            # Wrong contig
            else:
                if entry.seqid < seqid:
                    continue
                else:
                    break

    def _readline(self):
        return self._handle.readline()

    def __iter__(self):
        return self

    def __next__(self):
        for line in self._handle:
            if type(line) == bytes:
                line = line.decode('ascii')
            if line.startswith('#'):
                continue
            else:
                return GffEntry(line)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self._handle.close()

    def close(self):
        self._handle.close()

    def fileno(self):
        return self._handle.fileno

    def name(self):
        return self._handle.name
