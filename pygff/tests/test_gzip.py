import os
import pygff.core as gfp
import pytest

data_path = os.path.join(os.path.dirname(gfp.__file__), 'tests', 'data')

def test_not_gzip():
    assert not gfp._is_zipped(os.path.join(data_path, 'not_zipped.gff'))

def test_is_gzip():
    assert gfp._is_zipped(os.path.join(data_path, 'zipped.gff.gz'))
