import os
import pygff.core as gfp
import pytest

data_path = os.path.join(os.path.dirname(gfp.__file__), 'tests', 'data')

def test_not_version_3():
    with pytest.raises(TypeError):
        gfp.GffFile(os.path.join(data_path, 'not_zipped_v2.gff'))

def test_is_version_3():
    assert gfp.GffFile(os.path.join(data_path, 'not_zipped.gff'))
