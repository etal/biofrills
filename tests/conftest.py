import pytest
from alignments import ALIGNED_INPUTS, read_alignment


@pytest.fixture(params=ALIGNED_INPUTS)
def aligned_input(request):
    """Each aligned input of tests/data/inputs."""
    return read_alignment(request.param)
