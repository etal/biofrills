import pytest
from alignments import ALIGNED_INPUTS, read_alignment


@pytest.fixture(params=ALIGNED_INPUTS)
def aligned_input(request):
    """Each aligned input on which every library function succeeds."""
    return read_alignment(request.param)
