import os
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO not in sys.path:
    sys.path.insert(0, REPO)


def pytest_addoption(parser):
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Regenerate golden files in tests/golden/ instead of comparing.",
    )


@pytest.fixture
def update_golden(request):
    return request.config.getoption("--update-golden")
