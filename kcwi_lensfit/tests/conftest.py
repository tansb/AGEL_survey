import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: end-to-end tests that run the full pipeline (~1 min)")


@pytest.fixture(scope="session")
def demo(tmp_path_factory):
    from kcwilensfit.synthetic import make_demo_cubes
    d = tmp_path_factory.mktemp("demo")
    return make_demo_cubes(str(d))
