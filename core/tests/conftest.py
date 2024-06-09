import pytest


def pytest_addoption(parser):
    """Add custom command line arguments."""

    # flag for determining CI environment
    parser.addoption(
        "--env_type",
        action="store",
        default="latest",
        help="Testing environment type label",
        type=str,
    )


def pytest_configure(config):
    """PyTest session attributes, methods, etc."""

    env_type = config.getoption("env_type")

    # set ``ENV_TYPE``
    valid_env_types = ["oldest", "latest", "dev"]
    pytest.ENV_TYPE = env_type.split("-")[-1].split(".")[0]
    assert pytest.ENV_TYPE in valid_env_types
