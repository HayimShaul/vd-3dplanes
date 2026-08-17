"""Pytest hooks for the vertical-decomposition suite."""


def pytest_addoption(parser):
    parser.addoption(
        "--n-planes",
        action="store",
        type=int,
        default=6,
        help="number of random input planes (default: 6)",
    )
