def pytest_addoption(parser):
    parser.addoption("--dir", action="store", default="test_out")
