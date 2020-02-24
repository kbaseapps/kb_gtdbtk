def assert_exception_correct(got: Exception, expected: Exception):
    assert got.args == expected.args
    assert type(got) == type(expected)