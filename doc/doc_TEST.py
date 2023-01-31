# vim: set sts=4 tw=120 :

import unittest
import inspect
import os
import eos

class TestFailedError(BaseException):
    def __init__(self, msg):
        self.msg = msg

class DocTests:
    def check_000_Parameters(self):
        """Check the latex representation of all the parameters"""
        from matplotlib.texmanager import TexManager

        parameters = eos.Parameters.Defaults()
        texmanager = TexManager()

        for section in eos.Parameters.Defaults().sections():
            for group in section:
                for parameter in group:
                    try:
                        texmanager.get_text_width_height_descent(parameter.latex(), fontsize=12)
                    except:
                        raise TestFailedError(f"Cannot compile latex representation of parameter {parameter.name()}")

# Run legacy test cases
tests = DocTests()
for (name, testcase) in inspect.getmembers(tests, predicate=inspect.ismethod):
    print("%s: " % name)
    try:
        testcase()
    except TestFailedError as e:
        print("    failed with error: %s" % e.msg)
        exit(1)
    print("    passed")

# Run new tests
unittest.main()
