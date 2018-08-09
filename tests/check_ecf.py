import os
import sys
import shutil

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(ESPRESO_TESTS)
sys.path.append(os.path.join(ESPRESO_TESTS, "utils"))

from testing import *
import unittest

class ESPRESOECFChecker(unittest.TestCase):

    espreso = Espreso()

    def ecfchecker(self, path, file):
        self.espreso.ecfchecker(path, { "config": os.path.join(path, file) })


if __name__ == '__main__':

    print "ecf checker is not working"

#     benchmarks = TestCaseCreator.select(os.path.join(ROOT, "benchmarks"))
#     solver = TestCaseCreator.select(os.path.join(ROOT, "tests", "examples", "solver"))
#     input = TestCaseCreator.select(os.path.join(ROOT, "tests", "examples", "input"))
# 
#     for subdirectory in benchmarks:
#         for name, path, file in TestCaseCreator.gather(subdirectory, ".ecf"):
#             TestCaseCreator.create_test(ESPRESOECFChecker, ESPRESOECFChecker.ecfchecker, "benchmarks_" + name, path, file)
# 
#     for subdirectory in solver:
#         for name, path, file in TestCaseCreator.gather(subdirectory, ".ecf"):
#             TestCaseCreator.create_test(ESPRESOECFChecker, ESPRESOECFChecker.ecfchecker, "solver_" + name, path, file)
# 
#     for subdirectory in input:
#         for name, path, file in TestCaseCreator.gather(subdirectory, ".ecf"):
#             TestCaseCreator.create_test(ESPRESOECFChecker, ESPRESOECFChecker.ecfchecker, "input_" + name, path, file)
# 
#     unittest.main()