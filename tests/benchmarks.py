
import os
import sys
import shutil

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(ESPRESO_TESTS)
sys.path.append(os.path.join(ESPRESO_TESTS, "utils"))

from testing import *
import unittest

class ESPRESOBenchmarks(unittest.TestCase):

    espreso = Espreso()

    def benchmark(self, path, file, args, report):
        arguments = args.values()
        for key in args:
            arguments[int(key[3:])] = args[key]
        self.espreso.run(
            self.espreso.get_processes(os.path.join(path, file)),
            path,
            { "config": os.path.join(path, file), "ENV::VERBOSE_LEVEL": 0, "ENV::MEASURE_LEVEL": 0, "OUTPUT::RESULTS_STORE_FREQUENCY": "NEVER" },
            arguments
        )

        self.espreso.compare_monitors(
            os.path.join(path, report),
            os.path.join(path, "results", "last", file.replace(".ecf", ".emr"))
        )

        shutil.rmtree(os.path.join(path, "results"))

if __name__ == '__main__':

    benchmarks = TestCaseCreator.select(os.path.join(ROOT, "benchmarks"))

    def run(vars, range, name, path, file):
        def call(args):
            report = ""
            for var in vars:
                report += args[var] + "."
            report = ".".join([ args[var] for var in vars ])
            TestCaseCreator.create_test(ESPRESOBenchmarks, ESPRESOBenchmarks.benchmark, name + "." + report, path, file, args, report + ".emr")

        TestCaseCreator.iterate(call, range)

    for subdirectory in benchmarks:
        for name, path, file in TestCaseCreator.gather(subdirectory, ".ecf"):
            vars, range = Espreso.get_instances(os.path.join(path, file))
            run(vars, range, name, path, file)

    unittest.main()
