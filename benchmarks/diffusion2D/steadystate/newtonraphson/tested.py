
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 2, 2, 3, 2, 20, 30, "solver", "method" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "SQUARE4", "TRIANGLE3" ]:
        yield run, etype, "HYPRE", "0"
        for method in [ "TOTAL_FETI", "HYBRID_FETI" ]:
            yield run, etype, "FETI", method

def run(etype, solver, method):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[7] = solver
    ESPRESOTest.args[8] = method
    ESPRESOTest.run()
    ESPRESOTest.compare(".".join([etype, "emr"]))
    ESPRESOTest.report("espreso.time.xml")