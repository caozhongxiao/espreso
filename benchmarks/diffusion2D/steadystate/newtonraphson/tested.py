
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 2, 2, 3, 2, 20, 30, "method" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "SQUARE4", "TRIANGLE3" ]:
        for method in [ "TOTAL_FETI", "HYBRID_FETI" ]:
            yield run, etype, method

def run(etype, method):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[7] = method
    ESPRESOTest.run()
    ESPRESOTest.compare(".".join([etype, method, "emr"]))
    ESPRESOTest.report("espreso.time.xml")
