
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.processes = 2
    ESPRESOTest.args = [ "etype", 2, 1, 8, 2, 10, 10 ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "SQUARE8", "TRIANGLE6" ]:
        yield run, etype

def run(etype):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.run()
    ESPRESOTest.compare(".".join([etype, "emr"]))
    ESPRESOTest.report("espreso.time.xml")
