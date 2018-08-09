
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 4, 1, 1, 3, 1, 1, 4, 4, 4 ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "HEXA20", "TETRA10", "PRISMA15", "PYRAMID13" ]:
        yield run, etype

def run(etype):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.run()
    ESPRESOTest.compare(".".join([etype, "emr"]))
    ESPRESOTest.report("espreso.time.xml")
