
import os
from nose.tools import istest
from unittest.case import skip

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 2, 2, 1, 1, 2, 2, 5, 5, 5, "solver" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "HEXA8", "TETRA4", "TETRA10", "PRISMA6", "PRISMA15", "PYRAMID5", "PYRAMID13" ]:
        for solver in [ "FETI", "HYPRE" ]:
            yield run, etype, solver

@istest
@skip("Invalid HEXA20 base functions?")
def BY():
    pass

def run(etype, solver):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[10] = solver
    ESPRESOTest.run()
    ESPRESOTest.compare(".".join([etype, "emr"]))
    ESPRESOTest.report("espreso.time.xml")
