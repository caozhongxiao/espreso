
import os
from nose.tools import istest

from estest import ESPRESOTest
from unittest.case import skip

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = []

def teardown():
    ESPRESOTest.clean()

@istest
@skip("Not working for more MPIs")
def defaults():
    ESPRESOTest.run()
    ESPRESOTest.compare("espreso.emr")
    ESPRESOTest.report("espreso.time.xml")
