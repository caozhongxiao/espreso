
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = []

def teardown():
    ESPRESOTest.clean()

@istest
def defaults():
    ESPRESOTest.run()
    ESPRESOTest.compare("espreso.emr")
    ESPRESOTest.report("espreso.time.xml")
