import os
import pytest
import numpy as np

import murraylab_tools.echo as mt_echo

@pytest.fixture()
def test_dir():
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")

def gen_plate():
    fname = 'dplate.dat'
    if os.path.exists(fname):
        os.rm(fname)
    dplate = mt_echo.DestinationPlate(filename=fname)
    return dplate

@pytest.fixture(scope="session")
def unused_plate():
    return gen_plate()

@pytest.fixture(scope="session")
def used_plate():
    dplate = gen_plate()
    dplate.request_wells(10)
    return dplate

def test_request_wells(unused_plate):
    wells = unused_plate.request_wells(5)
    assert np.all(wells == np.array(['A01', 'A02', 'A03', 'A04', 'A05']))

def test_request_wells_from_used_plate(used_plate):
    wells = used_plate.request_wells(5)
    assert np.all(wells == np.array(['A11', 'A12', 'A13', 'A14', 'A15']))

def test_request_too_many_wells(unused_plate):
    with pytest.raises(Exception):
        wells = unused_plate.request_wells(500)


# def test_make_simple_picklist(test_dir):
#     # TODO Expand on assertions
#     dplate = mt_echo.DestinationPlate()
#     splate = mt_echo.SourcePlate()
#     splate.load_well_definitions(os.path.join(test_dir,
#                                 'test_def_good_column_names.csv'))
#     rxns = [
#                 [
#                     ['chem', 5, 10],
#                     ['h2o', 5]
#                 ],
#                 [
#                     ['chem', 5, 100],
#                     ['h2o', 5]
#                 ]
#             ]
#     picklist, _ = dplate.make_picklist(splate, rxns)
#     assert picklist[0][0] == 'A3'
#     assert picklist[2][0] == 'A5'
