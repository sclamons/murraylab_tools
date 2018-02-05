import os
import pytest
import numpy as np
import sys

import murraylab_tools.echo as mt_echo

class TestSourcePlate():

    test_plate = mt_echo.SourcePlate(filename='testplate.dat')

    test_dir = os.path.dirname(os.path.realpath(__file__))

    @pytest.mark.skip(reason="test not yet implemented")
    def test_sanity(self):
        assert 0

    @pytest.mark.skip(reason="test not yet implemented")
    def test_load_from_file(self):
        assert 0

    @pytest.mark.skip(reason="test not yet implemented")
    def test_write_to_file(self):
        assert 0

    @pytest.mark.skip(reason="test not yet implemented")
    def test_request_wells(self):
        assert 0

    @pytest.mark.skip(reason="test not yet implemented")
    def test_increment_position(self):
        assert 0

    def test_load_csv_fails_for_missing_column_names(self):
        with pytest.raises(AssertionError):
            self.test_plate.load_well_definitions(os.path.join(self.test_dir,
                                             'test_def_bad_column_names.csv'))

    def test_load_csv_passes_for_necessary_column_names(self):
        try:
            self.test_plate.load_well_definitions(os.path.join(self.test_dir,
                                             'test_def_good_column_names.csv'))
        except AssertionError:
            assert 0
        assert 1

    def test_get_chemical_location(self):
        loc = self.test_plate.get_location('chemical')
        assert loc == 'A1'

    def test_get_chemical_and_conc_location(self):
        loc = self.test_plate.get_location('chem', 10)
        assert loc == 'A3'
        loc = self.test_plate.get_location('chem', 100)
        assert loc == 'A5'

    def test_get_default_location_of_repeated_material(self):
        loc = self.test_plate.get_location('h20')
        assert loc == 'A6'

    def test_get_one_location_of_repeated_material(self):
        loc = self.test_plate.get_location('h20', i=1)
        assert loc == 'A7'

    def test_get_location_of_repeated_material_saturate_bounds(self):
        loc = self.test_plate.get_location('h20', i=5)
        assert loc == 'A8'
