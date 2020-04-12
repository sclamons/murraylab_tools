import os
import pytest
import numpy as np
import sys
import string

import murraylab_tools.echo as mt_echo

class TestSourcePlate():

    test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
    test_plate = mt_echo.SourcePlate(filename=os.path.join(test_dir,
                                     'testplate.dat'))

    def test_io(self):
        '''
        Checks that a plate file can be loaded and saved correctly.
        '''
        plate_in_file  = os.path.join(self.test_dir,
                                      "half_full_source_plate.dat")
        plate_out_file = os.path.join(self.test_dir,
                                      "half_full_source_plate_output.dat")
        plate = mt_echo.SourcePlate(filename = plate_in_file)

        reference_contents = dict()
        used_locations = ["A1", "A2", "A3", "B4", "B5", "B6"]
        for well in used_locations:
            reference_contents[well] = True
        well_used = plate.wells_used
        assert reference_contents == plate.wells_used

        plate.used_well_file = plate_out_file
        plate.write_to_file()
        with open(plate_in_file, 'r') as infile:
            with open(plate_out_file, 'r') as outfile:
                for in_line in infile:
                    out_line = next(outfile)
                    assert in_line.strip() == out_line.strip()

    def test_file_creation(self, tmpdir):
        '''
        Tests the case where no file exists yet for the plate.
        '''
        plate_file_name = tmpdir.mkdir("sub").join("hello.txt")
        plate = mt_echo.SourcePlate(filename = plate_file_name)
        plate.request_wells(5)
        plate.write_to_file()

        contents = ["A1", "A2", "A3", "A4", "A5", "A6"]
        with open(plate_file_name, 'r') as plate_file:
            for i, line in enumerate(plate_file):
                assert contents[i] == line.strip()

    @pytest.mark.skip(reason="test not yet implemented")
    def test_request_wells(self):
        assert 0

    @pytest.mark.skip(reason="test not yet implemented")
    def test_request_wells_skip_line(self):
        assert 0

    def test_increment_negative(self):
        '''
        Checks that exception is properly thrown when negative increments used.
        '''
        plate_in_file = os.path.join(self.test_dir, "deleteme.dat")
        plate = mt_echo.SourcePlate(filename = plate_in_file)
        with pytest.raises(ValueError):
            plate.increment_position(-1)

    def test_increment_position_1(self):
        '''
        Checks that incrementing by 1 works properly.
        '''
        plate_in_file = os.path.join(self.test_dir, "deleteme.dat")

        plate_types = ["384PP", "6"]
        target_locations = [
            [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(0,7),(0,8),(0,9),(0,10),
             (0,11),(0,12),(0,13),(0,14),(0,15),(0,16),(0,17),(0,18),(0,19),
             (0,20),(0,21),(0,22),(0,23),(1,0),(1,1),(1,2),(1,3),(1,4),(1,5),
             (1,6),(1,7),(1,8),(1,9),(1,10),(1,11),(1,12),(1,13),(1,14),(1,15),
             (1,16),(1,17),(1,18),(1,19),(1,20),(1,21),(1,22),(1,23),(2,0),
             (2,1),(2,2),(2,3),(2,4),(2,5),(2,6),(2,7),(2,8),(2,9),(2,10),
             (2,11),(2,12),(2,13),(2,14),(2,15),(2,16),(2,17),(2,18),(2,19),
             (2,20),(2,21),(2,22),(2,23),(3,0),(3,1),(3,2),(3,3),(3,4),(3,5),
             (3,6),(3,7),(3,8),(3,9),(3,10),(3,11),(3,12),(3,13),(3,14),(3,15),
             (3,16),(3,17),(3,18),(3,19),(3,20),(3,21),(3,22),(3,23),(4,0),
             (4,1),(4,2),(4,3),(4,4),(4,5),(4,6),(4,7),(4,8),(4,9),(4,10),
             (4,11),(4,12),(4,13),(4,14),(4,15),(4,16),(4,17),(4,18),(4,19),
             (4,20),(4,21),(4,22),(4,23),(5,0),(5,1),(5,2),(5,3),(5,4),(5,5),
             (5,6),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),(5,13),(5,14),(5,15),
             (5,16),(5,17),(5,18),(5,19),(5,20),(5,21),(5,22),(5,23),(6,0),
             (6,1),(6,2),(6,3),(6,4),(6,5),(6,6),(6,7),(6,8),(6,9),(6,10),
             (6,11),(6,12),(6,13),(6,14),(6,15),(6,16),(6,17),(6,18),(6,19),
             (6,20),(6,21),(6,22),(6,23),(7,0),(7,1),(7,2),(7,3),(7,4),(7,5),
             (7,6),(7,7),(7,8),(7,9),(7,10),(7,11),(7,12),(7,13),(7,14),(7,15),
             (7,16),(7,17),(7,18),(7,19),(7,20),(7,21),(7,22),(7,23),(8,0),(8,1),
             (8,2),(8,3),(8,4),(8,5),(8,6),(8,7),(8,8),(8,9),(8,10),(8,11),
             (8,12),(8,13),(8,14),(8,15),(8,16),(8,17),(8,18),(8,19),(8,20),
             (8,21),(8,22),(8,23),(9,0),(9,1),(9,2),(9,3),(9,4),(9,5),(9,6),
             (9,7),(9,8),(9,9),(9,10),(9,11),(9,12),(9,13),(9,14),(9,15),(9,16),
             (9,17),(9,18),(9,19),(9,20),(9,21),(9,22),(9,23),(10,0),(10,1),
             (10,2),(10,3),(10,4),(10,5),(10,6),(10,7),(10,8),(10,9),(10,10),
             (10,11),(10,12),(10,13),(10,14),(10,15),(10,16),(10,17),(10,18),
             (10,19),(10,20),(10,21),(10,22),(10,23),(11,0),(11,1),(11,2),
             (11,3),(11,4),(11,5),(11,6),(11,7),(11,8),(11,9),(11,10),(11,11),
             (11,12),(11,13),(11,14),(11,15),(11,16),(11,17),(11,18),(11,19),
             (11,20),(11,21),(11,22),(11,23),(12,0),(12,1),(12,2),(12,3),(12,4),
             (12,5),(12,6),(12,7),(12,8),(12,9),(12,10),(12,11),(12,12),(12,13),
             (12,14),(12,15),(12,16),(12,17),(12,18),(12,19),(12,20),(12,21),
             (12,22),(12,23),(13,0),(13,1),(13,2),(13,3),(13,4),(13,5),(13,6),
             (13,7),(13,8),(13,9),(13,10),(13,11),(13,12),(13,13),(13,14),
             (13,15),(13,16),(13,17),(13,18),(13,19),(13,20),(13,21),(13,22),
             (13,23),(14,0),(14,1),(14,2),(14,3),(14,4),(14,5),(14,6),(14,7),
             (14,8),(14,9),(14,10),(14,11),(14,12),(14,13),(14,14),(14,15),
             (14,16),(14,17),(14,18),(14,19),(14,20),(14,21),(14,22),(14,23),
             (15,0),(15,1),(15,2),(15,3),(15,4),(15,5),(15,6),(15,7),(15,8),
             (15,9),(15,10),(15,11),(15,12),(15,13),(15,14),(15,15),(15,16),
             (15,17),(15,18),(15,19),(15,20),(15,21),(15,22),(15,23)],
            [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)]]
        for i, plate_type in enumerate(plate_types):
            plate = mt_echo.SourcePlate(filename = plate_in_file,
                                        SPtype = plate_type)
            locs = target_locations[i]
            assert plate.current_row == locs[0][0]
            assert plate.current_col == locs[0][1]
            for i in range(1, len(locs)):
                plate.increment_position(1)
                row, col = locs[i]
                assert plate.current_row == row
                assert plate.current_col == col
            with pytest.raises(Exception):
                plate.increment_position(1)

    def test_increment_position_5(self):
        '''
        Checks that incrementing by 1 works properly.
        '''
        plate_in_file = os.path.join(self.test_dir, "deleteme.dat")

        plate_types = ["384PP", "6"]
        target_locations = [
            [(0,0),(0,5),(0,10),(0,15),(0,20),(1,1),(1,6),(1,11),(1,16),(1,21),
            (2,2),(2,7),(2,12),(2,17),(2,22),(3,3),(3,8),(3,13),(3,18),(3,23),
            (4,4),(4,9),(4,14),(4,19),(5,0),(5,5),(5,10),(5,15),(5,20),(6,1),
            (6,6),(6,11),(6,16),(6,21),(7,2),(7,7),(7,12),(7,17),(7,22),(8,3),
            (8,8),(8,13),(8,18),(8,23),(9,4),(9,9),(9,14),(9,19),(10,0),(10,5),
            (10,10),(10,15),(10,20),(11,1),(11,6),(11,11),(11,16),(11,21),
            (12,2),(12,7),(12,12),(12,17),(12,22),(13,3),(13,8),(13,13),(13,18),
            (13,23),(14,4),(14,9),(14,14),(14,19),(15,0),(15,5),(15,10),(15,15),
            (15,20)],
            [(0,0), (1, 2)]]
        for i, plate_type in enumerate(plate_types):
            plate = mt_echo.SourcePlate(filename = plate_in_file,
                                        SPtype = plate_type)
            locs = target_locations[i]
            assert plate.current_row == locs[0][0]
            assert plate.current_col == locs[0][1]
            for i in range(1, len(locs)):
                plate.increment_position(5)
                row, col = locs[i]
                assert plate.current_row == row
                assert plate.current_col == col
            with pytest.raises(Exception):
                plate.increment_position(5)

    # def test_load_csv_fails_for_missing_column_names(self):
    #     with pytest.raises(AssertionError):
    #         self.test_plate.load_well_definitions(os.path.join(self.test_dir,
    #                                          'test_def_bad_column_names.csv'))

    # def test_load_csv_passes_for_necessary_column_names(self):
    #     try:
    #         self.test_plate.load_well_definitions(os.path.join(self.test_dir,
    #                                          'test_def_good_column_names.csv'))
    #     except AssertionError:
    #         assert 0
    #     assert 1

    # def test_get_chemical_location(self):
    #     loc = self.test_plate.get_location('chemical')
    #     assert loc == 'A1'

    # def test_get_chemical_and_conc_location(self):
    #     loc = self.test_plate.get_location('chem', 10)
    #     assert loc == 'A3'
    #     loc = self.test_plate.get_location('chem', 100)
    #     assert loc == 'A5'

    # def test_get_default_location_of_repeated_material(self):
    #     loc = self.test_plate.get_location('h20')
    #     assert loc == 'A6'

    # def test_get_one_location_of_repeated_material(self):
    #     loc = self.test_plate.get_location('h20', i=1)
    #     assert loc == 'A7'

    # def test_get_location_of_repeated_material_saturate_bounds(self):
    #     loc = self.test_plate.get_location('h20', i=5)
    #     assert loc == 'A8'
