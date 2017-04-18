import numpy as np
import pandas as pd
import string
import math
import csv
import collections
import os
import warnings
from echo_functions import *

class SourcePlate():
    '''
    One Echo source plate. Responsible for allocating wells for
    EchoSourceMaterials.
    '''
    def __init__(self, SPname = None, SPtype = None, filename = None):
        '''
        SPname -- A(n arbitrary) string representing this plate. Default is
                    "Plate[1]"
        SPtype -- A string denoting the type of plate. Plate types are
                    restricted and determine the size of the plate. Default is
                    "384PP_AQ_BP"
        filename -- Name of a file holding a list of used wells for this plate.
                    If that file doesn't exist yet, it will be created and
                    populated when write_to_file is called.
        '''
        if SPname == None:
            self.name = "Source[1]"
        else:
            name_is_number = True
            for c in SPname:
                if c not in string.digits:
                    name_is_number = False
                    break
            if name_is_number:
                self.name = "Source[%s]" % SPname
            else:
                self.name = SPname
        if SPtype == None:
            self.type = "384PP_AQ_BP"
            self.rows = 16
            self.cols = 24
        elif SPtype.startswith("384PP"):
            self.rows = 16
            self.cols = 24
            self.type = SPtype
        elif SPtype.startswith("6"):
            self.rows = 2
            self.cols = 3
            self.type = SPtype
        else:
            raise ValueError("'%s' is not a recognized plate type." % SPtype)
        self.wells_used   = np.zeros((self.rows, self.cols), dtype=bool)
        self.current_row  = 0
        self.current_col  = 0

        self.used_well_file = filename
        if self.used_well_file:
            if os.path.isfile(filename):
                self.load_from_file(filename)

    def load_well_definitions(self, filename):
        """
        Import source plate definitions from CSV.

        Arguments:
            self: object
            filename: filname of CSV
        Returns:
            Nothing
        Raises:
            AssertionError when a necessary column is missing
        """
        in_plate = pd.read_csv(filename)
        necessary_cols = ['Location', 'Name', 'Concentration', 'Plate']
        for col in necessary_cols:
            assert col in in_plate.columns
        self.plate = in_plate

    def get_location(self, name, conc=None, i=0):
        """
        Get well location of material.

        Arguments:
            self: object
            name: name of material desired
            conc: select for concentration when available. Default is to select
                    the first entry.
        Returns:
            Well location, as a string
        Raises:
            AttributeError when named material is not found in source plate
        """
        def get_from_values(values):
            if type(values) is np.ndarray:
                n = min(len(values)-1, i)
                return values[n]
            return values

        names = self.plate.Name.values
        if name not in names:
            raise AttributeError(('No material named {} in source ' + \
                                  'plate').format(name))
        grouped_by_name = self.plate.groupby('Name')
        group = grouped_by_name.get_group(name)
        if conc is None:
            return get_from_values(group.Location.values)
        else:
            grouped_by_conc = group.groupby('Concentration')
            if conc in grouped_by_conc.groups:
                c_group = grouped_by_conc.get_group(conc)
                return get_from_values(c_group.Location.values)

    def load_from_file(self, filename):
        '''
        Reads which wells have been used from a data file. The well-use file
        lists wells that have been used, with one line per well used.
        '''
        with open(filename, 'r') as infile:
            for line in infile:
                line = line.strip()
                if line == "":
                    continue
                col_name = ""
                while line[-1] in string.digits:
                    col_name = line[-1] + col_name
                    line = line[:-1]
                col_num = int(col_name) - 1
                row_num = string.ascii_uppercase.find(line.upper())
                self.wells_used[row_num, col_num] = True

    def write_to_file(self):
        '''
        Writes which wells have been used to a data file. The well-use file
        lists wells that have been used, with one line per well used.

        Note that if no filename is given to the source plate, no well-use file
        will be written and this function will return silently.
        '''
        if not self.used_well_file:
            return
        used_well_indices = self.wells_used.nonzero()
        used_well_rows = used_well_indices[0]
        used_well_cols = used_well_indices[1]
        with open(self.used_well_file, 'w+') as outfile:
            for row_num, col_num in zip(used_well_rows, used_well_cols):
                row_string = string.ascii_uppercase[row_num]
                col_string = str(col_num + 1)
                outfile.write(row_string + col_string + "\n")

    def request_wells(self, n_wells, name = "None"):
        '''
        Called when an EchoSourceMaterial wants to get some wells. Returns a
        list of wells for the material, which are marked internally as used.

        Current logic is as follows: Scan from the top-left corner to the
        bottom-right corner, left to right, top to bottom, for a consecutive set
        of wells separated from used wells by a buffer well on the right and a
        buffer welll on the right. Assign the first block run across. If the
        number of wells requested is smaller than the number of wells per row,
        also require that the entire block be able to fit in one row.

        Alternatively, if the name of a material is passed, and the
        SourcePlate object knows where that material is stored, it can assign
        that well instead.
        '''
        if n_wells == 0:
            return []
        return_wells = np.empty(shape=(n_wells,), dtype=object)
        while True:
            # Check if this position will work
            if n_wells > self.cols or self.current_col + n_wells <= self.cols:
                flat_idx = self.current_row*self.cols + self.current_col
                if flat_idx + n_wells > self.rows * self.cols:
                    raise Exception("Source plate %s is out of available wells."\
                                   % self.name)
                block = self.wells_used.ravel()[flat_idx:flat_idx + n_wells]
                # If it will, return that block and mark it used
                if True not in block:
                    block[:]    = True
                    return_list = np.empty(n_wells, dtype=object)
                    for i in range(flat_idx, flat_idx + n_wells):
                        row = int(i / self.cols)
                        col = i % self.cols
                        return_list[i-flat_idx] = string.ascii_uppercase[row] +\
                                         ("%02d" % (col+1))
                    # Mark wells as used.
                    self.increment_position(n_wells)
                    # Move on, leaving an extra well as a buffer (unless it
                    # just crossed to a new line)
                    if self.current_col != 0:
                        self.wells_used[self.current_row, self.current_col] = \
                                                                        True
                        self.increment_position(1)
                    return return_list
            # Otherwise, increment
            self.increment_position(1)

    def increment_position(self, n):
        '''
        For internal use. Increments the plate's current position by n.
        '''
        self.current_col += n
        while self.current_col >= self.cols:
            self.current_col -= self.cols
            self.current_row += 1
        if self.current_row >= self.rows:
            raise Exception("Source plate %s is out of available wells." % \
                            self.name)