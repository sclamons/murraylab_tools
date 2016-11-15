# TODO: --Clean up example, publish to github.
#       --Add optional user-set source wells
#       --Positive control reaction on 2D dilution series
#

# Written by Victoria Hsiao, August 2016
# Modified by Samuel Clamons, October 2016

# EchoCSVMaker encapsulates several algorithms for writing .csv files that can
# be directly read by the Echo liquid handler, along with .txt files describing
# what needs to go in the source plate.
# Current features include:
#   -- Conversion from a pair of CSV templates
#   -- Automatic setup of an NxN dilution matrix with two inputs.
#   -- Picks from plates defined by CSVs, to wells defined by other CSVs


import numpy as np
import string
import math
import csv
import collections
import os
#import openpyxl as pyxl  # Required for reading excel files

__all__ = ["dna2nM_convert", "echo_round", "MasterMixMaterial",
           "SourcePlate", "EchoSourceMaterial", "Pick", "EchoRun"]

dead_volume = 15000 + 5000 # Dead volume per well in an Echo source plate, in nL
max_volume  = 65000 # Maximum Echoable volume in an Echo source plate, in nL
usable_volume = max_volume - dead_volume # Maximum pipettable volume per well
                                         # in an Echo source plate, in nL

rows = 16           # Number of rows in a 384 well plate
cols = 24           # Number of columns in a 384 well plate

def dna2nM_convert(dnaconc, dnalength):
    '''
    Converts DNA ng/ul to nM

    DNA_conc(ng/ul) * 1e6 (uL/L) * (bp*mol)/660g * 1/dna_length(bp = DNA (nM)
    Double stranded DNA is 660g/(bp*mol), single stranded DNA is 330g/(bp*mol)
    '''
    return (dnaconc*1e6)/(660*dnalength)

def echo_round(x, base = 25):
    '''Round a volume to the nearest 25 nl)'''
    return int(base *round(float(x)/base))

def floatify(element):
    '''
    Convert a string to a float, if possible; otherwise, return the string
    right back. Code ripped off code from Jacob Gabrielson's answer on
    StackOverflow thread 736043.

    Empty strings are converted to "0" to make vector addition work on relevant
    ranges.
    '''
    if not type(element) == str:
        raise ValueError("Why are you trying to floatify a " + type(element) + \
                         "?")
    if element == "":
        return 0.0
    partition=element.partition('.')
    if element.isdigit() or \
       (partition[0].isdigit() and partition[1]=='.' and \
        partition[2].isdigit()) \
        or \
        (partition[0]=='' and partition[1]=='.' and partition[2].isdigit()) \
        or \
        (partition[0].isdigit() and partition[1]=='.' and partition[2]==''):
        return float(element)
    else:
        return element

def process_column_argument(col):
    '''
    Convenience function for processing a column-name argument. Converts
    from either a string ("C") or a zero-indexed int (2) to a zero-indexed
    int. Nonetype arguments are returned as None.
    '''
    if col == None:
        return None
    if type(col) == str:
        if col == "":
            raise ValueError("Column argument can't be an empty string.")
        upper_col = col.upper()
        col_num = 0
        while len(upper_col) > 0:
            col_num *= 26
            col_num += string.ascii_uppercase.find(upper_col[-1]) + 1
            upper_col = upper_col[:-1]
        #Remember, it's zero-indexed!
        col_num -= 1
        return col_num
    elif type(col) == int:
        return col
    else:
        raise TypeError("Column argument must be a string ('C') or a " +\
                        "zero-indexed int (2)")


MasterMixMaterial = collections.namedtuple('MasterMixMaterial',
                                           ['name', 'stock', 'final'])


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
                    populated.
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
            self.type = "384_PP_AP_BP"
            self.rows = 16
            self.cols = 24
        elif SPtype.startswith("384PP"):
            self.rows = 16
            self.cols = 24
            self.type = SPtype
        else:
            raise ValueError("'%s' is not a recognized plate type." % SPtype)
        self.wells_used   = np.zeros((rows, cols), dtype=bool)
        self.current_row  = 0
        self.current_col  = 0

        self.used_well_file = filename
        if self.used_well_file:
            if os.path.isfile(filename):
                self.load_from_file(filename)

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

    def request_wells(self, n_wells):
        '''
        Called when an EchoSourceMaterial wants to get some wells. Returns a
        list of wells for the material, which are marked internally as used.

        Current logic is as follows: Scan from the top-left corner to the
        bottom-right corner, left to right, top to bottom, for a consecutive set
        of wells separated from used wells by a buffer well on the right and a
        buffer welll on the right. Assign the first block run across. If the
        number of wells requested is smaller than the number of wells per row,
        also require that the entire block be able to fit in one row.
        '''
        if n_wells == 0:
            return []
        return_wells = np.empty(shape=(n_wells,), dtype=object)
        while True:
            # Check if this position will work
            if n_wells > self.cols or self.current_col + n_wells <= self.cols:
                flat_idx = self.current_row*self.cols + self.current_col
                if flat_idx + n_wells > self.rows * self.cols:
                    raise Exeption("Source plate %s is out of available wells."\
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
                           self.name)
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


class EchoSourceMaterial():
    '''
    Container class for a single material on a source plate. Keeps track of how
    much material has been used and where they need to go. When requested,
    will commit to wells (assigned by a SourcePlate) and send a list of picks.
    '''
    def __init__(self, name, concentration, length, plate):
        '''
        name -- A string to associate with the name of this material.
        concentration -- the concentration of this material, in ng/uL if this is
                         a double-stranded DNA material, or in nM if this is
                         anything else
        length -- Length of the DNA in bp, if material is double-stranded DNA.
                  Otherwise, length should be 0.
        wells -- A list of source wells for the material.
        plate -- Name of the plate from which this material comes.
        '''
        self.name   = name
        self.length = length
        self.plate  = plate

        # DNA concentration in ng/uL, or other reagent concentration in nM.
        self.concentration = concentration
        if length > 0:
            self.nM = dna2nM_convert(concentration, length)
        else:
            self.nM = self.concentration

        self.wells    = None
        self.picklist = []
        self.total_volume_requested = 0
        self.well_volumes = None

    def request_material(self, destination_well, volume):
        '''
        Request some material be made available for picking.

        volume is the volume of material requested, in nL. If that volume is
        greater than the total amount of material that could possibly be in the
        wells designated for this material, throws an error.

        Updates this material's total requested volume and list of
        destionations.
        '''
        actual_volume = echo_round(volume)
        self.total_volume_requested += actual_volume
        self.picklist.append(Pick(self, None, destination_well, actual_volume))


    def request_picklist(self):
        '''
        Commit to wells, and return a list of finalized picks from this
        material.
        '''
        usable_volume  = max_volume - dead_volume
        n_source_wells = math.ceil(float(self.total_volume_requested) \
                                         / usable_volume)
        self.wells           = self.plate.request_wells(int(n_source_wells))
        self.well_volumes    = np.zeros(len(self.wells))
        self.well_volumes[0] = dead_volume
        self.current_well    = 0
        self.total_volume_requested += dead_volume


        for pick in self.picklist:
            volume_requested = pick.volume # For error-writing purposes.
            available_volume = max_volume - self.well_volumes[self.current_well]

            while pick.volume > available_volume:
                yield Pick(self, self.wells[self.current_well],
                           pick.destination_well, available_volume)
                self.well_volumes[self.current_well] = max_volume
                self.current_well += 1
                self.total_volume_requested += dead_volume
                if self.current_well >= len(self.wells):
                    raise ValueError(("Material %s has been asked to donate " +\
                                      "too much material with a call for " +\
                                      "%d nL") % \
                                    (self.name, volume_requested))
                self.well_volumes[self.current_well] = dead_volume
                pick.volume -= available_volume
                available_volume = usable_volume

            self.well_volumes[self.current_well] += pick.volume
            pick.source_well = self.wells[self.current_well]
            yield pick


class Pick():
    '''
    Container class for making a single Echo pick from some source well to
    a destination well.
    '''
    def __init__(self, source_material, source_well, destination_well, volume):
        self.source_material   = source_material
        self.source_well       = source_well
        self.destination_well  = destination_well
        self.volume            = volume


class EchoRun():
    '''
    Defines and prints an Echo picklist from one of several algorithms.
    '''
    def __init__(self, rxn_vol = 5, DPtype = None, plate = None):
        # The user gives the reaction volume in microliters; store it in nL
        self.rxn_vol = rxn_vol * 1e3

        if DPtype == None:
            self.DPtype = "Nunc_384_black_glassbottom"
        else:
            self.DPtype = DPtype
        if plate == None:
            self.plates = [SourcePlate("Source[1]", "384PP_AQ_BP")]
        else:
            self.plates = [plate]

        self.material_list = dict()
        self.picklist      = []
        self.make_master_mix = False
        self.extract_per_aliquot = 30000 # Volume of extract per aliquot tube
        self.buffer_per_aliquot  = 37000 # Volume of buffer per aliquot tube

    def define_plate(self, SPname, SPtype, DPtype):
        '''
        Convenience function for setting plate parameters.
        '''
        self.SPname = SPname
        self.SPtype = SPtype
        self.DPtype = DPtype

    def build_picklist_from_txtl_setup_excel(self, input_filename):
        '''
        CURRENTLY NONFUNCTIONAL DO NOT USE

        Build an Echo picklist based on a TX-TL setup spreadsheet (v2.1 or
        newer).
        '''
        # Open the workbook and identify all the important sheets
        workbook = pyxl.load_workbook(input_filename)
        recipe_sheet = None
        stock_sheet = None
        layout_sheet = None
        for sheet in workbook:
            if sheet.title == "Recipe":
                recipe = sheet
            elif sheet.title == "Stocks":
                stocks = sheet
            elif sheet.title == "Layout":
                layout = sheet

    def build_picklist_from_txtl_setup_csvs(self, stock_filename,
                                            recipe_filename):
        '''
        Build an Echo picklist based on a pair of CSV documents produced from
        a TX-TL setup spreadsheet (v2.1).

        The stock sheet is a CSV describing materials on the source plate
        (plasmids, usually).

        The recipe sheet is a CSV describing the master mix and what materials
        from the stock sheet go in what destination wells, in what quantity.
        '''
        self.make_master_mix = True
        self.master_mix_materials = []

        ####################
        # Read Input Files #
        ####################

        # Read in stock file
        stock_sheet = np.empty(shape = (12,5), dtype = object)
        with open(stock_filename, 'rU') as stock_file:
            stock_reader = csv.reader(stock_file)
            rownum = -1
            for row in stock_reader:
                rownum += 1
                for colnum in range(len(row)):
                    element = floatify(row[colnum])
                    stock_sheet[rownum, colnum] = element

        # Read in recipe file
        recipe_sheet = np.zeros(shape = (116, 16), dtype = object)
        with open(recipe_filename, 'rU') as recipe_file:
            recipe_reader = csv.reader(recipe_file)
            rownum = -1
            for row in recipe_reader:
                rownum += 1
                for colnum in range(len(row)):
                    element = floatify(row[colnum])
                    if element:
                        recipe_sheet[rownum, colnum] = element

        # Set some magic numbers based on the recipe file
        self.rxn_vol = float(recipe_sheet[11,10])
        self.extract_fraction = float(recipe_sheet[12,2])
        self.buffer_fraction  = 0.75 - self.extract_fraction
        self.mm_excess = float(recipe_sheet[11,12])


        ######################
        # Register Materials #
        ######################

        # Magic numbers here are spreadsheet positions of volumes to add.
        material_total_vols  = [0]*10
        for i in range(len(material_total_vols)):
            for j in range(20,116):
                if recipe_sheet[j, 5+i]:
                    material_total_vols[i] += recipe_sheet[j, 5+i] * 1e3

        # Assign source wells and register materials
        # Register TX-TL master mix
        if not "txtl_mm" in self.material_list:
            self.material_list["txtl_mm"] = EchoSourceMaterial("txtl_mm", 1, 0,
                                                               self.plates[0])
        txtl = self.material_list["txtl_mm"]

        # Register Water
        if not "water" in self.material_list:
            self.material_list['water'] = EchoSourceMaterial("water", 0, 0,
                                                              self.plates[0])
        water = self.material_list["water"]

        # Register other materials
        stocks = []
        for i in range(len(material_total_vols)):
            if material_total_vols[i] == 0 or stock_sheet[i+2,1] == "":
                continue
            material_name          = stock_sheet[i+2,1]
            material_concentration = stock_sheet[i+2, 2]
            material_length        = stock_sheet[i+2, 3]
            if material_name in self.material_list:
                existing_mat = self.material_list[material_name]
                if existing_mat.concentration == material_concentration and \
                   existing_mat.length == material_length:
                   continue
            self.material_list[material_name] = EchoSourceMaterial(\
                                                         material_name,
                                                         material_concentration,
                                                         material_length,
                                                         self.plates[0])
            stocks.append(self.material_list[material_name])

        ##################
        # Register picks #
        ##################

        first_row = 20
        last_row  = 115
        for rownum in range(first_row, last_row + 1):
            # Check to see if there's a name in this row; if not, skip it.
            if recipe_sheet[rownum, 2] == 0:
                continue
            destination_well = recipe_sheet[rownum, 1]
            if destination_well == 0:
                raise ValueError(("Error on row for ID #%d of recipe sheet: " +\
                                 "Must have a destination well.") % \
                                 (rownum - 21))

            # Material picks (magic number warning -- magic numbers define
            # positions of relevant blocks in the recipe sheet)
            for mat_num in range(len(material_total_vols)):
                colnum = mat_num + 5
                volume = recipe_sheet[rownum, colnum] * 1e3
                if volume > 0:
                    if mat_num >= len(stocks):
                        raise ValueError("Tried to pick a material that " + \
                                         "wasn't named in the stock sheet.")
                    source_material = stocks[mat_num]
                    source_material.request_material(destination_well, volume)

            # Water picks (magic number warning -- magic number defines
            # positions of relevant blocks in the recipe sheet)
            volume = recipe_sheet[rownum, 3] * 1e3
            water.request_material(destination_well, volume)

            # Master Mix picks (magic number warning -- magic number defines
            # positions of relevant blocks in the recipe sheet)
            volume = recipe_sheet[rownum, 4] * 1e3
            txtl.request_material(destination_well, volume)

        # Now we know enough to determine the recipe for the master mix
        for i in range(11,17):
            if recipe_sheet[i,4] == None or recipe_sheet[i,4] == 0:
                continue
            name  = recipe_sheet[i,0]
            stock = recipe_sheet[i,1]
            # Final concentration IN THE MASTER MIX, which has to be adjusted
            # to account for the fraction of master mix in the final solution.
            # Assuming they all have the same master mix amount.
            final = recipe_sheet[i,2] / (recipe_sheet[20,4] / self.rxn_vol)
            new_mm_material = MasterMixMaterial(name, stock, final)
            self.master_mix_materials.append(new_mm_material)


    def load_source_plate(self, input_filename, name_col, conc_col, len_col,
                          well_col, plate_col, header = True):
        '''
        Enter new materials from a CSV spreadsheet.

        Args:
            input_filename -- name of the CSV.
            name_col -- Name of the column containing the name of each material,
                        either as a string ("C") or a 0-indexed int (2)
            conc_col -- Name of the column containing the concentration of
                        each material, in ng/uL if the material is dsDNA or in
                        nM or relative concentration otherwise, either as a
                        string ("C") or a 0-indexed int (2)
            len_col  -- Name of the column containing the length of any dsDNA
                        material, either as a string ("C") or a 0-indexed int
                        (2), or None if no such column exists (for sheets
                        containing only non-dsDNA materials)
            well_col -- Name of the column containing the well location of each
                        material, either as a string ("C") or a 0-indexed int
                        (2)
            plate_col -- Name of the column containing the name of the plate
                         the material can be found on, either as a string ("C")
                         or a 0-indexed int (2), or None if no such column
                         exists. Plate will default to Plate[#], where # will
                         increment with each source plate file (without a
                         plate_col) loaded.
            header -- True iff there is a header row. Decides whether or not to
                        skip the first line of each file
        '''
        #####################
        # Process arguments #
        #####################
        name_idx  = process_column_argument(name_col)
        conc_idx  = process_column_argument(conc_col)
        len_idx   = process_column_argument(len_col)
        well_idx  = process_column_argument(well_col)
        plate_idx = process_column_argument(plate_col)

        #############
        # Read file #
        #############
        with open(input_filename, 'rU') as input_file:
            reader = csv.reader(input_file)
            # Skip first row if it's a header
            if header:
                reader.next()
            for row in reader:
                name          = row[name_idx]
                concentration = floatify(row[conc_idx])
                well          = row[well_idx]
                if len_idx != None:
                    length = int(floatify(row[len_idx]))
                else:
                    length = 0
                plate_name = row[plate_idx]
                plate = None
                for p in self.plates:
                    if p.name == plate_name:
                        plate = p
                        break
                if not plate:
                    plate = SourcePlate(SPname = plate_name)
                new_material = EchoSourceMaterial(name, concentration,
                                                  length, plate)
                self.material_list[name] = new_material

    def build_picklist_from_association_spreadsheet(self, input_filename,
                                                    well_column, header = True):
        '''
        Make an Echo picklist based on an association spreadsheet, which is a
        CSV file where 1) each line is a reaction and 2) lines contains
        alternating columns of material names and final concentrations (in nM),
        starting with material name.

        Args:
            input_filename -- name of the CSV.
            well_col -- Name of the column containing the destination well of
                        each reaction, either as a string ("C") or a 0-indexed
                        int (2)
            header -- True iff there is a header row. Decides whether or not to
                        skip the first line of each file
        '''
        well_idx = process_column_argument(well_column)
        with open(input_filename, 'rU') as input_file:
            reader = csv.reader(input_file)
            # Skip the first row if it's a header
            if header:
                reader.next()
            for row in reader:
                if well_idx >= len(row):
                    raise ValueError("Well column out of bounds for row '%s'" %\
                                     str(row))
                well        = row[well_idx]
                i           = 0
                is_name_col = True
                while i < len(row):
                    # Ignore it if it's the well column.
                    if i == well_idx:
                        i += 1
                        continue
                    # First column of each pair is the material
                    if is_name_col:
                        source_material = self.material_list[row[i]]
                        is_name_col     = False
                    # Second column of each pair is the final concentration of
                    # material
                    else:
                        final_conc = float(row[i])
                        well       = row[well_idx]
                        volume = self.rxn_vol * final_conc / source_material.nM
                        source_material.request_material(well, volume)
                        is_name_col = True
                    i += 1


    def build_dilution_series(self, dna1, dna2, dna1_final, dna2_final,
                              first_well, extract_fraction = 0.24):
        '''Calculate TXTL volumes and automatically generate a picklist for
        an NxN dilution matrix with two inputs.

        Args:
            dna1, dna2 -- First and second materials to serialy dilute
                          (EchoSourceMaterials)
            dna1_final, dna2_final -- Lists (or numpy array) defining the final
                                      concentrations of each material.
            first_well -- The upper-left-most corner of the block on the
                          destination plate you want to spit into.
        '''
        dna1.plate = self.plates[0]
        dna2.plate = self.plates[0]
        self.make_master_mix      = True
        self.master_mix_materials = []
        self.extract_fraction     = extract_fraction
        self.buffer_fraction      = 0.75 - self.extract_fraction
        self.master_mix_materials.append(MasterMixMaterial("Extract", 1,
                                                    self.extract_fraction/0.75))
        self.master_mix_materials.append(MasterMixMaterial("Buffer", 1,
                                                    self.buffer_fraction/0.75))
        self.mm_excess = 1.1

        self.material_list["dna1"] = dna1
        self.material_list["dna2"] = dna2
        '''if dna1.name == "pos_ctrl":
            self.material_list["pos_ctrl"] = dna1
        else:
            self.material_list = EchoSourceMaterial("pos_ctrl", 19, 3202, "")'''

        final_volume = self.rxn_vol
        # matrix size
        n_dna1 = len(dna1_final)
        n_dna2 = len(dna2_final)

        # Final volume of mastermix should be 75% of total volume
        # Buffer and extract are calculated separately so that it can do a final
        # calculation on how many tubes you need.
        txtlMM   = 0.75 * final_volume
        num_rxns = n_dna1 * n_dna2 + 2
        if not "txtl_mm" in self.material_list:
            self.material_list["txtl_mm"] = EchoSourceMaterial("txtl_mm", 1, 0,
                                                               self.plates[0])
        txtl = self.material_list["txtl_mm"]

        # Add water as a material (if it's not already there). Assume that a row
        # of water will do fine.
        if not "water" in self.material_list:
            self.material_list["water"] = EchoSourceMaterial("water", 0, 0,
                                                             self.plates[0])
        water = self.material_list["water"]

        # Fill in matrix with picks.
        first_row = string.ascii_uppercase.find(first_well[0])
        first_col = int(first_well[1:])
        for i in range(n_dna1):
            for j in range(n_dna2):
                # Well name
                destination = string.ascii_uppercase[first_row + i] + \
                              str(first_col + j)
                # DNA picks
                dna1_pick_vol = dna1_final[i]*(final_volume/dna1.nM)
                dna1.request_material(destination, dna1_pick_vol)

                dna2_pick_vol = dna2_final[j]*(final_volume/dna2.nM)
                dna2.request_material(destination, dna2_pick_vol)
                # TX-TL Master Mix pick
                txtl.request_material(destination, txtlMM)
                # Water pick
                water_vol = self.rxn_vol - dna1_pick_vol - dna2_pick_vol - \
                            txtlMM
                if water_vol < 0:
                    raise ValueError(("Well %s is overloaded! %s contains "+\
                                     "%.2f uL of %s, %.2f uL of %s, and " +\
                                     "%.2f uL of Master Mix.") % \
                                     (destination, destination, dna1_pick_vol,
                                      dna1.name, dna2_pick_vol, dna2.name,
                                      txtlMM))
                water.request_material(destination, water_vol)

        # Add positive control....

        # and negative control.
        neg_ctrl_well = string.ascii_uppercase[first_row + n_dna1] \
                        + str(first_col)
        txtl.request_material(neg_ctrl_well, txtlMM)
        water.request_material(neg_ctrl_well, self.rxn_vol - txtlMM)

    def generate_picklist(self):
        for mat in self.material_list.values():
            picks = mat.request_picklist()
            for pick in picks:
                yield pick

    def write_picklist(self, outputname):
        '''
        Write this EchoCSVMaker's protocol to an Echo picklist, and print any
        other necessary directions for the user.
        '''
        # Write picklist.
        with open((outputname + '_EchoInput.csv'), 'w') as outcsv:
            writer = csv.writer(outcsv)

            # Write header
            writer.writerow(["Source Plate","Source Plate Type","Source Well",
                             "Sample ID","Destination Plate Name",
                             "Destination Well","Transfer Volume",
                             "Sample Comment"])

            # Write picks
            for pick in self.generate_picklist():
                if pick.source_material.name == "txtl_mm" \
                    or pick.source_material.name == "water":
                    comment = ""
                else:
                    comment = "Actual concentration: %.2f nM" % \
                              (pick.source_material.nM*pick.volume/self.rxn_vol)
                plate = pick.source_material.plate
                row = [plate.name, plate.type, pick.source_well,
                       pick.source_material.name, self.DPtype,
                       pick.destination_well, pick.volume, comment]
                writer.writerow(row)

        # Write comment file
        with open((outputname+'_experiment_overview.txt'), 'w') as text_file:
            all_destination_wells = map(lambda pick:pick.destination_well,
                                        self.picklist)
            rxn_num = len(set(all_destination_wells))

            text_file.write("Materials used:")
            for material in self.material_list.values():
                text_file.write("\n%s:" % material.name)
                if not material.name == "txtl_mm":
                    text_file.write("\n\tstock concentration: %.2fnM" % \
                                    material.nM)
                if material.length > 0:
                    text_file.write(" (%.2f ng/uL)" % material.concentration)
                text_file.write("\n\ttotal volume: %.2f uL" % \
                                (material.total_volume_requested / 1000.0))
                if material.name == "txtl_mm":
                    total_mm = material.total_volume_requested
                    extract_tubes = math.ceil(total_mm * \
                                          (self.extract_fraction/0.75) / \
                                          self.extract_per_aliquot)
                    buffer_tubes = math.ceil(total_mm * \
                                          (self.buffer_fraction/0.75) / \
                                          self.extract_per_aliquot)
                    text_file.write("\n\tTubes of extract needed: %d" % \
                                    extract_tubes)
                    text_file.write("\n\tTubes of buffer needed: %d" % \
                                    buffer_tubes)
                    if self.make_master_mix:
                        text_file.write("\n\tMaster Mix:")
                        for mm_material in self.master_mix_materials:
                            vol_mm_material = total_mm * self.mm_excess * \
                                              mm_material.final / \
                                              mm_material.stock / 1000.0
                            text_file.write("\n\t\t%.2fuL %s" % \
                                            (vol_mm_material, mm_material.name))

            # Explicit loading instructions
            text_file.write("\n\nInstructions:")

            for material in self.material_list.values():
                vol_list = material.well_volumes
                if material.current_well != 0:
                    text_file.write("\n\t%.2f uL of %s in wells " % \
                                    (max_volume/1000.0, material.name))
                    first_entry = True
                    for full_well in material.wells[:material.current_well]:
                        if first_entry:
                            first_entry = False
                        else:
                            text_file.write(", ")
                        text_file.write(full_well)
                text_file.write("\n\t%.2f uL of %s in well %s" % \
                        (material.well_volumes[material.current_well] / 1000.0,
                         material.name,
                         material.wells[material.current_well]))

            # Make the plates write out their usage.
            for plate in self.plates:
                plate.write_to_file()

