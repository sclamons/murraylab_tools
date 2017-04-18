# coding=utf-8

class EchoRun():
    '''
    Defines and prints an Echo picklist from one of several algorithms.

    Parameters:
        -- rxn_vol: Volume of a single TX-TL reaction, in nL. Default 5000.
        -- DPtype: Destination plate type. Should be a string recognizable by
                    the Echo Plate Reformat software, though that is not
                    enforced in this code. Default "Nunc_384_black_glassbottom"
        --plate: Source plate. Should be a SourcePlate object. Required for
                    TX-TL setup experiments, but not for association spreadsheet
                    experiments.
    '''
    def __init__(self, rxn_vol = 5000, DPtype = None, plate = None,
                 master_mix = None):
        # The user gives the reaction volume in microliters; store it in nL
        self.rxn_vol = rxn_vol

        if DPtype == None:
            self.DPtype = "Nunc_384_black_glassbottom"
        else:
            self.DPtype = DPtype
        if '384' in self.DPtype:
            self.rows = 16
            self.cols = 24
        elif '96' in self.DPtype:
            self.rows = 8
            self.cols = 12
        else:
            raise ValueError(("Unrecognized plate type %s. If this is a " + \
                              "valid plate type, contact a murraylab_tools " +\
                              "developer to have it added.") % self.DPtype)
        if plate == None:
            self.plates = [SourcePlate("Source[1]", "384PP_AQ_BP")]
        else:
            self.plates = [plate]


        self.material_list   = dict()
        self.reactions       = dict()
        self.picklist        = []
        self.add_master_mix(master_mix) # This will set make_master_mix

    def define_plate(self, SPname, SPtype, DPtype):
        '''
        Convenience function for setting plate parameters.
        '''
        self.SPname = SPname
        self.SPtype = SPtype
        self.DPtype = DPtype

    def add_master_mix(self, master_mix):
        '''
        Add a master mix to all reactions.

        master_mix: A MasterMix object describing the new master mix.
        '''
        self.make_master_mix = master_mix != None
        if self.make_master_mix:
            self.material_list['txtl_mm'] = master_mix

    def remove_master_mix(self):
        '''
        Blanks the current master mix.
        '''
        self.material_list['txtl_mm'] = None
        self.make_master_mix          = False

    def add_material(self, material):
        '''
        Add a material to the materials list for this EchoRun, if an identical
        material is not already in the list. Attempting to add a material with
        the same name but different properties as another material already in
        this object's material list will raise a ValueError.

        Doesn't handle objects of class MasterMix or other subclasses of
        EchoSourceMaterial. Add these to the master mix manually, i.e.,

        self.master_mix.append(material)

        after whatever check is required to avoid duplications.

        Returns 0 if the material was added successfully; returns 1 if an
        identical material was already in this EchoRun object's material list,
        so the new material was not added.
        '''
        if material.name in self.material_list.keys():
            prior_mat = self.material_list[material.name]
            if prior_mat.name == material.name \
               and prior_mat.concentration == material.concentration \
               and prior_mat.length == material.length \
               and prior_mat.plate  == material.plate:
                return 1
            else:
                raise ValueError("Tried to add material " + material.name +    \
                    " with concentration " + str(material.concentration) +     \
                    ", length " + str(material.length) + ", and plate " +      \
                     str(material.plate) + "; that material already exists " + \
                     "with concentration " + str(prior_mat.concentration) +    \
                     ", length " + str(prior_mat.length) + ", and plate " +    \
                     str(prior_mat.plate) + ".")
        else:
            self.material_list[material.name] = material
            return 0

    def build_picklist_from_txtl_setup_excel(self, input_filename):
        '''
        CURRENTLY NONFUNCTIONAL DO NOT USE

        Build an Echo picklist based on a TX-TL setup spreadsheet (v2.1 or
        newer).
        '''
        raise NotImplementedError("'build_picklist_from_txtl_setup_excel' " + \
                                  "hasn't been implemented yet. Use " + \
                                  "'build_picklist_from_txtl_setup_csvs'.")
        # Open the workbook and identify all the important sheets
        # workbook = pyxl.load_workbook(input_filename)
        # recipe_sheet = None
        # stock_sheet = None
        # layout_sheet = None
        # for sheet in workbook:
        #     if sheet.title == "Recipe":
        #         recipe = sheet
        #     elif sheet.title == "Stocks":
        #         stocks = sheet
        #     elif sheet.title == "Layout":
        #         layout = sheet

    def build_picklist_from_txtl_setup_csvs(self, stock_filename,
                                            recipe_filename):
        '''
        Build an Echo picklist based on a pair of CSV documents produced from
        a TX-TL setup spreadsheet (v2.1).

        The stock sheet is a CSV describing materials on the source plate
        (plasmids, usually).

        The recipe sheet is a CSV describing the master mix and what materials
        from the stock sheet go in what destination wells, in what quantity.

        This function will overwrite any previous master mix defined for this
        EchoRun object, since the master mix is fully defined in the recipe
        sheet.
        '''
        self.make_master_mix = True

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
        recipe_sheet = np.zeros(shape = (384+2, 16), dtype = object)
        with open(recipe_filename, 'rU') as recipe_file:
            recipe_reader = csv.reader(recipe_file)
            rownum = -1
            for row in recipe_reader:
                rownum += 1
                if rownum >= recipe_sheet.shape[0]:
                    print("Warning -- You are trying to add to more than 384 "+\
                          "wells in the destination plate. Extra wells will " +\
                          "be clipped.")
                    break
                for colnum in range(len(row)):
                    element = floatify(row[colnum])
                    if element:
                        recipe_sheet[rownum, colnum] = element

        # Set some magic numbers based on the recipe file
        self.rxn_vol = float(recipe_sheet[11,10]) * 1e3
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
            self.material_list['txtl_mm'] = MasterMix(self.plates[0],
                                extract_fraction = self.extract_fraction,
                                mm_excess = self.mm_excess,
                                add_txtl = False,
                                rxn_vol = self.rxn_vol)
        txtl = self.material_list["txtl_mm"]

        # Register Water
        self.add_material(EchoSourceMaterial("water", 0, 0, self.plates[0]))
        water = self.material_list["water"]

        # Register other materials
        stocks = []
        for i in range(len(material_total_vols)):
            if material_total_vols[i] == 0 or stock_sheet[i+2,1] == "":
                continue
            material_name          = stock_sheet[i+2,1]
            material_concentration = stock_sheet[i+2, 2]
            material_length        = stock_sheet[i+2, 3]
            new_material           = EchoSourceMaterial(material_name,
                                                        material_concentration,
                                                        material_length,
                                                        self.plates[0])
            is_duplicate_material = self.add_material(new_material)
            if not is_duplicate_material:
                stocks.append(new_material)

        ##################
        # Register picks #
        ##################

        first_row = 20
        last_row  = 115
        n_rxns    = 0
        for rownum in range(first_row, last_row + 1):
            # Check to see if there's a name in this row; if not, skip it.
            if recipe_sheet[rownum, 2] == 0:
                continue
            n_rxns += 1
            well = recipe_sheet[rownum, 1]
            if well == 0:
                raise ValueError(("Error on row for ID #%d of recipe sheet: " +\
                                 "Must have a destination well.") % \
                                 (rownum - 21))
            if self.reactions[well]:
                raise ValueError("Well %s already has a reaction!" \
                                 % well)
            self.reactions[well] = WellReaction(self.rxn_vol, well)

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
                    self.reactions[well].add_volume_of_material(source_material,
                                                                volume)

            # Water picks (magic number warning -- magic number defines
            # positions of relevant blocks in the recipe sheet)
            volume = recipe_sheet[rownum, 3] * 1e3
            self.reactions[well].add_volume_of_material(water, volume)

            # Master Mix picks (magic number warning -- magic number defines
            # positions of relevant blocks in the recipe sheet)
            volume = recipe_sheet[rownum, 4] * 1e3
            self.reactions[well].add_volume_of_material(txtl, volume)

        for i in range(11,17):
            if recipe_sheet[i,4] == None or recipe_sheet[i,4] == 0:
                continue
            name  = recipe_sheet[i,0]
            stock = recipe_sheet[i,1]
            final = recipe_sheet[i,2]
            txtl.add_material(name, stock, final)


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
                next(reader)
            for row in reader:
                name          = row[name_idx]
                concentration = floatify(row[conc_idx])
                well          = row[well_idx]
                if len_idx != None:
                    length = int(floatify(row[len_idx]))
                else:
                    length = 0
                if plate_idx == None:
                    plate_name = "1"
                else:
                    plate_name = row[plate_idx]
                plate = None
                for p in self.plates:
                    if p.name == plate_name:
                        plate = p
                        break
                if not plate:
                    plate = SourcePlate(SPname = plate_name)
                material = EchoSourceMaterial(name, concentration,
                                              length, plate)
                self.add_material(material)
                if self.material_list[name].wells == None:
                    self.wells = [well]
                else:
                    self.material_list[name].wells.append(well)


    def build_picklist_from_association_spreadsheet(self, input_filename,
                                                    well_column, header = True,
                                                    fill_with_water = False,
                                                    water_name = None):
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
            fill_with_water -- If true, will fill all wells to the reaction size
                                with water. Requires water_name argument.
            water_name -- Determines the name of wells containing water. Must
                            match the name given in an association spreadsheet,
                            or otherwise defined.
        '''
        well_idx = process_column_argument(well_column)
        if fill_with_water:
            if not water_name:
                raise Exception("If 'Fill with water' option selected, must " +\
                                "set the name of wells contianing water")

        with open(input_filename, 'rU') as input_file:
            reader = csv.reader(input_file)
            # Skip the first row if it's a header
            if header:
                next(reader)
            for row in reader:
                if well_idx >= len(row):
                    raise ValueError("Well column out of bounds for row '%s'" %\
                                     str(row))
                well = row[well_idx]
                if not self.reactions[well]:
                    self.reactions[well] = WellReaction(self.rxn_vol, well)
                i = 0
                is_name_col = True
                volume_left = self.rxn_vol
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
                        self.reactions[well].add_material(source_material,
                                                          final_conc)
                        is_name_col = True
                    i += 1
                if fill_with_water:
                    water = self.material_list[water_name]
                    self.reactions[well].fill_with(water)

    def build_dilution_series(self, material_1, material_2,
                              material_1_final, material_2_final,
                              first_well, fill_with_water = True):
        '''Calculate TXTL volumes and automatically generate a picklist for
        an NxN dilution matrix with two inputs. If this EchoSource object has a
        MasterMix, then add that master mix as well. Fill the rest with water.

        Args:
            material_1, material_2 -- First and second materials to serialy
                                        dilute (EchoSourceMaterials)
            material_1_final, material_1_final -- Lists (or numpy array)
                                                    defining the final
                                                    concentrations of each
                                                    material.
            first_well -- The upper-left-most corner of the block on the
                          destination plate you want to spit into.
            fill_with_water -- Iff true, fills out any space left in the
                                reaction after adding master mix and materials
                                with water. Set to "False" if you want to do
                                other things to these reactions before using
                                them.
        '''
        material_1.plate = self.plates[0]
        material_2.plate = self.plates[0]

        self.add_material(material_1)
        self.add_material(material_2)

        # matrix size -- check to make sure it isn't going off the plate.
        n_material_1 = len(material_1_final)
        n_material_2 = len(material_2_final)
        first_row = string.ascii_uppercase.find(first_well[0])
        first_col = int(first_well[1:])
        if first_row + n_material_1 + 1 > self.rows \
            or first_col + n_material_2 > self.cols:
            raise ValueError(("Dilution series of size %dx%d starting in " +   \
                              "well %s runs off the edge of plate of type %s") \
                             % (n_material_1, n_material_2, first_well,
                                self.DPtype))


        # Add TX-TL master mix as a material, if applicable.
        if self.make_master_mix:
            if not "txtl_mm" in self.material_list:
                self.material_list["txtl_mm"] = MasterMix(self.plates[0],
                                                         rxn_vol = self.rxn_vol)
            txtl = self.material_list["txtl_mm"]
            txtl_mm_vol = txtl.current_vol_per_rxn()

        # Add water as a material (if it's not already there).
        self.add_material(EchoSourceMaterial("water", 0, 0, self.plates[0]))
        water = self.material_list["water"]

        # Fill in matrix with picks.
        for i in range(n_material_1):
            for j in range(n_material_2):
                # Initialize the reaction.
                well = string.ascii_uppercase[first_row + i] + \
                              str(first_col + j)
                if not self.reactions[well]:
                    self.reactions[well] = WellReaction(self.rxn_vol, well)
                # Diluted material picks
                self.reactions[well].add_material(material_1,
                                                  material_1_final[i])
                self.reactions[well].add_material(material_2,
                                                  material_2_final[j])
                if self.make_master_mix:
                    self.reactions[well].add_volume_of_material(txtl,
                                                                txtl_mm_vol)
                # Water pick
                if fill_with_water:
                    self.reactions[well].fill_with(water)

        # Add positive control....

        # and negative control.
        neg_ctrl_well = string.ascii_uppercase[first_row + n_material_1] \
                        + str(first_col)
        if self.make_master_mix:
            self.reactions[neg_ctrl_well].add_volume_of_material(txtl,
                                                                 txtl_mm_vol)
        if fill_with_water:
            self.reactions[neg_ctrl_well].fill_with(water)

    def add_material_to_well(self, material, final_conc, well):
        '''
        Add a single material, at a single concentration, to a single well.

        Parameters:
            material - An EchoSourceMaterial object representing the material
                        to add.
            final_conc - The final concentration of material, in nM (or the
                            same units as the material)
            well - Name of the well to add to.
        '''
        self.add_material_to_block(material, final_conc, well, well)

    def add_material_to_block(self, material, final_conc,
                              top_left, bottom_right):
        '''
        Add a single material, at a single concentration, to every well in a
        block on the destination plate.

        Parameters:
            material - An EchoSourceMaterial object representing the material
                        to add.
            final_conc - The final concentration of material, in nM (or the
                            same units as the material)
            top_left - top-left-most well of the block to add to.
            bottom_right - bottom-right-most well of the block to add to.
        '''
        self.add_material(material)
        material.plate = self.plates[0]

        start_row = string.ascii_uppercase.find(top_left[0])
        end_row   = string.ascii_uppercase.find(bottom_right[0])
        start_col = int(top_left[1:])-1
        end_col   = int(bottom_right[1:])-1

        for row in range(start_row, end_row+1):
            for col in range(start_col, end_col+1):
                destination = string.ascii_uppercase[row] + str(col+1)
                vol = final_conc * (self.rxn_vol / material.nM)
                self.reactions[well].add_material(material, final_conc)

    def generate_picklist(self):
        for mat in self.material_list.values():
            if mat:
                picks = mat.request_picklist()
                for pick in picks:
                    yield pick

    def write_picklist(self, outputname):
        '''
        Write this EchoCSVMaker's protocol to an Echo picklist, and print any
        other necessary directions for the user.
        '''
        # Finalize all of the reactions.
        for reaction in self.reactions.values():
            reaction.finalize_reaction()

        # Write picklist.
        # NOTE! This MUST come before writing the comment file; comments require
        # accurate count of total_volume_requested of each material, which is
        # only calculated once all of the picks are finalized and wells are
        # committed (which happens in this block).
        with open((outputname + '_EchoInput.csv'), 'w') as outcsv:
            writer = csv.writer(outcsv, lineterminator = "\n")

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
                              (pick.source_material.nM * pick.volume \
                              /self.rxn_vol)
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
                is_master_mix = (material.name == "txtl_mm" or \
                                 material.name == "master_mix")
                text_file.write("\n%s:" % material.name)
                if not is_master_mix:
                    text_file.write("\n\tstock concentration: %.2fnM" % \
                                    material.nM)
                if material.length > 0:
                    text_file.write(" (%.2f ng/uL)" % material.concentration)
                text_file.write("\n\ttotal volume: %.2f uL" % \
                                (material.total_volume_requested / 1000.0))
                # Rewrite with new MasterMixMaterial definitions (final concs
                # now in terms of final reaction)
                if is_master_mix and self.make_master_mix:
                    master_mix = self.material_list["txtl_mm"]
                    text_file.write("\n\tTubes of extract needed: %d" % \
                                    math.ceil(master_mix.n_extract_aliquots()))
                    text_file.write("\n\tTubes of buffer needed: %d" % \
                                    math.ceil(master_mix.n_buffer_aliquots()))
                    text_file.write("\n\tMaster Mix (including %d%% excess):"\
                                    %((master_mix.mm_excess-1) * 100))
                    for name, vol in master_mix.recipe():
                        text_file.write("\n\t\t%.2f uL %s" % \
                                        (vol / 1000, name))
            # Explicit loading instructions
            text_file.write("\n\nInstructions:")

            for material in self.material_list.values():
                vol_list = material.well_volumes
                if material.current_well < 0:
                    text_file.write("\n\t%s not used!" % material.name)
                    continue
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