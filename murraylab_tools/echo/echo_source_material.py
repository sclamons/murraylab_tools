# coding=utf-8

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

        # Check for commas in the name.
        if "," in self.name:
            warnings.warn(("Material %s has comma in its name; this may bug " +
                          "the echo when you run it.") % self.name, Warning)

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
        self.current_well = -1

    def __str__(self):
        if self.length > 0:
            conc_string = "ng/µL"
        else:
            conc_string = "nM"
        return "%s (.3f %s)" % (self.name, self.nM, conc_string)

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
        if actual_volume == 0:
            warnings.warn("Requesting 0 volume from material " + self.name + \
                          " into well " + destination_well + "; are you sure "+\
                          "you want to do this?")
        else:
            self.total_volume_requested += actual_volume
            self.picklist.append(Pick(self, None, destination_well,
                                      actual_volume))

    def request_picklist(self):
        '''
        Commit to wells, and return a list of finalized picks from this
        material.
        '''
        usable_volume  = max_volume - dead_volume
        n_source_wells = math.ceil(float(self.total_volume_requested) \
                                         / usable_volume)
        if n_source_wells == 0:
            print(("Warning: Material %s is requesting 0 wells in its " +\
                  "source plate to give %f total volume") % \
                  (self.name, self.total_volume_requested))
        if self.wells == None:
            self.wells = self.plate.request_wells(int(n_source_wells),self.name)
        if len(self.wells) < 1:
            warnings.warn(("Material %s has requested no wells. Are you sure "+\
                          "this is correct?") % self.name)
            raise StopIteration
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
                    print("self.wells: " + str(self.wells))
                    print("self.current_well: " + str(self.current_well))
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