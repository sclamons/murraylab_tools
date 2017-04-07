class Reaction(object):
    '''
    Container class for mixes of liquids.
    '''
    def __init__(self, rxn_vol, well = None):
        self.rxn_vol   = rxn_vol
        self.well      = well
        self.finalized = False
        self.materials = []

    def add_material(self, material, final_conc):
        '''
        Add a material at a known final concentration. Final concentrations are
        assumed to use the same units as the material (usually nM). Does NOT
        check final reaction volume -- that will not be checked until
        finalize_reaction is called (happens automatically when recipe is
        called).

        Rounding to Echo-compatible volumes occurs at this step.
        '''
        target_vol  = self.rxn_vol * final_conc / material.nM
        actual_vol  = echo_round(target_vol)
        actual_conc = actual_vol * material.nM / self.rxn_vol

        self.materials.append(material, actual_conc)

    def current_vol(self):
        '''
        Calculates the total volume of all of the materials currently in the
        reaction.
        '''
        return sum([vol for name, vol in self.recipe(finalize = False)])

    def fill_with(self, material):
        '''
        Fill all unfilled volume in the reaction with some material (usually
        water).
        '''
        fill_volume         = self.rxn_vol - self.current_vol()
        dilution_factor     = fill_volume / self.rxn_vol
        material_final_conc = material.nM * dilution_factor
        self.add_material(material, material_final_conc)

    def finalize_reaction(self):
        '''
        Checks reaction for consistency, throwing a ValueError if the reaction
        is overfilled or otherwise in obvious error, and raising a Warning if
        the reaction is underfull. Then, if everything checks out, volume is
        requested from the reaction's EchoSourceMaterials.
        '''
        current_vol = self.current_vol()
        if current_vol > self.rxn_vol:
            error_string = "Reaction "
            if self.well:
                error_string += "in well %s " % self.well
            error_string += "has %d nL volume but contains %.2f nL of " \
                            % (self.rxn_vol, current_vol)
            error_string += "ingredients:"
            for material, conc in self.materials:
                material_vol = conc * self.rxn_vol / material.nM
                error_string += "\n\t%d nL of %s" % (material_vol, material)
            raise ValueError(error_string)
        if current_vol < self.rxn_vol:
            warn_string = "Reaction "
            if self.well:
                warn_string += "in well %s " % self.well
            warn_string += "has %d nL volume but only contains %.2f nL of " \
                            % (self.rxn_vol, current_vol)
            warn_string += "ingredients. Are you sure you want to underfill " \
                            + "this reaction?"

            warnings.warn(warn_string, Warning)

        self.finalized = True

    def recipe(self, finalize = True):
        '''
        Iterator returning descriptors of what goes in the reaction. If the
        reaction hasn't been finalized, does so.

        Arguments:
            finalize -- Iff True (default), and if the reaction hasn't already
                            been finalized, makes sure that it gets finalized.

        Yields -- pairs of the form (material, vol), where 'material' is an
                    EchoSourceMaterial and 'vol' is the volume of that material
                    to add to the reaction, in nL.
        '''
        # Make sure everything's ready to go and materials have been requested.
        if finalize and not self.finalized:
            self.finalize_reaction()

        for material, final_conc in self.materials:
            name = str(material)
            vol  = final_conc * self.rxn_vol / material.nM
            yield (name, vol)



class MasterMix(EchoSourceMaterial):
    '''
    Container class for a list of materials to make up a master mix. Note that
    this class assumes that 75% of the final reaction volume will be extract +
    buffer.
    '''
    def __init__(self, plate, extract_fraction = 0.33, mm_excess = 1.1,
                 add_txtl = True, extract_per_aliquot = 30000,
                 buffer_per_aliquot = 37000):
        '''
        extract_fraction: If TX-TL is added, this is the fraction of the final
                            mix made up of TX-TL extract. Default 0.33 (lowest
                            protein concentration).
        mm_excess: The ratio of master-mix-to-make to total-mix-needed, i.e.,
                        mm_excess=1.1 => Make 10% excess, to account for
                        pipetting loss.
        add_txtl: If true, buffer and extract will automatically be added to
                    the master mix, using an extract percentage set by
                    extract_fraction. Default True.
        extract_per_aliquot: Volume of TX-TL extract in one aliquot, in nL.
                                Default 30000.
        buffer_per_aliquot: Volume of TX-TL buffer in one aliquot, in nL.
                                Default 37000.
        '''
        if add_txtl:
            self.name = "txtl_mm"
        else:
            self.name = "master_mix"
        self.length = 0
        self.plate = plate
        self.nM = 1   # Proxy nanomolar value

        self.wells = None
        self.picklist = []
        self.total_volume_requested = 0
        self.well_volumes = None

        self.rxn_vol = None
        self.mm_excess = mm_excess
        self.extract_fraction = extract_fraction
        self.extract_per_aliquot = extract_per_aliquot
        self.buffer_per_aliquot = buffer_per_aliquot
        self.materials = []
        if add_txtl:
            self.buffer_fraction = 0.75 - self.extract_fraction
            self.materials.append(MasterMixMaterial("Extract", 1,
                                                    self.extract_fraction))
            self.materials.append(MasterMixMaterial("Buffer", 1,
                                                    self.buffer_fraction))

    def add_material(self, name, init_conc, final_conc):
        '''
        Adds a new material to the master mix.

        name: String describing the material.
        init_conc: Stock concentration of the material.
        final_conc: Concentration of the material in the *final reaction*
                        (not in the master mix).
        '''
        self.materials.append(MasterMixMaterial(name, init_conc, final_conc))

    def set_rxn_vol(self, rxn_vol):
        '''
        Sets the (total) volume of each reaction that will use this master mix,
        in nL.
        '''
        self.rxn_vol = rxn_vol

    def one_rxn_recipe(self):
        '''
        Iterator returning descriptors of what goes in the master mix for a
        single reaction.

        self.rxn_vol must be set before calling this function(preferably with
        set_rxn_vol). Throws an AttributeError otherwise.
        '''
        if self.rxn_vol == None:
            raise AttributeError("rxn_vol must be set using set_rxn_vol " +
                                 "before requesting a recipe.")

        for material in self.materials:
            name = material.name
            vol  = material.final * self.rxn_vol / material.stock
            yield (name, vol)

    def recipe(self):
        '''
        Iterator returning descriptors of what goes in the total master mix.
        This is where the excess fraction is added in.

        Yields -- pairs of the form (name, vol), where 'name' is a string
                    describing a material in the master mix and 'vol' is the
                    volume of that material to add to the master mix, in nL.
        '''
        if self.total_volume_requested != 0:
            ingredients = self.one_rxn_recipe()
            one_rxn_vol = self.vol_per_rxn()
            for (name, vol) in ingredients:
                ingredient_fraction = vol / one_rxn_vol
                yield (name, self.mm_excess * ingredient_fraction \
                             * self.total_volume_requested)

        # if self.total_volume_requested != 0:
        #     for material in self.materials:
        #         name = material.name
        #         vol  = material.final * self.total_volume_requested \
        #                / material.stock / 0.75 * self.mm_excess
        #         yield (name, vol)


    def vol_per_rxn(self):
        '''
        Returns the calculated total volume of master mix in each reaction. NOT
        the same as self.rxn_vol (which is the *total* volume of each reaction,
        master mix plus everything else).

        self.rxn_vol must be set before calling this function
        (preferably with set_rxn_vol).Throws an AttributeError otherwise.
        '''
        if self.rxn_vol == None:
            raise AttributeError("rxn_vol must be set using set_rxn_vol " +
                                 "before calling vol_per_rxn.")

        return sum([vol for name, vol in self.one_rxn_recipe()])

    def n_extract_aliquots(self):
        '''
        Returns the number of extract aliquots required to make this master mix.

        self.total_volume_requested should be set before calling this function.
        Throws an AttributeError otherwise. total_volume_requested is set during
        a call to request_picklist, when picks are finalized. If
        total_volume_requested is not set, will return 0.
        '''
        for material in self.materials:
            if material.name == "Extract":
                extract_vol = self.total_volume_requested * material.final \
                                / material.stock / 0.75
                return extract_vol / self.extract_per_aliquot \
                        * self.mm_excess
        return 0

    def n_buffer_aliquots(self):
        '''
        Returns the number of buffer aliquots required to make this master mix.

        self.total_volume_requested must be set before calling this function.
        Throws an AttributeError otherwise. total_volume_requested is set during
        a call to request_picklist, when picks are finalized. If
        total_volume_requested is not set, will return 0.
        '''
        for material in self.materials:
            if material.name == "Buffer":
                buffer_vol = self.total_volume_requested * material.final \
                                / material.stock / 0.75
                return buffer_vol / self.buffer_per_aliquot \
                        * self.mm_excess

        return 0

    def fill_with(self, name, vol_per_rxn):
        '''
        Adds a new material (usually going to be water) such that it fills out
        the reaction to a specified volume.

        Don't try changing the reaction volume after calling this, or you'll
        get some weird answers.

        name: String describing the material.
        vol_per_rxn: The final volume you want of master mix, in nL, per
                        reaction.
        '''
        new_material_volume = vol_per_rxn - self.vol_per_rxn()
        dilution_factor     = new_material_volume / self.rxn_vol
        self.add_material(name, 1, dilution_factor)