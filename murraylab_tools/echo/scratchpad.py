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
        '''
        self.materials.append(material, final_conc)

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
        the reaction is underfull.
        '''
        current_vol = self.current_vol()
        if current_vol > self.rxn_vol:
            error_string = "Reaction "
            error_string += "has %d nL volume but contains %.2f nL of " \
                            % (self.rxn_vol, current_vol)
            error_string += "ingredients:"
            for material, conc in self.materials:
                material_vol = conc * self.rxn_vol / material.nM
                error_string += "\n\t%d nL of %s" % (material_vol, material)
            raise ValueError(error_string)

        if current_vol < self.rxn_vol:
            warn_string = "Reaction "
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


class WellReaction(Reaction):
    '''
    A reaction in a well on an Echo destination plate. Has a well, and has
    volumes that are rounded to Echo-compatible numbers.
    '''
    def __init__(self, rxn_vol, well):
        super(WellReaction, self).__init__(rxn_vol)
        self.well = well

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

        for material, conc in self.materials:
            vol = conc * self.rxn_vol / material.nM
            material.request_material(self.well, vol)

        self.finalized = True



class MasterMix(Reaction, EchoSourceMaterial):
    '''
    Container class for a list of materials that make up a master mix. This
    is any mix of materials that are combined into one single material that
    is in turn put into an Echo source well.
    '''
    def __init__(self, plate, extract_fraction = 0.33, mm_excess = 1.1,
                 rxn_vol = 10, add_txtl = True, extract_per_aliquot = 30000,
                 buffer_per_aliquot = 37000):
        '''
        extract_fraction: If TX-TL is added, this is the fraction of the final
                            mix made up of TX-TL extract. Default 0.33 (lowest
                            protein concentration).
        mm_excess: The ratio of master-mix-to-make to total-mix-needed, i.e.,
                        mm_excess=1.1 => Make 10% excess, to account for
                        pipetting loss.
        rxn_vol: Total volume of a single reaction using this master mix.
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

        self.wells    = None
        self.picklist = []
        self.total_volume_requested = 0
        self.well_volumes = None

        self.rxn_vol   = rxn_vol
        self.mm_excess = mm_excess
        self.extract_fraction = extract_fraction
        self.extract_per_aliquot = extract_per_aliquot
        self.buffer_per_aliquot = buffer_per_aliquot
        self.txtl_fraction = 0.75
        self.materials = []
        if add_txtl:
            self.buffer_fraction = self.txtl_fraction - self.extract_fraction
            txtl_extract = EchoSourceMaterial("Extract", 1, 0, None)
            txtl_buffer  = EchoSourceMaterial("Buffer",  1, 0, None)
            self.materials.append(txtl_extract, self.extract_fraction))
            self.materials.append(txtl_buffer, self.buffer_fraction))

    def one_rxn_recipe(self, finalize = True):
        '''
        Iterator returning descriptors of what goes in a SINGLE REACTION.

        Arguments:
            finalize -- Iff True (default), and if the reaction hasn't already
                            been finalized, makes sure that it gets finalized.

        Yields -- pairs of the form (material, vol), where 'material' is an
                    EchoSourceMaterial and 'vol' is the volume of that material
                    to add to the reaction, in nL.
        '''
        for material, vol in super(MasterMix, self).recipe(finalize):
            yield material, vol

    def current_vol_per_rxn(self):
        '''
        Returns the current PER REACTION volume based on the materials
        currently in the master mix.
        '''
        return sum([vol for name, vol in self.one_rxn_recipe(finalize = False)])

    def current_vol(self):
        '''
        Wraps current_vol_per_rxn
        '''
        return self.current_vol_per_rxn()

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
        if self.total_volume_requested != 0:
            ingredients = self.one_rxn_recipe()
            one_rxn_vol = self.current_vol_per_rxn()
            for (name, vol) in ingredients:
                ingredient_fraction = vol / one_rxn_vol
                yield (name, self.mm_excess * ingredient_fraction \
                             * self.total_volume_requested)

    def n_extract_aliquots(self, ):
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
                                / material.stock / self.txtl_fraction
                return extract_vol / self.extract_per_aliquot \
                        * self.mm_excess
        return 0

    def n_buffer_aliquots(self, ):
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
                                / material.stock / self.txtl_fraction
                return buffer_vol / self.buffer_per_aliquot \
                        * self.mm_excess

        return 0
