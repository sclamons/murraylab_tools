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