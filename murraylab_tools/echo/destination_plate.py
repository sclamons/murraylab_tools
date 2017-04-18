class DestinationPlate():
    '''
    One Echo destination plate. Responsible for allocating wells for
    reactions.
    '''
    def __init__(self, DPname = None, DPtype = None, filename = None):
        '''
        DPname -- A(n arbitrary) string representing this plate. Default is
                    "Destination[1]"
        DPtype -- A string denoting the type of plate. Plate types are
                    restricted and determine the size of the plate. Default is
                    "Nunc 384"
        filename -- Name of a file holding a list of used wells for this plate.
                    If that file doesn't exist yet, it will be created and
                    populated.
        '''
        if DPname == None:
            self.name = "Destination[1]"
        else:
            name_is_number = True
            for c in DPname:
                if c not in string.digits:
                    name_is_number = False
                    break
            if name_is_number:
                self.name = "Destination[%s]" % DPname
            else:
                self.name = DPname
        if DPtype == None or '384' in Dtype:
            self.type = "Nunc 384"
            self.rows = 16
            self.cols = 24
        elif '96' in DPtype:
            self.rows = 8
            self.cols = 12
            self.type = DPtype
        else:
            raise ValueError("'%s' is not a recognized plate type." % DPtype)

        # Construct use history object and helpers
        self.n_wells = self.rows * self.cols
        lets = string.ascii_uppercase[:self.rows]
        nums = range(1,1+self.cols)
        self.wells  = np.array(["{}{:02}".format(let, num) for let in lets for num in nums])
        self.indices = np.arange(self.n_wells)
        self.wells_used   = np.zeros((self.rows * self.cols,), dtype=bool)
        self.current_idx  = 0

        self.used_well_file = filename
        if self.used_well_file:
            if os.path.isfile(filename):
                self.load_from_file(filename)

    def make_picklist(self, source, rxns):
        n_rxns = len(rxns)
        d_wells = self.request_wells(n_rxns)
        picklist = []
        for rxn, well in zip(rxns, d_wells):
            for component in rxn:
                c_name = component[0]
                c_vol = component[1]
                if len(component)>2:
                    c_conc = component[2]
                else:
                    c_conc = None
                c_well = source.get_location(c_name, c_conc)
                picklist.append([c_well, c_name, well, c_vol, source.name,
                                 source.type, self.name])
        header = ["Source Well", "Sample Name", "Destination Well",
                  "Transfer Volume", "Source Plate Name", "Source Plate Type",
                  "Destination Plate Name"]
        return picklist, header

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
        if sum(self.wells_used) + n_wells > len(self.wells):
            raise Exception("Source plate %s is out of available wells." % self.name)
        unused_indices = self.indices[self.wells_used == False]
        return_indices = unused_indices[:n_wells]
        return_wells = self.wells[return_indices]
        self.wells_used[return_indices] = True
        return return_wells