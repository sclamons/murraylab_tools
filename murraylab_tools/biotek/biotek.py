import csv
import sys
import collections

####################
# Plate Reader IDs #
####################
plate_reader_ids = {"268449":'b1',
                    "271275":'b2',
                    "1402031D":'b3'}

####################
# Calibration Data #
####################
calibration_data = dict()
calibration_data["GFP"] = {'b3': {61:2261, 100:80850},
                           'b2': {61:1762, 100:62732},
                           'b1': {61:1517, 100:53258}}
calibration_data["Citrine"] = {'b3': {61:2379, 100:82680},
                               'b2': {61:1849, 100:64985},
                               'b1': {61:1487, 100:51725}}
calibration_data["RFP"] = {'b3': {61:70.38, 100:2541},
                           'b2': {61:67.64, 100:2392},
                           'b1': {61:47.96, 100:1689}}
calibration_data["CFP"] = {'b3': {61:578, 100:20722},
                           'b2': {61:513, 100:18136},
                           'b1': {61:387, 100:13809}}
calibration_data["Venus"] = {'b3': {61:2246, 100:79955},
                             'b2': {61:1742, 100:61130},
                             'b1': {61:1455, 100:50405}}
calibration_data["Cherry"] = {'b3': {61:79.39, 100:2850},
                              'b2': {61:80.84, 100:2822},
                              'b1': {61:51.07, 100:1782}}

def standard_channel_name(fp_name):
    upper_name = fp_name.upper()
    for std_name in ["GFP", "Citrine", "RFP", "CFP", "Venus", "Cherry"]:
        if std_name in upper_name:
            return std_name
    raise ValueError("Can't convert channel name '%s' to standard name." \
                     % fp_name)


def raw_to_uM(raw, protein, biotek, gain):
    if not protein in calibration_data or \
       not biotek in calibration_data[protein] or \
       not gain in calibration_data[protein][biotek]:
       return None
    return float(raw) / calibration_data[protein][biotek][gain]

ReadSet = collections.namedtuple('ReadSet', ['name', 'excitation', 'emission',
                                             'gain'])

def read_supplementary_info(input_filename):
    info = dict()
    with open(input_filename, 'rU') as infile:
        reader = csv.reader(infile)
        title_line = reader.next()
        for i in range(1, len(title_line)):
            info[title_line[i]] = dict()
        for line in reader:
            if line[0].strip() == "":
                continue
            for i in range(1, len(title_line)):
                info[title_line[i]][line[0]] = line[i]
    return info


def tidy_biotek_data(input_filename, supplementary_filename = None):
    supplementary_data = dict()
    if supplementary_filename:
        supplementary_data = read_supplementary_info(supplementary_filename)
    filename_base   = input_filename.rsplit('.')[0]
    output_filename = filename_base + "_tidy.csv"

    with open(input_filename, 'rU') as infile:
        with open(output_filename, 'w') as outfile:
            reader = csv.reader(infile)
            writer = csv.writer(outfile, delimiter = ',')
            title_row = ['Channel', 'Gain', 'Time (sec)', 'Well', 'AFU', 'uM',
                             'Excitation', 'Emission']
            for name in supplementary_data.keys():
                title_row.append(name)
            writer.writerow(title_row)

            # Read plate information
            read_sets = dict()
            for line in reader:
                if line[0].strip() == "Reader Serial Number:":
                    if line[1] in plate_reader_ids:
                        plate_reader_id = plate_reader_ids[line[1]]
                    else:
                        raise ValueError("Unknown plate reader id '%s'" % \
                                         line[1])
                    continue
                if line[0] == "Read":
                    read_name = line[1]
                    for line in reader:
                        if line[1].startswith("Read Height"):
                            break
                        if line[1].startswith("Filter Set"):
                            line = reader.next()
                            lineparts = line[1].split(",")
                            excitation = int(lineparts[0].split(":")[-1].strip())
                            emission = int(lineparts[1].split(":")[-1].strip())
                            line = reader.next()
                            gain = int(line[1].split(",")[-1].split(":")[-1].strip())
                            if not read_name in read_sets:
                                read_sets[read_name] = []
                            read_sets[read_name].append(ReadSet(read_name,
                                                                excitation,
                                                                emission, gain))
                if line[0] == "End Kinetic":
                    break

            # Read data blocks
            # Find a data block
            for line in reader:
                if line[0].strip() == "":
                    continue
                info = line[0].strip()
                if info in ["Layout", "Results"]:
                    continue
                read_name = info.split(":")[0]
                if not info.endswith(']'):
                    read_idx = 0
                else:
                    read_idx = int(info.split('[')[-1][:-1]) - 1
                read_properties = read_sets[read_name][read_idx]
                gain            = read_properties.gain
                excitation      = read_properties.excitation
                emission        = read_properties.emission

                line = reader.next() # Skip a line
                line = reader.next() # Chart title line
                well_names = line
                # Data lines
                for line in reader:
                    if line[1] == "":
                        break
                    raw_time = line[1]
                    time_parts = raw_time.split(':')
                    time = int(time_parts[2]) + 60*int(time_parts[1]) \
                           + 3600*int(time_parts[0])
                    temp = line[2]
                    for i in range(3,len(line)):
                        if line[i].strip() == "":
                            break
                        well_name = well_names[i]
                        afu       = line[i]
                        uM_conc   = raw_to_uM(line[i],
                                              standard_channel_name(read_name),
                                              plate_reader_id, gain)
                        row = [read_name, gain, time, well_name, afu, uM_conc,
                                         excitation, emission]
                        for name in supplementary_data.keys():
                            row.append(supplementary_data[name][well_name])
                        writer.writerow(row)

def main():
    input_filename = sys.argv[1]
    if len(sys.argv) > 2:
        supplementary_filename = sys.argv[2]
    else:
        supplementary_filename = None
    tidy_biotek_data(input_filename, supplementary_filename)

if __name__ == "__main__":
    main()