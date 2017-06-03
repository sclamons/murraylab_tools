# TODO:
#   -- Move calibration data out of source code
#   -- Add temperature records?

import csv
import sys
import collections
import pandas as pd
import numpy as np
import warnings
import scipy.interpolate

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

def standard_channel_name(fp_name, suppress_name_warning = False):
    upper_name = fp_name.upper()
    for std_name in ["GFP", "Citrine", "RFP", "CFP", "Venus", "Cherry"]:
        if std_name in upper_name:
            return std_name
    if not suppress_name_warning:
        warnings.warn(("Unable to convert channel %s into standard channel " + \
                       "name. Are you sure this is the right name?") % fp_name)
    return fp_name


def raw_to_uM(raw, protein, biotek, gain, volume):
    if not protein in calibration_data or \
       not biotek in calibration_data[protein] or \
       not gain in calibration_data[protein][biotek]:
       return None
    # Note that volume is in uL!
    if raw == "OVRFLW":
        raw = np.infty
    return float(raw) * 10.0 / calibration_data[protein][biotek][gain] / volume

ReadSet = collections.namedtuple('ReadSet', ['name', 'excitation', 'emission',
                                             'gain'])

def read_supplementary_info(input_filename):
    info = dict()
    with open(input_filename, 'rU') as infile:
        reader = csv.reader(infile)
        title_line = next(reader)
        for i in range(1, len(title_line)):
            info[title_line[i]] = dict()
        for line in reader:
            if line[0].strip() == "":
                continue
            for i in range(1, len(title_line)):
                info[title_line[i]][line[0]] = line[i]
    return info


def tidy_biotek_data(input_filename, supplementary_filename = None,
                     volume = None):
    '''
    Convert the raw output from a Biotek plate reader into tidy data.
    Optionally, also adds columns of metadata specified by a "supplementary
    file", which is a CSV spreadsheet mapping well numbers to metadata.

    Arguments:
        --input_filename: Name of a Biotek output file. Data file should be
                            standard excel output files, saved as a CSV.
        --supplementary_filename: Name of a supplementary file. Supplementary
                                    file must be a CSV wit a header, where the
                                    first column is the name of the well,
                                    additional columns define additional
                                    metadata, and each row after the header is a
                                    single well's metadata. Defaults to None
                                    (no metadata other than what can be mined
                                    from the data file).
        --volume: Volume of the TX-TL reactions. Note that default is 10 uL!
    Returns: None
    Side Effects: Creates a new CSV with the same name as the data file with
                    "_tidy" appended to the end. This new file is in tidy
                    format, with each row representing a single channel read
                    from a single well at a single time.

    '''
    if volume == None:
        print("Assuming default volume 10 uL. Make sure this is what you want!")
        volume = 10.0

    supplementary_data = dict()
    if supplementary_filename:
        supplementary_data = read_supplementary_info(supplementary_filename)
    filename_base   = input_filename.rsplit('.')[0]
    output_filename = filename_base + "_tidy.csv"

    with open(input_filename, 'rU') as infile:
        with open(output_filename, 'w') as outfile:
            reader = csv.reader(infile)
            writer = csv.writer(outfile, delimiter = ',')
            title_row = ['Channel', 'Gain', 'Time (sec)', 'Time (hr)', 'Well',
                         'AFU', 'uM', 'Excitation', 'Emission']
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
                        warnings.warn(("Unknown plate reader id '%s'; will " + \
                                      "not attempt to calculate molarity " + \
                                      "concentrations.") % line[1])
                        plate_reader_id = None
                    continue
                if line[0].strip() == "Read":
                    if line[1].strip() == "Fluorescence Endpoint":
                        read_name = ""
                    else:
                        read_name = line[1]
                    entered_layout = False
                    for line in reader:
                        print("processing line: " + line)
                        if line[0].strip() == "Layout":
                            entered_layout = True
                            break
                        if line[0].strip() == "Read":
                            if line[1].strip() == "Fluorescence Endpoint":
                                read_name = ""
                            else:
                                read_name = line[1].strip()
                            continue
                        if line[1].startswith("Filter Set"):
                            line = next(reader)
                            lineparts = line[1].split(",")
                            excitation = int(lineparts[0].split(":")[-1].split("/")[0].strip())
                            emission = int(lineparts[1].split(":")[-1].split("/")[0].strip())
                            line = next(reader)
                            gain = line[1].split(",")[-1].split(":")[-1].strip()
                            if gain != "AutoScale":
                                gain = int(gain)
                            if not read_name in read_sets:
                                read_sets[read_name] = []
                            read_sets[read_name].append(ReadSet(read_name,
                                                                excitation,
                                                                emission, gain))
                    if entered_layout:
                        break

            # Read data blocks
            # Find a data block
            for line in reader:
                if line[0].strip() == "":
                    continue
                info = line[0].strip()
                if info in ["Layout", "Results"]:
                    continue
                if info.strip() == "OD600":
                    reading_OD = True
                else:
                    reading_OD = False
                if reading_OD:
                    read_name = "OD600"
                    excitation = 600
                    emission   = -1
                    gain       = -1
                else:
                    read_name = info.split(":")[0]
                    if not info.endswith(']'):
                        read_idx = 0
                    else:
                        read_idx = int(info.split('[')[-1][:-1]) - 1
                    read_properties = read_sets[read_name][read_idx]
                    gain            = read_properties.gain
                    excitation      = read_properties.excitation
                    emission        = read_properties.emission

                line = next(reader) # Skip a line
                line = next(reader) # Chart title line
                well_names = line
                # Data lines
                for line in reader:
                    if line[1] == "":
                        break
                    raw_time = line[1]
                    time_parts = raw_time.split(':')
                    time_secs = int(time_parts[2]) + 60*int(time_parts[1]) \
                                + 3600*int(time_parts[0])
                    time_hrs = time_secs / 3600.0
                    temp = line[2]
                    for i in range(3,len(line)):
                        if line[i].strip() == "":
                            break
                        well_name = well_names[i]
                        # Check to see if there's any supplementary information
                        # on this well.
                        if supplementary_filename and \
                          not well_name in list(supplementary_data.values())[0]:
                            warnings.warn("No supplementary data for well " + \
                                          "%s; throwing out data for that well."\
                                          % well_name)
                            continue
                        afu       = line[i]
                        # Check for overflow
                        if afu.upper() == "OVRFLW":
                            afu = np.infty
                        if reading_OD:
                            uM_conc = afu
                        else:
                            uM_conc = raw_to_uM(line[i],
                                            standard_channel_name(read_name),
                                            plate_reader_id, gain, volume)
                        if uM_conc == None:
                            uM_conc = -1
                        row = [read_name, gain, time_secs, time_hrs, well_name,
                               afu, uM_conc, excitation, emission]
                        for name in supplementary_data.keys():
                            row.append(supplementary_data[name][well_name])
                        writer.writerow(row)


def background_subtract(df, negative_control_wells):
    '''
    Create a new version of a dataframe with background removed. Background is
    inferred from one or more negative control wells. If more than one negative
    control is specified, the average of the wells is used as a background
    value.

    Note that this function assumes that every measurement has a corresponding
    negative control measurement in each of the negative control wells (same
    channel and gain).

    Arguments:
        df -- DataFrame of Biotek data, pulled from a tidy dataset of the form
                produced by tidy_biotek_data.
        negative_control_wells -- String or iterable of Strings specifying one
                                    or more negative control wells.
    Returns: A new DataFrame with background subtracted out.
    '''
    if type(negative_control_wells) == str:
        negative_control_wells = [negative_control_wells]
    return_df = pd.DataFrame()
    # Split the dataframe by channel and gain
    for channel in df.Channel.unique():
        channel_df = df[df.Channel == channel]
        for gain in channel_df.Gain.unique():
            condition_df = channel_df[channel_df.Gain == gain]
            neg_ctrl_df  = pd.DataFrame()
            for well in negative_control_wells:
                well_df = condition_df[condition_df.Well == well]
                neg_ctrl_df = neg_ctrl_df.append(well_df)
            grouped_neg_ctrl = neg_ctrl_df.groupby(["Time (sec)"])
            avg_neg_ctrl = grouped_neg_ctrl.aggregate(np.average)
            avg_neg_ctrl.sort_index(inplace = True)
            avg_neg_ctrl.reset_index(inplace = True)
            # Easiest thing to do is to apply the background subtraction to each
            # well separately
            for well in condition_df.Well.unique():
                well_df = condition_df[condition_df.Well == well].copy()
                well_df.sort_values("Time (sec)", inplace = True)
                well_df.reset_index(inplace = True)
                well_df.AFU = well_df.AFU - avg_neg_ctrl.AFU
                if channel in calibration_data and \
                   gain in calibration_data[channel]['b1']:
                    well_df.uM = well_df.uM - avg_neg_ctrl.uM
                return_df = return_df.append(well_df)
    return return_df


def window_averages(df, start, end, units = "seconds"):
    '''
    Converts a dataframe of fluorescence data to a dataframe of average
    fluorescences from a specified time window. The time window can be specified
    in seconds, hours, or index.

    Args:
        df -- Dataframe of fluorescence data
        start -- First frame to be included (inclusively)
        end -- Last frame to be included (also inclusively)
        units -- Either "seconds" (default), "hours", or "index". If "seconds"
                    or "hours", takes data within a time window (using the
                    obvious columns). If "index", takes a slice using an index,
                    where the first time is 0, etc.
    '''
    # Find times within the given window.
    if units.lower() == "index":
        all_times = df["Time (sec)"].unique()
        all_times.sort()
        window_times = all_times[start:end+1]
        window_df = df[df["Time (sec)"].isin(window_times)]
    else:
        if units.lower() == "seconds":
            col = "Time (sec)"
        elif units.lower() == "hours":
            col = "Time (hr)"
        else:
            raise ValueError(('Unknown unit "{0}"; units must be "seconds", ' \
                            +'"hours", or "index"').format(units))
        window_df = df[(df[col] >= start) & (df[col] <= end)]

    # Group by channel, gain, and well
    grouped_df = window_df.groupby(["Channel", "Gain", "Well"])

    # Figure out which columns are numeric and which have strings.
    column_names = df.columns.values.tolist()
    functions    = dict()
    for col in column_names:
        if col in ["Channel", "Gain", "Well"]:
            continue
        # Check to see if the column is numerically typed
        if np.issubdtype(window_df[col].dtype, np.number):
            # Numbers get averaged
            functions[col] = np.average
        else:
            # Non-numbers get a copy of the last value
            functions[col] = lambda x:x.iloc[-1]

    # Calculate windowed average
    averages_df = grouped_df.agg(functions)

    averages_df.reset_index(inplace = True)

    return averages_df


def endpoint_averages(df, window_size = 10):
    '''
    Converts a dataframe of fluorescence data to a dataframe of endpoint
    average fluorescence. Averages taken over the last window_size points.
    '''
    all_times = df["Time (sec)"].unique()
    start = len(all_times) - window_size
    end = len(all_times)-1
    return window_averages(df, start, end, "index")


def spline_fit(df, column = "uM", smoothing_factor = None):
    '''
    Adds a spline fit of the uM traces of a dataframe of the type made by
    tidy_biotek_data.

    Params:
        df - DataFrame of time traces, of the kind produced by tidy_biotek_data.
        column - Column to find spline fit over. Defaults to "uM".
        smoothing_factor - Parameter determining the tightness of the fit.
                            Default is the number of time points. Smaller
                            smoothing factor produces tighter fit; 0 smoothing
                            factor interpolates every point. See parameter 's'
                            in scipy.interpolate.UnivariateSpline.
    Returns:
        A DataFrame of df augmented with columns for a spline fit.
    '''
    # Fit 3rd order spline
    grouped_df = df.groupby(["Channel", "Gain", "Well"])
    splined_df = pd.DataFrame()
    for name, group in grouped_df:
        spline = scipy.interpolate.UnivariateSpline(group["Time (sec)"],
                                                    group[column],
                                                    s = smoothing_factor)
        group[column + " spline fit"] = spline(group["Time (sec)"])
        splined_df = splined_df.append(group)
    return splined_df

def smoothed_derivatives(df, column = "uM", smoothing_factor = None):
    '''
    Calculates a smoothed derivative of the time traces in a dataframe. First
    fits a spline, then adds the derivatives of the spline to a copy of the
    DataFrame, which is returned.

    Args:
        df - DataFrame of time traces, of the kind produced by tidy_biotek_data.
        column - Column to fit derivatives to. Defaults to "uM"
        smoothing_factor - Parameter determining the tightness of the spline
                            fit made before derivative calculation.
                            Default is the number of time points. Smaller
                            smoothing factor produces tighter fit; 0 smoothing
                            factor interpolates every point. See parameter 's'
                            in scipy.interpolate.UnivariateSpline.
    Returns:
        A DataFrame of df augmented with columns for a spline fit and a
        derivative
    '''
    splined_df = spline_fit(df, column, smoothing_factor)
    grouped_df = splined_df.groupby(["Channel", "Gain", "Well"])
    deriv_df   = pd.DataFrame()
    for name, group in grouped_df:
        group[column + "/sec"] = np.gradient(group[column + " spline fit"])
        deriv_df = deriv_df.append(group)
    return deriv_df

def normalize(df, norm_channel = "OD600", norm_channel_gain = -1):
    '''
    Normalize expression measurements by dividing each measurement by the value
    of a reference channel at that time (default OD600).

    Args:
        df - DataFrame of time traces, of the kind produced by tidy_biotek_data.
        norm_channel - Name of a channel to normalize by. Default "OD600"
        norm_channel_gain - Gain of the channel you want to normalize by.
                            Default -1 (for OD600)
    Returns:
        A DataFrame of df augmented with columns for normalized AFU
        ("AFU/<norm_channel>").
    '''
    # Do some kind of check to make sure the norm channel exists with the given
    # channel...
    if not norm_channel in df.Channel:
        raise ValueError("No data for channel '%s' in dataframe." % \
                         norm_channel)
    if not norm_channel_gain in df[df.Channel == norm_channel].Gain:
        raise ValueError("Channel %s does not use gain %d." % \
                         (norm_channel, norm_channel_gain))

    # Iterate over channels/gains, applying normalization
    grouped_df = df.groupby(["Channel", "Gain", "Well"])
    normalized_df = pd.DataFrame()
    for name, group in grouped_df:
        group["AFU/" + norm_channel] = group["AFU"] / group[norm_channel]
        normalized_df = normalized_df.append(group)
    return normalized_df

