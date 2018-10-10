# TODO:
#   -- Move calibration data out of source code
#       Have a string identifying the calibration data file, use that.
#       Default to most recent calibration data.
#   -- Add temperature records?

import sys
import collections
import pandas as pd
import numpy as np
import warnings
import scipy.interpolate
import csv
import os
import matplotlib.pyplot as plt
from collections import namedtuple
from ..utils import *
import pkg_resources
from datetime import datetime

####################
# Plate Reader IDs #
####################
plate_reader_ids = {"268449":'b1',
                    "271275":'b2',
                    "1402031D":'b3',
                    "18060417":'b4'}

def calibration_data_df(filename = None):
    '''
    Loads all available calibration data into a pandas DataFrame. By default,
    loads from a data file distributed with the package; optionally, you can
    instead read from your own file by setting filename. If you do, you'll need
    the following columns:
        Fluorophore
        Date
        Biotek
        Gain
        AFU per uM

    Params:
        filename: The name of the (CSV) file to read calibration data from.
                    Default None, in which case data is loaded from a
                    package-distributed data file.
    Returns: A dataframe containing all available calibration data.
    '''
    if filename is None:
        filename = pkg_resources.resource_filename('murraylab_tools',
                                   os.path.join('data', 'calibration_data.csv'))
    return pd.read_csv(filename)

def calibration_dates(filename = None):
    '''
    Returns information about when each channel was calibrated.

    Params:
        filename: The name of the (CSV) file to read calibration data from.
                    Default None, in which case data is loaded from a
                    package-distributed data file.
    Returns: dictionary of the form {Fluorophore -> [dates]}.
    '''
    df = calibration_data_df(filename)
    date_dict = dict()
    for fluor in df.Fluorophore:
        date_dict[fluor] = df[df.Fluorophore == fluor].Date.unique()
    return date_dict

def calibration_data(date = None, filename = None):
    '''
    Returns calibration data for the bioteks in the form of a nested dictionary:
        fluorophore -> {biotek -> {gain -> AFU/uM}}

    Params:
        date: Date of calibration data you would like to use, as a string of the
                form "MM/DD/YY" (e.g. 06/23/18). Defaults to None, in which case
                the latest date of calibration for each fluorophore is used.
        filename: The name of the (CSV) file to read calibration data from.
                    Default None, in which case data is loaded from a
                    package-distributed data file.
    Returns: A nested dictionary containing calibration data from a single date,
                or the most recent calibration data.
    '''
    df = calibration_data_df(filename)
    if date is not None:
        data_df = df[df.Date == date]
        if len(df) == 0:
            raise Warning("No calibration data found on date %s." % date)
    else:
        data_list = []
        for fluor in df.Fluorophore.unique():
            fluor_df = df[df.Fluorophore == fluor]
            dates = list(map(lambda d: datetime.strptime(d, "%m/%d/%y").date(),
                             fluor_df.Date))
            most_recent_date = max(dates).strftime("%m/%d/%y")
            data_list.append(fluor_df[fluor_df.Date == most_recent_date])
        data_df = pd.concat(data_list)

    # Convert data to dictionary form.
    calibration_dict = dict()
    for index, row in data_df.iterrows():
        fluor = row.Fluorophore
        date  = row.Date
        bt    = 'b' + str(row.Biotek)
        gain  = row.Gain
        AFU   = row["AFU per uM"]
        if fluor not in calibration_dict:
            calibration_dict[fluor] = dict()
            calibration_dict[fluor]["date"] = date
        if bt not in calibration_dict[fluor]:
            calibration_dict[fluor][bt] = dict()
        calibration_dict[fluor][bt][gain] = AFU
    return calibration_dict


def standard_channel_name(fp_name, calibration_dict,
                          suppress_name_warning = False):
    upper_name = fp_name.upper()
    all_names  = calibration_dict.keys()
    for name in all_names:
        if name.upper() == upper_name:
            return name
    if not suppress_name_warning:
        warnings.warn(("Unable to convert channel %s into standard channel " + \
                       "name. Are you sure this is the right name?") % fp_name)
    return fp_name

def raw_to_uM(calibration_dict, raw, protein, biotek, gain, volume):
    '''
    Convert an AFU measurement (in TX-TL) to a uM measurement, if possible.

    Params:
        calibration_dict: A nested dictionary of calibration data, as returned
                            by calibration_data().
        raw: An AFU fluorescence reading.
        protein: Name of the fluorescent protein or channel. Must match a
                    channel name in calibration_data, but isn't case-sensitive.
        biotek: Number of the biotek used, e.g. 3.
        gain: Gain used. Usually should be 61 or 100.
        volume: Volume of the TX-TL reaction.
    '''
    protein = standard_channel_name(protein, calibration_dict)
    if not protein in calibration_dict or \
       not biotek in calibration_dict[protein] or \
       not gain in calibration_dict[protein][biotek]:
       return None
    # Note that volume is in uL!
    if raw == "OVRFLW":
        raw = np.infty
    return float(raw) * 10.0 / calibration_dict[protein][biotek][gain] / volume

ReadSet = collections.namedtuple('ReadSet', ['name', 'excitation', 'emission',
                                             'gain'])

def read_supplementary_info(input_filename):
    info = dict()
    with mt_open(input_filename, 'rU') as infile:
        reader = csv.reader(infile)
        title_line = next(reader)
        title_line = list(map(lambda s:s.strip(), title_line))
        for i in range(1, len(title_line)):
            info[title_line[i]] = dict()
        for line in reader:
            line = list(map(lambda s:s.strip(), line))
            if line[0].strip() == "":
                continue
            for i in range(1, len(title_line)):
                info[title_line[i]][line[0]] = line[i]
    return info


def tidy_biotek_data(input_filename, supplementary_filename = None,
                     volume = None, convert_to_uM = True,
                     calibration_dict = None, override_plate_reader_id=None):
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
                    if you don't care about the volume set convert_to_uM to
                    False.
        --convert_to_uM: Flag that decides whether or not to calculate
                            micromolar concentrations from biotek data.
                            Default True. Note that OD data is not, by default,
                            normalized, but other channels will be even if you
                            are using cells, unless you set this flag to False.
        --calibration_dict: Dictionary of calibrations you want to use to
                                convert AFU readings to uM measurements, as
                                returned by calibration_data. Default none,
                                in which case the most recent calibration will
                                be used for each channel.
        --override_plate_reader_id: If not None, the plate reader ID will be
                                        set to this. Default None.
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
    filename_base   = input_filename.rsplit('.', 1)[0]
    output_filename = filename_base + "_tidy.csv"

    if calibration_dict is None:
        calibration_dict = calibration_data()

    # If the user gave you an excel file, convert it to a CSV so we can read
    # it properly.
    file_extension = input_filename.rpartition(".")[2]
    if file_extension.startswith("xls"):
        excel_filename = input_filename
        input_filename = excel_filename.rpartition(".")[0] + ".csv"
        pd.read_excel(excel_filename).to_csv(input_filename, index = False)

    # Open data file and tidy output file at once, so that we can stream data
    # directly from one to the other without having to store much.
    with mt_open(input_filename, 'rU') as infile:
        with mt_open(output_filename, 'w') as outfile:
            # Write a header to the tidy output file.
            reader = csv.reader(infile)
            writer = csv.writer(outfile, delimiter = ',')
            title_row = ['Channel', 'Gain', 'Time (sec)', 'Time (hr)', 'Well',
                         'Measurement', 'Units', 'Excitation', 'Emission']
            for name in supplementary_data.keys():
                title_row.append(name)
            title_row.append('ChanStr')
            writer.writerow(title_row)

            # Read plate information
            # Basic reading flow looks like:
            #   1) Read lines until a line that reads "Read", recording
            #       information about plate reader ID.
            #   2) For each line until the next line that reads "Layout":
            #       2.1) Look for a line starting with "Filter Set:"
            #       2.2) Get read set information from next two lines, store it.
            #   3) Read lines until the line that reads "Layout", looking for
            #       information about read settings.
            #   4) For each line:
            #       4.1) If line contains information about this read setting,
            #               store it.
            #       4.2) If line is the start of data, then for each line until
            #             empty line:
            #           4.2.1) Rewrite data on that line to tidy data file,
            #                   converting to uM if possible.
            read_sets = dict()
            for line in reader:
                if len(line) == 0:
                    continue
                if line[0].strip() == "Reader Serial Number:":
                    if override_plate_reader_id != None:
                            warnings.warn(("Plate reader id overridden to be '%s'") \
                                    % overrride_plate_reader_id)
                            plate_reader_id= overrride_plate_reader_id
                    elif line[1] in plate_reader_ids:
                        plate_reader_id = plate_reader_ids[line[1]]
                    else:
                        warnings.warn(("Unknown plate reader id '%s'; will " + \
                                      "not attempt to calculate molarity "   + \
                                      "concentrations.") % line[1])
                        plate_reader_id = None
                    continue
                if line[0].strip() == "Read":
                    if line[1].strip() == "Fluorescence Endpoint":
                        read_name = ""
                    else:
                        read_name = line[1].strip()
                    entered_layout = False
                    hit_data       = False
                    for line in reader:
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
                            # Sometimes excitation and emission get split to
                            # different cells; check for this.
                            if len(lineparts) == 1:
                                emission_cell = line[2]
                            else:
                                emission_cell = lineparts[1]
                            emission = int(emission_cell.split(":")[-1].split("/")[0].strip())
                            line = next(reader)
                            gain = line[1].split(",")[-1].split(":")[-1].strip()
                            if gain != "AutoScale":
                                gain = int(gain)
                            if not read_name in read_sets:
                                read_sets[read_name] = []
                            read_sets[read_name].append(ReadSet(read_name,
                                                                excitation,
                                                                emission, gain))
                        maybe_read_name = line[0].split(":")[0].strip()
                        if maybe_read_name in read_sets.keys() or \
                           maybe_read_name.startswith("OD"):
                            hit_data = True
                            break
                    if entered_layout or hit_data:
                        break

            # Read data blocks
            # Find a data block
            while line != None:
                info = line[0].strip()
                if info == "":
                    line = next(reader, None)
                    continue
                if info in ["Layout", "Results"]:
                    line = next(reader, None)
                    continue
                if info.upper().startswith("OD"):
                    reading_OD = True
                else:
                    reading_OD = False
                if reading_OD:
                    read_name = info.split(":")[0].strip()
                    excitation = int(info.split(":")[1])
                    emission   = -1
                    gain       = -1
                else:
                    if ":" in info:
                        info_parts = info.split(":")
                        read_name  = info_parts[0]
                    else:
                        info_parts = [info]
                        read_name = ""
                    if not info.endswith(']'):
                        read_idx = 0
                    else:
                        read_idx = int(info.split('[')[-1][:-1]) - 1
                    read_properties = read_sets[read_name][read_idx]
                    gain            = read_properties.gain
                    if len(info_parts) > 1:
                        excitation = info_parts[1].split("[")[0].split(",")[0]
                        excitation = int(excitation)
                        emission   = info_parts[1].split("[")[0].split(",")[1]
                        emission   = int(emission)
                    else:
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
                            continue
                        well_name = well_names[i]
                        # Check to see if there's any supplementary information
                        # on this well.
                        if supplementary_filename and \
                          not well_name in list(supplementary_data.values())[0]:
                            warnings.warn("No supplementary data for well " + \
                                        "%s; throwing out data for that well."\
                                          % well_name)
                            continue
                        afu = line[i]
                        # Check for overflow
                        if afu.upper() == "OVRFLW":
                            afu = np.infty
                        if reading_OD:
                            measurement = afu
                            units = "absorbance"
                        else:
                            if(convert_to_uM):
                                uM = raw_to_uM(calibration_dict, line[i],
                                               read_name, plate_reader_id, gain,
                                               volume)
                            else:
                                uM = None
                            if uM != None:
                                measurement = uM
                                units = "uM"
                            else:
                                measurement = afu
                                units = "AFU"
                        row = [read_name, gain, time_secs, time_hrs, well_name,
                               measurement, units, str(excitation),
                               str(emission)]

                        for name in supplementary_data.keys():
                            row.append(supplementary_data[name][well_name])
                        row.append(read_name + str(gain) + str(excitation) +\
                                   str(emission))
                        try:
                            writer.writerow(row)
                        except TypeError as e:
                            print("Error writing line: " + str(row))
                            raise e
                line = next(reader, None)


def extract_trajectories_only(df):
    '''
    Given a DataFrame that has been read in from a tidied piece of BioTek data,
    return a DataFrame that collapses each channel's AFU value as a column in a
    new DataFrame whose only columns are WellID, Time, and each measured
    channel.

    Assumptions:
        - The time is in a channel called 'Time (hr)', which becomes Time
        - There is at least 1 measured channel

    Arguments:
        df -- DataFrame of Biotek data, pulled from a tidy dataset of the form
                produced by tidy_biotek_data.
    Returns: A new DataFrame with the only columns being Well, Time, and each
                measured channel, in whatever units they were in the original
                DataFrame.
    '''
    # Create dictionary to make new data-frame
    master_dict = {}
    master_dict['Time'] = []
    master_dict['Well'] = []

    # get all channel names and initialize columns
    all_channels = df.Channel.unique()
    for channel in all_channels:
        master_dict[channel] = []

    all_wells = df.Well.unique()

    # go through and build the trajectory for each well.
    for well in all_wells:
        well_df = df[df.Well == well]
        time_series = well_df[well_df.Channel == all_channels[0]]['Time (hr)']
        series_length = len(time_series)
        well_series = [well] * series_length

        master_dict['Time'].extend(time_series)
        master_dict['Well'].extend(well_series)
        for channel in all_channels:
            channel_series = well_df[well_df.Channel == channel].Measurement
            assert (len(channel_series) == series_length)
            master_dict[channel].extend(channel_series)

    return_df = pd.DataFrame(master_dict)
    return return_df


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
                well_df.Measurement = well_df.Measurement - \
                                      avg_neg_ctrl.Measurement
                return_df = return_df.append(well_df)
    return return_df


def window_averages(df, start, end, units = "seconds",
                    grouping_variables = None):
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
        grouping_variables - Optional list of column names on which to group.
                                Use this option primarily to separate multiple
                                plates' worth of data with overlapping wells.
    '''
    group_cols = ["Channel", "Gain","Excitation","Emission", "Well"]
    if grouping_variables:
        group_cols += grouping_variables

    # Start by screening out everything outside the desired time window.
    def pick_out_window(df):
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
                raise ValueError(('Unknown unit "{0}"; units must be ' \
                                + '"seconds", "hours", or ' \
                                + '"index"').format(units))
            window_df = df[(df[col] >= start) & (df[col] <= end)]
        return window_df

    grouped_df = df.groupby(group_cols)
    window_df = grouped_df.apply(pick_out_window)

    # Figure out which columns are numeric and which have strings.
    column_names = window_df.columns.values.tolist()
    functions    = dict()
    for col in column_names:
        if col in group_cols:
            continue
        # Check to see if the column is numerically typed
        if col == "Measurement":
            #only average the measurement!! Why would you averaged
            #anything else
            #np.issubdtype(window_df[col].dtype, np.number):
            # Numbers get averaged
            functions[col] = np.average
        else:
            # Non-numbers get a copy of the last value
            functions[col] = lambda x:x.iloc[-1]

    # Calculate windowed average
    averages_df = window_df.groupby(group_cols).agg(functions)

    averages_df.reset_index(inplace = True)

    return averages_df


def endpoint_averages(df, window_size = 10, grouping_variables = None):
    '''
    Converts a dataframe of fluorescence data to a dataframe of endpoint
    average fluorescence.

    Params:
        window_size - Averages are taken over the last window_size points.
        grouping_variables - Optional list of column names on which to group.
                                Use this option primarily to separate multiple
                                plates' worth of data with overlapping wells.
    '''
    group_cols = ["Channel", "Gain", "Well"]
    if grouping_variables:
        group_cols += grouping_variables
    grouped_df = df.groupby(group_cols)
    last_time_dfs = []
    for name, group in grouped_df:
        all_times = group["Time (hr)"].unique()
        first_last_time = np.sort(all_times)[-window_size]
        last_time_dfs.append(group[group["Time (hr)"] >= first_last_time])
    end_time_dfs = pd.concat(last_time_dfs)
    return window_averages(end_time_dfs, 0, window_size, "index",
                           grouping_variables)


def spline_fit(df, column = "Measurement", smoothing_factor = None):
    '''
    Adds a spline fit of the uM traces of a dataframe of the type made by
    tidy_biotek_data.

    Params:
        df - DataFrame of time traces, of the kind produced by tidy_biotek_data.
        column - Column to find spline fit over. Defaults to "Measurement".
                    Should probably always be this.
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
        group["spline fit"] = spline(group["Time (sec)"])
        splined_df = splined_df.append(group)
    return splined_df

def smoothed_derivatives(df, column = "Measurement", smoothing_factor = None):
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
        deriv_name = "%s (%s/sec)" % (column, group.Units.unique()[0])
        group[deriv_name] = np.gradient(group["spline fit"])
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
                            Default -1 (for OD600).
    Returns:
        A DataFrame of df augmented with columns for normalized AFU/uM
        ("Normalized Measurement"). Will have units of
        "<measurement units>/<normalization units>", or "<measurement units>/OD"
        if normalizing with an OD.
    '''
    # Do some kind of check to make sure the norm channel exists with the given
    # channel...
    if not norm_channel in df.Channel.unique():
        raise ValueError("No data for channel '%s' in dataframe." % \
                         norm_channel)
    if not norm_channel_gain in df[df.Channel == norm_channel].Gain.unique():
        raise ValueError("Channel %s does not use gain %d." % \
                         (norm_channel, norm_channel_gain))

    # Iterate over channels/gains, applying normalization
    grouped_df = df.groupby(["ChanStr","Well"])
    norm_channel_df = df[(df.Channel == norm_channel) & \
                         (df.Gain == norm_channel_gain)]
    # normalized_df = pd.DataFrame()
    normalized_df_list = [None] * len(grouped_df)
    i = 0
    # Set this Pandas flag to force Pandas to NOT run garbage collection
    # all the time. Improves speed by ~3x.
    pd.set_option('mode.chained_assignment', None)
    for names, group in grouped_df:
        group.reset_index(inplace = True)
        chanstr,well = names
        norm_data = norm_channel_df[norm_channel_df.Well == well]
        norm_data.reset_index(inplace = True)
        #print("DFS!")
        #print(group.head())
        #print(norm_data.head())
        """
        print("before")
        print(group[group.Excitation != 600].Measurement.iloc[0])
        print("norm")
        print(norm_data.Measurement.iloc[0])
        print("after")
        print((group[group.Excitation != 600].Measurement / norm_data.Measurement).iloc[0])
        """
        group["Measurement"] = group.Measurement \
                               / norm_data.Measurement
        orig_units = group.Units.unique()[0]
        norm_units = "OD" if norm_channel.startswith("OD") \
                          else norm_data.Units.unique()[0]
        group.Units = "%s/%s" % (orig_units, norm_units)

        normalized_df_list[i] = group
        i += 1
    normalized_df = pd.concat(normalized_df_list)
    # Undo Pandas flag change
    pd.set_option('mode.chained_assignment', 'warn')
    return normalized_df.reset_index()


CellSpec = namedtuple("CellSpec", ['well_name', 'color', 'label'])

class BiotekCellPlotter(object):
    '''
    Class in charge of plotting fluorescence data from multiple wells, with a
    background plot of each well's OD curve in the background.

    Each BiotekCellPlotter plots data from a single dataframe (of the form
    produced by tidy_biotek_data), in a single channel, with a single gain.
    Channel data will, by default, be normalized by OD.
    '''

    def __init__(self, df, channel, gain, normalize_by_od = True,
                 od_channel = None):
        self.df        = df
        self.channel   = channel
        self.gain      = gain
        self.well_list = [] # Stores CellSpec objects, which hold a few pieces
                            # of information about the well and how it should be
                            # plotted.
        self.normalize_by_od = normalize_by_od
        if not od_channel:
            for c in self.df.Channel:
                if c.upper().startswith("OD"):
                    self.od_channel = c
                    break
        else:
            self.od_channel = od_channel
        if not self.od_channel:
            raise Exception("No OD channel specified or detected")

    def add_well(self, well, color = None, label = None):
        self.well_list.append(CellSpec(well, color, label))

    def add_condition(self, idxs, color = None, label = None):
        cond_df = self.df[idxs]
        if len(cond_df) == 0:
            raise Warning("No data for specified condition.")
        for well in cond_df.Well.unique():
            self.add_well(well, color, label)
            label = ""   # If multiple wells pulled out, only the first one gets
                         # the label. Should read more cleanly this way.

    def plot(self, title = None, split_plots = False, filename = None,
             show = True, figsize = (8, 4)):
        '''
        Plot/show/save the figure.

        Arguments:
            title -- String that will go at the top of the figure. Default "".
            split_plots -- Boolean controlling how data will be presented. If
                            True, will produce two plots -- one with OD data,
                            one with normalized fluorescence data. If False
                            (default), will show both sets of data on the same
                            plot.
            filename -- If set, will save this figure in addition to showing it.
            show -- Boolean that sets whether the plot will actually be shown.
                    Default True; set to False if you want to save without
                    viewing it.
            figsize -- 2-tuple of the size of the figure. Is passed directly
                        to figure creation.
        '''
        plt.clf()

        # Slice out and normalize all relevant data
        well_names = [ws.well_name for ws in self.well_list]
        df = self.df[self.df.Well.isin(well_names)]
        if self.normalize_by_od:
            norm_df = normalize(df, norm_channel = self.od_channel)
        else:
            norm_df = df
        norm_df = norm_df[(norm_df.Channel == self.channel) & \
                          (norm_df.Gain == self.gain)]

        # Plot out all of the fluorescence data, keeping track of the largest
        # plotted measurement.
        fig, ax1 = plt.subplots(figsize = figsize)
        for well_spec in self.well_list:
            well_df = norm_df[norm_df.Well == well_spec.well_name]
            ax1.plot(well_df["Time (hr)"], well_df.Measurement,
                     color = well_spec.color, label = well_spec.label)
        ax1.set_xlabel("Time (hr)")
        if self.normalize_by_od:
            ax1.set_ylabel("%s/OD (gain %d)" % (self.channel, self.gain))
        else:
            ax1.set_ylabel("%s (gain %d)" % (self.channel, self.gain))

        if split_plots:
            handles, labels = ax1.get_legend_handles_labels()
            ax1.legend(handles, labels)
            if title:
                plt.title(title, y = 1.08)
            fig.tight_layout()
            if filename:
                plt.savefig("fluor_" + filename, dpi = 400)
            if show:
                plt.show()
            plt.clf()

        # Plot out all OD data, scaling it to roughly the same size as the other
        # data.
        if split_plots:
            fig, ax2 = plt.subplots(figsize = figsize)
        else:
            ax2 = ax1.twinx()
        all_od_df = self.df[(self.df.Channel == self.od_channel) & \
                            (self.df.Well.isin(well_names))]
        for well_spec in self.well_list:
            od_df     = all_od_df[all_od_df.Well == well_spec.well_name]
            linestyle = "-" if split_plots else ":"
            ax2.plot(od_df["Time (hr)"], od_df.Measurement,
                     color = well_spec.color, linewidth = 1,
                     linestyle = linestyle,
                     label = well_spec.label if split_plots else "")
        ax2.set_ylabel(self.od_channel)

        handles, labels = ax1.get_legend_handles_labels()
        ax2.legend(handles, labels)
        if title:
            plt.title(title, y = 1.08)
        fig.tight_layout()
        if filename:
            if split_plots:
                filename = "OD_" + filename
            plt.savefig(filename, dpi = 400)
        if show:
            plt.show()

        return plt.gca()
