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
import math
import matplotlib.pyplot as plt
import seaborn as sns
from collections import namedtuple
from ..utils import *
import pkg_resources
from datetime import datetime
from copy import deepcopy as dc

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
                     volume = None, convert_to_uM = False,
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
            read_set_idxs = dict()
            next_line = ""
            while True:
                if next_line != "":
                    line = next_line
                    next_line = ""
                else:
                    try:
                        line = next(reader)
                    except StopIteration:
                        break
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

                    # Process all the information for this read set into one
                    # solid text block, which we will search for relevant
                    # information
                    read_information_block = ""
                    for line in reader:
                        if len(line) == 0:
                            continue
                        if line[0].strip() == "Layout":
                            entered_layout = True
                            break
                        maybe_read_name = line[0].split(":")[0].strip()
                        if maybe_read_name in read_sets.keys() or \
                           maybe_read_name.startswith("OD"):
                            hit_data = True
                            break
                        if line[0].strip() == "Read":
                            next_line = line
                            break
                        line[-1] = line[-1].strip()
                        read_information_block += ",".join(line)

                    # Now go through the block to figure out information for the
                    # read set.
                    info_parts = read_information_block.split(",")
                    for idx in range(len(info_parts)):
                        block = info_parts[idx]
                        block = block.strip()
                        if block.startswith("Wavelengths"):
                            emission = int(block.split(":")[-1].strip()\
                                           .split("/")[0].strip())
                            excitation = emission
                        if block.startswith("Filter Set"):
                            for i in range(idx+1, len(info_parts)):
                                sub_block = info_parts[i]
                                sub_block = sub_block.strip()
                                if sub_block.startswith("Excitation"):
                                    excitation = int(sub_block.split(":")[-1]\
                                                     .split("/")[0].strip())
                                if sub_block.startswith("Emission"):
                                    emission = int(sub_block.split(":")[-1]\
                                                     .split("/")[0].strip())
                                if sub_block.startswith("Gain"):
                                    gain = sub_block.split(",")[-1]\
                                           .split(":")[-1].strip()
                                    if gain != "AutoScale":
                                        gain = int(gain)
                                if sub_block.startswith("Filter Set"):
                                    break
                            if not read_name in read_sets:
                                read_sets[read_name] = []
                                read_set_idxs[read_name] = 0
                            read_sets[read_name].append(ReadSet(read_name,
                                                            excitation,
                                                            emission, gain))
                    if entered_layout or hit_data:
                        break
            # Read data blocks
            # Find a data block
            while line != None:
                if len(line) == 0:
                    line = next(reader, None)
                    continue
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
                    if read_name not in read_sets:
                        line = next(reader, None)
                        continue
                    if not info.endswith(']'):
                        read_idx = read_set_idxs[read_name]
                        read_set_idxs[read_name] += 1
                    else:
                        read_idx = int(info.split('[')[-1][:-1]) - 1

                    read_channel    = read_sets[read_name]
                    read_properties = read_channel[read_idx]
                    gain            = read_properties.gain
                    print(line)
                    if len(info_parts) > 1:
                        if line[1] != "":
                            excitation = info_parts[1].split("[")[0]
                            emission = line[1].split("[")[0]
                        else:
                            excitation = info_parts[1].split("[")[0].split(",")[0]
                            emission   = info_parts[1].split("[")[0].split(",")[1]
                        excitation = int(excitation)
                        emission   = int(emission)
                    else:
                        excitation      = read_properties.excitation
                        emission        = read_properties.emission

                line = next(reader) # Skip a line
                line = next(reader) # Chart title line
                well_names = line
                # Data lines
                for line in reader:
                    idx = 0
                    while idx < len(line):
                        if line[idx] != "":
                            raw_time = line[idx]
                            break
                        idx += 1
                    if idx == len(line):
                        break
                    time_parts = raw_time.split(':')
                    days=0
                    minutes=int(time_parts[1])
                    seconds=int(time_parts[2])

                    if("1900" in raw_time):
                      tstamp=pd.to_datetime(raw_time)
                      days=tstamp.day
                      hours=int(tstamp.hour)+24*days
                    else:
                      hours=int(time_parts[0])
                    time_secs = int(seconds) + 60*int(minutes) \
                                + 3600*int(hours)

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


def logistic_growth(t, rate, cap, floor, init):
    '''
    Model function for logistic growth with a noise floor.

    Params:
        t -- Time.
        rate -- Growth rate parameter.
        init -- Initial population.
        cap -- Maximum population size.
        floor -- Noise floor (i.e., OD reading for zero cells).
    Returns: Model OD reading at the given time, for the given parameters.
    '''

    rate = np.abs(rate)
    init = np.abs(init)
    cap = np.abs(cap)
    floor = np.abs(floor)
    return floor + cap * init * np.exp(rate * t) \
            / (cap + init * (np.exp(rate * t) - 1))

def summarize_single_well_growth(well_df, growth_threshold = None,
                                 fixed_init = None, verbose = False):
    '''
    Summarizes the growth characteristics of a single well's worth of dataframe,
    returning the results as a dictionary describing a single dataframe line.
    See summarize_growth for measurement details.

    This function is intended as a helper function for summarize_growth; it uses
    the helper function logistic_growth as a growth model.

    Params:
        well_df -- A DataFrame of Biotek data from a single well and a single
                channel (presumably an OD channel).
        growth_threshold -- If set, determines the fraction of of total
                                population to find the time of, i.e., if
                                growth_threshold = 0.25, this function will
                                report the time that each well crossed 25%
                                of the total population size for that well.
        fixed_init -- Sets a fixed value for the initial population parameter.
                        If None, this value is optimized with the rest of the
                        parameters.
        verbose -- Iff True, print some hints about how it's progressing.
    Returns: A dictionary containing the well name, growth characteristics, and
                any supplementary data from df.
    '''
    well_df.reset_index(inplace = True)
    if verbose:
        print("Summarizing from well %s" % well_df.Well[0])

    # Some empirically-reasonable guesses for most growth experiments.
    if fixed_init == None:
        param_guess = (1.3, 1, 0.05, 0)
        opt_func = logistic_growth
    else:
        param_guess = (1.3, 1, 0.05)
        opt_func = lambda t, r, c, f: logistic_growth(t, r, c, f, fixed_init)

    times = well_df["Time (hr)"]

    opt_params = scipy.optimize.curve_fit(opt_func, times,
                                          well_df["Measurement"],
                                          p0 = param_guess,
                                          maxfev = int(1e4))[0]

    # To keep parameters positive, logistic_growth uses the absolute value
    # of whatever parameters it gets, so optimization will sometimes return
    # negative parameter values; have to correct these.
    opt_params = np.abs(opt_params)

    return_dict = dict()
    return_dict["Rate"]  = opt_params[0]
    return_dict["Cap"]   = opt_params[1]
    return_dict["Floor"] = opt_params[2]
    if fixed_init == None:
        return_dict["Init"] = opt_params[3]
    else:
        return_dict["Init"] = fixed_init

    # Calculate threshold time, if it is specified
    if growth_threshold:
        raise NotImplementedError()

    # Add supplemental data.
    for column in well_df.columns.values:
        if column in ["Channel", "Gain", "Time (sec)", "Time (hr)",
                      "Measurement", "Units", "Excitation", "Emission"]:
            continue
        return_dict[column] = well_df[column][0]

    return return_dict

def summarize_growth(df, channel, fixed_init = None, growth_threshold = None,
                     verbose = False):
    '''
    Summarizes the growth characteristics of OD curves from a dataframe of
    Biotek data. Performs the following summaries on each well:
        * Fits OD curve to a logistic-plus-floor, finding an initial value, a
            noise floor, a rate constant (R, not to be interpreted directly),
            and a population maximum. The rate parameter has time units of
            hours.
        * Optionally finds the time when the population crosses some fraction of
            maximum population. By default, does not calculate this -- set
            the growth_threshold parameter to add this calculation.

    Params:
        df -- A DataFrame of Biotek data with at least one channel of growth
                data.
        channel -- The name of the channel with growth data. Should be a channel
                    with only one Gain, or this will do weird things.
        growth_threshold -- If set, determines the fraction of of total
                                population to find the time of, i.e., if
                                growth_threshold = 0.25, this function will
                                report the time that each well crossed 25%
                                of the total population size for that well.
        fixed_init -- Sets a fixed value for the initial population parameter.
                        If None, this value is optimized with the rest of the
                        parameters.
        verbose -- Iff True, prints some hints about what it's doing. Use if it's taking
                    a while and you want to make sure it's making progress. Default False.
    Returns: A new dataframe where each row summarizes the growth
                characteristics of one from the original dataframe.
    '''
    channel_df = df[df.Channel == channel]

    # Split dataframe into a list of dataframes for individual wells.
    well_dfs = [channel_df[channel_df.Well == w] \
                for w in channel_df.Well.unique()]
    measurement_summary_rows = \
        list(map(lambda df: summarize_single_well_growth(df, growth_threshold,
                                                         fixed_init, verbose),
                 well_dfs))
    return pd.DataFrame(measurement_summary_rows)


def window_averages(df, start, end, units = "seconds",
                    grouping_variables = None):
    '''
    Converts a dataframe of fluorescence data to a dataframe of average
    fluorescences from a specified time window. The time window can be specified
    in seconds, hours, or index.

    Params:
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

def moving_average_fit(df, column = "Measurement", window_size = 1,
                       units = "hours", grouping_variables = None):
    '''
    Returns a DataFrame in which measurements (or values from some other column)
    have been smoothed with a moving average filter. The size of the moving
    average window canbe specified in frames, seconds, hours. Window sizes are
    fixed across the time series in *frames*, but not necessarily fixed in
    amount of *time* if the times in the DataFrame are not evenly spaced.

    Whichever column is smooothed will be replaced by the smoothed data.

    Note that smoothing will clip floor(n/2) frames from each end of the data,
    where n is the size of the moving average window in frames.

    Assumes that any data outside of the column to be smoothed is identical
    across time within each well. If this is not true for some column, then the
    values of that column will be overwritten by the value for the earliest
    value of that column in that well.

    Params:
        df -- DataFrame of fluorescence data
        column -- The column to smooth. Default "Measurement", which smooths OD
                    and fluorescence data.
        window_size -- Specifies the number of frames to average over. Can be
                        in units of frames ("index") or time ("hours" or
                        "seconds"). Window size will be rounded up to the
                        nearest odd number.
        units -- Either "hours" (default), "seconds", or "index". If "seconds"
                    or "hours", all measurements must be equally spaced in time.
        grouping_variables - Optional list of column names on which to group.
                                Use this option primarily to separate multiple
                                plates' worth of data with overlapping wells.
    Returns: A DataFrame in which measurements (or some other column) are
                replaced by a moving average of those measurements.
    '''
    group_cols = ["Channel", "Gain","Excitation","Emission", "Well"]
    if grouping_variables:
        group_cols += grouping_variables
    grouped_df = df.groupby(group_cols)
    smoothed_groups = []
    for name, group in grouped_df:
        group = group.sort_values("Time (sec)", ascending = True,
                                  inplace = False)

        # Figure out the window size, in frames.
        if units == "index":
            n_frames = window_size
        else:
            if units == "hours":
                times = group["Time (hr)"]
            elif units == "seconds":
                times = group["Time (sec)"]
            # Time-per-frame can vary substantially between frames. We'll
            # use the median time-difference as an estimate for
            # the time-per-frame.
            time_per_frame = np.median(np.diff(times))
            n_frames = math.ceil(window_size / time_per_frame)
        if n_frames%2 == 0:
            n_frames += 1

        # Calculate smoothed data. Code taken from StackOverflow user Jamie.
        column_data = group[column].to_numpy()
        smoothed_data = np.cumsum(column_data)
        smoothed_data[n_frames:] = \
                            smoothed_data[n_frames:] - smoothed_data[:-n_frames]
        smoothed_data = smoothed_data[n_frames - 1:] / n_frames

        # Clip out data from either end of the group and replace the column of
        # interest with the smoothed version.
        group = group.iloc[n_frames//2 - 1 : -n_frames//2]
        group[column] = smoothed_data

        # This will get concatenated into a total dataframe later.
        smoothed_groups.append(group)
    return pd.concat(smoothed_groups)


def smoothed_derivatives(df, column = "Measurement", window_size = 1,
                         units = "hours", grouping_variables = None):
    '''
    Calculates a smoothed derivative of the time traces in a dataframe. Returns
    a new DataFrame with the measurements in a column replaced by a derivative
    of a moving window average of those measurements.

    Args:
        df - DataFrame of time traces, of the kind produced by tidy_biotek_data.
        column - Column to fit derivatives to. Defaults to "Measurement".
        window_size -- Specifies the number of frames to average over. Can be
                        in units of frames ("index") or time ("hours" or
                        "seconds"). Window size will be rounded up to the
                        nearest odd number.
        units -- Either "hours" (default), "seconds", or "index". If "seconds"
                    or "hours", all measurements must be equally spaced in time.
        grouping_variables - Optional list of column names on which to group.
                                Use this option primarily to separate multiple
                                plates' worth of data with overlapping wells.
    Returns:
        A DataFrame of df augmented with columns for a spline fit and a
        derivative
    '''
    smoothed_df = moving_average_fit(df, column, window_size = window_size,
                             units = units,
                             grouping_variables = grouping_variables)
    group_cols = ["Channel", "Gain","Excitation","Emission", "Well"]
    if grouping_variables:
        group_cols += grouping_variables
    grouped_df = smoothed_df.groupby(group_cols)
    deriv_df   = pd.DataFrame()
    for name, group in grouped_df:
        if column == "Measurement":
            deriv_name = "%s (%s/sec)" % (column, group.Units.unique()[0])
            group.Units = deriv_name
        group = group.copy()
        group[column] = np.gradient(group[column].to_numpy(),
                                    group["Time (sec)"].to_numpy())
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
    OD_channel_string = df[(df.Channel == norm_channel) & \
                 (df.Gain== norm_channel_gain)].ChanStr.unique()[0]
    od_df = df[df.ChanStr == OD_channel_string].reset_index()
    #dflst = []
    channel_list = df.ChanStr.unique().tolist()
    del channel_list[channel_list.index(OD_channel_string)]
    dflist = [df[df.ChanStr == a].reset_index() for a in channel_list]
    normalized_df = od_df.copy()
    for channel_df in dflist:
        channel_df.Measurement = channel_df.Measurement/od_df.Measurement
        orig_units = channel_df.Units.unique()[0]
        norm_units = "OD" if norm_channel.startswith("OD") \
                          else norm_data.Units.unique()[0]
        channel_df.Units = "%s/%s" % (orig_units, norm_units)
        normalized_df = normalized_df.append(channel_df,ignore_index=True)
    normalized_df.reset_index()

    return normalized_df

def apply_by_well(df, summary_function, split_channels = True):
    # README:
    # May be able to rewrite a bunch of other functions using this!!!!!
    '''
    Applies a function over the wells in a DataFrame of Biotek measurements.

    The summary function should be a function that takes one well's worth of
    data and returns a list representing a new row in the summarized DataFrame.

    Args:
        df - The dataframe to be summarized.
        summary_function - A function that summarizes data from a single well,
                            e.g., finds the maximum expression or the median
                            measurement. The summary function should take one
                            well's worth of data, as a DataFrame, and return a
                            list representing one or more new rows in the
                            summary DataFrame that will be returned.
        split_channels - If True, the summary funtion will be applied separately
                            for each channel (for example, when producing a
                            smoothed version of each time trace). If False,
                            the summary function will receive a DataFrame with
                            ALL of the channels for each well (for example,
                            when normalizing measurements against OD). Default
                            True.
    '''
    wells = df.Well.unique()
    if split_channels:
        channels = df.Channel.unique()
        gains    = df.Gain.unique()
    else:
        channels = [None]
        gains    = [None]

    groups = ["Well"]
    if split_channels:
        groups.append("Channel")
        groups.append("Gain")
    grouped_df = df.groupby(groups)
    summarized_list = []
    for name, group in grouped_df:
        summarized_list.append(summary_function(group))
    return pd.DataFrame(summarized_list)


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

    def plot(self, title = None, column = "Measurement", split_plots = False,
             filename = None, show = True, figsize = (8, 4), linewidth = 1,
             ymax = None, od_ymax = None, show_legend = True):
        '''
        Plot/show/save the figure.

        Arguments:
            title -- String that will go at the top of the figure. Default "".
            column -- Column containing the data you want plotted. Defaults to
                        "Measurement", which is standard for fluorescence and
                        OD raw data. You may need to use other columns if you're
                        plotting a DataFrame of derivatives or fits.
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
            linewidth -- Default 1. Passed to plot commands.
            ymax -- Optional flag to set a maximum y limit on the fluorescence
                    axis.
            od_ymax -- Optional flag to set a maximum y limit on the
                        OD axis.
            show_legend -- Boolean that determines whether the legend is
                            displayed.
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
        highest_observed_y = 0
        for well_spec in self.well_list:
            well_df = norm_df[norm_df.Well == well_spec.well_name]
            ax1.plot(well_df["Time (hr)"], well_df[column],
                     color = well_spec.color, label = well_spec.label,
                     linewidth = linewidth)
            highest_observed_y = max(highest_observed_y,
                                     np.nanmax(well_df[column]))
        ax1.set_xlabel("Time (hr)")
        if self.normalize_by_od:
            ax1.set_ylabel("%s/OD (gain %d)" % (self.channel, self.gain))
        else:
            ax1.set_ylabel("%s (gain %d)" % (self.channel, self.gain))

        if split_plots:
            if show_legend:
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
                     color = well_spec.color, linewidth = linewidth,
                     linestyle = linestyle,
                     label = well_spec.label if split_plots else "")
        ax2.set_ylabel(self.od_channel)
        ax2.set_xlabel("Time (hr)")

        if ymax:
            ax1.set_ylim(0, ymax)
        else:
            ax1.set_ylim(0, highest_observed_y * 1.1)
        if od_ymax:
            ax2.set_ylim(0, od_ymax)
        else:
            ax2.set_ylim(0, np.nanmax(all_od_df.Measurement) * 1.1)

        if show_legend:
            if split_plots:
                handles, labels = ax2.get_legend_handles_labels()
            else:
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

def applyFunc(df, inputs, dofunc, output="Calculation", newunits="unknown"):
    '''
    Apply an arbitrary function to over measurements from each channel in a
    dataframe of Biotek data.

    Concatenates the results of that function (another dataframe) as a new set
    of measurements.

    Args:
        df - DataFrame of time traces, of the kind produced by tidy_biotek_data.
        inputs - Name of channels to use as inputs. List of strings.
        dofunc - A function that takes a Series of measurements and returns a
                    dataframe based on those measurements.
        output - Name of the output channel. Can specify new channel or
                    existing.
        newunits - The units of the output channel. Defaults to "unknown"
    Returns:
        A DataFrame of df augmented with columns for whatever your calculation
        was.
    '''
    nout = output
    i = 1
    # the following makes sure that we aren't making duplicate calculation
    # channels
    while(nout in df.Channel.unique()):
        nout = output+str(i)
        i+=1
    output = nout
    #extract out dataframes corresponding to the interesting channels
    indfs = [df[df.Channel == a].reset_index() for a in inputs]
    #make sure they are all the same length!! it won't work otherwise
    testl = len(indfs[0])
    for chi in range(len(inputs)):
        dfl = len(indfs[chi])
        if(dfl != testl):
            raise ValueError("channel '%s' is %s members long, which is " +\
                             "different from %s" % \
                             inputs[chi],dfl,testl)
    # this next part is for building the output. Pretty much copy one of the
    # inputs
    calcdf = indfs[0].copy()
    #now we just take the measurements....
    inmeasures = [a.Measurement for a in indfs]
    #then we apply the function
    calcdf.Measurement = dofunc(inmeasures)
    #set the units and the new channel name
    calcdf.Units = newunits
    calcdf.Channel = output
    outDF = df.append(calcdf,ignore_index=True)
    outDF.drop("index",1)
    return outDF

def hmap_plt(indf,yaxis,xaxis,fixedinds=[],fixconcs=[],construct=None,\
            chan="RFP",axes=None,labels = (1,1), annot=True,vmin=0.6,vmax=0.9,
            cmap="RdBu"):
    '''
    2D heatmap using seaborn heatmap and matplotlib pivot_table. If you don't
    specify any columns, then it will by default average them!!
    Args:
        indf - DataFrame of endpoint data, of the kind produced by
                tidy_biotek_data.
        yaxis - column name for y axis
        xaxis - column name for x axis
        fixedinds - list of fixed inducer names.
        fixconcs - list of fixed concentrations, corresponding to fixedinds.
        construct - name of construct to use
        chan - channel name to plot!
        axes - matplotlib axes to plot onto
        labels - this list determines whether the x or y axis labels are
                displayed
                default: both axis labels are displayed
        hmset - settings for the heatmap!
                values: annot, vmin, vmax, cmap
    Returns:
        Nothing. Creates a matplotlib plot
    '''
    logiclist = (indf.Channel == chan)
    if((construct == None) and not ("Construct" in fixedinds)):

        conpick = indf.Construct.unique()[0]
        logiclist = logiclist&(indf.Construct==conpick)
        print("you didn't specify a construct so we are showing "+conpick)
    for fixedstf in zip(fixedinds,fixconcs):
        logiclist = logiclist&(indf[fixedstf[0]]==fixedstf[1])
    slicedf = indf[logiclist]
    #print(slicedf.head())
    pivtabl = pd.pivot_table(slicedf, index=yaxis, columns=xaxis, values='Measurement')
    #print(pivtabl)

    outax = sns.heatmap(pivtabl,\
                annot=annot,\
                vmin = vmin,\
                vmax = vmax,\
                cmap= cmap, \
                cbar = False, \
                ax = axes)
    if(not labels[0]):
        outax.get_xaxis().set_visible(False)
    if(not labels[1]):
        outax.get_yaxis().set_visible(False)


def multiPlot(dims_in,plotdf,fixedinds_in,fixconcs_in,constructs,FPchan,\
                                annot=False,vmin=None,vmax=None,cmap="RdBu"):
    '''
    3D or 4D heatmap plot of dataframe
    Args:
        dims - List of inducer names to vary. Max of four!
        plotdf - DataFrame of endpoint data, of the kind produced by tidy_biotek_data.
        fixedinds - list of fixed inducer names.
        fixconcs - list of fixed concentrations, corresponding to fixedinds.
        constructs - list of constructs to plot. You can put multiple constructs
                    to plot here, but then dims should be only three elements long,
                    since constructs constitutes another dimension
        FPchan - channel name to plot!
    Returns:
        Nothing. Creates a matplotlib plot
    '''
    dims = dc(dims_in)
    fixedinds = dc(fixedinds_in)
    fixconcs = dc(fixconcs_in)
    def allcomb(listoflists):
        """creates all paths through a list"""
        if(len(listoflists)==1):
            return([[a] for a in listoflists[0]])
        outlist = []
        for element in listoflists[0]:
            for lists in allcomb(listoflists[1:]):
                outlist+=[[element]+lists]
        return outlist
    dimlist = [sorted(plotdf[a].unique()) for a in dims]
    logiclist = (plotdf.Channel == FPchan)
    if(len(fixedinds)>0):
        for fixedval in zip(fixedinds,fixconcs):
            logiclist=logiclist&(plotdf[fixedval[0]]==fixedval[1])
    loglist2 = plotdf["Construct"]==constructs[0]
    if(len(constructs)>1):
        for cons in constructs[1:]:
            loglist2=loglist2|(plotdf["Construct"]==cons)
    logiclist = logiclist&loglist2
    subdf = plotdf[logiclist]
    #print(subdf.head())
    maxval = max(subdf.Measurement)
    minval = min(subdf.Measurement)
    if(vmin == None):
        vmin = minval
    if(vmax == None):
        vmax = maxval
    if(len(constructs)>1):
        dimlist += [constructs]
        dims+= ["Construct"]
    else:
        fixedinds += ["Construct"]#,"ATC"]
        fixconcs += constructs#,250]
    enddims = range(len(dims[2:])) #cut out the first two dimensions
    axiscombs = allcomb(dimlist[2:]) #all combinations of conditions
    plotpos = allcomb([range(len(a)) for a in dimlist[2:]]) #positions on the graph grid that above will go
    rows = len(dimlist[2])
    if(len(dimlist)>=4):
        cols = len(dimlist[3])
    else:
        cols = 1
    fig, axes = plt.subplots(cols, rows,figsize=(rows*2,cols*2))
    #four dimensions is about the best we can do
    #this next part populates the axes list with blanks so that the
    #plotting code still runs like a 4x4 figure
    if(cols == 1):
        axes = [axes]
    if(len(plotpos[0])<2):
        plotpos = [[a[0],0] for a in plotpos]


    fig.subplots_adjust(hspace=0.05,wspace=0.05)
    #this next part creates the title
    titstr = FPchan +" " #what channel we are plotting
    #this next part determines what inducers there are
    #and which ones weren't specified.
    #the pivottable function actually averages those channels
    #so it would be nice for the user to know that!!
    defaultColumns = ["Channel","Gain","Excitation",\
                      "Emission","Well","ChanStr",\
                      "Construct","Measurement",\
                      "Time (hr)","Time (sec)",\
                      "Units","index","level_0",\
                      "level_1","level_2","level_3"]
    inducerCols = []
    for c in list(subdf.loc[:,(subdf != 0).any(axis=0)].columns):
        if(not (c in defaultColumns)):
            inducerCols+=[c]
    notspecified = []
    for c in inducerCols:
        if(not(c in dims+fixedinds)):
            notspecified+=[c]
    for ind in zip(fixedinds,fixconcs):
        titstr+= "; "+ind[0]+" is "+str(ind[1])
    if(len(notspecified)>0):
        for nsind in notspecified:
            titstr+= "; aggregate of " +nsind
    fig.suptitle(titstr, fontsize=16)

    rowcount = 0
    for colax in axes:
        colcount = 0
        for curax in colax:
            #iterate through each plot basically

            #curax = fourthaxis[2]
            plotlocation = [colcount,rowcount]
            #print(plotlocation)

            if(plotlocation in plotpos):#if we are populating this plot:
                plotind = plotpos.index(plotlocation)
                curconcs = axiscombs[plotind]

                hmap_plt(plotdf,dims[0],dims[1],\
                         fixedinds+dims[2:],\
                         fixconcs+curconcs,\
                         None,\
                         FPchan,curax,\
                        (rowcount==len(axes)-1,colcount==0),\
                        annot=annot,\
                        vmin=vmin,\
                        vmax=vmax,\
                        cmap=cmap)
                if(len(dims) == 4):
                    if(dims[2:][1]=="Construct"):
                        curax.set_ylabel(str(curconcs[1])+"\n"+dims[0])
                    else:
                        curax.set_ylabel(str(curconcs[1])+" of "+dims[2:][1]+"\n"+dims[0])
                if(dims[2:][0]=="Construct"):
                    curax.set_xlabel(dims[1]+"\n"+str(curconcs[0]))
                else:
                    curax.set_xlabel(dims[1]+"\n"+str(curconcs[0])+" of "+dims[2:][0])
            colcount+=1
        rowcount+=1
CellSpec = namedtuple("CellSpec", ['well_name', 'color', 'label'])
