# coding=utf-8

# TODO:
#       --Add optional user-set source wells
#       --Positive control reaction on 2D dilution series
#           * Alternatively, add manual positive control reaction function?
#

# Written by Victoria Hsiao, August 2016
# Modified by Samuel Clamons, October 2016

# EchoCSVMaker encapsulates several algorithms for writing .csv files that can
# be directly read by the Echo liquid handler, along with .txt files describing
# what needs to go in the source plate.
# Current features include:
#   -- Conversion from a pair of CSV templates
#   -- Automatic setup of an NxN dilution matrix with two inputs.
#   -- Picks from plates defined by CSVs, to wells defined by other CSVs



#import openpyxl as pyxl  # Required for reading excel files

from echo_source_material import EchoSourceMaterial
from reaction import Reaction, WellReaction, MasterMix
from destination_plate import DestinationPlate
from source_plate import SourcePlate
from pick import Pick
from echo_run import EchoRun
from echo_functions import dna2nM_convert, echo_round, floatify, \
                           process_column_argument














