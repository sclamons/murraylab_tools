from echo_functions import *
import os

example_loc = 'example_projects'

#####################
# Association Lists #
#####################
stock_filename = os.path.join(example_loc, 'association_list',
                              'source_sheet.csv')
association_filename = os.path.join(example_loc, 'association_list',
                                    'association_sheet.csv')
aa_picklist_filename = os.path.join(example_loc, 'association_list_test',
                                    'assocation_test'))
aa_echo_calculator = EchoRun()
aa_echo_calculator.rxn_vol = 50000
name_col  = 'B'
conc_col  = 'C'
len_col   = 'D'
well_col  = 'A'
plate_col = 'E'
aa_echo_calculator.load_source_plate(stock_filename, name_col, conc_col,
                                     len_col, well_col, plate_col)
aa_echo_calculator.build_picklist_from_association_spreadsheet(
                                                    association_filename,
                                                    'A')
aa_echo_calculator.write_picklist(aa_picklist_filename)


###########################
# TX-TL Setup Spreadsheet #
###########################
plate_file = "tx_tl_test/test_plate.dat"
txtl_plate = SourcePlate(filename = plate_file)
txtl_echo_calculator = EchoRun(plate = txtl_plate)

stock_filename  = os.path.join(example_loc, 'tx_tl_test',
                               'TXTL_eSC3_template_2.1_stocks.csv')
recipe_filename = os.path.join(example_loc, 'tx_tl_test',
                               'TXTL_eSC3_template_2.1_recipe.csv')
tx_tl_picklist_filename = os.path.join(example_loc, 'tx_tl_test', 'example')
txtl_echo_calculator.build_picklist_from_txtl_setup_csvs(stock_filename,
                                                    recipe_filename)
txtl_echo_calculator.write_picklist(tx_tl_picklist_filename)

######################
# 2D Dilution Series #
######################
dilution_plate_file = os.path.join(example_loc, "dilution_series",
                                   "test_plate.dat")
dilution_plate = SourcePlate(filename = dilution_plate_file)
dilution_echo_calculator = EchoRun(plate = dilution_plate)

dna1_final = range(0,6,1) # in nM
dna2_final = range(0,6,1) # in nM

rxn_vol = 5
# Integrase plasmid
dna1_conc = 300 # 1704 ng/ul stock
dna1_len = 4524 # bp
dna1 = EchoSourceMaterial('Integrase', dna1_conc, dna1_len, ["C1"])

# Reporter plasmid
dna2_conc = 300 #564 ng/ul stock
dna2_len = 4656 # bp
dna2 = EchoSourceMaterial("Reporter", dna2_conc, dna2_len, ["C2"])

dilution_echo_calculator.build_dilution_series(dna1, dna2, dna1_final,
                                                dna2_final, "D2")
dilution_picklist_filename = os.path.join(example_loc, 'dilution_series',
                                          'bxb1')
dilution_echo_calculator.write_picklist(dilution_picklist_filename)

# TP901
tp901_conc = 300 #1017 ng/ul stock
tp901_len = 4521 # bp
tp901 = EchoSourceMaterial("TP901", tp901_conc, tp901_len, ["D1"])

tp901r_conc = 300 #492 ng/ul stock
tp901r_len = 4538 # bp
tp901r = EchoSourceMaterial("TP901r", tp901r_conc, tp901r_len, ["D2"])

dilution_plate_2 = SourcePlate(filename = dilution_plate_file)
dilution_echo_calculator_2 = EchoRun(plate = dilution_plate_2)
dilution_echo_calculator_2.build_dilution_series(tp901, tp901r, dna1_final,
                                                 dna2_final, "D2")
dilution_picklist_filename_2 = os.path.join(example_loc, 'dilution_series',
                                            'tp901')
dilution_echo_calculator_2.write_picklist()
