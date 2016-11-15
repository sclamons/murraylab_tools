from biotek import *
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    filename = "161103_reporter_titration_results_tidy.csv"
    df = pd.read_csv(filename)
    df = background_subtract(df, ["E10", "E15", "K5"])