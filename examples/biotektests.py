import murraylab_tools.biotek as mt_biotek
import os
gitexamplepath = "C:\\Users\\Andrey\\Documents\\GitHub\\"+\
                "murraylab_tools\\examples\\biotek_examples\\"
data_filename = gitexamplepath+\
                "cell_data_traces.csv"
#os.path.join("biotek_examples", "RFP_GFP_traces.csv")
mt_biotek.tidy_biotek_data(data_filename, volume = 5,convert_to_uM=False)

import pandas as pd

tidy_filename = gitexamplepath+"cell_data_traces_tidy.csv"
df = pd.read_csv(tidy_filename)
df.head()
normdf = mt_biotek.normalize(df)
normdf.head()
