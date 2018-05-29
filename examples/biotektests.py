import murraylab_tools.biotek as mt_biotek
import os
import pandas as pd
import seaborn as sns

gitexamplepath =".\\examples\\biotek_examples\\"
data_filename = gitexamplepath+\
                "180515_big384wellplate.csv"
supplementary_filename = gitexamplepath+\
                "supp_inductiongrid.csv"

#mt_biotek.tidy_biotek_data(data_filename, supplementary_filename, convert_to_uM = False)



tidy_filename = gitexamplepath+"180515_big384wellplate_tidy.csv"
df = pd.read_csv(tidy_filename)
df.loc[df.Excitation==580,"Channel"] = "RFP"
df.loc[df.Excitation==485,"Channel"] = "GFP"

normdf = mt_biotek.normalize(df,norm_channel= "OD")
normdf = normdf.drop("index",axis=1)
calcdf = mt_biotek.applyFunc(normdf,("GFP","RFP"),lambda x:x[0]/(x[0]+x[1]))
calcdf = calcdf.drop("index",axis=1)
#normdf[normdf.Gain==100].head()
end_df = mt_biotek.window_averages(normdf,15,17,"hours")
end_df.aTC.unique()
dims = ["Ara","IPTG"]#,"ATC"]
fixedinds = ["aTC"]#,"Construct"]
fixconcs = [250]
plotdf = end_df
FPchan = "RFP"
constructs = ["pQi41","pQi42","pQi51","pQi52"]

mt_biotek.multiPlot(dims,end_df,fixedinds,fixconcs,constructs,FPchan,annot=True,vmin=None,vmax=None)


end_df.Excitation.unique()
slicedf = end_df[(end_df.Gain == 100 )&(end_df.Construct=="pQi41")&(end_df.aTC==250)]
end_df[(end_df.Gain == 100 )&(end_df.Construct=="pQi41")&(end_df.aTC==250)].head()
