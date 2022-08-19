import numpy as np
import pandas as pd
from pathlib import Path

# # TODO: this should not just be about running plots, but should also output all the data into tables too
# # 1) get thermalisation estimates from autocorr time * 10
# # 2) plot action and other histories of important observables, indicating thermalisation cutoff, perhaps

# TODO:
# rename python files etc
# remove superfluous stuff, maybe, in common.py
# think about how to deal with thermalisation data... 1) read from file 2) calculate autocorrelation data and then strip...
# what do I need to do to make the tables, and where should I put it?
# 

filename=Path("data-files/N16/S16/M8/T029/P02/summary.csv")
res=pd.read_csv(filename,header=0,index_col=0)

