import numpy as np
import gvar
import csv
from pathlib import Path
from typing import List
from common import Point, Points_Same_Type

def output_T025_WL_table_csv(display_Ls):
    T=0.25
    S = 24
    P_FILLER = 0.5
    ns = 16, 24, 32, 48, 64, 96, 128
    W_TYPE = ("dec")

    with open(Path("products", "tables", "T025.csv"), "w+") as f:
        header = ["N"] + [f"W({l})" for l in display_Ls]
        writer = csv.writer(f)
        writer.writerow(header)
        for i, n in enumerate(ns):
            pt = Point(n, n, P_FILLER, T, S, W_TYPE)
            line = [n] + [str(gvar.gvar(pt.WL_data_av[l-1], pt.WL_errs[l-1])) for l in display_Ls]
            writer.writerow(line)
        

def main():
    display_Ls = 1,2,3,4,5
    output_T025_WL_table_csv(display_Ls)

if __name__ == "__main__":
    main()