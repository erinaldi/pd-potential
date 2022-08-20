import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import gvar
from typing import List
from common import Points_Same_Type, Point, Points_Combine_MixCon
import plots_common as p_common
import csv


def get_fit_coeffs(pointset_mixcon):
    Ls = pointset_mixcon.Ls
    USE_POINTS = [0,1,2,3]

    gW = gvar.gvar(pointset_mixcon.WL_continuum_largeN_together, pointset_mixcon.WL_cln_together_err)
    y_for_fit = np.log(gW[USE_POINTS])
    (c, d), cov = np.polyfit(Ls[USE_POINTS], gvar.mean(y_for_fit), w=1/gvar.sdev(y_for_fit), deg=1, cov='unscaled')

    c_err, d_err = np.sqrt(np.diag(cov))
    g_c = -gvar.gvar(c, c_err)
    g_d = gvar.gvar(d, d_err)

    return g_c, g_d
    

def main():
    T = 0.29
    points:List[Point] = []
    w_types = ("con", "mix")
    nmp = ((16,8,0.2), (16,10,0.25), (24,11,0.2), (24,15,0.25), (32,15,0.2), (32,20,0.25), (64,29,0.2), (64,40,0.25)) 
    Ss = (16,24,32)
    for i, (n, m, p) in enumerate(nmp):
        for j, s in enumerate(Ss):
            for k, w_type in enumerate(w_types):
                pt = Point(n, m, p, T, s, w_type)
                points.append(pt)

    remove_Ns = (None, 16, 24, 32, 64)
    
    c_cons_p02 = {}
    c_cons_p025 = {}

    for rmvN in remove_Ns:
        WL_con_P020 = Points_Same_Type([pt for pt in points if pt.p==0.2 and pt.WL_type=="con" and pt.n != rmvN])
        WL_mix_P020 = Points_Same_Type([pt for pt in points if pt.p==0.2 and pt.WL_type=="mix" and pt.n != rmvN])
        WL_con_P025 = Points_Same_Type([pt for pt in points if pt.p==0.25 and pt.WL_type=="con" and pt.n != rmvN])
        WL_mix_P025 = Points_Same_Type([pt for pt in points if pt.p==0.25 and pt.WL_type=="mix" and pt.n != rmvN])

        WL_mix_con_P020 = Points_Combine_MixCon(WL_mix_P020.all_points, WL_con_P020.all_points)
        WL_mix_con_P025 = Points_Combine_MixCon(WL_mix_P025.all_points, WL_con_P025.all_points)

        c_con_p02, _ = get_fit_coeffs(WL_mix_con_P020)
        c_con_p025, _ = get_fit_coeffs(WL_mix_con_P025)
        
        c_cons_p02[rmvN] = c_con_p02
        c_cons_p025[rmvN] = c_con_p025

    header = ["P"] + [f"rmv. $N={rmvN}$" if rmvN is not None else "all $N$" for rmvN in remove_Ns]
    with open(Path("products", "tables", "c_con_removeN.csv"), "w+") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        p02_line = ["$0.2$"] + [c_cons_p02[rmvN] for rmvN in remove_Ns]
        p025_line = ["$0.25$"] + [c_cons_p025[rmvN] for rmvN in remove_Ns]
        writer.writerow(p02_line)
        writer.writerow(p025_line)

if __name__ == "__main__":
    main()