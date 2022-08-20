import numpy as np
import gvar
import csv
from pathlib import Path
from typing import List
from common import Point, Points_Same_Type

# we will export as a csv and then use a tool within LaTeX to put this in the right format

def output_WL_table_csv(all_points, display_ls):
    # create table for each wl_type
    # sort as follows: P-> N -> S

    WL_types = set([pt.WL_type for pt in all_points])
    points_by_wl_type = {wlt : [pt for pt in all_points if pt.WL_type == wlt] for wlt in WL_types}

    for wlt, wlt_points in points_by_wl_type.items():
        # maybe I should just do for all points, and then sort afterwards...
        sorted_points = sorted(wlt_points, key = lambda pt: (pt.p, pt.n, pt.n_t))

        header = ["P", "N", "S"] + [f"W({l})" for l in display_ls]
        last_p, last_n = None, None
        with open(Path("products", "tables", f"{wlt}.csv"), "w+") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            for pt in sorted_points:
                line = [pt.p if pt.p != last_p else "", pt.n if pt.n != last_n else "", pt.n_t] + [str(gvar.gvar(pt.WL_data_av[l-1], pt.WL_errs[l-1])) for l in display_ls]
                writer.writerow(line)
                last_p, last_n = pt.p, pt.n


def main():
    T=0.29
    points:List[Point] = []
    w_types = ("con", "mix", "dec")
    nmp = ((16,8,0.2), (16,10,0.25), (24,11,0.2), (24,15,0.25), (32,15,0.2), (32,20,0.25), (64,29,0.2), (64,40,0.25)) 
    Ss = (16,24,32)
    for i, (n, m, p) in enumerate(nmp):
        for j, s in enumerate(Ss):
            for k, w_type in enumerate(w_types):
                pt = Point(n, m, p, T, s, w_type)
                points.append(pt)
    display_Ls = (1,2,3,4,5)
    output_WL_table_csv(points, display_Ls)

if __name__ == "__main__":
    main()










# def create_latex_tables(nmp, w_types, w_avs, errs, valid, T):
#     # could make it adaptable for arbitrary l
#     tables = {}
#     for k, w_type in enumerate(w_types):
#         header = r"""\begin{table}
# \centering
# \begin{tabular}{|c|c|c||l|l|l|l|l|}
#     \hline
#     $N$ &\ $P$ &\ $S$ &\ \quad$W(1)$ & \ \quad$W(2)$& \ \quad$W(3)$& \ \quad$W(4)$& \ \quad$W(5)$\\
# \hline
# \hline
#     """
        
#         last_n = None
#         last_p = None 
#         table_core = ""

#         for i, (n, m, p) in enumerate(nmp):
#             for j, s in enumerate(Ss):
#                 if n == last_n:
#                     table_core += " & "
#                 else:
#                     table_core += f" ${n}$ & "
#                 if p == last_p:
#                     table_core += " & "
#                 else:
#                     table_core += f" ${p}$ & "
#                 last_n = n
#                 last_p = p
#                 table_core += f" ${s}$ "
#                 for l in Ls:
#                     wav = w_avs[i,j,k, l-1]
#                     err = errs[i,j,k,l-1]
#                     table_core += f" & ${wav:#.3} \pm {err:#.2}$"
#                     if not valid[i,j,k, l-1]:
#                         table_core+= "**"
#                 table_core += r" \\" + "\n  "
#             footer = r"""  \hline
#     \hline
# \end{tabular}
# \caption{""" + f"$W_{{\\rm {w_type}}}$ for $L=1,2,3,4,5$, constrained, $T={T}$." + "} \n\end{table}"
#         tables[w_type] = header + table_core + footer
#     return tables