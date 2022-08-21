import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import gvar
from typing import List
from common import Points_Same_Type, Point
import plots_common as p_common

plt.rcParams["figure.figsize"] = (10,10 / 1.618)

def plot_2D_interpolation(pointset):
    # T = 0.29
    # BETA = 1 / T
    fig, ax = plt.subplots()
    ax.scatter(pointset.Ls, pointset.WL_continuum_largeN_together, marker='x', s=80, c='#0000a7')
    ax.errorbar(pointset.Ls, pointset.WL_continuum_largeN_together, yerr=pointset.WL_cln_together_err, linestyle='none', capsize=4, c='#0000a7')

    # ypred_xs = np.arange(0, pointset.Ls[-1])
    # ypred = np.exp(- BETA * ypred_xs / 2)
    # ax.plot(ypred_xs, ypred, c=COLOURS['pred'])
    ax.set_yscale('log')
    
    USE_POINTS = [0,1,2,3]
    (m, c), cov = np.polyfit(pointset.Ls[USE_POINTS], np.log(pointset.WL_continuum_largeN_together[USE_POINTS]), w=pointset.WL_continuum_largeN_together[USE_POINTS]/pointset.WL_cln_together_err[USE_POINTS], deg=1, cov='unscaled')

    m_err, c_err = np.sqrt(np.diag(cov))
    g_m = gvar.gvar(m, m_err)
    g_c = gvar.gvar(c, c_err)
    with open(Path(p_common.SAVE_NUMBERS_FOLDER, f"P{str(pointset.p).replace('.','')}_{pointset.WL_type}_fit_2Dinterp.txt"),"w+") as f:
        f.write(f"c={-g_m}\nd={g_c}")

    if pointset.p > 0:
        if pointset.WL_type == "con":
            ylabel = "$\\frac{{1}}{{N(1-\\frac{{M}}{{N}})^2}} W(L)_{{\\rm con}}$"
        elif pointset.WL_type == "mix":
            ylabel = "$\\frac{{1}}{{M(1-\\frac{{M}}{{N}})}} W(L)_{{\\rm mix}}$"
        elif pointset.WL_type == "dec": # maybe will have to except for P=0
            ylabel = "$\\frac{{N}}{{M^2}} W(L)_{{\\rm dec}}$"
        else:
            ylabel = "?"
    else:
        ylabel="TODO label"

    xlim = np.array((0.5, 6.5))
    ylim = np.array((10**-6, 0.5))

    xl_points = np.linspace(xlim[0], xlim[1], 100)
    if pointset.WL_type != "dec":
        yfit = np.exp(g_c + g_m * xl_points)
        # ax.fill_between(xl_points, gvar.mean(yfit-gvar.sdev(yfit)), gvar.mean(yfit+gvar.sdev(yfit)), alpha=0.3, color='#424287')#, c=COLOURS['fit'])
        ax.fill_between(xl_points, gvar.mean(yfit)-gvar.sdev(yfit), gvar.mean(yfit)+gvar.sdev(yfit), alpha=0.3, color='#424287')#, c=COLOURS['fit'])
        ax.plot(xl_points, gvar.mean(yfit), c='#424287')#, linestyle='--')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel("$L$", fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.grid(linestyle='--')

    return fig


def main():
    T = 0.29
    points:List[Point] = []
    w_types = ("con", "mix", "dec")
    nmp = ((24,11,0.2), (24,15,0.25), (32,15,0.2), (32,20,0.25), (64,29,0.2), (64,40,0.25)) 
    Ss = (16,24,32)
    for i, (n, m, p) in enumerate(nmp):
        for j, s in enumerate(Ss):
            for k, w_type in enumerate(w_types):
                pt = Point(n, m, p, T, s, w_type)
                points.append(pt)
    WL_con_P020 = Points_Same_Type([pt for pt in points if pt.p==0.2 and pt.WL_type=="con"])
    WL_mix_P020 = Points_Same_Type([pt for pt in points if pt.p==0.2 and pt.WL_type=="mix"])
    WL_dec_P020 = Points_Same_Type([pt for pt in points if pt.p==0.2 and pt.WL_type=="dec"])
    WL_con_P025 = Points_Same_Type([pt for pt in points if pt.p==0.25 and pt.WL_type=="con"])
    WL_mix_P025 = Points_Same_Type([pt for pt in points if pt.p==0.25 and pt.WL_type=="mix"])
    WL_dec_P025 = Points_Same_Type([pt for pt in points if pt.p==0.25 and pt.WL_type=="dec"])

    all_point_sets = [WL_con_P020, WL_mix_P020, WL_dec_P020, WL_con_P025, WL_mix_P025, WL_dec_P025]

    for ps in all_point_sets:
        fig = plot_2D_interpolation(ps)
        fig.savefig(Path(p_common.SAVE_FOLDER, f'largeN_continuum_2D_P{str(ps.p).replace(".","")}_{ps.WL_type}.pdf'))

if __name__ == "__main__":
    main()