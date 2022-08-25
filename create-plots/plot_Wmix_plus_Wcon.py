import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import gvar
from typing import List
from common import Points_Same_Type, Point, Points_Combine_MixCon
import plots_common as p_common

plt.rcParams["figure.figsize"] = (10,10 / 1.618)

def plot_Wmix_plus_Wcon(pointset_mixcon):
    T = 0.29
    BETA = 1 / T
    fig, ax = plt.subplots()

    Ls = pointset_mixcon.Ls

    ax.scatter(Ls, pointset_mixcon.WL_continuum_largeN_together, marker='x', s=80, c='#0000a7')
    ax.errorbar(Ls, pointset_mixcon.WL_continuum_largeN_together, pointset_mixcon.WL_cln_together_err, linestyle='none', capsize=4, c='#0000a7')

    ypred_xs = np.arange(0, Ls[-1])
    ypred = np.exp(- BETA * ypred_xs / 2)
    ax.plot(ypred_xs, ypred, c=p_common.COLOURS['pred'])
    ax.set_yscale('log')
    
    USE_POINTS = [0,1,2,3]
    gW = gvar.gvar(pointset_mixcon.WL_continuum_largeN_together, pointset_mixcon.WL_cln_together_err)
    y_for_fit = np.log(gW[USE_POINTS])
    (m, c), cov = np.polyfit(Ls[USE_POINTS], gvar.mean(y_for_fit), w=1/gvar.sdev(y_for_fit), deg=1, cov='unscaled')

    m_err, c_err = np.sqrt(np.diag(cov))
    g_m = gvar.gvar(m, m_err)
    g_c = gvar.gvar(c, c_err)
    with open(Path(p_common.SAVE_NUMBERS_FOLDER, f"P{str(pointset_mixcon.p).replace('.','')}_WconPlusMix_2D.txt"),"w+") as f:
        f.write(f"c={-g_m}\nd={g_c}")

    ylabel = "$\\frac{{1}}{{N-M}} (W(L)_{{\\rm con}} + W(L)_{{\\rm mix}})$" # TODO: check it is actually this...

    xlim = np.array((0.5, 6.5))
    ylim = np.array((10**-6, 0.5))

    xl_points = np.linspace(xlim[0], xlim[1], 100)
    yfit = np.exp(g_c + g_m * xl_points)
    ax.fill_between(xl_points, gvar.mean(yfit)-gvar.sdev(yfit), gvar.mean(yfit)+gvar.sdev(yfit), alpha=0.3, color='#424287')#, c=COLOURS['fit'])
    ax.plot(xl_points, gvar.mean(yfit), c='#424287')#, linestyle=':')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel("$L$", fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.grid(linestyle='--')

    return fig

def main():
    # TODO: DRY with plot_2D_interpolation.py ??
    T = 0.29
    points:List[Point] = []
    w_types = ("con", "mix", "dec")
    nmp = ((24,11,0.2), (24,15,0.25), (32,15,0.2), (32,20,0.25), (64,29,0.2), (64,40,0.25)) 
    # Ss = (16,24,32)
    Ss = (24,32)
    for i, (n, m, p) in enumerate(nmp):
        for j, s in enumerate(Ss):
            for k, w_type in enumerate(w_types):
                pt = Point(n, m, p, T, s, w_type)
                points.append(pt)
    WL_con_P020 = Points_Same_Type([pt for pt in points if pt.p==0.2 and pt.WL_type=="con"])
    WL_mix_P020 = Points_Same_Type([pt for pt in points if pt.p==0.2 and pt.WL_type=="mix"])
    WL_con_P025 = Points_Same_Type([pt for pt in points if pt.p==0.25 and pt.WL_type=="con"])
    WL_mix_P025 = Points_Same_Type([pt for pt in points if pt.p==0.25 and pt.WL_type=="mix"])

    WL_mix_con_P020 = Points_Combine_MixCon(WL_mix_P020.all_points, WL_con_P020.all_points)
    WL_mix_con_P025 = Points_Combine_MixCon(WL_mix_P025.all_points, WL_con_P025.all_points)
    fig = plot_Wmix_plus_Wcon(WL_mix_con_P020)
    fig.savefig(Path(p_common.SAVE_FOLDER,f'WL_mixPlusCon_P020.pdf'))
    fig = plot_Wmix_plus_Wcon(WL_mix_con_P025)
    fig.savefig(Path(p_common.SAVE_FOLDER,f'WL_mixPlusCon_P025.pdf'))

if __name__ == "__main__":
    main()