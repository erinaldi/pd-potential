import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import gvar
from typing import List
from common import Points_Same_Type, Point, Points_Combine_MixCon
import plots_common as p_common

plt.rcParams["figure.figsize"] = (10,10 / 1.618)

def plot_why_remove_N16(pointset, show_L):
    large_nt_WL = {}
    large_nt_WL_err = {}
    ms = {}
    ns = []

    for n, pts in pointset.points_by_N.items():
        nts = [pt.n_t for pt in pts]
        WL_avs = np.array([pt.WL_data_av for pt in pts])
        WL_errs = np.array([pt.WL_errs for pt in pts]) # TODO: WL_errs is not built-in to the object! But it is important enough to be
        cs_this_N = np.zeros(pointset.num_Ls)
        errs_this_N = np.zeros(pointset.num_Ls)
        for i, l in enumerate(pointset.Ls):
            c, err = pointset._take_large_nt(nts, WL_avs[:,i], WL_errs[:,i])
            cs_this_N[i] = c
            errs_this_N[i] = err
        large_nt_WL[n] = cs_this_N
        large_nt_WL_err[n] = errs_this_N # indexed by n_t; gives W[i]
        ms[n] = pts[0].m

    ns = pointset.points_by_N.keys()

    xs = []
    ys = []
    y_errs= []
    for n in ns:
        m = ms[n]
        xs.append(1 / n / (n-m))
        ys.append(large_nt_WL[n])
        y_errs.append(large_nt_WL_err[n])
    
    xs_reduced = []
    ys_reduced = []
    y_errs_reduced = []
    ns_reduced = [n for n in ns if n!=16]
    for n in ns_reduced:
        m = ms[n]
        xs_reduced.append(1 / n / (n-m))
        ys_reduced.append(large_nt_WL[n])
        y_errs_reduced.append(large_nt_WL_err[n])

    xs = np.array(xs)
    ys = np.array(ys)
    y_errs = np.array(y_errs)
    xs_reduced = np.array(xs_reduced)
    ys_reduced = np.array(ys_reduced)
    y_errs_reduced = np.array(y_errs_reduced)

    show_L = 1
    show_ys = ys[:,show_L-1]
    show_y_errs = y_errs[:,show_L-1]

    show_ys_r = ys_reduced[:,show_L-1]
    show_y_errs_r = y_errs_reduced[:,show_L-1]

    fig, ax = plt.subplots()
    ax.scatter(xs, show_ys, marker='x', s=80, zorder=10)
    ax.errorbar(xs, show_ys, yerr=show_y_errs, linestyle='none', capsize=4, zorder=10)

    USE_POINTS = [0,1,2,3]

    (c2_r, c1_r), ccov_r = np.polyfit(xs_reduced, show_ys_r, w=1/show_y_errs_r, deg=1, cov='unscaled')
    c2_err_r, c1_err_r = np.sqrt(np.diag(ccov_r))

    gc2r = gvar.gvar(c2_r, c2_err_r)
    gc1r = gvar.gvar(c1_r, c1_err_r)

    (c2, c1), ccov = np.polyfit(xs, show_ys, w=1/show_y_errs, deg=1, cov='unscaled')
    c2_err, c1_err = np.sqrt(np.diag(ccov))

    gc2 = gvar.gvar(c2, c2_err)
    gc1 = gvar.gvar(c1, c1_err)

    (quad3, quad2, quad1), qcov = np.polyfit(xs, show_ys, w=1/show_y_errs, deg=2, cov='unscaled')
    q3_err, q2_err, q1_err = np.sqrt(np.diag(qcov))

    gq3 = gvar.gvar(quad3, q3_err)
    gq2 = gvar.gvar(quad2, q2_err)
    gq1 = gvar.gvar(quad1, q1_err)

    x_lin = np.linspace(0,np.amax(xs))
    plt.plot(x_lin, c2 * x_lin + c1, label="linear all $N$", color='green')
    plt.plot(x_lin, c2_r * x_lin + c1_r, label="linear no $N=16$", color='red')
    plt.plot(x_lin, quad3 * x_lin ** 2 + quad2 * x_lin + quad1, label="quadratic", color='orange')

    yq = gq3 * x_lin ** 2 + gq2 * x_lin + gq1
    ycr = gc2r * x_lin + gc1r
    yc = gc2 * x_lin + gc1

    ax.fill_between(x_lin, gvar.mean(yc)-gvar.sdev(yc), gvar.mean(yc)+gvar.sdev(yc), alpha=0.3, color='green')#, color='#424287')#, c=COLOURS['fit'])
    ax.fill_between(x_lin, gvar.mean(ycr)-gvar.sdev(ycr), gvar.mean(ycr)+gvar.sdev(ycr), alpha=0.3, color='red')#, color='#424287')#, c=COLOURS['fit'])
    ax.fill_between(x_lin, gvar.mean(yq)-gvar.sdev(yq), gvar.mean(yq)+gvar.sdev(yq), alpha=0.3, color='orange')#, color='#424287')#, c=COLOURS['fit'])

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

    ax.set_xlim(left=0)
    ax.legend(fontsize=20)
    ax.set_xlabel("$\\frac{{1}}{{N(N-M)}}$", fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    
    X_TICKER_SPACING = 0.002
    ax.set_xticks(np.arange(0, np.amax(xs)+X_TICKER_SPACING/2, X_TICKER_SPACING))
    ax.tick_params(axis='both', which='major', labelsize=16)
    fig.tight_layout()


    return fig


def main():
    # TODO: DRY with plot_Wmix_plus_Wcon.py and plot_2D_interpolation.py ?
    T=0.29
    points_with_N16:List[Point] = []
    w_types = ("con", "mix", "dec")
    nmp = ((16,8,0.2), (16,10,0.25), (24,11,0.2), (24,15,0.25), (32,15,0.2), (32,20,0.25), (64,29,0.2), (64,39,0.25)) 
    Ss = (16,24,32)
    for i, (n, m, p) in enumerate(nmp):
        for j, s in enumerate(Ss):
            for k, w_type in enumerate(w_types):
                pt = Point(n, m, p, T, s, w_type)
                points_with_N16.append(pt)

    WL_con_P020 = Points_Same_Type([pt for pt in points_with_N16 if pt.p==0.2 and pt.WL_type=="con"])
    WL_con_P025 = Points_Same_Type([pt for pt in points_with_N16 if pt.p==0.25 and pt.WL_type=="con"])

    show_L=1
    fig = plot_why_remove_N16(WL_con_P020, show_L=show_L)
    fig.savefig(Path(p_common.SAVE_FOLDER,f'compare_N16removed_con_P020_L{show_L}.pdf'))
    fig = plot_why_remove_N16(WL_con_P025, show_L=show_L)
    fig.savefig(Path(p_common.SAVE_FOLDER,f'compare_N16removed_con_P025_L{show_L}.pdf'))

if __name__ == "__main__":
    main()