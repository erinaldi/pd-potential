import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import gvar
from typing import List
from common import Points_Same_Type, Point
import plots_common as p_common

plt.rcParams["figure.figsize"] = (10,10 / 1.618)

# move into common module?
def take_large_N_unconstrained(ns, WLs, WL_errs):
    ns = np.array(ns)
    WLs = np.array(WLs)
    WL_errs = np.array(WL_errs)
    (m, c), cov = np.polyfit(1/ns/ns, WLs, deg=1, w=1/WL_errs, cov='unscaled')
    merr, cerr = np.sqrt(np.diag(cov))
    largeN_wl = gvar.gvar(c, cerr)
    grad = gvar.gvar(m, merr)
    return largeN_wl, grad

def plot_T025():
    # ** COLLECT POINTS WE WILL USE **
    UNCONSTRAINED_P_FILLER = "0.5"
    CONSTANT_NT = 24
    T = 0.25
    WL_TYPE = "dec"
    Ns = (16, 24, 32, 48, 64, 96, 128)
    # Ns = (48, 64, 96, 128)
    points = {}
    for n in Ns:
        points[n] = Point(n, n, UNCONSTRAINED_P_FILLER, T, CONSTANT_NT, WL_TYPE)
    
    # ** PLOT L vs WL **
    fig_Lfit, ax = plt.subplots()
    xlim = np.array((0.5, 5.5))
    ylim =  np.array((10**-5, 1))

    # Large N extrapolation
    ## " For the extrapolation to N = ∞, we used data points at N = 24, 32, 48, 64, 96, 128 for L = 1, 2, 3, 4 and N = 48, 64, 96, 128 for L = 5. "
    largeN_wl = []
    largeN_extrap_grad = []
    largeN_by_L = {}
    largeN_extrap_grad_by_L = {}
    Ls = (1,2,3,4,5)
    low_Ns = (24, 32, 48, 64, 96, 128)
    high_Ns = (48, 64, 96, 128)
    for l in Ls:
        if l<5:
            use_Ns = low_Ns
        else:
            use_Ns = high_Ns
        WLs = [points[n].WL_data_av[l-1] for n in use_Ns]
        errs = [points[n].WL_errs[l-1] for n in use_Ns]
        lNwl, grad = take_large_N_unconstrained(use_Ns, WLs, errs)
        largeN_wl.append(lNwl)
        largeN_extrap_grad.append(grad)
        largeN_by_L[l] = lNwl
        largeN_extrap_grad_by_L[l] = grad
        
    for n in Ns:
        ws, err = points[n].WL_data_av, points[n].WL_errs
        ls = np.arange(1, len(ws)+1)
        n_label = str(n)
        ax.scatter(ls, ws, label=f"$N={n_label}$", marker=p_common.MARKERS_BY_N[n], s=70, linewidths=1, c=p_common.COLOURS_BY_N[n])
        ax.errorbar(ls, ws, yerr=err, linestyle='none', capsize=4, linewidth=2, c=p_common.COLOURS_BY_N[n])
    
    ax.scatter(Ls, gvar.mean(largeN_wl), label=f"$N=\\infty$", marker=p_common.MARKERS_BY_N['large'], s=70, linewidths=1, c=p_common.COLOURS_BY_N['large'])
    ax.errorbar(Ls, gvar.mean(largeN_wl), yerr=gvar.sdev(largeN_wl), linestyle='none', capsize=4, linewidth=2, c=p_common.COLOURS_BY_N['large'])
    yy = np.log(largeN_wl)
    (m, c), cov = np.polyfit(Ls, gvar.mean(yy), w=1/gvar.sdev(yy), deg=1, cov='unscaled')
    merr, cerr = np.sqrt(np.diag(cov))
    grad = gvar.gvar(m, merr)
    intcpt = gvar.gvar(c, cerr)
    yfit = np.exp(grad * np.array(xlim) + intcpt)
    ax.fill_between(xlim, gvar.mean(yfit)-gvar.sdev(yfit), gvar.mean(yfit)+gvar.sdev(yfit), alpha=0.3, color=p_common.COLOURS_BY_N['large'])#, c=COLOURS['fit'])
    ax.plot(xlim, gvar.mean(yfit), c=p_common.COLOURS_BY_N['large'])

    predicted_y = np.exp(-1.0/2.0/0.25*xlim)
    ax.plot(xlim, predicted_y, color='red')

    ax.set_xlabel("$L$", fontsize=20)
    ax.set_ylabel("$W(L)/N$", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.grid(linestyle='--')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_yscale('log')
    fig_Lfit.legend(loc=7, fontsize=20, frameon=False)
    fig_Lfit.tight_layout()
    fig_Lfit.subplots_adjust(right=0.8)  

    # write fit data to file
    with open(Path(p_common.SAVE_NUMBERS_FOLDER, f"T025_fit_largeNonly.txt"),"w+") as f:
        f.write(f"c={-grad}\nd={intcpt}")

    # ** PLOT 1/N/N vs WL and fits **
    fig_Nfit, ax = plt.subplots()
    
    show_Ls = (2,3,4)
    xlim = np.array((0, 0.002))
    ylim = (0, 0.025)
    for l in show_Ls:
        xx, wl, wl_errs = zip(*[(1/pt.n/pt.n, pt.WL_data_av[l-1], pt.WL_errs[l-1]) for pt in points.values()])
        grad = largeN_extrap_grad_by_L[l]
        intcpt = largeN_by_L[l]
        
        ax.scatter(xx, wl, label=f"$L={l}$", marker=p_common.MARKERS_BY_L[l], s=70, linewidths=1, c=p_common.COLOURS_BY_L[l])
        ax.errorbar(xx, wl, yerr=wl_errs, linestyle='none', capsize=4, linewidth=2, c=p_common.COLOURS_BY_L[l])
    
        yfit = xlim*grad + intcpt
        ax.plot(xlim, gvar.mean(yfit), c=p_common.COLOURS_BY_L[l])
        ax.fill_between(xlim, gvar.mean(yfit)-gvar.sdev(yfit), gvar.mean(yfit)+gvar.sdev(yfit), alpha=0.3, color=p_common.COLOURS_BY_L[l])#, c=COLOURS['fit'])

    ax.set_xlabel("$1/N^2$", fontsize=20)
    ax.set_ylabel("$W(L)/N$", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)

    ax.grid(linestyle='--')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    X_TICKER_SPACING = 0.0005
    ax.set_xticks(np.arange(xlim[0], xlim[1]+X_TICKER_SPACING, X_TICKER_SPACING))
    fig_Nfit.legend(loc=7, fontsize=20, frameon=False)
    fig_Nfit.tight_layout()
    fig_Nfit.subplots_adjust(right=0.8)  

    fig_Lfit.savefig(Path(p_common.SAVE_FOLDER,'T025_L_fit.pdf'))
    fig_Nfit.savefig(Path(p_common.SAVE_FOLDER,'T025_N_fit.pdf'))


def main():
    plot_T025()
    # plot_MH_T025()

if __name__ == "__main__":
    main()


# TODO:
# 2) get rid of surplus code here
# 3) ... why do decays not match well??
#   we do not take continuum limitt
# 4) output the numbers