import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import gvar
from typing import List
from common import Points_Same_Type, Point
import plots_common as p_common

plt.rcParams["figure.figsize"] = (10,10 / 1.618)

def plot_MH_T025():
    figs = []

    input_filename = "MH_plot_data/T025_summary_editedForPython.txt"

    res = np.loadtxt(input_filename)
    N16 = res[:5]
    N24 = res[5:10]
    N32 = res[10:15]
    N48 = res[15:20]
    N64 = res[20:25]
    N96 = res[25:30]
    N128 = res[30:35]

    L2 = res[35:42]
    L3 = res[42:49]
    L4 = res[49:56]
    L1 = res[56:63]
    L5 = res[63:70]

    largeN = res[70:]

    wl_by_N = {16 : N16, 24 : N24, 32: N32, 64 : N64, 96 : N96, 128 : N128, 'large': largeN}
    wl_by_L = {2 : L2, 3 : L3, 4 : L4, 1 : L1, 5 : L5}
    
    # marker = itertools.cycle(('+', 'x', '1', '2', '3')) 

    fig, ax = plt.subplots()
    for n, data in wl_by_N.items():
        ls, ws, err = data.T
        n_label = str(n)
        if n == "large":
            n_label="\\infty"
        ax.scatter(ls, ws, label=f"$N={n_label}$", marker=p_common.MARKERS_BY_N[n], s=70, linewidths=1, c=p_common.COLOURS_BY_N[n])
        ax.errorbar(ls, ws, yerr=err, linestyle='none', capsize=4, linewidth=2, c=p_common.COLOURS_BY_N[n])

    xlim = (0.5, 5.5)
    ylim = (10**-5, 1)

    xs = np.linspace(xlim[0], xlim[1], 5)
    extrap1 = np.exp(-1.0/2.0/0.25*xs)
    extrap2 = np.exp(-(2.09421+0.00698366)*xs)
    ax.plot(xs, extrap1)#, label="exp(-1.0/2.0/0.25*xs)")
    ax.plot(xs, extrap2)#, label="exp(-2.09421*xs+0.00698366)")

    ax.set_xlabel("$L$", fontsize=20)
    ax.set_ylabel("$W(L)/N$", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)
    # ax.tick_params(axis='both', which='minor', labelsize=8)

    ax.grid(linestyle='--')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_yscale('log')
    fig.legend(loc=7, fontsize=20, frameon=False)
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)  
    figs.append(fig)

    fig, ax = plt.subplots()
    for l, data in wl_by_L.items():
        if l == 1 or l == 5:
            continue
        ns, ws, err = data.T
        ax.scatter(1 / ns ** 2, ws, label=f"$L={l}$", marker=p_common.MARKERS_BY_L[l], s=70, linewidths=1, c=p_common.COLOURS_BY_L[l])
        ax.errorbar(1 / ns ** 2, ws, yerr=err, linestyle='none', capsize=4, linewidth=2, c=p_common.COLOURS_BY_L[l])
    
    ax.set_xlabel("$1/N^2$", fontsize=20)
    ax.set_ylabel("$W(L)/N$", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)
    # ax.tick_params(axis='both', which='minor', labelsize=8)

    xlim = np.array((0, 0.002))
    ylim = (0, 0.025)
    yl2 = 0.0152932+4.7833*xlim
    yl3 = 0.00182453+5.97845*xlim
    yl4 = 0.000226379+6.96711*xlim

    ax.plot(xlim, yl2, c=p_common.COLOURS_BY_L[2])
    ax.plot(xlim, yl3, c=p_common.COLOURS_BY_L[3])
    ax.plot(xlim, yl4, c=p_common.COLOURS_BY_L[4])

    ax.grid(linestyle='--')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    # ax.set_yscale('log')
    X_TICKER_SPACING = 0.0005
    ax.set_xticks(np.arange(xlim[0], xlim[1]+X_TICKER_SPACING, X_TICKER_SPACING))
    fig.legend(loc=7, fontsize=20, frameon=False)
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)  
    figs.append(fig)

    return figs

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
    # ax.fill_between(xlim, gvar.mean(yfit)-gvar.sdev(yfit), gvar.mean(yfit)+gvar.sdev(yfit), alpha=0.3, color='#8B008B')#p_common.COLOURS_BY_N['large'])#, c=COLOURS['fit'])
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

    # ** PLOT 1/N/N vs WL and fits **
    fig_Nfit, ax = plt.subplots()
    
    show_Ls = (2,3,4)
    xlim = np.array((0, 0.002))
    ylim = (0, 0.025)
    # grad = {}
    # intcpt = {}
    for l in show_Ls:
        xx, wl, wl_errs = zip(*[(1/pt.n/pt.n, pt.WL_data_av[l-1], pt.WL_errs[l-1]) for pt in points.values()])
        # (m, c), cov = np.polyfit(xx, wl, deg=1, w=1/np.array(wl_errs), cov='unscaled')
        # merr, cerr = np.sqrt(np.diag(cov))
        # grad = gvar.gvar(m, merr)
        # intcpt = gvar.gvar(c, cerr)
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
    # ax.tick_params(axis='both', which='minor', labelsize=8)

    ax.grid(linestyle='--')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    # ax.set_yscale('log')
    X_TICKER_SPACING = 0.0005
    ax.set_xticks(np.arange(xlim[0], xlim[1]+X_TICKER_SPACING, X_TICKER_SPACING))
    fig_Nfit.legend(loc=7, fontsize=20, frameon=False)
    fig_Nfit.tight_layout()
    fig_Nfit.subplots_adjust(right=0.8)  

    fig_Lfit.savefig(Path(p_common.SAVE_FOLDER,'T029_L_fit.pdf'))
    fig_Nfit.savefig(Path(p_common.SAVE_FOLDER,'T029_N_fit.pdf'))


def main():
    plot_T025()
    # plot_MH_T025()

if __name__ == "__main__":
    main()


# TODO:
# 1) change so that only do interpolation once; use results of grad and intercept in same plot
# 2) get rid of surplus code here
# 3) ... why do decays not match well??
# 4) output the numbers