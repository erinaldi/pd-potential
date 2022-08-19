import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import gvar
from typing import List
from common import Points_Same_Type, Point
import plots_common as p_common

def plot_subtraction_method():
    t=0.29
    p=0.25
    m_n_ratio = 20/32#+0.001
    n_t = 32
    wl_type = "dec"
    # # folder_path = "first_kind_data"
    p16 = Point(16, 16, p, t, n_t, "dec")#, folder_path) #TODO: update
    p24 = Point(24, 24, p, t, n_t, "dec")#, folder_path)
    p32 = Point(32, 32, p, t, n_t, "dec")#, folder_path)

    p16_dec = Point(16, 10, p, t, n_t, "dec")
    p24_dec = Point(24, 15, p, t, n_t, "dec")
    p32_dec = Point(32, 20, p, t, n_t, "dec")

    p16_mix = Point(16, 10, p, t, n_t, "mix")
    p24_mix = Point(24, 15, p, t, n_t, "mix")
    p32_mix = Point(32, 20, p, t, n_t, "mix")

    # p16.n_t = n_t # pretend for now
    point_set_firstKind = Points_Same_Type([p16, p24, p32], first_kind=True)
    point_set_dec = Points_Same_Type([p16_dec, p24_dec, p32_dec], first_kind=False)
    point_set_mix = Points_Same_Type([p16_mix, p24_mix, p32_mix], first_kind=False)

    Ls = point_set_firstKind.Ls

    fig, ax = plt.subplots()
    # ax.scatter(Ls, point_set_firstKind.large_N_WL[n_t], marker='x', label=f"1st-kind large $N$, $n_t={n_t}$", s=80)
    ax.scatter(Ls, point_set_firstKind.large_N_WL[n_t], marker='x', label=f"$W_{{\\rm full, 1st}}$", s=80, c='orange')
    ax.errorbar(Ls, point_set_firstKind.large_N_WL[n_t], yerr=point_set_firstKind.large_N_WL_err[n_t], linestyle='none', capsize=4, c='orange')

    # ax.scatter(Ls, point_set_dec.large_N_WL[n_t], marker='x', label=f"deconfined large $N$, $n_t={n_t}$", s=80)
    ax.scatter(Ls, point_set_dec.large_N_WL[n_t], marker='x', label=f"$W_{{\\rm dec, 2nd}}$", s=80, c='#0000a7')
    ax.errorbar(Ls, point_set_dec.large_N_WL[n_t], yerr=point_set_dec.large_N_WL_err[n_t], linestyle='none', capsize=4, c='#0000a7')

    # sustracted_ys = 1 / (1 - m_n_ratio ** 2) * (point_set_firstKind.large_N_WL[n_t] - m_n_ratio ** 2 * point_set_dec.large_N_WL[n_t])
    # sustracted_ys_err = np.sqrt(1 / (1 - m_n_ratio ** 2) * (point_set_firstKind.large_N_WL_err[n_t] ** 2 + (m_n_ratio ** 2 * point_set_dec.large_N_WL_err[n_t]) ** 2))

    g_first = gvar.gvar(point_set_firstKind.large_N_WL[n_t], point_set_firstKind.large_N_WL_err[n_t])
    g_dec = gvar.gvar(point_set_dec.large_N_WL[n_t], point_set_dec.large_N_WL_err[n_t])
    g_mix = gvar.gvar(point_set_mix.large_N_WL[n_t], point_set_mix.large_N_WL_err[n_t])

    # sustracted_ys = 1 / (1 - m_n_ratio ** 2) * (point_set_firstKind.large_N_WL[n_t] - m_n_ratio ** 2 * point_set_dec.large_N_WL[n_t] - )
    
    sustracted_ys =  ( g_first - m_n_ratio ** 2 * g_dec - (m_n_ratio) * (1 - m_n_ratio) * g_mix) / (1 - m_n_ratio)

    # sustracted_ys_err = np.sqrt(1 / (1 - m_n_ratio ** 2) * (point_set_firstKind.large_N_WL_err[n_t] ** 2 + (m_n_ratio ** 2 * point_set_dec.large_N_WL_err[n_t]) ** 2))

    # ax.scatter(Ls, gvar.mean(sustracted_ys), marker='x', label=f"subtracted", s=80, c='#008176')
    ax.scatter(Ls, gvar.mean(sustracted_ys), marker='x', label=f"$W_{{\\rm subtracted}}$", s=80, c='#008176')
    ax.errorbar(Ls, gvar.mean(sustracted_ys), yerr=gvar.sdev(sustracted_ys), linestyle='none', capsize=4, c='#008176')

    beta = 1 / t
    xlim = np.array((0.5,5.5))
    ypred = np.exp(- beta * xlim / 2)
    ax.plot(xlim, ypred, c=p_common.COLOURS['pred'])
    
    # ax.set_xlim(xlim)
    # ax.set_ylim(ylim)
    # ax.legend(fontsize=16)

    ax.set_xlabel("$L$", fontsize=20)
    # ax.set_ylabel("$\\frac{{1}}{{N-M}} \\left( W_{\\rm 1st} - W_{\\rm dec} - W_{\\rm mix} \\right)$", fontsize=20)
    # ax.set_ylabel("$\\frac{{1}}{{N-M}} W_{{\\rm subtracted}}$", fontsize=20)
    ax.set_ylabel("Normalised $W(L)$", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)
    # ax.tick_params(axis='both', which='minor', labelsize=8)
    # ax.set_xticks(np.arange(1,4+1,1))
    ax.grid(linestyle='--')
    ax.set_xlim(xlim)
    # ax.set_ylim(ylim)
    ax.set_yscale('log')
    # fig.legend(loc=7, fontsize=20, frameon=False)

    fig.legend(loc=7, fontsize=20, frameon=False)
    fig.tight_layout()
    # fig.legend(fontsize=20, frameon=True,loc='lower left')

    fig.subplots_adjust(right=0.75) 

    return fig

def main():
    fig = plot_subtraction_method()
    fig.savefig(Path(p_common.SAVE_FOLDER,"subtraction_method_nt32.pdf"))

if __name__ == "__main__":
    main()