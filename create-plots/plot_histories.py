import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from common import Point

SAVE_DIR_ROOT = Path("products", "plots_debug")
FIG_SIZE = 8

all_con_files = sorted(Path('data-files').rglob('*cor_con.csv'))

for con_path in all_con_files:
    _, nstr, sstr, mstr, tstr, pstr, _ = con_path.parts

    n = int(nstr[1:]) 
    n_t = int(sstr[1:])
    m = int(mstr[1:])
    t = float("0."+tstr[2:]) # assumes 0<T<1
    p = float("0."+pstr[2:]) # assumes 0<P<1

    if m == n:
        continue

    summary_filename = con_path.parent / "summary.csv"
    summary_data = pd.read_csv(summary_filename,header=0,index_col=0)
    
    pt = Point(n, m, p, t, n_t, "con")
    tcut = pt.therm_cut

    plot_jobs = []
    plot_jobs.append(("action", summary_data.action, 200))
    # plot_jobs.append(("$dH$", summary_data.dH))
    plot_jobs.append(("perm count", summary_data.permcount, 0))
    plot_jobs.append(("$W_{{\\rm con}}$ $L=1$", pt.WL_raw_data_unthermalised_head[:,0], 200))
    plot_jobs.append(("$W_{{\\rm con}}$ $L=2$", pt.WL_raw_data_unthermalised_head[:,1], 200))

    fig, axs = plt.subplots(2,2,figsize=(FIG_SIZE*1.618, FIG_SIZE))
    for ax, (y_label, y_data, run_av_num) in zip(axs.flatten(), plot_jobs):
        x_right = len(y_data)-1
        ax.set_xlim(left=0, right=x_right)
        ax.plot(y_data)
        y_lower, y_upper = ax.get_ylim()
        ax.set_ylim(auto=False)
        ax.plot((tcut, tcut), (y_lower, y_upper), c='r')
        ax.set_ylabel(y_label, fontsize=16)
        if run_av_num > 0:
            rav = np.convolve(y_data, np.ones(run_av_num)/run_av_num, mode='valid')
            rav_x = np.arange(run_av_num // 2, len(y_data)-run_av_num // 2+1)
            ax.plot(rav_x, rav)

    save_dir = Path(SAVE_DIR_ROOT, con_path.parent)
    save_dir.mkdir(parents=True, exist_ok=True)
    fig.suptitle(f"$L=1$ $\\tau={pt.autocor_taus[0]:.3f}$\n$L=2$ $\\tau={pt.autocor_taus[1]:.3f}$\nthermal cut $={tcut}$")
    fig.tight_layout()
    fig.savefig(Path(save_dir, "history.png"))
    plt.close(fig)
