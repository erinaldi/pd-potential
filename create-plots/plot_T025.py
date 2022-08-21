import numpy as np
import matplotlib.pyplot as plt
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

def main():
    plot_MH_T025()

if __name__ == "__main__":
    main()