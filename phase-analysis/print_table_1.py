import os
import pandas as pd


def mvalue(N,NT,T,P):
    T = f"{T:.2f}"  # deal with T=0.30
    P = f"{P}"
    data_folder = f"./m-values/N{N}/S{NT}/M{N}/T{T.replace('.','')}/P{P.replace('.','')}"
    try:
        df_2 = pd.read_csv(os.path.join(data_folder,"NF2/results.csv"), header=0,index_col=0)
        df_3 = pd.read_csv(os.path.join(data_folder,"NF3/results.csv"), header=0,index_col=0)
        df_4 = pd.read_csv(os.path.join(data_folder,"NF4/results.csv"), header=0,index_col=0)
        df_5 = pd.read_csv(os.path.join(data_folder,"NF5/results.csv"), header=0,index_col=0)
    except FileNotFoundError:
        return None
    return pd.concat([df_2,df_3,df_4,df_5], axis=0)


for p in [0.2, 0.25]:
    j = 0
    for n in [16, 24, 32, 64]:
        for i, nt in enumerate([16, 24, 32]):
            df = mvalue(n, nt, 0.29, p)
            try:
                ms = df.query("NF==3")[["M", "plusM", "minusM"]].values[0]
                if i == 0 and j == 0:
                    print(
                        f"{p} & {n} & {nt} & ${ms[0]:.2f}^{{+{ms[1]:.2f}}}_{{-{ms[2]:.2f}}}$\\\\"
                    )
                elif i == 0 and j > 0:
                    print(
                        f"    & {n} & {nt} & ${ms[0]:.2f}^{{+{ms[1]:.2f}}}_{{-{ms[2]:.2f}}}$\\\\"
                    )
                else:
                    print(
                        f"    &     & {nt} & ${ms[0]:.2f}^{{+{ms[1]:.2f}}}_{{-{ms[2]:.2f}}}$\\\\"
                    )
            except AttributeError:
                pass
        j += 1
    for n in [96]:
        for i, nt in enumerate([16, 24]):
            df = mvalue(n, nt, 0.29, p)
            try:
                ms = df.query("NF==3")[["M", "plusM", "minusM"]].values[0]
                if i == 0 and j == 0:
                    print(
                        f"{p} & {n} & {nt} & ${ms[0]:.2f}^{{+{ms[1]:.2f}}}_{{-{ms[2]:.2f}}}$\\\\"
                    )
                elif i == 0 and j > 0:
                    print(
                        f"    & {n} & {nt} & ${ms[0]:.2f}^{{+{ms[1]:.2f}}}_{{-{ms[2]:.2f}}}$\\\\"
                    )
                else:
                    print(
                        f"    &     & {nt} & ${ms[0]:.2f}^{{+{ms[1]:.2f}}}_{{-{ms[2]:.2f}}}$\\\\"
                    )
            except AttributeError:
                pass
        j += 1
    print("\\hline")