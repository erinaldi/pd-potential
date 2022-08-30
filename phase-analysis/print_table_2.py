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
    for i,n in enumerate([16, 24, 32]):
        for j,t in enumerate([0.29, 0.30, 0.31]):
            df = mvalue(n, 24, t, p)
            try:
                ms = df.query("NF==3")[["M", "plusM", "minusM"]].values[0]
                if i == 0:
                    print(
                        f"{p} & {n}", end=' & ')
                else:
                    print(
                        f"    & {n}", end=' & ')
                if j<2:
                    print(
                        f"${ms[0]:.2f}^{{+{ms[1]:.2f}}}_{{-{ms[2]:.2f}}}$", end=' & '
                    )
                else:
                    print(
                        f"${ms[0]:.2f}^{{+{ms[1]:.2f}}}_{{-{ms[2]:.2f}}}$\\\\"
                    )
                    print(
                        f"    & {n} & ${ms[0]:.2f}^{{+{ms[1]:.2f}}}_{{-{ms[2]:.2f}}}$\\\\"
                    )
            except AttributeError:
                pass
