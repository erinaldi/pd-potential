import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np


def mvalue(data_folder):
    try:
        df_2 = pd.read_csv(
            os.path.join(data_folder, "NF2/results.csv"), header=0, index_col=0
        )
        df_3 = pd.read_csv(
            os.path.join(data_folder, "NF3/results.csv"), header=0, index_col=0
        )
        df_4 = pd.read_csv(
            os.path.join(data_folder, "NF4/results.csv"), header=0, index_col=0
        )
        df_5 = pd.read_csv(
            os.path.join(data_folder, "NF5/results.csv"), header=0, index_col=0
        )
    except FileNotFoundError:
        return None
    return pd.concat([df_2, df_3, df_4, df_5], axis=0)


parser = argparse.ArgumentParser()
parser.add_argument(
    "--N", type=int, default=16, help="Size of matrices/number of phases (N)"
)
parser.add_argument("--T", type=float, default=0.29, help="Temperature (T)")
parser.add_argument("--P", type=float, default=0.2, help="Polyakov Loop (P)")
parser.add_argument(
    "--data_dir", type=str, default="./m-values", help="Folder with result files"
)
args = parser.parse_args()

N = args.N
T = args.T
P = args.P
data_folder = args.data_dir

all_nt = os.listdir(f"{data_folder}/N{N}/")

df_dict = {}
for nts in np.sort(all_nt):
    nt = nts.split("S")[-1]
    assert os.path.isdir(os.path.join(data_folder, f"N{N}", nts))
    m_folder = f"{data_folder}/N{N}/S{nt}/M{N}/T{str(T).replace('.','')}/P{str(P).replace('.','')}"
    df_dict[nt] = mvalue(m_folder)


for nt,df in df_dict.items():
    plt.errorbar(df.NF,df.M,[df.minusM,df.plusM],marker='s',label=f"$n_t=${nt}")
plt.xticks([2,3,4,5])
plt.xlabel(r"$\Lambda$")
plt.ylabel(r"$M$",rotation=0)
plt.legend(loc="upper left")
plt.title(f"$N=${N} $T=${T} $P=${P}")
plt.tight_layout()
plt.savefig(os.path.join(data_folder, f"N{N}","MvsL.png"))
