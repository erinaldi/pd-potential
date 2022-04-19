# %% [markdown]
# # Nested sampling for bayesian analysis of phases distributions

# %%
import numpy as np
import pandas as pd
from emcee import autocorr
from corner import corner
import arviz as az
import xarray as xr
import ultranest
from ultranest.plot import PredictionBand
import matplotlib.pyplot as plt
import argparse


# %%
# define command line arguments:
parser = argparse.ArgumentParser()
parser.add_argument("--nf", type=int, default=3, help="Number of Fourier coefficients")
parser.add_argument(
    "--N", type=int, default=32, help="Size of matrices/number of phases"
)
parser.add_argument(
    "--freq", type=int, default=2, help="Frequency of measurements"
)
parser.add_argument(
    "--thermcut", type=int, default=0, help="Thermalization cut"
)
parser.add_argument("--data_dir", type=str, default="./", help="Folder with data files")
parser.add_argument(
    "--s",
    nargs="+",
    required=True,
    help="Summary file. Add more to build a list.",
)
parser.add_argument(
    "--p",
    nargs="+",
    required=True,
    help="Phases file. Add more to build a list",
)
args = parser.parse_args()

cols = ["tj", "pc", "dH", "pdec", "pcon", "e", "w", "acc"]
N = args.N
Nf = args.nf
parameters = [f"k{i+1}" for i in np.arange(Nf)]
data_folder = args.data_dir
summary_files = args.s
phases_files = args.p
freq = args.freq
thermcut = args.thermcut
print(f"Summary files: {summary_files}")
print(f"Phases files: {phases_files}")


# %%
def row(lst: list, n: int) -> np.ndarray:
    """Grab successive n-sized chunks from list of lines in file and return a numerical array of elements"""
    rows = []
    for i in np.arange(0, len(lst), n):
        one_line = " ".join(lst[i : i + n]).strip().split()
        rows.append(list(map(float, one_line)))
    return np.asarray(rows)


def get_summary(filenames: list) -> pd.DataFrame:
    all_data_from_summary = []
    for summary_name in filenames:
        with open(data_folder + summary_name) as fp:
            lines = fp.readlines()
        data_from_summary = row(lines[15:], 3)
        all_data_from_summary.append(data_from_summary)
    summary = pd.DataFrame(data=np.asarray(all_data_from_summary).reshape(-1,8), columns=cols)
    print(summary.info())
    return summary


def get_phases(filenames: list) -> pd.DataFrame:
    all_data_from_phases = []
    for phases_name in filenames:
        with open(data_folder + phases_name) as fp:
            lines = fp.readlines()
        data_from_phases = row(lines, int(np.ceil(N / 3)))
        all_data_from_phases.append(data_from_phases)
    phases = pd.DataFrame(data=np.asarray(all_data_from_phases).reshape(-1,N))
    phases.columns = [f"theta{i}" for i in phases.columns]
    print(phases.info())
    return phases


def autocorrelations(summary: pd.DataFrame) -> None:
    tau = autocorr.integrated_time(summary.e, tol=0)
    print(f"Action autocorrelation integrated time: {tau[0]:.2f}")
    tau = autocorr.integrated_time(summary.w, tol=0)
    print(f"|Wilson| autocorrelation integrated time: {tau[0]:.2f}")

    _, ax = plt.subplots(figsize=(12, 10))
    summary[["e", "w"]].plot(subplots=True, ax=ax)
    plt.savefig("energy.png")


def select_data(phases: pd.DataFrame, freq: int = 2, thermcut: int = 0) -> np.ndarray:
    alphas = phases.iloc[thermcut::freq].values.flatten()
    print(f"We have a total of {alphas.shape[0]} phases")
    return alphas


# %%
def prior_transform(cube):
    # the argument, cube, consists of values from 0 to 1
    # we have to convert them to physical scales
    params = cube.copy()
    # all parameters should around 0: change the limits if needed
    lo = -0.1
    hi = 0.1
    params[:] = cube[:] * (hi - lo) + lo
    return params


def prob_model(alpha, params):
    fourier_terms = np.array(
        [params[k] * np.cos((k + 1) * alpha) for k in np.arange(Nf)]
    )
    prob = 1.0 / (2.0 * np.pi) + fourier_terms.sum(axis=0)  # the k=0 term is 1/2pi
    if (prob < 0).any():  # if some probs are negative, the total prob is 0
        return np.zeros_like(prob)
    return prob


def log_likelihood(params):
    # compute the probability for each alpha point
    probs_alphas = prob_model(alphas, params)
    assert probs_alphas.shape[0] == alphas.shape[0]
    # the total probability is the product of the individual ones, we assume independent
    # for numerical stability, we work in log and avoid zeros
    loglike = np.log(probs_alphas + 1e-100).sum()
    return loglike


# get summary data and plot autocorrelations
summary = get_summary(summary_files)
autocorrelations(summary)
# get phases data
phases = get_phases(phases_files)
alphas = select_data(phases, freq, thermcut)
alphas_folded = np.fabs(alphas)

sampler = ultranest.ReactiveNestedSampler(
    parameters,
    log_likelihood,
    prior_transform,
    log_dir="sample_rho",
    resume="overwrite",
)
results = sampler.run(min_num_live_points=400, viz_callback=False)
sampler.print_results()
sampler.plot()

# %%
df = pd.DataFrame(sampler.results["samples"], columns=parameters)
df["chain"] = 0
df["draw"] = np.arange(len(df), dtype=int)
df = df.set_index(["chain", "draw"])
xdata = xr.Dataset.from_dataframe(df)

dataset = az.InferenceData(posterior=xdata)

# %%
labels = [r"$\tilde{\rho}$" + f"{k+1}" for k in range(Nf)]
fig = corner(
    data=dataset,
    labels=labels,
    show_titles=True,
    title_fmt=".5f",
    quantiles=[0.16, 0.5, 0.84],
    title_kwargs={"fontsize": 14},
    label_kwargs={"fontsize": 16},
)
plt.savefig("corner.png")

# %%
plt.figure(figsize=(20, 10))
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\rho$")
bins = np.linspace(0, np.pi, 100)
plt.hist(alphas_folded, bins, density=True, histtype="step")

alpha_grid = np.linspace(0, np.pi, 500)
band = PredictionBand(alpha_grid)

# go through the solutions
for params in sampler.results["samples"]:
    # compute for each time the probability from the model
    band.add(2 * prob_model(alpha_grid, params))

# add central values
band.line(color="k")
# add 1 sigma quantile
band.shade(color="k", alpha=0.3)
# add wider quantile (0.01 .. 0.99)
band.shade(q=0.49, color="gray", alpha=0.2)
plt.savefig("rho_alpha.png")

rhopi = np.asarray([rho[-1] / 2 for rho in band.ys])
errs = np.diff(np.quantile(rhopi, [0.16, 0.5, 0.84]))
print(f"Value of rho(pi) = {rhopi.mean():.5f} + {errs[1]:.5f} - {errs[0]:.5f}")
M = N * (1.0 - 2 * np.pi * rhopi)
errs = np.diff(np.quantile(M, [0.16, 0.5, 0.84]))
print(f"Value of M = {M.mean():.5f} + {errs[1]:.5f} - {errs[0]:.5f}")
