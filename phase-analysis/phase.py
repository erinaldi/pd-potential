# %% [markdown]
# # Nested sampling for bayesian analysis of phases distributions

# %%
import os
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
parser.add_argument("--NF", type=int, default=3, help="Number of Fourier coefficients")
parser.add_argument(
    "--N", type=int, default=32, help="Size of matrices/number of phases (N)"
)
parser.add_argument("--NT", type=int, default=16, help="Size of lattice (Nt)")
parser.add_argument("--M", type=int, default=16, help="Size of deconfined sector (M)")
parser.add_argument("--T", type=float, default=0.29, help="Temperature (T)")
parser.add_argument("--P", type=float, default=0.2, help="Polyakov Loop (P)")
parser.add_argument("--freq", type=int, default=2, help="Frequency of measurements")
parser.add_argument("--thermcut", type=int, default=100, help="Thermalization cut")
parser.add_argument(
    "--data_dir", type=str, default="./data-files", help="Folder with data files"
)
parser.add_argument(
    "--out_dir", type=str, default="./phase-analysis/m-values", help="Folder with output files"
)
args = parser.parse_args()

NF = args.NF
N = args.N
M = args.M
NT = args.NT
T = args.T
P = args.P
data_folder = args.data_dir
out_folder = args.out_dir
file_folder = f"N{N}/S{NT}/M{M}/T{str(T).replace('.','')}/P{str(P).replace('.','')}"
parameters = [f"k{i+1}" for i in np.arange(NF)]
phases_file = os.path.join(data_folder, file_folder, "phase.csv")
assert os.path.isfile(phases_file)
out_dir = os.path.join(out_folder,file_folder,f"NF{NF}")
os.makedirs(out_dir, exist_ok = True)
freq = args.freq
thermcut = args.thermcut

print(f"Phases files: {phases_file}")
print(f"Output written to: {out_dir}")


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
        [params[k] * np.cos((k + 1) * alpha) for k in np.arange(NF)]
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


# get phases data
phases = pd.read_csv(phases_file, header=0, index_col=0)
alphas = select_data(phases, freq, thermcut)
alphas_folded = np.fabs(alphas)

sampler = ultranest.ReactiveNestedSampler(
    parameters,
    log_likelihood,
    prior_transform,
    log_dir=os.path.join(out_dir,"sample_rho"),
    resume="overwrite",
)
sampler.run(min_num_live_points=400, viz_callback=False)
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
labels = [r"$\tilde{\rho}$" + f"{k+1}" for k in range(NF)]
fig = corner(
    data=dataset,
    labels=labels,
    show_titles=True,
    title_fmt=".5f",
    quantiles=[0.16, 0.5, 0.84],
    title_kwargs={"fontsize": 14},
    label_kwargs={"fontsize": 16},
)
plt.savefig(os.path.join(out_dir,"corner.png"))

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
plt.savefig(os.path.join(out_dir,"rho_alpha.png"))

rhopi = np.asarray([rho[-1] / 2 for rho in band.ys])
errs = np.diff(np.quantile(rhopi, [0.16, 0.5, 0.84]))
print(f"Value of rho(pi) = {rhopi.mean():.5f} + {errs[1]:.5f} - {errs[0]:.5f}")
M = N * (1.0 - 2 * np.pi * rhopi)
mcmc = np.percentile(M, [16, 50, 84])
errs = np.diff(mcmc)
print(f"Value of M = {mcmc[1]:.5f} + {errs[1]:.5f} - {errs[0]:.5f}")
M_df = pd.DataFrame(
    data={
        "M": mcmc[1],
        "plusM": errs[1],
        "minusM": errs[0],
        "NF": NF,
        "logZ": sampler.results["logz"],
        "logZerr": sampler.results["logzerr"],
    },
    index=[0]
)
M_df.to_csv(os.path.join(out_dir,"results.csv"))
