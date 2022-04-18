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
import fire


# %%
def row(lst, n):
    """Grab successive n-sized chunks from list of lines in file and return a numerical array of elements"""
    rows = []
    for i in np.arange(0, len(lst), n):
        one_line = " ".join(lst[i : i + n]).strip().split()
        rows.append(list(map(float, one_line)))
    return np.asarray(rows)


# %% [markdown]
# ## Import phases from file on disk
#
# Th file has rows with 16 different columns, one for each measured phase.
# The number of rows is the number of lattice configurations used for the measurements.
cols = ["tj", "pc", "dH", "pdec", "pcon", "e", "w", "acc"]

# %%
N = 64
S = 16

# %%
data_folder = f"/Users/enrythebest/Dropbox (University of Michigan)/Projects/PD/EK_v14/N{N}M{N}/S{S}/T029/P02/"
summary_name = f"EK_D3N{N}M{N}S{S}T0290P020_4.txt"
phase_name = f"EK_D3N{N}M{N}S{S}T0290P020_phase_4.txt"
# %% [markdown]
# ### Summary file

# %%
with open(data_folder + summary_name) as fp:
    lines = fp.readlines()
data_from_summary = row(lines[15:], 3)
summary = pd.DataFrame(data=data_from_summary, columns=cols)
print(summary.info())

# %% [markdown]
# We can look at some of the autocorrelation times

# %%
tau = autocorr.integrated_time(summary.e, tol=0)
print(f"Action autocorrelation integrated time: {tau[0]:.2f}")
tau = autocorr.integrated_time(summary.w, tol=0)
print(f"|Wilson| autocorrelation integrated time: {tau[0]:.2f}")

# %%
fig, ax = plt.subplots(figsize=(12, 10))
summary[["e", "w"]].plot(subplots=True, ax=ax)
plt.savefig("energy.png")

# %% [markdown]
# ### Phase file

# %%
with open(data_folder + phase_name) as fp:
    lines = fp.readlines()
print(f"Chunk size: {np.ceil(N/3)}")
data_from_phases = row(lines, int(np.ceil(N / 3)))
phases = pd.DataFrame(data=data_from_phases)
phases.columns = [f"theta{i}" for i in phases.columns]

# %% [markdown]
# * We can concatenate multiple files

# %%
# phases_more = pd.read_csv(data_folder+phase_name.replace('_10','_9'), sep="\s+", header=None)

# %%
# phases_more.columns = [f"theta{i}" for i in phases_more.columns]

# %%
# concatenate and reindex from 0 to the total number of measurements
# phases = pd.concat((phases,phases_more),ignore_index=True)

# %%
print(phases.info())


# %% [markdown]
# Create an histogram for entire dataset of phases. Here we can also choose to apply cuts for thermalization or select measurements with a certain frequency

# %%
freq = 2
thermcut = 0

# %%
alphas = phases.iloc[thermcut::freq].values.flatten()

# %%
print(f"We have a total of {alphas.shape[0]} phases")

# %% [markdown]
# Fold the distribution by taking the absolute value

# %%
alphas_folded = np.fabs(alphas)


# %% [markdown]
# ## Fit the probability distribution based on Fourier expansion

# %% [markdown]
# We use this as our model for the probability distribution of $\alpha$ with unknown coefficients $\tilde{\rho}_k$:
#
# $$\rho(\alpha)=\frac{1}{2 \pi}+\sum_{k=1}^{\infty} \tilde{\rho}_{k} \cos (k \alpha)$$

# %% [markdown]
# Choose the number of parameters and their names based on how many Fourier coefficients to keep:

# %%
Nf = 3  # number of Fourier coefficient
parameters = [f"k{i+1}" for i in np.arange(Nf)]


def prior_transform(cube):
    # the argument, cube, consists of values from 0 to 1
    # we have to convert them to physical scales
    params = cube.copy()
    # all parameters should around 0: change the limits if needed
    lo = -0.1
    hi = 0.1
    params[:] = cube[:] * (hi - lo) + lo
    return params


# %% [markdown]
# Define the model above

# %% [markdown]
# Make sure to only keep positive probabilities

# %%
def prob_model(alpha, params):
    fourier_terms = np.array(
        [params[k] * np.cos((k + 1) * alpha) for k in np.arange(Nf)]
    )
    prob = 1.0 / (2.0 * np.pi) + fourier_terms.sum(axis=0)  # the k=0 term is 1/2pi
    if (prob < 0).any():  # if some probs are negative, the total prob is 0
        return np.zeros_like(prob)
    return prob


# %% [markdown]
# This likelihood ignores the correlations between differen $\alpha$ values. We define it as:
#
# $$ \mathcal{L} \propto \prod_{n=1}^{n_{\text {config }}} \prod_{i=1}^{N} \rho\left(\alpha_{i}^{(n)}\right) $$

# %%
def log_likelihood(params):
    # compute the probability for each alpha point
    probs_alphas = prob_model(alphas, params)
    assert probs_alphas.shape[0] == alphas.shape[0]
    # the total probability is the product of the individual ones, we assume independent
    # for numerical stability, we work in log and avoid zeros
    loglike = np.log(probs_alphas + 1e-100).sum()
    return loglike


# %% [markdown]
# The log-likelihood is often ridiculously small (and negative) because the probabilities are zero hence the sum of logs is -$\infty$

# %% [markdown]
# ### Run the sampler

# %% [markdown]
# Define the reactive nested sampler

# %%
sampler = ultranest.ReactiveNestedSampler(
    parameters,
    log_likelihood,
    prior_transform,
    log_dir="sample_rho",
    resume="overwrite",
)

# %% [markdown]
# If needed you can add a slice sampler for a more efficient sampling in high dimensions

# %%
# import ultranest.stepsampler

# # have to choose the number of steps the slice sampler should take
# # after first results, this should be increased and checked for consistency.

# nsteps = 2 * len(parameters)
# # create step sampler:
# sampler.stepsampler = ultranest.stepsampler.RegionSliceSampler(nsteps=nsteps)

# %%
results = sampler.run(min_num_live_points=400, viz_callback=False)

# %% [markdown]
# ### Plot results

# %%
sampler.print_results()

# %% [markdown]
# The result includes an estimate of the evidence `logZ` (or marginal likelihood) which tells us how probable the data is under this model. This is the strength of nested sampling methods. It allows us to discriminate between models. For example, we can run the code again with the same data but a different number of Fourier coefficients. The difference in `logZ` will tell us which model is favored by the data! (e.g. see [this tutotial](https://johannesbuchner.github.io/UltraNest/example-sine-modelcomparison.html))
#
# Table of results:
#
#  Model  | Evidence |
# | ----- | ----------- |
# | Nf = 5  |  logZ = -70674.178 +- 0.413 |
# | Nf = 4  |  logZ = -70670.095 +- 0.315 |
# | Nf = 3  |  logZ = -70667.004 +- 0.259 |
# | Nf = 2  |  logZ = -70672.906 +- 0.273 |
# | Nf = 1  |  logZ = -70781.606 +- 0.116 |

# %%
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

# %% [markdown]
# ## Extract the value of $\rho$ at $\pi$


# %% [markdown]
# Take all the samples at the value of $\pi$ (remember that we saved 2$\rho$ in the prediction band)

# %%
rhopi = np.asarray([rho[-1] / 2 for rho in band.ys])

# %%
errs = np.diff(np.quantile(rhopi, [0.16, 0.5, 0.84]))
print(f"Value of rho(pi) = {rhopi.mean():.5f} + {errs[1]:.5f} - {errs[0]:.5f}")


# %% [markdown]
# Value of $M$: derived from the equation
# $$
# \rho(\pi) = \frac{1}{2\pi} \left( 1 - \frac{M}{N} \right)
# $$

# %%
M = N * (1.0 - 2 * np.pi * rhopi)
errs = np.diff(np.quantile(M, [0.16, 0.5, 0.84]))
print(f"Value of M = {M.mean():.5f} + {errs[1]:.5f} - {errs[0]:.5f}")
