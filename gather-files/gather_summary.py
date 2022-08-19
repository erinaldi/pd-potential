# this script will create the CSV files for each lattice parameter set
# given N, NT, M, T we collect all the observables along the monte carlo trajectory
# including the parameters of the monte carlo integration
# We have three main output files: summary, phases, correlators
import numpy as np
import pandas as pd
from pathlib import Path
from linecache import getline
import re
import argparse


def get_chunk_and_skip_size(pfile: str):
    """Automatically get the number of lines the data is chunked with and the number of lines for the header of the summary file

    Args:
        pfile (str): The file with the summary data

    Returns:
        chunk, skip: two elements, the chunk size and the number of rows to skip
    """
    # a regexp that counts the blocks of spaces in a line
    pattern = r"\s+"
    # read the 16th line of the file (one with data, after the header)
    line = getline(pfile, 16)
    # count the blocks of spaces
    num_blocks = len(re.findall(pattern, line))
    # if the number of blocks is smaller than 9, the file has been created with ifort
    # and the chunk size is 3 and the number of rows to skip is 15
    if num_blocks < 9:
        chunk = 3
        skips = 15
    else:
        chunk = 1
        skips = 14
    return chunk, skips


# %%
# extract the MCMC parameters for a run (they are contained in the header)
# they will be used as new columns in the dataframe so we can plot them for each run
def get_params_summary(pfile: str) -> dict:
    """Extract the parameters of a MC chain by looking at the header of the output file.txt
    Args:
        pfile (str): a MC output file passed as a string
    Returns:
        dict: the parameters of the MC run returned as a dictionary
    """
    pattern = r"[:=]\s+(\d.*)$"
    line_nmat = getline(pfile, 1)
    nmat = re.findall(pattern, line_nmat)[0]
    line_ndec = getline(pfile, 2)
    ndec = re.findall(pattern, line_ndec)[0]
    line_ntime = getline(pfile, 3)
    ntime = re.findall(pattern, line_ntime)[0]
    line_temp = getline(pfile, 5)
    temp = re.findall(pattern, line_temp)[0]
    line_ntau = getline(pfile, 7)
    ntau = re.findall(pattern, line_ntau)[0]
    line_udtau = getline(pfile, 8)
    udtau = re.findall(pattern, line_udtau)[0]
    line_adtau = getline(pfile, 9)
    adtau = re.findall(pattern, line_adtau)[0]
    line_pdec = getline(pfile, 10)
    pdec = re.findall(pattern, line_pdec)[0]
    line_pcon = getline(pfile, 11)
    pcon = re.findall(pattern, line_pcon)[0]
    line_pcoeff = getline(pfile, 12)
    pcoeff = re.findall(pattern, line_pcoeff)[0]
    return {
        "nmat": int(nmat),
        "ndec": int(ndec),
        "ntime": int(ntime),
        "ntau": int(ntau),
        "udtau": float(udtau),
        "adtau": float(adtau),
        "temperature": float(temp),
        "pdec": float(pdec.split("+/-")[0]),
        "pdec_var": float(pdec.split("+/-")[1]),
        "pcon": float(pcon),
        "pcoeff": float(pcoeff),
    }


def row(lst: list, n: int) -> np.ndarray:
    """Grab successive n-sized chunks from list of lines 
    in file and return a numerical array of elements
    Args:
        lst (list): A list of lines from a file
        n (int): the size of each chunk of lines to consider
    Returns:
        np.ndarray: an array with all the entries in all the lines in the chunk flattened
    """
    rows = []
    for i in np.arange(0, len(lst), n):
        one_line = " ".join(lst[i : i + n]).strip().split()
        rows.append(list(map(float, one_line)))
    return np.asarray(rows)


# %%
# read the observables from a txt file and create a dataframe
# also add the MCMC params as new columns, using the trajectory number as index
def create_dataframe_summary(pfile: str) -> pd.DataFrame:
    """Create the dataframe containing the observables for each trajectory.
    Then add the MCMC parameters from the header as additional columns.
    This is advantageous because these parameters can be used to group specific trajectories together.
    Args:
        pfile (str): The name of the txt file where the observables are saved
    Returns:
        pd.DataFrame: a `pandas` dataframe containing the observables and the MCMC parameters of the run
    """
    try:
        mc_params = get_params_summary(pfile)
        chunk, skips = get_chunk_and_skip_size(pfile)
    except:
        print("Problem getting the mcmc parameters")
    # column names
    cols = [
        "sweep",
        "permcount",
        "dH",
        "pdeconfined",
        "pconfined",
        "action",
        "wilson",
        "acceptance",
    ]
    # read lines
    with open(pfile) as fp:
        lines = fp.readlines()
    data_array = row(lines[skips:], chunk)
    data = pd.DataFrame(data=data_array, columns=cols)
    # add mcmc params as columns
    for k, v in mc_params.items():
        data[k] = v
    # add traj_length and mdtu and save frequency columns
    data["utau"] = data.udtau * data.ntau
    data["atau"] = data.adtau * data.ntau
    data["mdtu"] = data.utau * data.sweep
    data["freq"] = data.mdtu.diff()
    return data


# %%
# main function to gather the data from a folder or many folders
def gather_data_summary(
    data_folder: str = "../runs/hokusai",
    run_folder: str = "N16/S16/M16/T029/P02",
    out_folder: str = "data-files",
    do_all: bool = False,
):
    """Collect all the data for the observables in different output files for the same set of parameters
    Args:
        data_folder (str, optional): The main data folder where all the different parameters were run. Defaults to "../runs/hokusai".
        run_folder (str, optional): The specific run folder with N/S/M/T/P. Defaults to "N16/S16/M16/T029/P02/".
        put_folder (str, optional): The folder where to save the N/S/M/T/P/summary.csv file. Defaults to "data-files".
        do_all (bool, optional): If we should collect the data from all runs or just the one specified. Defaults to False
    """
    pdata = Path(data_folder)
    pout = Path(out_folder)
    assert pdata.is_dir()
    assert pout.is_dir()
    if do_all:
        all_runs = [x for x in pdata.glob("N*/S*/M*/T*/P*")]
    else:
        assert (pdata / run_folder).is_dir()
        all_runs = [pdata / run_folder]
    # loop over runs: they are Path objects
    print(f"We have a total of {len(all_runs)} runs to gather...")
    for run in all_runs:
        # get the output files: should end with a number before the extension
        pfiles = [x for x in run.glob("EK_*[0-9]_[0-9].txt") if x.is_file()]
        pfiles.sort()
        pfiles1 = [x for x in run.glob("EK_*[0-9]_[0-9][0-9].txt") if x.is_file()]
        pfiles1.sort()
        pfiles += pfiles1
        if len(pfiles) > 0:
            # create dataframes for each file
            print(f"- We have a total of {len(pfiles)} files in run {run}")
            frames = [create_dataframe_summary(str(f)) for f in pfiles]
            # concatenate in single dataframe
            result = pd.concat(frames, ignore_index=True)
            print(f"-- total data size: {result.shape}")
            # make output dir for this run
            pnew = pout.joinpath('/'.join(run.parts[-5:]))
            pnew.mkdir(parents=True,exist_ok=True)
            # save to output file
            outputfile = pnew / "summary.csv"
            result.sort_index().to_csv(outputfile, header=True)
            print(f"-- file saved in {outputfile.as_posix()}")


if __name__ == "__main__":
    # define command line arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=16, help="Size of matrices (N)")
    parser.add_argument("--NT", type=int, default=16, help="Size of lattice (Nt)")
    parser.add_argument(
        "--M", type=int, default=16, help="Size of deconfined sector (M)"
    )
    parser.add_argument("--T", type=float, default=0.29, help="Temperature (T)")
    parser.add_argument("--P", type=float, default=0.2, help="Polyakov Loop (P)")
    parser.add_argument(
        "--data",
        type=str,
        default="/home/enrico/Projects/LP/runs/hokusai",
        help="Folder with data files",
    )
    parser.add_argument(
        "--out",
        type=str,
        default="/home/enrico/Projects/LP/pd-potential/data-files",
        help="Folder where to store output files",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="If we need to run this gather scritp an ALL sub-folders in `data`",
    )
    args = parser.parse_args()

    # build folder path
    N = args.N
    M = args.M
    NT = args.NT
    T = args.T
    P = args.P
    data_folder = args.data
    out_folder = args.out
    run_folder = f"N{N}/S{NT}/M{M}/T{str(T).replace('.','')}/P{str(P).replace('.','')}"
    if args.all:
        gather_data_summary(data_folder, run_folder, out_folder, True)
    else:
        gather_data_summary(data_folder, run_folder, out_folder)

