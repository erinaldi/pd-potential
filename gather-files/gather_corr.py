# this script will create the CSV files for each lattice parameter set
# given N, NT, M, T we collect all the cors along the monte carlo trajectory
import numpy as np
import pandas as pd
from pathlib import Path
from linecache import getline
import re
import argparse


def get_chunk_and_skip_size(pfile: str, N: int = 15):
    """Automatically get the number of lines the data is chunked with and the number of lines for the header of the cor file.

    Args:
        pfile (str): The file with the cor data
        N (int, optional): Number of expected corr distances. Defaults to 15

    Returns:
        chunk, skip: two elements, the chunk size and the number of rows to skip
    """
    # a regexp that counts the blocks of spaces in a line
    pattern = r"\s+"
    # read the 1st line of the file
    line = getline(pfile, 1)
    # count the blocks of spaces
    num_blocks = len(re.findall(pattern, line))
    # if the number of blocks is smaller than (N+1), the file has been created with ifort
    # and the chunk size is 3 and the number of rows to skip is 0
    if num_blocks < N:
        chunk = int(np.ceil(N / 3))
        skips = 0
    else:
        chunk = int(1)
        skips = 0
    return chunk, skips


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
def create_dataframe_cor(pfile: str, N: int = 15) -> pd.DataFrame:
    """Create the dataframe containing the cors for each measured configuration.
    Args:
        pfile (str): The name of the txt file where the cors are saved
        N (int): The number of correlation distances. Default to 15.
    Returns:
        pd.DataFrame: a `pandas` dataframe containing the cors of the run
    """
    try:
        chunk, skips = get_chunk_and_skip_size(pfile)
    except:
        print("Problem getting the chunk and skip sizes")
    # read lines
    with open(pfile) as fp:
        lines = fp.readlines()
    data_array = row(lines[skips:], chunk)
    # there are always N=15 saved values of L (because the minimal NT is 16)
    data = pd.DataFrame(data=data_array.reshape(-1,N))
    data.columns = [f"L{i}" for i in data.columns]
    return data


# %%
# main function to gather the data from a folder or many folders
def gather_data_cor(
    corr_type: str = "con",
    data_folder: str = "../runs/hokusai",
    run_folder: str = "N16/S16/M16/T029/P02",
    out_folder: str = "data-files",
    do_all: bool = False,
):
    """Collect all the data for the observables in different output files for the same set of parameters
    Args:
        corr_type (str,optional): The type of correlator, betwen `con`, `mix`, and `dec`. Defaults to "con".
        data_folder (str, optional): The main data folder where all the different parameters were run. Defaults to "../runs/hokusai".
        run_folder (str, optional): The specific run folder with N/S/M/T/P. Defaults to "N16/S16/M16/T029/P02/".
        put_folder (str, optional): The folder where to save the N/S/M/T/P/cor.csv file. Defaults to "data-files".
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
        # get the cors files: should end with a number before the extension after corr_type
        pfiles = [x for x in run.glob(f"EK_*_{corr_type}_[0-9].txt") if x.is_file()]
        pfiles.sort()
        pfiles1 = [x for x in run.glob(f"EK_*_{corr_type}_[0-9][0-9].txt") if x.is_file()]
        pfiles1.sort()
        pfiles += pfiles1
        if len(pfiles) > 0:
            # create dataframes for each file
            print(f"- We have a total of {len(pfiles)} files in run {run}")
            # N=15 is the number of correlation distances
            frames = [create_dataframe_cor(str(f),15) for f in pfiles]
            # concatenate in single dataframe
            result = pd.concat(frames, ignore_index=True)
            print(f"-- total data size: {result.shape}")
            # make output dir for this run
            pnew = pout.joinpath('/'.join(run.parts[-5:]))
            pnew.mkdir(parents=True,exist_ok=True)
            # save to output file
            outputfile = pnew / f"cor_{corr_type}.csv"
            result.sort_index().to_csv(outputfile, header=True)
            print(f"-- file saved in {outputfile.as_posix()}")
        else:
            print("- no files... skipping.")


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
        gather_data_cor("con", data_folder, run_folder, out_folder, True)
        gather_data_cor("dec", data_folder, run_folder, out_folder, True)
        gather_data_cor("mix", data_folder, run_folder, out_folder, True)
    else:
        gather_data_cor("con", data_folder, run_folder, out_folder)
        gather_data_cor("dec", data_folder, run_folder, out_folder)
        gather_data_cor("mix", data_folder, run_folder, out_folder)

