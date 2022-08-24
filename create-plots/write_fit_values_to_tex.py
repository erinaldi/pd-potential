# Latex macros to filename
from pathlib import Path

NUMBERS_FOLDER = Path("products", "numbers")
SAVE_FILENAME = Path(NUMBERS_FOLDER, "fit_numbers_tex.txt")

files_to_command = {
    "P0_comp-conf_fit_2Dinterp" : "CompletelyConfinedPZero",
    "P02_con_fit_2Dinterp" : "ConPpointTwo",
    "P02_mix_fit_2Dinterp" : "MixPpointTwo",
    "P02_WconPlusMix_2D" : "ConPlusMixPpointTwo",
    "P025_con_fit_2Dinterp" : "ConPpointTwoFive",
    "P025_mix_fit_2Dinterp" : "MixPpointTwoFive",
    "P025_WconPlusMix_2D" : "ConPlusMixPpointTwoFive",
    "T025_fit_largeNonly" : "TZeroTwoFive"
}

def main():
    lines = []
    for file_stem, command_stem in files_to_command.items():
        filename = Path(NUMBERS_FOLDER, file_stem+'.txt')
        with open(filename,'r') as f:
            c = f.readline()[2:].rstrip('\n')
            d = f.readline()[2:].rstrip('\n')
        lines.append(f"\\newcommand\\c{command_stem}{{{c}}}\n")
        lines.append(f"\\newcommand\\d{command_stem}{{{d}}}\n")
    with open(SAVE_FILENAME, "w+") as f:
        f.writelines(lines)

if __name__ == "__main__":
    main()