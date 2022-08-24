import plot_2D_interpolation
import plot_P0
import plot_subtraction_method
import plot_T025
import plot_why_remove_N16
import plot_Wmix_plus_Wcon
import plot_histories
import make_WL_tables
import make_rmvN16_fit_table
import write_fit_values_to_tex

# # TODO: this should not just be about running plots, but should also output all the data into tables too
# # 1) get thermalisation estimates from autocorr time * 10
# # 2) plot action and other histories of important observables, indicating thermalisation cutoff, perhaps

# TODO:
# rename python files etc
# remove superfluous stuff, maybe, in common.py
# think about how to deal with thermalisation data... 1) read from file 2) calculate autocorrelation data and then strip...
# what do I need to do to make the tables, and where should I put it?
# # ISN'T THERE A DIFFERENT TAU FOR EACH L?? SO HOW DO I DEAL WITH THIS?
# make sure that all the points are included - no data left out

# TODO: tables
# table 4 - why remove N=16
# table 6 - first kind...
# table 7 - T-0.29, P=0.2, 0.25, data...
# table 8 ?
# table 9 T=0.25
# table 10,11,12 : Wdec, Wmix, Wcon with errors

# DONE - In draft: update all occurrences of c, d

# TODO:
# this plot: T = 0.29 with the constraint P = 0
# update captions

# TODO:
# T=0.25 data
# New N=64 data
# subtraction data...

# TODO
# 1. Do T=0.25 with error bounds and new data
# 2. Run all plots - perhaps with the central file
# 3. pull request on GIT (check for surplus code and comments...)
# 4. subtraction plot - i) gather the data ii) feed into the plot iii) decide on which set of N=24 data to use...

def main():
    print("Making plots")
    plot_2D_interpolation.main()
    plot_P0.main()
    plot_subtraction_method.main()
    plot_T025.main()
    plot_why_remove_N16.main()
    plot_Wmix_plus_Wcon.main()

    print("Making debug plots")
    plot_histories.main()

    print("Making tables")
    make_WL_tables.main()
    make_rmvN16_fit_table.main()

    print("Writing TeX definitions")
    write_fit_values_to_tex.main()

    print("Done")

if __name__ == "__main__":
    main()