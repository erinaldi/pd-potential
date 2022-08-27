import plot_2D_interpolation
import plot_P0
import plot_subtraction_method
import plot_T025
import plot_why_remove_N16
import plot_Wmix_plus_Wcon
import plot_histories
import make_WL_tables
import make_T025_table
import make_rmvN16_fit_table
import write_fit_values_to_tex

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
    make_T025_table.main()
    make_rmvN16_fit_table.main()

    print("Writing TeX definitions")
    write_fit_values_to_tex.main()

    print("Done")

if __name__ == "__main__":
    main()