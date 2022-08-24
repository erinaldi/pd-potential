import numpy as np
import pandas as pd
import emcee
import copy
import gvar
from pathlib import Path
import math

# TODO:
# 1. replace filename derivation procedure
# 2. replace uft.get_wl_data
# 3. maybe how we compute averages...

class Point():
    STANDARD_PATH_FOLDER =  Path(Path.cwd(), "data-files") # depends on cwd definition!
    def __init__(self, n, m, p, t, n_t, WL_type, path_folder=STANDARD_PATH_FOLDER):
        self.n = n
        self.m = m
        self.p = p
        self.t = t
        self.n_t = n_t
        self.WL_type = WL_type
        # self.path_folder = path_folder
        # self.last_suffix = last_suffix
        self.WL_filename = self.__construct_filename(path_folder)
        self.WL_raw_data_unthermalised_head = np.array(pd.read_csv(self.WL_filename, header=0, index_col=0))
        if WL_type == "con":
            self.WL_raw_data_unthermalised_head = self.WL_raw_data_unthermalised_head / (self.n - self.m) * self.n
        elif WL_type == "dec":
            self.WL_raw_data_unthermalised_head = self.WL_raw_data_unthermalised_head / self.m * self.n
        elif WL_type == "mix":
            self.WL_raw_data_unthermalised_head = self.WL_raw_data_unthermalised_head /  self.m / (1 - self.m/self.n) # ??
        else:
            raise("WL_type is incorrect")

        self.WL_data_raw, self.therm_cut = self.__get_thermalisation_cut(self.WL_raw_data_unthermalised_head)

        self.WL_data_av = np.average(self.WL_data_raw,axis=0)
        self.WL_errs = np.zeros_like(self.WL_data_av)
        self.errs_validity = np.zeros_like(self.WL_data_av)
        self.WL_errs, self.errs_validity, self.autocor_taus = self.__get_errors()

    def __construct_filename(self, path_folder=STANDARD_PATH_FOLDER):
        p_str_few_dp = f"{self.p}".replace(".","")
        t_str_2dp = f"{self.t:.2f}".replace(".","")
        rel_filename=f"N{self.n}/S{self.n_t}/M{self.m}/T{t_str_2dp}/P{p_str_few_dp}/cor_{self.WL_type}.csv"
        filename = Path(path_folder, rel_filename)
        return filename

    def __get_errors(self):
        num_Ls = len(self.WL_data_av)
        Ls = np.arange(1, num_Ls+1)
        errs = np.zeros(num_Ls)
        valid = np.ones(num_Ls)
        taus = np.ones(num_Ls)
        for i, l in enumerate(Ls):
            w_data = self.WL_data_raw[:,l-1]
            try:
                tau = emcee.autocorr.integrated_time(w_data, c=3)[0]
                if tau > 1:
                    err = np.sqrt(tau / len(w_data)) * np.std(w_data)
                else:
                    err = np.sqrt(1 / len(w_data)) * np.std(w_data)
            except:
                tau = emcee.autocorr.integrated_time(w_data, c=3, quiet=True)[0]
                if tau > 1:
                    err = np.sqrt(tau / len(w_data)) * np.std(w_data)
                else:
                    err = np.sqrt(1 / len(w_data)) * np.std(w_data)
                valid[i] = 0
            errs[i] = err
            taus[i] = tau
        return errs, valid, taus

    def __get_thermalisation_cut(self, w_data):
        tau = emcee.autocorr.integrated_time(w_data[:,0], c=3, quiet=True)
        tau = max(tau, 1)
        cut = math.ceil(10 * tau)
        assert(cut < len(w_data))
        return w_data[cut:,], cut

# def combine_many_data_files(filename_base, last_letter):
#     suffices = [""] + list(string.ascii_uppercase[1:].split(last_letter)[0]) + [last_letter]
#     data = ()
#     for s in suffices:
#         pre, ext = filename_base.split(".")
#         new_data = np.loadtxt(pre+s+"."+ext)
#         if data != ():
#             data = np.concatenate((data, new_data))
#         else:
#             data = new_data
#     return data



class Points_Same_Type():
    def __init__(self, points, completely_confined=False, first_kind=False):
        assert(all(pt.p == points[0].p for pt in points))
        assert(all(pt.WL_type == points[0].WL_type for pt in points))
        self.p = points[0].p
        assert(self.p!=0 or completely_confined)
        self.completely_confined = completely_confined
        if completely_confined:
            self.WL_type = "comp-conf"
        else:
            self.WL_type = points[0].WL_type
        self.all_points = points

        self.num_Ls = len(points[0].WL_data_av)
        self.Ls = np.arange(1, self.num_Ls+1)

        # organise into by-n
        self.points_by_N = {}
        for pt in self.all_points:
            if pt.n in self.points_by_N:
                self.points_by_N[pt.n].append(pt)
            else:
                self.points_by_N[pt.n] = [pt]

        # organise into by-nt
        self.points_by_nt = {}
        for pt in self.all_points:
            if pt.n_t in self.points_by_nt:
                self.points_by_nt[pt.n_t].append(pt)
            else:
                self.points_by_nt[pt.n_t] = [pt]

        self.large_N_WL = {} # indexed by n_t; gives W[i]
        self.large_N_WL_err = {} # indexed by n_t; gives W[i]
        for nt, pts in self.points_by_nt.items():
            ns = [pt.n for pt in pts]
            ms = [pt.m for pt in pts]
            if self.completely_confined or first_kind:
                ms = [0 for pt in pts]
            WL_avs = np.array([pt.WL_data_av for pt in pts])
            WL_errs = np.array([pt.WL_errs for pt in pts]) # TODO: WL_errs is not built-in to the object! But it is important enough to be
            cs_this_nt = np.zeros(self.num_Ls)
            errs_this_nt = np.zeros(self.num_Ls)
            for i, l in enumerate(self.Ls):
                c, err = self._take_large_N(ns, ms, WL_avs[:,i], WL_errs[:,i])
                cs_this_nt[i] = c
                errs_this_nt[i] = err
            self.large_N_WL[nt] = cs_this_nt
            self.large_N_WL_err[nt] = errs_this_nt # indexed by n_t; gives W[i]

        self.WL_large_N_then_large_nt = np.zeros(self.num_Ls)
        self.WL_large_N_then_large_nt_err = np.zeros(self.num_Ls)
        if len(self.points_by_nt) > 1:
            for i,l in enumerate(self.Ls):
                c, cerr = self._take_large_nt(list(self.large_N_WL.keys()), np.array(list(self.large_N_WL.values()))[:,i], np.array(list(self.large_N_WL_err.values()))[:,i])
                self.WL_large_N_then_large_nt[i] = c
                self.WL_large_N_then_large_nt_err[i] = cerr
        else:
            self.WL_large_N_then_large_nt[:] = None
            self.WL_large_N_then_large_nt_err[:] = None

        if len(self.all_points) > 3:
            self.WL_continuum_largeN_together, self.WL_cln_together_err = self._take_large_N_large_nt_together()
        else:
            self.WL_continuum_largeN_together, self.WL_cln_together_err = None, None

    def _take_large_N(self, ns, ms, WLs, WL_errs):
        (m, c), cov = np.polyfit([1/n/(n-m) for n, m in zip(ns, ms)], WLs, deg=1, w=1/WL_errs, cov='unscaled')
        # (m, c), cov = np.polyfit([1/n/(n) for n, m in zip(ns, ms)], WLs, deg=1, w=1/WL_errs, cov='unscaled')
        cerr = np.sqrt(np.diag(cov))[1]
        return c, cerr

    def _take_large_nt(self, nts, WLs, WL_errs):
        (m, c), cov = np.polyfit([1/n_t for n_t in nts], WLs, deg=1, w=1/WL_errs, cov='unscaled')
        cerr = np.sqrt(np.diag(cov))[1]
        return c, cerr

    def _take_large_N_large_nt_together(self):
        import scipy.optimize as opt
        # for a given L

        interpolated_WLs = np.zeros(self.num_Ls)
        interpolated_errs = np.zeros(self.num_Ls)

        for li, l in enumerate(self.Ls):
            def ansatz(independent_vars, a, b, c, d):
                m, n, n_t = independent_vars
                return a + b / n / (n - m) + c / n_t + d / n / (n - m) / n_t
            # x0 = (0, 0, 0, 0)
            # build up array with (m, n, n_t, WL) for the erro function
            data_list = np.array(list(map(lambda p: (p.m, p.n, p.n_t, p.WL_data_av[li], p.WL_errs[li]), self.all_points))).T
            if self.WL_type == "comp-conf":
                data_list[0]=0
            popt, pcov = opt.curve_fit(ansatz, xdata=data_list[:3], ydata=data_list[3], sigma=data_list[4], absolute_sigma=True)
            interpolated_WLs[li] = popt[0]
            interpolated_errs[li] = np.sqrt(np.diag(pcov)[0])
        return interpolated_WLs, interpolated_errs




class Points_Combine_MixCon():
    # Maybe supply pointset instead??
    def __init__(self, points_mix, points_con):
        self.p = points_mix[0].p
        assert(self.p!=0)
        assert(all(pt.p == self.p for pt in points_mix))
        assert(all(pt.p == self.p for pt in points_con))
        self.WL_type = "con+mix"

        # self.all_points = points1
        self.num_Ls = len(points_mix[0].WL_data_av)
        self.Ls = np.arange(1, self.num_Ls+1)

        # for each N, n_t, we add the values W_con + W_mix to get the new W
        # we then do 2D extrapolation in this. We need the correct M values at each point...

        self.all_points = []
        for pt_mix in points_mix:
            pt_con = [ptc for ptc in points_con if ptc.n==pt_mix.n and ptc.n_t==pt_mix.n_t]
            assert(len(pt_con)==1)
            pt_con = pt_con[0]
            # just need averages and errors

            gmix = gvar.gvar(pt_mix.WL_data_av,pt_mix.WL_errs)
            gcon = gvar.gvar(pt_con.WL_data_av,pt_con.WL_errs)
            m_n_ratio = pt_mix.m / pt_mix.n
            new_wl = ((1 - m_n_ratio) * gcon + m_n_ratio * gmix) # ????
            
            new_pt = copy.deepcopy(pt_mix)
            new_pt.WL_type="mix+con"
            new_pt.WL_data_av = gvar.mean(new_wl)
            new_pt.WL_err = gvar.sdev(new_wl)
            new_pt.errs_validity = np.logical_and(pt_mix.errs_validity,pt_con.errs_validity)
            self.all_points.append(new_pt)

            # w_con_plus_mix_over_N = 1 / ((1 - m_n_ratio ) ** 2 + m_n_ratio * (1 - m_n_ratio ) ) * ((1 - m_n_ratio ) ** 2 * g_con + (m_n_ratio ) * (1 - m_n_ratio ) * g_mix)
        self.WL_continuum_largeN_together, self.WL_cln_together_err = self._take_large_N_large_nt_together()


    def _take_large_N_large_nt_together(self):
        import scipy.optimize as opt
        # for a given L

        interpolated_WLs = np.zeros(self.num_Ls)
        interpolated_errs = np.zeros(self.num_Ls)

        for li, l in enumerate(self.Ls):
            def ansatz(independent_vars, a, b, c, d):
                m, n, n_t = independent_vars
                return a + b / n / (n - m) + c / n_t + d / n / (n - m) / n_t
            # x0 = (0, 0, 0, 0)
            # build up array with (m, n, n_t, WL) for the erro function
            data_list = np.array(list(map(lambda p: (p.m, p.n, p.n_t, p.WL_data_av[li], p.WL_errs[li]), self.all_points))).T
            if self.WL_type == "comp-conf":
                data_list[0]=0
            popt, pcov = opt.curve_fit(ansatz, xdata=data_list[:3], ydata=data_list[3], sigma=data_list[4], absolute_sigma=True)
            interpolated_WLs[li] = popt[0]
            interpolated_errs[li] = np.sqrt(np.diag(pcov)[0])
        return interpolated_WLs, interpolated_errs


def main():
    pt = Point(16,8,0.2,0.29,16,"con")
    pass

if __name__ == "__main__":
    main()