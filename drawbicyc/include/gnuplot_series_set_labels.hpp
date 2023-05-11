#pragma once
#include <UWLCM_plotters/common.hpp>
#include "gnuplot.hpp"

void gnuplot_series_set_labels(Gnuplot &gp, std::string plt)
{
  gp << "set yrange[*:*]\n";
  gp << "set xrange[*:*]\n";

  if (plt == "cl_acnv25_rico" || plt == "cl_acnv25_dycoms")
  {
  //  res_pos *= 60.;
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set title 'autoconersion rate (r>25um) (cloudy) [g/(m3*s}]'\n";
  }
  if (plt == "cloud_cover_dycoms" || plt == "cloud_cover_rico")
  {
  //  res_pos *= 60.;
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set title 'cloud cover'\n";
  }
  else if (plt == "inversion_height_rico")
  {
    gp << "set ylabel 'inversion_height [m]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'inversion height (RICO)'\n";
  }
  else if (plt == "ract_com")
  {
    gp << "set ylabel 'q_c center of mass [km]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'center of mass of cloud water'\n";
  }
  else if (plt == "th_com")
  {
    gp << "set ylabel 'th prtrb center of mass [km]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'center of mass height'\n";
  }
  else if (plt == "ract_avg")
  {
    gp << "set ylabel '<q_c>  [g/kg]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'avg cloud water mixing ratio'\n";
  }
  else if (plt == "ract_std_dev")
  {
    gp << "set ylabel 'sigma(q_c) / <q_c>'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'relative std dev of q_c'\n";
  }
  else if (plt == "cloud_avg_act_conc")
  {
    gp << "set ylabel '<N_c>  [1/cm^3]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'avg concentration of activated droplets'\n";
  }
  else if (plt == "cloud_std_dev_act_conc")
  {
    gp << "set ylabel 'sigma(N_c) / <N_c>'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'relative std dev N_c'\n";
  }
  else if (plt == "cloud_avg_supersat")
  {
    gp << "set ylabel '<S> [%]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'avg supersaturation'\n";
  }
  else if (plt == "cloud_std_dev_supersat")
  {
    gp << "set ylabel 'sigma(S) / <S>'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'relative std dev S'\n";
  }
  else if (plt == "cloud_avg_act_rad")
  {
    gp << "set ylabel '<r_{mean}> [um]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'average mean radius of activated droplets'\n";
  }
  else if (plt == "cloud_avg_std_dev_act_rad")
  {
    gp << "set ylabel '<sigma(r)> [um]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'average std dev of radius of activated droplets'\n";
  }
  else if (plt == "sd_conc")
  {
    gp << "set ylabel '<N_{SD}>'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'average SD number'\n";
  }
  else if (plt == "cl_sd_conc")
  {
    gp << "set ylabel '<N_{SD}>'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'average SD number (cloudy cells)'\n";
  }
  else if (plt == "cl_sd_conc_act")
  {
    gp << "set ylabel '<N_{SD}^{act}>'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'average activated SD number (cloudy cells)'\n";
  }
  else if (plt == "tot_water")
  {
    gp << "set title 'mean water mass density'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'mean (hydrometeors+vapor) water mass density [g/m^3]'\n";
  }
  else if (plt == "com_vel")
  {
    gp << "set title 'vertical velocity of the COM\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'w [m/s]'\n";
  }
  else if (plt == "com_supersat")
  {
    gp << "set title 'supersaturation at the COM\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'S [%]'\n";
  }
  else if (plt == "com_mom0")
  {
    gp << "set title 'cloud drops concentration at the center of mass\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'N_c [cm^{-3}]'\n";
  }
  else if (plt == "com_mom1")
  {
    gp << "set title 'mean wet radius at the center of mass\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel '<r> [micron]'\n";
  }
  else if (plt == "com_mom2")
  {
    gp << "set title 'std dev of r at the center of mass\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'std dev [micron]'\n";
  }
  else if (plt == "com_sd_conc")
  {
    gp << "set title 'number of SDs at the center of mass\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'SD no'\n";
  }
  else if (plt == "nc")
    gp << "set title 'average cloud drop conc [1/cm^3]'\n";
  else if (plt == "ntot")
    gp << "set title 'average particle conc [1/cm^3]'\n";
  else if (plt == "cl_nc" || plt == "cl_nc_rico")
  {
    gp << "set title 'average cloud drop conc [1/cm^3] in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "nr")
  {
    gp << "set title 'average rain drop conc [1/cm^3]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_nr" || plt == "cl_nr_rico")
  {
    gp << "set title 'average rain drop conc [1/cm^3] in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "wvarmax")
  {
    gp << "set title 'max variance of w [m^2 / s^2]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "surf_precip")
  {
    gp << "set title 'surface precipitation [mm/d]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "acc_precip")
  {
    gp << "set title 'accumulated surface precipitation [mm]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cloud_top_height")
  {
    gp << "set title 'cloud top height [m]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'height [m]'\n";
  }
  else if (plt == "total_droplets_number")
  {
    gp << "set title 'Total droplets number [#]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "acc_vol_precip")
  {
    gp << "set title 'accumulated volume precipitation [mm^3]'\n";
    gp << "set xlabel 'time [h]'\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "mass_dry")
    gp << "set title 'total dry mass [g]'\n";
  else if (plt == "RH_max")
  {
    gp << "set title 'max RH in the domain'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'RH'\n";
  }
  else if (plt == "lwp")
  {
    gp << "set title 'liquid water path [g / m^2]'\n";
    gp << "set xlabel 'time[h]'\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "lwm")
  {
    gp << "set title 'liquid water mass [kg]'\n";
    gp << "set xlabel 'time[h]'\n";
    gp << "set ylabel ''\n";
  }else if (plt == "cwm")
  {
    gp << "set title 'cloud water mass [kg]'\n";
    gp << "set xlabel 'time[h]'\n";
    gp << "set ylabel ''\n";
  }else if (plt == "rwm")
  {
    gp << "set title 'rain water mass [kg]'\n";
    gp << "set xlabel 'time[h]'\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "rwp")
  {
    gp << "set title 'rain water path [g / m^2]'\n";
    gp << "set xlabel 'time [h]'\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cwp")
  {
    gp << "set title 'cloud water path [g / m^2]'\n";
    gp << "set xlabel 'time [h]'\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cloud_base_dycoms")
  {
    gp << "set title 'cloud base [m]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "min_cloud_base_rico")
  {
    gp << "set title 'lowest cloud base [m]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_gccn_conc")
  {
    gp << "set title 'average gccn conc [1/cm^3] in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_gccn_to_non_gccn_conc_ratio")
  {
    gp << "set title 'average ratio of gccn to non-gccn conc in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_non_gccn_conc")
  {
    gp << "set title 'average non gccn conc [1/cm^3] in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_gccn_meanr")
  {
    gp << "set title 'average wet radius [um] of GCCNs in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_meanr")
  {
    gp << "set title 'average wet radius [um] of cloud droplets in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_avg_cloud_rad")
  {
    gp << "set title 'average wet radius [um] of cloud droplets in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "gccn_conc")
  {
    gp << "set title 'average gccn conc [1/cm^3]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "non_gccn_conc")
  {
    gp << "set title 'average non gccn conc [1/cm^3]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_rd_geq_0.8um_conc")
  {
    gp << "set title 'average conc with r_d>=0.8um in cloudy cells [1/cm^3]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "rd_geq_0.8um_conc")
  {
    gp << "set title 'average conc with r_d>=0.8um [1/cm^3]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "surf_flux_sensible")
  {
    gp << "set title 'sensible surf flux [W/m^2]'\n";
  }
  else if (plt == "surf_flux_latent")
  {
    gp << "set title 'latent surf flux [W/m^2]'\n";
  }
  else if (plt == "er")
  {
    gp << "set title 'entrainment rate [cm / s]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "tot_tke")
  {
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set title 'TKE (resolved + sgs) [m^2 / s^2]'\n";
  }
  else if (plt == "res_tke")
  {
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set title 'resolved TKE [m^2 / s^2]'\n";
  }
  else if (plt == "supersat")
  {
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set yrange[0:*]\n";
    gp << "set title 'supersaturation, >12.5 cm from walls [%]'\n";
  }
  else if (plt == "tot_tke")
  {
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set title 'TKE (resolved + sgs), >12.5 cm from walls [m^2 / s^2]'\n";
  }
  else if (plt == "uw_tot_tke")
  {
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set title 'TKE from u and w components (resolved + 2/3 SGS) [m^2 / s^2]'\n";
  }
  else if (plt == "uw_tot_tke_running_avg")
  {
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set title 'TKE from u and w, turb vs running avg. (resolved + 2/3 SGS) [m^2 / s^2]'\n";
  }
  else if (plt == "sgs_tke")
  {
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set title 'SGS TKE (Smg) [m^2 / s^2]'\n";
  }
  else if (plt == "sgs_tke_sd")
  {
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set title 'SGS TKE (SD) [m^2 / s^2]'\n";
  }
  else if (plt == "clb_bigrain_mean_inclt")
  {
    gp << "set title 'big rain (r>40um) at clbase: <time since activation> [s]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "clb_bigrain_mean_rd")
  {
    gp << "set title 'big rain (r>40um) at clbase: <dry radius> [um]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "clb_bigrain_mean_kappa")
  {
    gp << "set title 'big rain (r>40um) at clbase: <kappa>'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "clb_bigrain_mean_conc")
  {
    gp << "set title 'big rain (r>40um) at clbase: <conc.> [1/cc]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "clb_bigrain_mean_gccn_fraction")
  {
    gp << "set title 'big rain (r>40um) at clbase: <fraction with kappa>0.61>'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "T")
  {
    gp << "set title 'temperature [K]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "Qv")
  {
    gp << "set title 'q_v [kg/kg]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "RH")
  {
    gp << "set title 'RH [1]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "LWC")
  {
    gp << "set title 'LWC [g/kg]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "LWC_gm-3")
  {
    gp << "set title 'LWC [g/m^3]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "N_drop")
  {
    gp << "set title 'N_{drop} [cm^{-3}]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "N_aerosol")
  {
    gp << "set title 'N_{aerosol} [cm^{-3}]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "r_mean1")
  {
    gp << "set title '<r>_{drop} [um]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "r_mean2")
  {
    gp << "set title 'r^{effective}_{drop} [um'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "disp_r")
  {
    gp << "set title 'relative dispersion of r_{drop} [1]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "N_removal")
  {
    gp << "set title 'N_{drop} removal rate [cm^{-3} min^{-1}]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "Sigma2_S")
  {
    gp << "set title 'supersat variance [1]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "Sigma2_T")
  {
    gp << "set title 'temperature variance [K^2]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "Sigma2_Qv")
  {
    gp << "set title 'q_v variance [(kg/kg)^2]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "epsilon")
  {
    gp << "set title 'SGS TKE dissipation rate [m^2 s^{-3}]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "S_drop")
  {
    gp << "set title 'supersat. at droplet location\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "Sigma2_S_drop")
  {
    gp << "set title 'supersat. at droplet location variance\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "H_flux_t")
  {
    gp << "set title 'heat flux at the top [W/m$^2$]\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set xrange[0.001:*]\n";
  }
  else if (plt == "H_flux_b")
  {
    gp << "set title 'heat flux at the bot [W/m$^2$]\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set xrange[0.001:*]\n";
  }
  else if (plt == "qv_flux_t")
  {
    gp << "set title 'moist. flux at the top [kg/kg * m/s]\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set xrange[0.001:*]\n";
  }
  else if (plt == "qv_flux_b")
  {
    gp << "set title 'moist. flux at the bottom [kg/kg * m/s]\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set xrange[0.001:*]\n";
  }
}
