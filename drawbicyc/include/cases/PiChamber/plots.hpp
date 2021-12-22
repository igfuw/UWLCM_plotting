#pragma once

const std::vector<std::string> series_PiChamber({
// "clb_bigrain_mean_rd",
// "clb_bigrain_mean_kappa",
// "clb_bigrain_mean_conc",
// "clb_bigrain_mean_inclt",
//, "clb_bigrain_mean_gccn_fraction"
// "clfrac", 
 "lwp"
, "rwp"
, "surf_precip"

// "acc_precip",
, "nc"
, "ntot"
, "nr"
//, "cl_nr"
 

// "cloud_base",
// "surf_flux_latent",
// "surf_flux_sensible"
 ,"sd_conc"
//"mass_dry",
// ,"cl_gccn_conc", "gccn_conc"
// ,"cl_non_gccn_conc", "non_gccn_conc", "cl_gccn_to_non_gccn_conc_ratio"
//, "cl_gccn_meanr"
//,"cl_avg_cloud_rad"
//, "sd_conc_std_dev"
// "tot_water"

,"tot_tke"
,"tot_tke_nowall"
,"uw_tot_tke"
,"sgs_tke"
,"supersat_nowall"
,"uw_tot_tke_running_avg"
});

std::vector<std::string> profs_PiChamber({
  /*
"00rtot"
//, "rliq", "thl", "wvar", 
 ,"prflux"
//,"clfrac"
//,"sd_conc"
,"cl_nc"
//,"cl_nc_up"
//,"w"
//,"u", "v"
,"base_prflux_vs_clhght"
//, "non_gccn_rw_cl"
//, "gccn_rw_cl"
//,"sat_RH_up"
//, "N_c", 
//,"vel_div"
//, "nc_up" 
//,"sat_RH_up"
//, "act_conc_up" 
//, "nc_down" 
*/
}); // rtot has to be first


std::vector<std::string> fields_PiChamber({
  /*
"rl", "nc",
 "rr", "nr",
//"ef", "na", 
"th", "rv",     
"u", "w", 
//"sd_conc",//, "r_dry", 
//"RH", "supersat",
//"lib_pres", "lib_temp"
*/
});
