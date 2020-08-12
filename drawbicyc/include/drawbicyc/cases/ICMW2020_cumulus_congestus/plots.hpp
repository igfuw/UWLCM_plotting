#pragma once


const std::vector<std::string> series_ICMW2020({
  "lwp",
  "rwp",
  "cwp",
  "lwm",
  "cwm",
  "rwm",
  "surf_precip",
  "acc_precip",
  "acc_precip_vol",
  "cloud_base",
  "RH_max",
  "cloud_top_height",
  "total_droplets_number"
});
 
std::vector<std::string> profs_ICMW2020({
  "00rtot",
  "rliq",
  "thl",
  "wvar",
  "prflux",
  "clfrac",
  "sd_conc",
  "cl_nc",
  "cl_nc_up",
  "w",
  "u",
  "v",
  "base_prflux_vs_clhght",
  "non_gccn_rw_cl",
  "gccn_rw_cl,"
  //, "N_c",
  //,"vel_div"
  //, "nc_up"
  //,"sat_RH_up"
  //, "act_conc_up"
  //, "nc_down"
});

std::vector<std::string> fields_ICMW2020({
  "rl",
  "nc",
  "rr",
  "nr",
  //"ef", "na",
  "th",
  "rv",
  "u",
  "w"
  //"sd_conc",//, "r_dry",
  //"RH", "supersat",
});
