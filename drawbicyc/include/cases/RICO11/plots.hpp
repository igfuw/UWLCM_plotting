#pragma once

const std::vector<std::string> series_rico({
"RH_max",
"cloud_cover",
"min_cloud_base",
"inversion_height_rico",
"tot_water",
"lwp",
"rwp",
"surf_precip",
"acc_precip",
"cl_nc",
"cl_nr",
"cl_avg_supersat",
"cl_avg_th",
"cl_avg_rv",
"cl_avg_cloud_meanr",
"cl_avg_cloud_stddevr",
"wvarmax",
 "surf_flux_latent",
 "surf_flux_sensible",
"sd_conc"


// "clb_bigrain_mean_rd",
// "clb_bigrain_mean_kappa",
// "clb_bigrain_mean_conc",
// "clb_bigrain_mean_inclt",
//, "clb_bigrain_mean_gccn_fraction"

//"cl_acnv25",
//"cl_accr25",
//,"nc"
//,"nr"
//TODO (po usprawnieniu cloud mask i ujednoliceniu tego:
/*
,"cl_avg_cloud_rad"
"cloud_avg_std_dev_act_rad"
 * */
//,"rd_lt_0.8um_conc"
//,"rd_geq_0.8um_conc"
//,"cl_rd_lt_0.8um_conc"
//,"cl_rd_geq_0.8um_conc"

/*
 "cloud_base",
 ,"cl_sd_conc"
//"mass_dry",
 ,"cl_gccn_conc", "gccn_conc"
 ,"cl_non_gccn_conc", "non_gccn_conc", "cl_gccn_to_non_gccn_conc_ratio"
//, "cl_gccn_meanr"
//,"cl_avg_cloud_rad"
// "cl_sd_conc_std_dev",
// "tot_water"
*/
});

std::vector<std::string> profs_rico({
"00rtot"
//,"rd_lt_0.8um_conc"
//,"rd_geq_0.8um_conc"
,"sd_conc"
, "rliq"
, "thl", "wvar" 
 ,"prflux"
,"clfrac"
,"cl_nc"
//,"cl_nc_up"
//,"w"
,"u", "v"
//,"base_prflux_vs_clhght"
//, "non_gccn_rw_cl"
//, "gccn_rw_cl"
//,"sat_RH_up"
//, "N_c", 
//,"vel_div"
//, "nc_up" 
//,"sat_RH_up"
//, "act_conc_up" 
//, "nc_down" 
}); // rtot has to be first


std::vector<std::string> fields_rico({
"rl", "nc",
 "rr", "nr",
//"ef", "na", 
"th", "rv",     
"u", "w", 
//"sd_conc",//, "r_dry", 
//"RH", "supersat",
//"lib_pres", "lib_temp"
});
