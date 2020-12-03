#pragma once

const std::vector<std::string> series_PiChamberICMW({
"T_nowall",
"Qv_nowall",
"RH_nowall",
"LWC_nowall",
"LWC_gm-3_nowall",
"N_drop_nowall",
"N_aerosol_nowall",
"N_removal",
"disp_r_nowall",
"r_mean1_nowall",
"r_mean2_nowall",
"Sigma2_S_nowall",
"Sigma2_T_nowall",
"Sigma2_Qv_nowall",
"epsilon_nowall"
});

std::vector<std::string> profs_PiChamberICMW({
"T",
"rv",
"RH"
}); // rtot has to be first


std::vector<std::string> fields_PiChamberICMW({
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
