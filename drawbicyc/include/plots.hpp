#pragma once

#include "cases/RICO11/plots.hpp"
#include "cases/Dycoms_RF02/plots.hpp"
#include "cases/moist_thermal/plots.hpp"
#include "cases/PiChamber/plots.hpp"
#include "cases/PiChamber/plots_icmw.hpp"
#include "cases/Lasher_Trapp/plots.hpp"
#include "cases/common_plots/GCCN_CCN_conc/plots.hpp"
#include "cases/Cumulus_Congestus/plots.hpp"

const std::vector<std::string> series_sgs({
// "tot_tke"
});

std::vector<std::string> profs_sgs({
/*
 "sgs_tke"
,"k_m"
,"sgs_tht_flux"
,"sgs_rv_flux"
//,"sgs_rc_flux"
,"sgs_u_flux"
*/
});


class Plots
{
  public:
    std::vector<std::string> series;
    std::vector<std::string> profs;
    std::vector<std::string> fields;

  Plots(const std::string &type, bool sgs)
  {
    if(type == "dycoms") { 
      profs.insert(profs.end(), profs_dycoms.begin(), profs_dycoms.end());
      series.insert(series.end(), series_dycoms.begin(), series_dycoms.end());
      fields.insert(fields.end(), fields_dycoms.begin(), fields_dycoms.end());
    }
    else if(type == "rico") { 
      profs.insert(profs.end(), profs_rico.begin(), profs_rico.end());
      series.insert(series.end(), series_rico.begin(), series_rico.end());
      fields.insert(fields.end(), fields_rico.begin(), fields_rico.end());
    }
    else if(type == "moist_thermal") { 
      profs.insert(profs.end(), profs_moist_thermal.begin(), profs_moist_thermal.end());
      series.insert(series.end(), series_moist_thermal.begin(), series_moist_thermal.end());
      fields.insert(fields.end(), fields_moist_thermal.begin(), fields_moist_thermal.end());
    }
    else if(type == "pi_chamber") { 
      profs.insert(profs.end(), profs_PiChamber.begin(), profs_PiChamber.end());
      series.insert(series.end(), series_PiChamber.begin(), series_PiChamber.end());
      fields.insert(fields.end(), fields_PiChamber.begin(), fields_PiChamber.end());
    }
    else if(type == "Lasher_Trapp") {
     // profs.insert(profs.end(), profs_Lasher_Trapp.begin(), profs_Lasher_Trapp.end());
      series.insert(series.end(), series_Lasher_Trapp.begin(), series_Lasher_Trapp.end());
      fields.insert(fields.end(), fields_Lasher_Trapp.begin(), fields_Lasher_Trapp.end());
    }
    else if(type == "cumulus_congestus") {
      profs.insert(profs.end(), profs_Cumulus_Congestus.begin(), profs_Cumulus_Congestus.end());
      series.insert(series.end(), series_Cumulus_Congestus.begin(), series_Cumulus_Congestus.end());
     // fields.insert(fields.end(), fields_Cumulus_Congestus.begin(), fields_Cumulus_Congestus.end());
    }
    else if(type == "pi_chamber_icmw") { 
      profs.insert(profs.end(), profs_PiChamberICMW.begin(), profs_PiChamberICMW.end());
      series.insert(series.end(), series_PiChamberICMW.begin(), series_PiChamberICMW.end());
      fields.insert(fields.end(), fields_PiChamberICMW.begin(), fields_PiChamberICMW.end());
    }
    else if(type == "gccn_ccn_conc") { 
      profs.insert(profs.end(), profs_gccn_ccn_conc.begin(), profs_gccn_ccn_conc.end());
      series.insert(series.end(), series_gccn_ccn_conc.begin(), series_gccn_ccn_conc.end());
      fields.insert(fields.end(), fields_gccn_ccn_conc.begin(), fields_gccn_ccn_conc.end());
    }
    /*
    else if(type == "base_prflux_vs_clhght") { 
      profs.insert(profs.end(), profs_base_prflux_vs_clhght.begin(), profs_base_prflux_vs_clhght.end());
    }*/
    else
      throw std::runtime_error("drawbicyc Plots.hpp: unknown 'type'.");
    
    if (sgs)// && type != "base_prflux_vs_clhght")
    {
      profs.insert(profs.end(), profs_sgs.begin(), profs_sgs.end());
      series.insert(series.end(), series_sgs.begin(), series_sgs.end());
    }
  }
};
