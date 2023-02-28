#pragma once
#include "Plotter2d.hpp"
#include "Plotter3d.hpp"
//#include <libcloudph++/common/moist_air.hpp>

template<int NDims>
class PlotterCommon : public Plotter_t<NDims> 
{
  protected:
  using parent_t = Plotter_t<NDims>;

  public:
  using arr_t = typename parent_t::arr_t;
  using parent_t::parent_t;

  protected:
  const double L_evap = 2264.76e3; // latent heat of evaporation [J/kg]

  public:

  // ---- functions for diagnosing statistics ----
  //
  arr_t& multiply_by_rhod(arr_t arr)
  {
    if(arr.extent(NDims-1) == this->map_prof["rhod"].extent(0))
      arr *= this->map_prof["rhod"](this->LastIndex);
    else if(arr.extent(NDims-1) == this->map_prof["refined rhod"].extent(0))
      arr *= this->map_prof["refined rhod"](this->LastIndex);
    else
      throw std::runtime_error("multiply_by_rhod: input array is neither normal grid size nor refined grid size");
    return arr;
  }

  arr_t& multiply_by_CellVol(arr_t arr)
  {
    if(arr.extent(NDims-1) == this->map_prof["rhod"].extent(0))
      arr *= this->CellVol;
    else if(arr.extent(NDims-1) == this->map_prof["refined rhod"].extent(0))
      arr *= this->CellVol_ref;
    else
      throw std::runtime_error("multiply_by_CellVol: input array is neither normal grid size nor refined grid size");
    return arr;
  }
  
  // helper function that calculates staistics (mean and std_dev) of a field
  std::pair<double, double> hlpr(arr_t arr, int at)
  {
    std::pair<double, double> res;
    res.first = blitz::mean(arr);
    arr = pow(arr - res.first, 2); 
    res.second = sqrt(blitz::mean(arr)); 
    return res;
  }

  // mean and std_dev of supersaturation in cells with positive supersaturation [%] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> positive_supersat_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    std::pair<double, double> res;

    // read RH 
    arr_t RH(this->h5load_timestep("RH", at));
    RH -= 1.; 
    arr_t sat(RH.copy());
    sat = iscloudy_sat(RH);
    RH *= sat; //apply the cloudiness mask
    RH *= 100; // to get %
    if(blitz::sum(sat) > 0)
      res.first = blitz::sum(RH) / blitz::sum(sat); 
    else
      res.first = 0;
  
    RH = pow(RH - res.first, 2); 
    RH *= sat; // apply filter
    if(res.first>0)
      res.second = sqrt(blitz::sum(RH) / blitz::sum(sat)); 
    else
      res.second = 0.; 

    return res;
  }

  // mean and std_dev of supersaturation at droplet locations (i.e. supersat weighted by the number of droplets)
  std::pair<double, double> drop_supersat_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    std::pair<double, double> res;

    // read supersat 
    arr_t ssat(this->h5load_timestep("RH", at));
    ssat -= 1.; 

    // read number of droplets per cell
    arr_t nc(this->h5load_timestep("cloud_rw_mom0", at));
    multiply_by_rhod(nc);
    multiply_by_CellVol(nc);


    if(blitz::sum(nc) > 0)
    {
      res.first = blitz::sum(ssat * nc) / blitz::sum(nc);
      ssat = pow(ssat - res.first, 2); 
      res.second = sqrt(blitz::sum(ssat * nc) / blitz::sum(nc));
    }
    else
    {
      res.first = 0;
      res.second = 0;
    }

    return res;
  }

  // mean and std_dev of number of SDs (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> sdconc_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    arr_t sdconc(this->h5load_timestep("sd_conc", at));
    return hlpr(sdconc, at);
  }

  // mean and std_dev of temperature [K] (characteristics of the spatial distribution at this timestep, without near-wall cells)
  std::pair<double, double> T_stats_timestep(int at)
  {   
    // read theta away from walls 
    //arr_t tht(this->nowall(arr_t(this->h5load_timestep("th", at)), distance_from_walls));
    arr_t tht(arr_t(this->h5load_timestep("th", at)));
    tht *= pow(this->map_prof["p_e"](this->LastIndex) / p_1000, R_d / c_pd); // tht -> T
    return hlpr(tht, at);
  }

  // mean and std_dev of r_v [1] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> rv_stats_timestep(int at)
  {   
    arr_t rv(arr_t(this->h5load_timestep("rv", at)));
    return hlpr(rv, at);
  }

  // mean and std_dev of RH [1] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> RH_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    std::pair<double, double> res;

    arr_t RH(this->h5load_timestep("RH", at));
    return hlpr(RH, at);
  }

  // surface precipitation since last output [mm/day]
  double calc_surf_precip(double prec_vol_diff)
  {
    if(this->micro == "lgrngn")
      return prec_vol_diff / this->DomainSurf / (double(this->map["outfreq"]) * this->map["dt"] / 3600. / 24.) * 1e3; // SDM
    if(this->micro == "blk_1m" || this->micro == "blk_2m")
      return prec_vol_diff / double(this->map["outfreq"]) // flux in [kg / m^3 / s] averaged over time since last output and over cells on the bottom
                     / (this->map["x"] * this->map["y"])
                     * 3600. * 24. // per day
                     * this->map["dz"]     // per m^2
                     / 1e3         // to m^3 of water
                     * 1e3;        // to mm
  }

  // accumulated surface precipitation [mm]
  double calc_acc_surf_precip(double prec_vol)
  {
    if(this->micro == "lgrngn")
      return prec_vol / this->DomainSurf * 1e3; 
    if(this->micro == "blk_1m" || this->micro == "blk_2m")
      return prec_vol * this->map["dt"] 
                     / (this->map["x"] * this->map["y"])
                     * this->map["dz"]     // per m^2
                     / 1e3         // to m^3 of water
                     * 1e3;        // to mm
  }

  // droplet removal rate (at boundaries) since last output [1/(cm^3 s)]
  double calc_prtcl_removal(double prtcl_removal_diff)
  {
    if(this->micro == "lgrngn")
      return prtcl_removal_diff / this->DomainVol / (double(this->map["outfreq"]) * this->map["dt"]) / 1e6;
    if(this->micro == "blk_1m")
      return 0;
  }

  // heat flux thru boundary since last output [W/m^2]
  double calc_heat_flux(double tot_th_diff, int z_idx) // input in [K]
  {
    if(this->micro == "lgrngn")
    {
      tot_th_diff *= pow(this->map_prof["p_e"](z_idx) / p_1000, R_d / c_pd); // tht -> T
      double ret = tot_th_diff * c_pd * this->map_prof["rhod"](z_idx)         // sum of th diff over boundary cells since last output (K) * c_pd * density 
                   * this->map["dz"] / ((this->map["x"]-1) * (this->map["y"]-1)) // multiply by cell volume and divide by domain surface area (without walls)
                   * (double(this->map["outfreq"]) * this->map["dt"]);    // divide by time since last output
      return ret;
    }
    if(this->micro == "blk_1m")
      return 0;
  }

  double calc_heat_flux_top(double mean_th_diff, bool errfix)
  {
    double tot_th_diff = mean_th_diff * this->map["x"] * this->map["y"]; 
    if(errfix)
    {
 //     th_diff += (this->map["x"] * this->map["y"] - 1) * 280; // to counter to error in tot_th_diff calculation in UWLCM
      tot_th_diff -= (2*this->map["x"] + 2*(this->map["y"] - 1)) * (280 - 285); // dont count side wall
    }
    return calc_heat_flux(tot_th_diff, this->map["z"]-1);
  }

  double calc_heat_flux_bot(double mean_th_diff, bool errfix)
  {
    double tot_th_diff = mean_th_diff * this->map["x"] * this->map["y"]; 
    if(errfix)
    {
 //     th_diff += (this->map["x"] * this->map["y"] - 1) * 299;
      tot_th_diff -= (2*this->map["x"] + 2*(this->map["y"] - 1)) * (299 - 285);
    }
    return calc_heat_flux(tot_th_diff, 0);
  }

  // kinematic rv flux thru boundary since last output [kg/kg * m / s]
  double calc_moist_flux(double tot_rv_diff) // rv change summed over horizontal plane  [kg/kg]
  {
    if(this->micro == "lgrngn")
    {
      return tot_rv_diff * this->map["dz"]
             * (double(this->map["outfreq"]) * this->map["dt"]);    // divide by time since last output
    }
    if(this->micro == "blk_1m")
      return 0;
  }

  double calc_moist_flux_top(double mean_rv_diff, bool errfix)
  {
    // 3D assumed here!
    double tot_rv_diff = mean_rv_diff * this->map["x"] * this->map["y"]; 
    if(errfix)
    {
//      rv_diff += (this->map["x"] * this->map["y"] - 1) * 0.0062192674278; // to counter to error in tot_rv_diff calculation in UWLCM
      tot_rv_diff -= (2*this->map["x"] + 2*(this->map["y"] - 1)) * (0.0062192674278 - 0.00611718803008); // dont count side wall
    }
    return calc_moist_flux(tot_rv_diff);
  }

  double calc_moist_flux_bot(double mean_rv_diff, bool errfix)
  {
    // 3D assumed here!
    double tot_rv_diff = mean_rv_diff * this->map["x"] * this->map["y"]; 
    if(errfix)
    {
//      rv_diff += (this->map["x"] * this->map["y"] - 1) * 0.0213489271007; // to counter to error in tot_rv_diff calculation in UWLCM
      tot_rv_diff -= (2*this->map["x"] + 2*(this->map["y"] - 1)) * (0.0213489271007 - 0.00611718803008); // dont count side wall
    }
    return calc_moist_flux(tot_rv_diff);
  }

  //ctor
  /*
  PlotterCommon(const string &file, const string &micro):
    parent_t(file),
    micro(micro),
    res(this->tmp.shape()),
    rhod(this->h5load(file + "/const.h5", "G"))
  {
  }
    */
};

