#pragma once
#include "PlotterCommon.hpp"

// TODO: make two: plotterlgrngn and plotter blk1m
template<int NDims>
class PlotterMicro : public PlotterCommon<NDims> 
{
  protected:
  using parent_t = PlotterCommon<NDims>;

  public:
  using arr_t = typename parent_t::arr_t;

  protected:
  std::string micro;
  arr_t res;
  arr_t rhod;
  const double L_evap = 2264.76e3; // latent heat of evaporation [J/kg]

  public:

  // functions for diagnosing fields
  //
  // aerosol droplets mixing ratio
  auto load_ra_timestep(
    int at
  ) //-> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
      return arr_t(this->h5load_timestep("aerosol_rw_mom3", at) * 4./3. * 3.1416 * 1e3);
    
    else if(this->micro == "blk_1m")
    {
      res = 0;
      return res;
     // return blitz::safeToReturn(res + 0);
    }
  }

  // cloud droplets mixing ratio
  auto load_rc_timestep(
    int at
  ) //-> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
      return arr_t(this->h5load_timestep("cloud_rw_mom3", at) * 4./3. * 3.1416 * 1e3);
    else if(this->micro == "blk_1m")
      return arr_t(this->h5load_timestep("rc", at));
    else if(this->micro == "blk_2m")
      return arr_t(this->h5load_timestep("rc", at));
  }

  // rain droplets mixing ratio
  auto load_rr_timestep(
    int at
  ) //-> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
      return arr_t(this->h5load_timestep("rain_rw_mom3", at) * 4./3. * 3.1416 * 1e3);
    else if(this->micro == "blk_1m")
      return arr_t(this->h5load_timestep("rr", at));
    else if(this->micro == "blk_2m")
      return arr_t(this->h5load_timestep("rr", at));
  }

  // activated drops mixing ratio
  auto load_ract_timestep(
    int at
  ) 
  {
    if(this->micro == "lgrngn")
    {
      return arr_t(
	       arr_t(this->h5load_timestep("cloud_rw_mom3", at) * 4./3. * 3.1416 * 1e3) + 
               arr_t(this->h5load_timestep("rain_rw_mom3", at) * 4./3. * 3.1416 * 1e3)
	     );
    }
    else if(this->micro == "blk_1m")
    {
      res = this->h5load_timestep("rc", at);
      res += arr_t(this->h5load_timestep("rr", at));
    }
    else if(this->micro == "blk_2m")
    {
      res = this->h5load_timestep("rc", at);
      res += arr_t(this->h5load_timestep("rr", at));
    }
   // return blitz::safeToReturn(res + 0);
    return res;
  }

  // cloud droplets concentration [1/cm^3]
  auto load_nc_timestep(
    int at
  ) //-> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
      return arr_t(this->multiply_by_rhod(arr_t(this->h5load_timestep("cloud_rw_mom0", at) / 1e6)));
    else if(this->micro == "blk_1m")
    {
      res = 0;
      return res;
     // return blitz::safeToReturn(res + 0);
    }
    else if(this->micro == "blk_2m")
      return arr_t(this->multiply_by_rhod(arr_t(this->h5load_timestep("nc", at) / 1e6)));
  }

  // precipitation flux [W/m2]
  auto load_prflux_timestep(
    int at
  )// -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
    {
      return arr_t(this->h5load_timestep("precip_rate", at)
              *  4./3 * 3.14 * 1e3 // to get mass
              / this->CellVol    // averaged over cell volume, TODO: make precip rate return specific moment? wouldnt need the dx and dy
              * L_evap);
    }
    else if(this->micro == "blk_1m")
      try
      {
        res = this->h5load_timestep("precip_rate", at); // precip_rate is the difference between influx and outflux
        for(int z = this->map["z"] - 2; z>=0; --z)
        {
          res(this->hrzntl_slice(z)) = res(this->hrzntl_slice(z+1)) - res(this->hrzntl_slice(z)); 
        }
        res *= rhod * this->map["dz"] * L_evap;
      }
      catch(...)
      {
        res = 0;
      }
   // return blitz::safeToReturn(res + 0);
     return res;
  }

  // RH
  auto load_RH_timestep(
    int at
  ) //-> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
      return arr_t(this->h5load_timestep("RH", at));
    else if(this->micro == "blk_1m")
      res = 0;
   // return blitz::safeToReturn(res + 0);
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
    nc *= rhod; // [1/m^3]
    nc *= this->dv;


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
    return this->hlpr(sdconc, at);
  }

  // mean and std_dev of temperature [K] (characteristics of the spatial distribution at this timestep, without near-wall cells)
  std::pair<double, double> T_stats_timestep(int at)
  {   
    // read theta away from walls 
    //arr_t tht(this->nowall(arr_t(this->h5load_timestep("th", at)), distance_from_walls));
    arr_t tht(arr_t(this->h5load_timestep("th", at)));
    tht *= pow(this->map_prof["p_e"](this->LastIndex) / p_1000, R_d / c_pd); // tht -> T
    return this->hlpr(tht, at);
  }

  // mean and std_dev of r_v [1] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> rv_stats_timestep(int at)
  {   
    arr_t rv(arr_t(this->h5load_timestep("rv", at)));
    return this->hlpr(rv, at);
  }

  // mean and std_dev of RH [1] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> RH_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    std::pair<double, double> res;

    arr_t RH(this->h5load_timestep("RH", at));
    return this->hlpr(RH, at);
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
  // accumulated volume precipitation [m^3]
  double calc_acc_surf_precip_volume(double prec_vol)
  {
    if(this->micro == "lgrngn")
      return calc_acc_surf_precip(prec_vol) * this->DomainSurf / 1000.;
    if(this->micro == "blk_1m")
      return calc_acc_surf_precip(prec_vol) * this->DomainSurf / 1000.;// to m^3 of water
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
  PlotterMicro(const string &file, const string &micro):
    parent_t(file),
    micro(micro),
    res(this->tmp.shape())
    //,rhod(this->h5load(file + "/const.h5", "G"))
  {
  }
};

