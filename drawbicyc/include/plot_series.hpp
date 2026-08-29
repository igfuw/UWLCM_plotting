#include <UWLCM_plotters/PlotterMicro.hpp>
#include <boost/tuple/tuple.hpp>
#include "plots.hpp"
#include "gnuplot_series_set_labels.hpp"
#include "gnuplot.hpp"

template<class Plotter_t>
void plot_series(Plotter_t plotter, Plots plots, std::string type)
{
  auto& n = plotter.map;
  auto& n_prof = plotter.map_prof;
  for(auto elem : n)
  {
     std::cout << elem.first << " " << elem.second << std::endl;
  }
  Gnuplot gp;
  string file = plotter.file + "_" + type + "_series.svg";
  int hor = min<int>(plots.series.size(), 4);
  int ver = double(plots.series.size()) / 4. + 0.99999;
  init_prof(gp, file, ver, hor); 

  string prof_file = plotter.file + "_" + type + "_series.dat";
  std::ofstream oprof_file(prof_file);

  // read in density
  auto tmp = plotter.h5load(plotter.file + "/const.h5", "G");
  typename Plotter_t::arr_t rhod(tmp);
  typename Plotter_t::arr_t rtot(rhod.shape());

  typename Plotter_t::arr_t res_tmp(rhod.shape());
  
  // for calculating running averages of u and w, needed in TKE calc in Pi chamber LES
  std::vector<typename Plotter_t::arr_t> prev_u_vec, prev_w_vec;
  // container for the running sum
  typename Plotter_t::arr_t run_sum_u(rhod.shape()), run_sum_w(rhod.shape());
  run_sum_u = 0;
  run_sum_w = 0;
  // number of timesteps over which the running avg of u and w is calculated (1 min interval)
  int run_avg_n_step = 60. / (n["dt"] * n["outfreq"]) + 0.5;


  // read opts
  po::options_description opts("profile plotting options");
  opts.add_options()
    ("series_start", po::value<int>()->default_value(0) , "time in sec when we start drawin series")
    ("series_end", po::value<int>()->default_value(0) , "time in sec when we end drawing series")
  ;
  po::variables_map vm; 
  handle_opts(opts, vm);

  int first_timestep =  vm["series_start"].as<int>() / n["dt"] / n["outfreq"];
  int last_timestep =  vm["series_end"].as<int>() / n["dt"] / n["outfreq"];
  if(last_timestep == 0) last_timestep = n["t"]-1;

  Array<double, 1> res_pos(last_timestep - first_timestep + 1),
    com_N_c(last_timestep - first_timestep + 1), // particles concentration at the center of mass
    com_miu(last_timestep - first_timestep + 1); // to keep mean particle radius at the center of mass
  Array<int, 1> com_z_idx(last_timestep - first_timestep + 1), 
    com_x_idx(last_timestep - first_timestep + 1); // index of the center of mass cell

  // save time steps to the series file
  oprof_file << "position" << endl;
  oprof_file << plotter.timesteps;

  std::map<std::string, bool> data_found;
  std::map<std::string, Array<double, 1>> res_series, res_series_std_dev;

  for (auto &plt : plots.series)
  {
    data_found[plt] = true;
    res_series.emplace(std::make_pair(plt, Array<double, 1>(last_timestep - first_timestep + 1)));
    res_series_std_dev.emplace(std::make_pair(plt, Array<double, 1>(last_timestep - first_timestep + 1)));
    res_series[plt] = 0;
    res_series_std_dev[plt] = 0;
  }

  res_pos = 0;

  double prec_vol = 0.;
  double prec_vol_prev;
  double removed_droplets = 0.;
  double removed_droplets_prev;
  double removed_ice = 0.;
  double removed_ice_prev;

  double tot_acc_acnv_prev = 0;
  double tot_acc_accr_prev = 0;

  for (int at = first_timestep; at <= last_timestep; ++at) // TODO: mark what time does it actually mean!
  {
    // used in pi chamber
    double th_change_top = 0.;
    double th_change_bot = 0.;
    double rv_change_top = 0.;
    double rv_change_bot = 0.;

    res_pos(at) = at * n["outfreq"] * n["dt"] / 3600.;

    // store accumulated precip volume
    prec_vol_prev = prec_vol;
    try
    {
      prec_vol = plotter.h5load_attr_timestep(at * n["outfreq"], "liquid_volume", "puddle");
    }
    catch(...){;}

    // store accumulated number of removed droplets
    removed_droplets_prev = removed_droplets;
    removed_ice_prev = removed_ice;
    try
    {
      removed_droplets = plotter.h5load_attr_timestep(at * n["outfreq"], "water_number", "puddle");
      removed_ice = plotter.h5load_attr_timestep(at * n["outfreq"], "ice_number", "puddle");
    }
    catch(...){;}

    // th and rv flues thru top and bot
    try
    {
      th_change_top = plotter.h5load_attr_timestep(at * n["outfreq"], "acc_mean_th_change_top");
      th_change_bot = plotter.h5load_attr_timestep(at * n["outfreq"], "acc_mean_th_change_bot");
      rv_change_top = plotter.h5load_attr_timestep(at * n["outfreq"], "acc_mean_rv_change_top");
      rv_change_bot = plotter.h5load_attr_timestep(at * n["outfreq"], "acc_mean_rv_change_bot");
    }
    catch(...){;}

    for (auto &plt : plots.series)
    {
      if (plt == "cloud_cover_dycoms")
      {
        try
        {
          auto tmp = plotter.h5load_rc_timestep(at * n["outfreq"]) * 1e3; //g/kg
          typename Plotter_t::arr_t snap(tmp); 
          snap += plotter.h5load_rr_timestep(at * n["outfreq"]) * 1e3; //g/kg
          snap *= rhod; // water per cubic metre (should be wet density...)
          plotter.k_i = blitz::sum(snap, plotter.LastIndex) * n["dz"]; // LWP [g/m2] in the column 
          plotter.k_i = where(plotter.k_i > 20 , 1 , 0); // cloudiness as in Ackermann et al. 
          res_series[plt](at) = blitz::mean(plotter.k_i);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cloud_cover_rico")
      {
        try
        {
          // cloud fraction (fraction of columns with at least one cloudy cell, i.e. cell with  q_c > 0.01 g/kg)
          auto tmp = plotter.h5load_rc_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap = iscloudy_rc_rico(snap); // find cells with rc>1e-5
          snap(plotter.hrzntl_slice(0)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux
          plotter.k_i = blitz::first((snap == 1), plotter.LastIndex); 
          plotter.k_i = blitz::sum(snap, plotter.LastIndex);
          plotter.k_i = where(plotter.k_i > 0, 1, 0);
          res_series[plt](at) = blitz::mean(plotter.k_i); 
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // max RH in the domain
      else if (plt == "RH_max")
      {
        try
        {
          // read RH 
          auto tmp = plotter.h5load_timestep("RH", at * n["outfreq"]);

          typename Plotter_t::arr_t snap(tmp);
          res_series[plt](at) = blitz::max(snap);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // r_act averaged over cloudy cells
      else if (plt == "ract_avg")
      {
        try
        {
          auto stats = plotter.cloud_ract_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.first;
          res_series_std_dev[plt](at) = stats.second;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // cloud top height
      else if (plt =="cloud_top_height")
      {
        try
        {
	  // cloud fraction (cloudy if ql > 1e-5))
	  auto tmp = plotter.h5load_rc_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap = iscloudy_rc_rico(snap);
          plotter.k_i = blitz::last((snap == 1), plotter.LastIndex);
          auto cloudy_column = plotter.k_i.copy();
          cloudy_column = blitz::sum(snap, plotter.LastIndex);
          cloudy_column = where(cloudy_column > 0, 1, 0);
          plotter.k_i = where(cloudy_column == 0, 0, plotter.k_i);
          if(blitz::sum(cloudy_column) > 0)
           res_series[plt](at) = blitz::max(plotter.k_i)*n["dz"];
         else
           res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // average sd_conc
      else if (plt == "sd_conc")
      {
        try
        {
          auto stats = plotter.sdconc_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.first;
          res_series_std_dev[plt](at) = stats.second;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // average sd_conc in cloudy cells
      else if (plt == "cl_sd_conc")
      {
        try
        {
          auto stats = plotter.cloud_sdconc_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.first;
          res_series_std_dev[plt](at) = stats.second;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // average activated sd_conc in clloudy cells
      else if (plt == "cl_sd_conc_act")
      {
        try
        {
          auto stats = plotter.cloud_sdconc_act_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.first;
          res_series_std_dev[plt](at) = stats.second;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "tot_water")
      {
        try
        {
/*  
          {
            auto tmp = plotter.h5load_timestep("aerosol_rw_mom3", at * n["outfreq"]) * 4./3. * 3.1416 * 1e3;
            typename Plotter_t::arr_t snap(tmp);
            snap *= rhod;
            res_series[plt](at) = blitz::mean(snap);
          }
*/  
          {
            auto tmp = plotter.h5load_timestep("cloud_rw_mom3", at * n["outfreq"]) * 4./3. * 3.1416 * 1e3;
            typename Plotter_t::arr_t snap(tmp);
            snap *= rhod;
            res_series[plt](at) += blitz::mean(snap);
          }
          {
            auto tmp = plotter.h5load_timestep("rain_rw_mom3", at * n["outfreq"]) * 4./3. * 3.1416 * 1e3;
            typename Plotter_t::arr_t snap(tmp);
            snap *= rhod;
            res_series[plt](at) += blitz::mean(snap);
          }
          {
            auto tmp = plotter.h5load_timestep("rv", at * n["outfreq"]);
            typename Plotter_t::arr_t snap(tmp);
            snap *= rhod;
            res_series[plt](at) += blitz::mean(snap);
          } 
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "ract_com")
      {
        // center of mass of activated droplets
        try
        {
          res_series[plt](at) = plotter.act_com_z_timestep(at * n["outfreq"]);
        }        
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "com_vel")
      {
        // vertical velocity at the center of mass of activated droplets
        try
        {
          auto tmp = plotter.h5load_ract_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          typename Plotter_t::arr_t snap2(tmp);
          typename Plotter_t::arr_t snap3(tmp);
          
          snap2 = snap2 * plotter.LastIndex;
          snap3 = snap3 * blitz::tensor::i;
          if(blitz::sum(snap) > 1e-3)
          {
            int z_idx = blitz::sum(snap2) / blitz::sum(snap); 
            int x_idx = blitz::sum(snap3) / blitz::sum(snap); 
            auto tmp2 = plotter.h5load_timestep("w", at * n["outfreq"]);
            typename Plotter_t::arr_t snap_mom(tmp2);
            res_series[plt](at) = snap_mom(x_idx, z_idx);
          } 
          else 
            res_series[plt](at) = 0.;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "com_supersat")
      {
        // supersaturation at the center of mass of activated droplets
        try
        {
          auto tmp = plotter.h5load_ract_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          typename Plotter_t::arr_t snap2(tmp);
          typename Plotter_t::arr_t snap3(tmp);
          
          snap2 = snap2 * plotter.LastIndex;
          snap3 = snap3 * blitz::tensor::i;
          if(blitz::sum(snap) > 1e-3)
          {
            int z_idx = blitz::sum(snap2) / blitz::sum(snap); 
            int x_idx = blitz::sum(snap3) / blitz::sum(snap); 
            auto tmp2 = plotter.h5load_timestep("RH", at * n["outfreq"]);
            typename Plotter_t::arr_t snap_mom(tmp2);
            res_series[plt](at) = snap_mom(x_idx, z_idx) - 1;
          } 
          else 
            res_series[plt](at) = 0.;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "com_mom0")
      {
        // 0th moment of rw distribution at the center of mass of activated droplets (particles concentration), 2D only
        try
        {
          auto tmp = plotter.h5load_ract_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          typename Plotter_t::arr_t snap2(tmp);
          typename Plotter_t::arr_t snap3(tmp);
          
          snap2 = snap2 * plotter.LastIndex;
          snap3 = snap3 * blitz::tensor::i;
          if(blitz::sum(snap) > 1e-3)
          {
            com_z_idx(at) = blitz::sum(snap2) / blitz::sum(snap); 
            com_x_idx(at) = blitz::sum(snap3) / blitz::sum(snap); 
            std::cout << at << ": (" << com_x_idx(at) << "," << com_z_idx(at) << ")" << std::endl;
            auto tmp2 = plotter.h5load_timestep("actrw_rw_mom0", at * n["outfreq"]);
            typename Plotter_t::arr_t snap_mom(tmp2);
            com_N_c(at) = snap_mom(com_x_idx(at), com_z_idx(at)); // 0th raw moment / mass [1/kg]
            snap_mom *= rhod; // now per m^3
            res_series[plt](at) = snap_mom(com_x_idx(at), com_z_idx(at));
          }
          else 
          {
            com_N_c(at) = 0.;
            res_series[plt](at) = 0.;
          }
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "com_mom1")
      {
        // mean droplet radius at the center of mass
        try
        {
          auto tmp = plotter.h5load_timestep("actrw_rw_mom1", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp); // 1st raw moment / mass [m / kg]
          if(com_N_c(at) > 0)
            res_series[plt](at) = snap(com_x_idx(at), com_z_idx(at)) / com_N_c(at);
          else
            res_series[plt](at) = 0.;
          com_miu(at) = res_series[plt](at); // mean radius [m]
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "com_mom2")
      {
        // std deviation of distribution of radius at center of mass
        try
        {
          auto tmp = plotter.h5load_timestep("actrw_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t zeroth_raw_mom(tmp); // 0th raw moment / mass [1 / kg]
          tmp = plotter.h5load_timestep("actrw_rw_mom1", at * n["outfreq"]);
          typename Plotter_t::arr_t first_raw_mom(tmp); // 1st raw moment / mass [m / kg]
          tmp = plotter.h5load_timestep("actrw_rw_mom2", at * n["outfreq"]);
          typename Plotter_t::arr_t second_raw_mom(tmp); // 2nd raw moment / mass [m^2 / kg]
          tmp = plotter.h5load_timestep("sd_conc", at * n["outfreq"]);
          typename Plotter_t::arr_t sd_conc(tmp); // number of SDs
          if(com_N_c(at) > 0)
          {
            double SD_no = sd_conc(com_x_idx(at), com_z_idx(at));
            if(SD_no > 1 && com_miu(at) > 0)
            {
              res_series[plt](at) = ( 
                SD_no / (SD_no - 1) /
                com_N_c(at) * (
                  second_raw_mom(com_x_idx(at), com_z_idx(at)) - 
                  2. * com_miu(at) * first_raw_mom(com_x_idx(at), com_z_idx(at)) + 
                  com_miu(at) * com_miu(at) * zeroth_raw_mom(com_x_idx(at), com_z_idx(at))
                )
              );
              
              // could not be true due to numerics?
              if(res_series[plt](at) > 0.) 
                res_series[plt](at) = sqrt(res_series[plt](at));
              else 
                res_series[plt](at) = 0.;
            }
          }
          else
            res_series[plt](at) = 0.;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "com_sd_conc")
      {
        // number of SDs at the center of mass
        try
        {
          tmp = plotter.h5load_timestep("sd_conc", at * n["outfreq"]);
          typename Plotter_t::arr_t sd_conc(tmp); // number of SDs
          res_series[plt](at) = sd_conc(com_x_idx(at), com_z_idx(at));
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "th_com")
      {
        // center of mass of temp perturb
        try
        {
          auto tmp = plotter.h5load_timestep("th", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          
          res_tmp = is_th_prtrb(snap); // find cells with th>300.1
          snap *= res_tmp; // apply filter
          res_tmp = snap * plotter.LastIndex * n["dz"];
          
          if(blitz::sum(res_tmp) > 0.)
            res_series[plt](at) = blitz::sum(res_tmp) / blitz::sum(snap); 
          else
            res_series[plt](at) = 0.;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "inversion_height_rico")
      {
        // height of cells with largest gradient of theta
        try
        {
          typename Plotter_t::arr_t th(plotter.h5load_timestep("th", at * n["outfreq"]));
          auto grad = plotter.cent_diff_vert(th);
          auto max_index = blitz::maxIndex(grad, plotter.LastIndex);
          res_series[plt](at) = (blitz::mean(max_index) + 1) * n["dz"];
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "nc")
      {
        // cloud droplet (0.5um < r < 25 um) concentration
        try
        {
          auto tmp = plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap /= 1e6; // per cm^3
          snap *= rhod; // b4 it was per milligram
          res_series[plt](at) = blitz::mean(snap); 
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "ntot")
      {
        // concentration of all particles
        try
        {
          auto tmp = plotter.h5load_timestep("all_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap /= 1e6; // per cm^3
          snap *= rhod; // b4 it was per milligram
          res_series[plt](at) = blitz::mean(snap); 
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "nr")
      {
        // rain droplet ( r > 25 um) concentration
        try
        {
          auto tmp = plotter.h5load_timestep("rain_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap /= 1e6; // per cm^3
          snap *= rhod; // b4 it was per milligram
          res_series[plt](at) = blitz::mean(snap); 
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cl_nc")
      {
        // cloud droplet (0.5um < r < 25 um) concentration in cloudy grid cells
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          auto tmp = plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          typename Plotter_t::arr_t snap2;
          snap2.resize(snap.shape());
          snap2=snap;
          snap = iscloudy(snap); // cloudiness mask
          snap2 *= snap;
          if(blitz::sum(snap) > 0)
            res_series[plt](at) = blitz::sum(snap2) / blitz::sum(snap); 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cl_nc_rico")
      {
        // cloud droplet (0.5um < r < 25 um) concentration in cloudy grid cells
        try
        {
	  // rico cloud mask
          auto tmp = plotter.h5load_rc_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap = iscloudy_rc_rico(snap); // find cells with rc>1e-5
          snap(plotter.hrzntl_slice(0)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux

          typename Plotter_t::arr_t snap2(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          snap2 *= rhod; // b4 it was specific moment
          snap2 /= 1e6; // per cm^3

          snap2 *= snap;
          if(blitz::sum(snap) > 0)
            res_series[plt](at) = blitz::sum(snap2) / blitz::sum(snap); 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cl_nr")
      {
        // rain drop (25um < r) concentration in cloudy grid cells
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          auto tmp = plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          auto tmp2 = plotter.h5load_timestep("rain_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap2(tmp2);
          snap2 *= rhod; // b4 it was specific moment
          snap2 /= 1e6; // per cm^3
          snap2 *= snap;
          if(blitz::sum(snap) > 0)
            res_series[plt](at) = blitz::sum(snap2) / blitz::sum(snap); 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cl_nr_rico")
      {
        // rain drop (25um < r) concentration in cloudy grid cells
        try
        {
	  // rico cloud mask
          auto tmp = plotter.h5load_rc_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap = iscloudy_rc_rico(snap); // find cells with rc>1e-5
          snap(plotter.hrzntl_slice(0)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux

          auto tmp2 = plotter.h5load_timestep("rain_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap2(tmp2);
          snap2 *= rhod; // b4 it was specific moment
          snap2 /= 1e6; // per cm^3

          snap2 *= snap;
          if(blitz::sum(snap) > 0)
            res_series[plt](at) = blitz::sum(snap2) / blitz::sum(snap); 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "total_droplets_number")
      {
        // Total number of Cloud and Rain Droplets
        try
        {
          auto tmp = plotter.h5load_timestep("actrw_rw_mom0", at *n["outfreq"]) * rhod;
	        typename Plotter_t::arr_t snap(tmp);
          res_series[plt](at) = blitz::sum(snap)*plotter.CellVol;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;} 
      }
      else if (plt == "cloud_base_dycoms")
      {
        // average cloud base in the domain
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          auto tmp = plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          snap(plotter.hrzntl_slice(0)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux
          plotter.k_i = blitz::first((snap == 1), plotter.LastIndex); 
          auto cloudy_column = plotter.k_i.copy();
          cloudy_column = blitz::sum(snap, plotter.LastIndex);
          cloudy_column = where(cloudy_column > 0, 1, 0);
          plotter.k_i = where(cloudy_column == 0, 0, plotter.k_i);
          if(blitz::sum(cloudy_column) > 0)
            res_series[plt](at) = double(blitz::sum(plotter.k_i)) / blitz::sum(cloudy_column) * n["dz"];
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cloud_base_precip_dycoms")
      {
        // domain-averaged (all columns) precipitation at the average cloud base height [mm/day]
        try
        {
          // -- average cloud base, almost exactly as in "cloud_base_dycoms"... --
          // cloud fraction (cloudy if N_c > 20/cm^3)
          auto tmp = plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          snap(plotter.hrzntl_slice(0)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux
          plotter.k_i = blitz::first((snap == 1), plotter.LastIndex); 
          auto cloudy_column = plotter.k_i.copy();
          cloudy_column = blitz::sum(snap, plotter.LastIndex);
          cloudy_column = where(cloudy_column > 0, 1, 0);
          plotter.k_i = where(cloudy_column == 0, 0, plotter.k_i);
          int cloud_base_idx;
          if(blitz::sum(cloudy_column) > 0)
            cloud_base_idx = double(blitz::sum(plotter.k_i)) / blitz::sum(cloudy_column) + 0.5;
          else
            cloud_base_idx = 0;

          if(cloud_base_idx == 0)
            res_series[plt](at) = 0;
          else
          {
            // -- precipitation at this height averaged over all cells, cloudy or not -- 
            auto prflux = plotter.h5load_prflux_timestep(at * n["outfreq"]); // prflux in [W/m^2]
            res_series[plt](at) = blitz::mean(prflux(plotter.hrzntl_slice(cloud_base_idx))) / 2264.705 * 3.6 * 24; // convert to [mm/day]
          }
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "min_cloud_base_rico")
      {
        // lowest cloud base in the domain
        try
        {
          // cloud fraction (cloudy if r_c > 1e-5)
          typename Plotter_t::arr_t snap(plotter.h5load_rc_timestep(at * n["outfreq"]));
          snap = iscloudy_rc_rico(snap); 
          snap(plotter.hrzntl_slice(0)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux
          plotter.k_i = blitz::first((snap == 1), plotter.LastIndex); 
          auto cloudy_column = plotter.k_i.copy();
          cloudy_column = blitz::sum(snap, plotter.LastIndex);
          cloudy_column = where(cloudy_column > 0, 1, 0);
          plotter.k_i = where(cloudy_column == 0, 1e6, plotter.k_i); // 1e6 denotes no clouds in the column
          if(blitz::sum(cloudy_column) > 0)
            res_series[plt](at) = blitz::min(plotter.k_i) * n["dz"];
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // average concentration of activated droplets in cloudy cells
      else if (plt == "cloud_avg_act_conc")
      {
        try
        {
          auto stats = plotter.cloud_actconc_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.first;
          res_series_std_dev[plt](at) = stats.second;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // average supersaturation in cells with S>0
      else if (plt == "cloud_avg_supersat")
      {
        try
        {
          auto stats = plotter.positive_supersat_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.first;
          res_series_std_dev[plt](at) = stats.second;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // spatial average of mean radius of activated droplets in cloudy cells
      else if (plt == "cl_avg_cloud_rad")
      {
        try
        {
          auto stats = plotter.cloud_meanr_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.first;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      // spatial average of standard deviation of the acttivated droplets radius distribution in cloudy cells
      else if (plt == "cloud_avg_std_dev_act_rad")
      {
        try
        {
          auto stats = plotter.cloud_stddevr_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.first;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else if (plt == "mass_dry")
      {
        // total dry mass
        double rho_dry = 1769; //[kg/m^3] - density of ammonium sulfate from wikipedia
        try
        {
          auto tmp = plotter.h5load_timestep("rd_rng000_mom3", at * n["outfreq"]) * 4./3. * 3.14 * rho_dry * 1e3;
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod * plotter.CellVol; // turn mixing ratio in g/kg to total mass in g
          res_series[plt](at) = blitz::sum(snap); 
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "surf_precip")
      {
        // surface precipitation [mm/day]
        try
        {
          res_series[plt](at) = plotter.calc_surf_precip(prec_vol - prec_vol_prev);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "acc_precip")
      {
        // accumulated surface precipitation [mm]
        try
        {
          res_series[plt](at) = plotter.calc_acc_surf_precip(prec_vol);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "acc_vol_precip")
      {
        // accumulated surface precipitation [m^3]
        try
        {
          res_series[plt](at) = plotter.calc_acc_surf_precip_volume(prec_vol);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cl_acnv25_dycoms")
      {
        // autconversion rate with rain threshold r=25um (cloudy cells as in Dycoms) [g/(m3*s)]
        // rather coarse estimate, sum of acnv accumulated over ALL cells since the last output
        // is divided by the instantaneous volume of all cloudy cells
        // TODO: output instantaneous acnv rate in libcloud, not the accumulated one?
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask

          typename Plotter_t::arr_t acc_acnv(plotter.h5load_timestep("acc_acnv25", at * n["outfreq"]));
          auto tot_acc_acnv = blitz::sum(acc_acnv);

          if(blitz::sum(snap) > 0)
            res_series[plt](at) =  4./3. * 3.1416 * 1e6 * (tot_acc_acnv - tot_acc_acnv_prev) / ((blitz::sum(snap) * plotter.CellVol) * (n["outfreq"] * n["dt"])); 
          else
            res_series[plt](at) = 0;

          tot_acc_acnv_prev = tot_acc_acnv;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cl_accr25_dycoms")
      {
        // accretion rate with rain threshold r=25um (cloudy cells as in Dycoms) [g/(m3*s)]
        // rather coarse estimate, sum of accr accumulated over ALL cells since the last output
        // is divided by the instantaneous volume of all cloudy cells
        // TODO: output instantaneous accr rate in libcloud, not the accumulated one?
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask

          typename Plotter_t::arr_t acc_accr(plotter.h5load_timestep("acc_accr25", at * n["outfreq"]));
          double tot_acc_accr = blitz::sum(acc_accr); 

          if(blitz::sum(snap) > 0)
            res_series[plt](at) = 4./3. * 3.14166 * 1e6 * (tot_acc_accr - tot_acc_accr_prev) / ((blitz::sum(snap) * plotter.CellVol) * (n["outfreq"] * n["dt"])); 
          else
            res_series[plt](at) = 0;

          tot_acc_accr_prev = tot_acc_accr;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cl_acnv25_rico")
      {
        // autconversion rate with rain threshold r=25um (cloudy cells as in Rico) [g/(m3*s)]
        // rather coarse estimate, sum of acnv accumulated over ALL cells since the last output
        // is divided by the instantaneous volume of all cloudy cells
        // TODO: output instantaneous acnv rate in libcloud, not the accumulated one?
        try
        {
          // cloud fraction (cloudy if r_c > 1e-5)
          typename Plotter_t::arr_t snap(plotter.h5load_rc_timestep(at * n["outfreq"]));
          snap = iscloudy_rc_rico(snap); 

          typename Plotter_t::arr_t acc_acnv(plotter.h5load_timestep("acc_acnv25", at * n["outfreq"]));
          auto tot_acc_acnv = blitz::sum(acc_acnv);

          if(blitz::sum(snap) > 0)
            res_series[plt](at) =  4./3. * 3.1416 * 1e6 * (tot_acc_acnv - tot_acc_acnv_prev) / ((blitz::sum(snap) * plotter.CellVol) * (n["outfreq"] * n["dt"])); 
          else
            res_series[plt](at) = 0;

          tot_acc_acnv_prev = tot_acc_acnv;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cl_accr25_rico")
      {
        // accretion rate with rain threshold r=25um (cloudy cells as in Rico) [g/(m3*s)]
        // rather coarse estimate, sum of accr accumulated over ALL cells since the last output
        // is divided by the instantaneous volume of all cloudy cells
        // TODO: output instantaneous accr rate in libcloud, not the accumulated one?
        try
        {
          // cloud fraction (cloudy if r_c > 1e-5)
          typename Plotter_t::arr_t snap(plotter.h5load_rc_timestep(at * n["outfreq"]));
          snap = iscloudy_rc_rico(snap); 

          typename Plotter_t::arr_t acc_accr(plotter.h5load_timestep("acc_accr25", at * n["outfreq"]));
          double tot_acc_accr = blitz::sum(acc_accr); 

          if(blitz::sum(snap) > 0)
            res_series[plt](at) = 4./3. * 3.14166 * 1e6 * (tot_acc_accr - tot_acc_accr_prev) / ((blitz::sum(snap) * plotter.CellVol) * (n["outfreq"] * n["dt"])); 
          else
            res_series[plt](at) = 0;

          tot_acc_accr_prev = tot_acc_accr;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "lwp")
      {   
        // liquid water path
        try
        {
          {
            auto tmp = plotter.h5load_rc_timestep(at * n["outfreq"]) * 1e3; //g/kg
            typename Plotter_t::arr_t snap(tmp); 
            snap += plotter.h5load_rr_timestep(at * n["outfreq"]) * 1e3; //g/kg
            snap *= rhod; // water per cubic metre (should be wet density...)
            res_series[plt](at) = blitz::mean(snap); 
          }
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }   
      else if (plt == "rwp")
      {   
        // rain water path
        try
        {
          {
            typename Plotter_t::arr_t snap(plotter.h5load_rr_timestep(at * n["outfreq"]));
            snap *= rhod * 1e3; // water per cubic metre (should be wet density...) & g/kg
            res_series[plt](at) = blitz::mean(snap); 
          }
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }       
      else if (plt == "cwp")
      {   
        // cloud water path
        try
        {
          {
            typename Plotter_t::arr_t snap(plotter.h5load_rc_timestep(at * n["outfreq"]));
            snap *= rhod * 1e3; // water per cubic metre (should be wet density...) & g/kg
            res_series[plt](at) = blitz::mean(snap); 
          }
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }       
      else if (plt == "lwm")
      {
        //liquid water mass
        try
        {
          { 
            auto tmp = plotter.h5load_rc_timestep(at * n["outfreq"]) * rhod;
            typename Plotter_t::arr_t snap(tmp);
            snap += plotter.h5load_rr_timestep(at * n["outfreq"]) * rhod;
            snap *= plotter.CellVol;
            snap(plotter.hrzntl_slice(0)) = snap(plotter.hrzntl_slice(0))/2;
            snap(plotter.hrzntl_slice(-1)) = snap(plotter.hrzntl_slice(-1))/2;
            res_series[plt](at) = blitz::sum(snap);
          }
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cwm")
      {
	//cloud water mass
	try
	{
	  { 
            auto tmp = plotter.h5load_rc_timestep(at * n["outfreq"]);
            typename Plotter_t::arr_t snap(tmp);
            snap *= rhod * plotter.CellVol;
            snap(plotter.hrzntl_slice(0)) = snap(plotter.hrzntl_slice(0))/2;
            snap(plotter.hrzntl_slice(-1)) = snap(plotter.hrzntl_slice(-1))/2;
            res_series[plt](at) = blitz::sum(snap);
          }
	}
	catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "rwm")
      {
	//liquid water mass
	try
	{
	  {
	    auto tmp = plotter.h5load_rr_timestep(at * n["outfreq"]);
            typename Plotter_t::arr_t snap(tmp);
            snap *= rhod * plotter.CellVol;
            snap(plotter.hrzntl_slice(0)) = snap(plotter.hrzntl_slice(0))/2;
            snap(plotter.hrzntl_slice(-1)) = snap(plotter.hrzntl_slice(-1))/2;
            res_series[plt](at) = blitz::sum(snap);
	  }
	}
	catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "surf_flux_latent")
      {   
        try
        {
          {
            typename Plotter_t::arr_t snap(plotter.h5load_timestep("latent surface flux", at * n["outfreq"], true)); 
            res_series[plt](at) = blitz::mean(snap); 
          }
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }   
      else if (plt == "surf_flux_sensible")
      {   
        try
        {
          {
            typename Plotter_t::arr_t snap(plotter.h5load_timestep("sensible surface flux", at * n["outfreq"], true)); 
            res_series[plt](at) = blitz::mean(snap); 
          }
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }   
      else if (plt == "er")
      {   
        //entrainment rate as in the 2009 Ackerman paper
        // to store total mixingg ratio
        try
        {
          {
            auto tmp = plotter.h5load_rc_timestep(at * n["outfreq"]) * 1e3; //g/kg
            typename Plotter_t::arr_t snap(tmp); 
            snap += plotter.h5load_rr_timestep(at * n["outfreq"]) * 1e3; //g/kg
            rtot = snap;
          }
          {
            auto tmp = plotter.h5load_timestep("rv", at * n["outfreq"]) * 1e3;
            typename Plotter_t::arr_t snap(tmp); // vapor mixing ratio [g/kg]
            rtot += snap;
          }
          plotter.k_i = 0;
          plotter.k_i = blitz::first((rtot < 8.), plotter.LastIndex); 
          res_series[plt](at) = blitz::mean(plotter.k_i);
        }
        catch (...) {;}
      }
      else if (plt == "wvarmax")
      {
        // maximum variance of vertical velocity
        try
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          Array<double, 1> mean(n["z"]);
          snap = snap * snap; // 2nd power, w_mean = 0
          // mean variance of w in horizontal
//          mean = blitz::mean(snap(tensor::j, tensor::i), tensor::j); // mean over x and y
          mean = plotter.horizontal_mean(snap);
          res_series[plt](at) = blitz::max(mean); // the max value
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "tot_tke")
      {
        try
        {
          typename Plotter_t::arr_t u(plotter.h5load_timestep("u", at * n["outfreq"]));
          plotter.subtract_horizontal_mean(u);
          u = u * u;
          res_series[plt](at) = blitz::mean(plotter.horizontal_mean(u));

          typename Plotter_t::arr_t w(plotter.h5load_timestep("w", at * n["outfreq"]));
          plotter.subtract_horizontal_mean(w);
          w = w * w;
          res_series[plt](at) += blitz::mean(plotter.horizontal_mean(w));

          if (Plotter_t::n_dims > 2)
          {
            typename Plotter_t::arr_t v(plotter.h5load_timestep("v", at * n["outfreq"]));
            plotter.subtract_horizontal_mean(v);
            v = v * v;
            res_series[plt](at) += blitz::mean(plotter.horizontal_mean(v));
          }
          
          res_series[plt](at) *= 0.5; // * n["dz"];

          typename Plotter_t::arr_t tke(plotter.h5load_timestep("tke", at * n["outfreq"]));
          typename Plotter_t::arr_t snap;
          snap.reference(tke);

          res_series[plt](at) += blitz::mean(plotter.horizontal_mean(snap));
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "uw_tot_tke") // As in the Thomas et al. 2019 paper about Pi chamber LES, but without substractin the running average.
      {
        try
        {
          typename Plotter_t::arr_t u(plotter.h5load_timestep("u", at * n["outfreq"]));
          plotter.subtract_horizontal_mean(u);
          u = u * u;
          res_series[plt](at) = blitz::mean(plotter.horizontal_mean(u));

          typename Plotter_t::arr_t w(plotter.h5load_timestep("w", at * n["outfreq"]));
          plotter.subtract_horizontal_mean(w);
          w = w * w;
          res_series[plt](at) += blitz::mean(plotter.horizontal_mean(w));

          res_series[plt](at) *= 0.5;// * n["dz"];

          typename Plotter_t::arr_t tke(plotter.h5load_timestep("tke", at * n["outfreq"]));
          if (Plotter_t::n_dims == 3)
          {
            // assume that sgs tke is isotropic, hence 2/3 are in the uw plane
            res_series[plt](at) += 2./3. * blitz::mean(plotter.horizontal_mean(tke));
          }
          if (Plotter_t::n_dims == 2)
          {
            res_series[plt](at) += blitz::mean(plotter.horizontal_mean(tke));
          }
          
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "uw_tot_tke_running_avg") // As in the Thomas et al. 2019 paper about Pi chamber LES. 
      {
        try
        {
          typename Plotter_t::arr_t u(plotter.h5load_timestep("u", at * n["outfreq"]));
          typename Plotter_t::arr_t w(plotter.h5load_timestep("w", at * n["outfreq"]));
          int step_no = at - first_timestep;
          if(step_no == 0) // the first step
          {
            prev_u_vec.push_back(u);
            run_sum_u += u;

            prev_w_vec.push_back(w);
            run_sum_w += w;
          }
          else // not the first step
          {
            if(step_no < run_avg_n_step) // less steps so far than the averaging window
            {
              u -= run_sum_u / step_no;
              prev_u_vec.push_back(u); 
              run_sum_u += u;

              w -= run_sum_w / step_no;
              prev_w_vec.push_back(w); 
              run_sum_w += w;
            }
            else // full averaging window
            {
              u -= run_sum_u / run_avg_n_step;
              // substract oldest u from the running sum
              int oldest_position = step_no % run_avg_n_step;
              run_sum_u -= prev_u_vec.at(oldest_position);
              // add current u to the running sum
              run_sum_u += u;
              // replace the oldest u with current
              prev_u_vec.at(oldest_position) = u;

              w -= run_sum_w / run_avg_n_step;
              // substract oldest u from the running sum
              run_sum_w -= prev_w_vec.at(oldest_position);
              // add current u to the running sum
              run_sum_w += w;
              // replace the oldest u with current
              prev_w_vec.at(oldest_position) = w;
            }
          }

          u = u * u;
          w = w * w;
          res_series[plt](at) = blitz::mean(plotter.horizontal_mean(u));
          res_series[plt](at) = blitz::mean(plotter.horizontal_mean(w));

          res_series[plt](at) *= 0.5;// * n["dz"];

          typename Plotter_t::arr_t tke(plotter.h5load_timestep("tke", at * n["outfreq"]));
          if (Plotter_t::n_dims == 3)
          {
            // assume that sgs tke is isotropic, hence 2/3 are in the uw plane
            res_series[plt](at) += 2./3. * blitz::mean(plotter.horizontal_mean(tke));
          }
          if (Plotter_t::n_dims == 2)
          {
            res_series[plt](at) += blitz::mean(plotter.horizontal_mean(tke));
          }
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "sgs_tke") // TODO: make it nowall?
      {
        try
        {
          typename Plotter_t::arr_t tke(plotter.h5load_timestep("tke", at * n["outfreq"]));
          res_series[plt](at) = blitz::mean(plotter.horizontal_mean(tke));
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "sgs_tke_sd") // TKE of SGS motion of SD (turb_adve) 
      {
        try
        {
          typename Plotter_t::arr_t tot_m0(plotter.h5load_timestep("aerosol_rw_mom0", at * n["outfreq"]));
          tot_m0 += typename Plotter_t::arr_t(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          tot_m0 += typename Plotter_t::arr_t(plotter.h5load_timestep("rain_rw_mom0", at * n["outfreq"]));

          typename Plotter_t::arr_t tke(plotter.h5load_timestep("all_up_mom2", at * n["outfreq"]));
          tke += typename Plotter_t::arr_t(plotter.h5load_timestep("all_wp_mom2", at * n["outfreq"]));

          if (Plotter_t::n_dims > 2)
            tke += typename Plotter_t::arr_t(plotter.h5load_timestep("all_vp_mom2", at * n["outfreq"]));

          tke = blitz::where(tot_m0 > 0., 0.5 * tke / tot_m0, 0); // tke in each cell

          res_series[plt](at) = blitz::sum(tke) / blitz::count(tot_m0 > 0.);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cl_gccn_conc")
      {
        // gccn (r_d > 2 um) concentration in cloudy grid cells
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          typename Plotter_t::arr_t snap2(plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]));
          snap2 *= rhod; // b4 it was specific moment
          snap2 /= 1e6; // per cm^3
          snap2 *= snap;
          if(blitz::sum(snap) > 0)
            res_series[plt](at) = blitz::sum(snap2) / blitz::sum(snap); 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cl_non_gccn_conc")
      {
        // gccn (r_d < 2 um) concentration in cloudy grid cells
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          typename Plotter_t::arr_t snap2(plotter.h5load_timestep("non_gccn_rw_mom0", at * n["outfreq"]));
          snap2 *= rhod; // b4 it was specific moment
          snap2 /= 1e6; // per cm^3
          snap2 *= snap;
          if(blitz::sum(snap) > 0)
            res_series[plt](at) = blitz::sum(snap2) / blitz::sum(snap); 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cl_gccn_to_non_gccn_conc_ratio")
      {
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          typename Plotter_t::arr_t snap2(plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]));
          typename Plotter_t::arr_t snap3(plotter.h5load_timestep("non_gccn_rw_mom0", at * n["outfreq"]));

          snap2 = where(snap3 > 0, snap2 / snap3, 0); // even if snap3=0, they are noncloudy anyway
          snap2 *= snap;

          if(blitz::sum(snap) > 0)
            res_series[plt](at) = blitz::sum(snap2) / blitz::sum(snap); 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "cl_gccn_meanr")
      {
        // gccn (r_d > 2 um) mean radius in cloudy grid cells
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          typename Plotter_t::arr_t snap_m0(plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]));
          typename Plotter_t::arr_t snap_m1(plotter.h5load_timestep("gccn_rw_mom1", at * n["outfreq"]) * 1e6); // in microns
          snap_m0 *= snap;
          snap_m1 *= snap;
          auto tot_gccn_m0 = blitz::sum(snap_m0);
          if(tot_gccn_m0 > 0)
            res_series[plt](at) = blitz::sum(snap_m1) / tot_gccn_m0; 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "gccn_conc")
      {
        // gccn (r_d > 2 um) concentration
        try
        {
          typename Plotter_t::arr_t snap2(plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]));
          snap2 /= 1e6; // per cm^3
          snap2 *= rhod; // b4 it was per milligram
          res_series[plt](at) = blitz::mean(snap2); 
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "non_gccn_conc")
      {
        // gccn (r_d > 2 um) concentration
        try
        {
          typename Plotter_t::arr_t snap2(plotter.h5load_timestep("non_gccn_rw_mom0", at * n["outfreq"]));
          snap2 /= 1e6; // per cm^3
          snap2 *= rhod; // b4 it was per milligram
          res_series[plt](at) = blitz::mean(snap2); 
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }


      else if (plt == "cl_rd_geq_0.8um_conc")
      {
        // gccn (r_d >= 0.8 um) concentration in cloudy grid cells
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          typename Plotter_t::arr_t snap2(plotter.h5load_timestep("rd_geq_0.8um_rw_mom0", at * n["outfreq"]));
          snap2 *= rhod; // b4 it was specific moment
          snap2 /= 1e6; // per cm^3
          snap2 *= snap;
          if(blitz::sum(snap) > 0)
            res_series[plt](at) = blitz::sum(snap2) / blitz::sum(snap); 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "rd_geq_0.8um_conc")
      {
        // gccn (r_d >= 0.8 um) concentration
        try
        {
          typename Plotter_t::arr_t snap2(plotter.h5load_timestep("rd_geq_0.8um_rw_mom0", at * n["outfreq"]));
          snap2 /= 1e6; // per cm^3
          snap2 *= rhod; // b4 it was per milligram
          res_series[plt](at) = blitz::mean(snap2); 
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }


      else if (plt == "cl_rd_lt_0.8um_conc")
      {
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          typename Plotter_t::arr_t snap2(plotter.h5load_timestep("rd_lt_0.8um_rw_mom0", at * n["outfreq"]));
          snap2 *= rhod; // b4 it was specific moment
          snap2 /= 1e6; // per cm^3
          snap2 *= snap;
          if(blitz::sum(snap) > 0)
            res_series[plt](at) = blitz::sum(snap2) / blitz::sum(snap); 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "rd_lt_0.8um_conc")
      {
        try
        {
          typename Plotter_t::arr_t snap2(plotter.h5load_timestep("rd_lt_0.8um_rw_mom0", at * n["outfreq"]));
          snap2 /= 1e6; // per cm^3
          snap2 *= rhod; // b4 it was per milligram
          res_series[plt](at) = blitz::mean(snap2); 
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else if (plt == "cl_meanr")
      {
        // cloud droplets mean radius in cloudy grid cells
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("actrw_rw_mom0", at * n["outfreq"]));
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          typename Plotter_t::arr_t snap_m0(plotter.h5load_timestep("actrw_rw_mom0", at * n["outfreq"]));
          typename Plotter_t::arr_t snap_m1(plotter.h5load_timestep("actrw_rw_mom1", at * n["outfreq"])*1e6); // in microns
          snap_m0 *= snap;
          snap_m1 *= snap;
          auto tot_m0 = blitz::sum(snap_m0);
          if(tot_m0 > 0)
            res_series[plt](at) = blitz::sum(snap_m1) / tot_m0; 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // cloud base mean incloud time of bigrain (r>40um)
      else if (plt == "clb_bigrain_mean_inclt")
      {
        try
        {
          // find cloud base (cloudy if q_c > 0.1 g/kg)
          typename Plotter_t::arr_t snap(plotter.h5load_rc_timestep(at * n["outfreq"]));
          snap = iscloudy_rc(snap); // cloudiness mask
          //for(int i=0;i<10;++i)
          //  snap(plotter.hrzntl_slice(i)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux
          plotter.k_i = blitz::first((snap == 1), plotter.LastIndex); 

          // 0-th specific mom of bigrain cloud drops
          typename Plotter_t::arr_t bigrain_conc(plotter.h5load_timestep("bigrain_rw_mom0", at * n["outfreq"]));
          // concentration of bigrain cloud drops at cloud base
          plotter.tmp_float_hrzntl_slice = plotter.get_value_at_hgt(bigrain_conc, plotter.k_i) * plotter.get_value_at_hgt(rhod, plotter.k_i); // we need to multiply by rhod here, because different cloud bases can mean different rhod

          // 1st specific mom of incloud time of bigrain drops
          typename Plotter_t::arr_t bigrain_inclt_mom1(plotter.h5load_timestep("bigrain_incloud_time_mom1", at * n["outfreq"]));
          // 1st mom of incloud time at cloud base
          plotter.tmp_float_hrzntl_slice2 = plotter.get_value_at_hgt(bigrain_inclt_mom1, plotter.k_i) * plotter.get_value_at_hgt(rhod, plotter.k_i);  // same as above

          if(blitz::sum(plotter.tmp_float_hrzntl_slice) > 0) // if any bigrain drops in the domain
            res_series[plt](at) = double(blitz::sum(plotter.tmp_float_hrzntl_slice2)) / blitz::sum(plotter.tmp_float_hrzntl_slice);
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // cloud base mean rd of bigrain (r>40um)
      else if (plt == "clb_bigrain_mean_rd")
      {
        try
        {
          // find cloud base (cloudy if q_c > 0.1 g/kg)
          typename Plotter_t::arr_t snap(plotter.h5load_rc_timestep(at * n["outfreq"]));
          snap = iscloudy_rc(snap); // cloudiness mask
          //for(int i=0;i<10;++i)
          //  snap(plotter.hrzntl_slice(i)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux
          plotter.k_i = blitz::first((snap == 1), plotter.LastIndex); 

          // 0-th specific mom of bigrain cloud drops
          typename Plotter_t::arr_t bigrain_conc(plotter.h5load_timestep("bigrain_rw_mom0", at * n["outfreq"]));
          // concentration of bigrain cloud drops at cloud base
          plotter.tmp_float_hrzntl_slice = plotter.get_value_at_hgt(bigrain_conc, plotter.k_i) * plotter.get_value_at_hgt(rhod, plotter.k_i); // we need to multiply by rhod here, because different cloud bases can mean different rhod

          // 1st specific mom of rd of bigrain drops
          typename Plotter_t::arr_t bigrain_rd_mom1(plotter.h5load_timestep("bigrain_rd_mom1", at * n["outfreq"]));
          // 1st mom of rd at cloud base
          plotter.tmp_float_hrzntl_slice2 = plotter.get_value_at_hgt(bigrain_rd_mom1, plotter.k_i) * plotter.get_value_at_hgt(rhod, plotter.k_i);  // same as above

          if(blitz::sum(plotter.tmp_float_hrzntl_slice) > 0) // if any bigrain drops in the domain
            res_series[plt](at) = double(blitz::sum(plotter.tmp_float_hrzntl_slice2)) / blitz::sum(plotter.tmp_float_hrzntl_slice);
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // cloud base mean rkappa of bigrain (r>40um)
      else if (plt == "clb_bigrain_mean_kappa")
      {
        try
        {
          // find cloud base (cloudy if q_c > 0.1 g/kg)
          typename Plotter_t::arr_t snap(plotter.h5load_rc_timestep(at * n["outfreq"]));
          snap = iscloudy_rc(snap); // cloudiness mask
          //for(int i=0;i<10;++i)
          //  snap(plotter.hrzntl_slice(i)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux
          plotter.k_i = blitz::first((snap == 1), plotter.LastIndex); 

          // 0-th specific mom of bigrain cloud drops
          typename Plotter_t::arr_t bigrain_conc(plotter.h5load_timestep("bigrain_rw_mom0", at * n["outfreq"]));
          // concentration of bigrain cloud drops at cloud base
          plotter.tmp_float_hrzntl_slice = plotter.get_value_at_hgt(bigrain_conc, plotter.k_i) * plotter.get_value_at_hgt(rhod, plotter.k_i); // we need to multiply by rhod here, because different cloud bases can mean different rhod

          // 1st specific mom of rd of bigrain drops
          typename Plotter_t::arr_t bigrain_kappa_mom1(plotter.h5load_timestep("bigrain_kappa_mom1", at * n["outfreq"]));
          // 1st mom of rd at cloud base
          plotter.tmp_float_hrzntl_slice2 = plotter.get_value_at_hgt(bigrain_kappa_mom1, plotter.k_i) * plotter.get_value_at_hgt(rhod, plotter.k_i);  // same as above

          if(blitz::sum(plotter.tmp_float_hrzntl_slice) > 0) // if any bigrain drops in the domain
            res_series[plt](at) = double(blitz::sum(plotter.tmp_float_hrzntl_slice2)) / blitz::sum(plotter.tmp_float_hrzntl_slice);
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // cloud base mean concentration of bigrain
      else if (plt == "clb_bigrain_mean_conc")
      {
        try
        {
          // find cloud base (cloudy if q_c > 0.1 g/kg)
          typename Plotter_t::arr_t snap(plotter.h5load_rc_timestep(at * n["outfreq"]));
          snap = iscloudy_rc(snap); // cloudiness mask
          //for(int i=0;i<10;++i)
          //  snap(plotter.hrzntl_slice(i)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux
          plotter.k_i = blitz::first((snap == 1), plotter.LastIndex); 

          // 0-th specific mom of bigrain cloud drops
          typename Plotter_t::arr_t bigrain_conc(plotter.h5load_timestep("bigrain_rw_mom0", at * n["outfreq"]));
          // concentration of bigrain cloud drops at cloud base [1/m^3]
          plotter.tmp_float_hrzntl_slice = plotter.get_value_at_hgt(bigrain_conc, plotter.k_i) * plotter.get_value_at_hgt(rhod, plotter.k_i); // we need to multiply by rhod here, because different cloud bases can mean different rhod

          // number of cloudy columns
          plotter.k_i = where(plotter.k_i > 0, 1, 0);
          if(blitz::sum(plotter.k_i) > 0)
            res_series[plt](at) = double(blitz::sum(plotter.tmp_float_hrzntl_slice)) / blitz::sum(plotter.k_i) / 1e6; // [1/cm^3]
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // cloud base mean fraction of bigrain formed on gccn
      else if (plt == "clb_bigrain_mean_gccn_fraction")
      {
        try
        {
          // find cloud base (cloudy if q_c > 0.1 g/kg)
          typename Plotter_t::arr_t snap(plotter.h5load_rc_timestep(at * n["outfreq"]));
          snap = iscloudy_rc(snap); // cloudiness mask
          //for(int i=0;i<10;++i)
          //  snap(plotter.hrzntl_slice(i)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux
          plotter.k_i = blitz::first((snap == 1), plotter.LastIndex); 

          // 0-th specific mom of bigrain cloud drops
          typename Plotter_t::arr_t bigrain_conc(plotter.h5load_timestep("bigrain_rw_mom0", at * n["outfreq"]));
          // concentration of bigrain cloud drops at cloud base [1/m^3]
          plotter.tmp_float_hrzntl_slice = plotter.get_value_at_hgt(bigrain_conc, plotter.k_i) * plotter.get_value_at_hgt(rhod, plotter.k_i); // we need to multiply by rhod here, because different cloud bases can mean different rhod

          // 0-th specific mom of bigrain cloud drops formed on gccn
          typename Plotter_t::arr_t bigrain_gccn_conc(plotter.h5load_timestep("bigrain_gccn_rw_mom0", at * n["outfreq"]));
          // concentration of bigrain cloud drops at cloud base [1/m^3]
          plotter.tmp_float_hrzntl_slice2 = plotter.get_value_at_hgt(bigrain_gccn_conc, plotter.k_i) * plotter.get_value_at_hgt(rhod, plotter.k_i); // we need to multiply by rhod here, because different cloud bases can mean different rhod

          if(blitz::sum(plotter.tmp_float_hrzntl_slice) > 0) // if any bigrain drops in the domain
            res_series[plt](at) = double(blitz::sum(plotter.tmp_float_hrzntl_slice2)) / blitz::sum(plotter.tmp_float_hrzntl_slice);
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      // ------ plots specific to the Pi Chamber ICMW case, averaged over the domain with the exception of near-wall cells ------

      else if (plt == "Qv")
      {
        try
        {
          auto stats = plotter.rv_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.first;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else if (plt == "T")
      {
        try
        {
          auto stats = plotter.T_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.first;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else if (plt == "S_drop")
      {
        try
        {
          auto stats = plotter.drop_supersat_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.first;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else if (plt == "LWC")
      {
        try
        {
          typename Plotter_t::arr_t lwc(plotter.h5load_rc_timestep(at * n["outfreq"])); // cloud water, no rain water in pi chamber icmw
          //res_series[plt](at) = blitz::mean(typename Plotter_t::arr_t(plotter.nowall(lwc, distance_from_walls)));
          res_series[plt](at) = blitz::mean(lwc);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else if (plt == "LWC_gm-3")
      {
        try
        {
          typename Plotter_t::arr_t lwc(plotter.h5load_rc_timestep(at * n["outfreq"])); // cloud water, no rain water in pi chamber icmw
          lwc *= rhod;
          res_series[plt](at) = blitz::mean(lwc);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else if (plt == "IWC_gm-3")
      {
        try
        {
          typename Plotter_t::arr_t iwc(plotter.h5load_ri_timestep(at * n["outfreq"])); // cloud water, no rain water in pi chamber icmw
          iwc *= rhod;
          res_series[plt](at) = blitz::mean(iwc);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else if (plt == "RH")
      {
        try
        {
          auto stats = plotter.RH_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.first;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else if (plt == "N_drop")
      {
        try
        {
          typename Plotter_t::arr_t nc(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          nc *= rhod; // 1/kg -> 1/m^3
          res_series[plt](at) = blitz::mean(nc);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else if (plt == "N_ice")
      {
        try
        {
          typename Plotter_t::arr_t ni(plotter.h5load_timestep("ice_rw_mom0", at * n["outfreq"]));
          ni *= rhod; // 1/kg -> 1/m^3
          res_series[plt](at) = blitz::mean(ni);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else if (plt == "N_aerosol")
      {
        try
        {
          typename Plotter_t::arr_t na(plotter.h5load_timestep("aerosol_rw_mom0", at * n["outfreq"]));
          na *= rhod; // 1/kg -> 1/m^3
          res_series[plt](at) = blitz::mean(na);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else if (plt == "res_tke") // resolved tke
      {
        try
        {
          typename Plotter_t::arr_t u(plotter.h5load_timestep("u", at * n["outfreq"]));
          plotter.subtract_horizontal_mean(u);
          u = u * u;
          res_series[plt](at) = blitz::mean(u);

          typename Plotter_t::arr_t w(plotter.h5load_timestep("w", at * n["outfreq"]));
          plotter.subtract_horizontal_mean(w);
          w = w * w;
          res_series[plt](at) += blitz::mean(w);

          if (Plotter_t::n_dims > 2)
          {
            typename Plotter_t::arr_t v(plotter.h5load_timestep("v", at * n["outfreq"]));
            plotter.subtract_horizontal_mean(v);
            v = v * v;
            res_series[plt](at) += blitz::mean(v);
          }
          
          res_series[plt](at) *= 0.5; // * n["dz"];
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "tot_tke")
      {
        try
        {
          typename Plotter_t::arr_t u(plotter.h5load_timestep("u", at * n["outfreq"]));
          plotter.subtract_horizontal_mean(u);
          u = u * u;
          res_series[plt](at) = blitz::mean(u);

          typename Plotter_t::arr_t w(plotter.h5load_timestep("w", at * n["outfreq"]));
          plotter.subtract_horizontal_mean(w);
          w = w * w;
          res_series[plt](at) += blitz::mean(w);

          if (Plotter_t::n_dims > 2)
          {
            typename Plotter_t::arr_t v(plotter.h5load_timestep("v", at * n["outfreq"]));
            plotter.subtract_horizontal_mean(v);
            v = v * v;
            res_series[plt](at) += blitz::mean(v);
          }
          
          res_series[plt](at) *= 0.5; // * n["dz"];

          typename Plotter_t::arr_t tke(plotter.h5load_timestep("tke", at * n["outfreq"]));
          typename Plotter_t::arr_t snap;
          snap.reference(tke);

          res_series[plt](at) += blitz::mean(plotter.horizontal_mean(snap));
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "r_dmean")
      {
        // droplets mean radius away from walls
        try
        {
          typename Plotter_t::arr_t m0(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          typename Plotter_t::arr_t m1(plotter.h5load_timestep("cloud_rw_mom1", at * n["outfreq"]));
          //auto tot_m0 = blitz::sum(typename Plotter_t::arr_t(plotter.nowall(m0, distance_from_walls)));
          auto tot_m0 = blitz::sum(m0);
          if(tot_m0 > 0)
            res_series[plt](at) = blitz::sum(m1) / tot_m0; 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "r_imean")
      {
        // ice mean radius away from walls
        try
        {
          typename Plotter_t::arr_t m0(plotter.h5load_timestep("ice_rw_mom0", at * n["outfreq"]));
          typename Plotter_t::arr_t m1(plotter.h5load_timestep("ice_rw_mom1", at * n["outfreq"]));
          //auto tot_m0 = blitz::sum(typename Plotter_t::arr_t(plotter.nowall(m0, distance_from_walls)));
          auto tot_m0 = blitz::sum(m0);
          if(tot_m0 > 0)
            res_series[plt](at) = blitz::sum(m1) / tot_m0; 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "r_mean2")
      {
        // droplets effective radius away from walls
        try
        {
          typename Plotter_t::arr_t m2(plotter.h5load_timestep("cloud_rw_mom2", at * n["outfreq"]));
          typename Plotter_t::arr_t m3(plotter.h5load_timestep("cloud_rw_mom3", at * n["outfreq"]));
          auto tot_m2 = blitz::sum(m2);
          if(tot_m2 > 0)
            res_series[plt](at) = blitz::sum(m3) / tot_m2; 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // spatial variance of supersaturation  [1]
      else if (plt == "Sigma2_S")
      {
        try
        {
          auto stats = plotter.RH_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.second * stats.second; // std_dev -> variance; sigma(RH) = sigma(S)
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // spatial variance of supersaturation weighted by droplet number [1]
      else if (plt == "Sigma2_S_drop")
      {
        try
        {
          auto stats = plotter.drop_supersat_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.second * stats.second; // std_dev -> variance
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // spatial variance of T [K^2]
      else if (plt == "Sigma2_T")
      {
        try
        {
          auto stats = plotter.T_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.second * stats.second; // std_dev -> variance
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      // spatial variance of rv  [1]
      else if (plt == "Sigma2_Qv")
      {
        try
        {
          auto stats = plotter.rv_stats_timestep(at * n["outfreq"]);
          res_series[plt](at) = stats.second * stats.second; // std_dev -> variance
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "disp_r")
      {
        // relative dispersion (std dev / mean) of droplet radius distribution averaged over cells away from walls 
        try
        {
          //typename Plotter_t::arr_t m0(plotter.nowall(typename Plotter_t::arr_t(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"])), distance_from_walls));
          typename Plotter_t::arr_t m0(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          typename Plotter_t::arr_t m1(plotter.h5load_timestep("cloud_rw_mom1", at * n["outfreq"]));
          typename Plotter_t::arr_t m2(plotter.h5load_timestep("cloud_rw_mom2", at * n["outfreq"]));
          // calculate stddev of radius, store in m2
          m2 = where(m0 > 0,
            m2 / m0 - m1 / m0 * m1 / m0, 0.);
          // might be slightly negative due to numerical errors
          m2 = where(m2 < 0, 0, m2);
          m2 = sqrt(m2); // sqrt(variance)

          // calculate mean radius, store in m1
          m1 = where(m0 > 0,
            m1 / m0, 0.);

          auto tot_m1 = blitz::sum(m1);
          if(tot_m1 > 0)
            res_series[plt](at) = blitz::sum(m2) / tot_m1; 
          else
            res_series[plt](at) = 0;
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "epsilon")
      {
        // SGS TKE dissipation rate away from walls [m^2/s^3]
        try
        {
          const float C_E = 0.845;
          typename Plotter_t::arr_t tke(plotter.h5load_timestep("tke", at * n["outfreq"]));
          tke = pow(tke, 3./2.); 
          tke /= n_prof["mix_len"](plotter.LastIndex); // divide by SGS mixing length
//          res_series[plt](at) = blitz::mean(plotter.nowall(tke, distance_from_walls)) * C_E; 
          res_series[plt](at) = blitz::mean(tke) * C_E; 
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "N_drop_removal")
      {
        // droplet removal rate [1/(cm^3 * min)]
        try
        {
          res_series[plt](at) = plotter.calc_prtcl_removal(removed_droplets - removed_droplets_prev);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "N_ice_removal")
      {
        // ice removal rate [1/(cm^3 * min)]
        try
        {
          res_series[plt](at) = plotter.calc_prtcl_removal(removed_ice - removed_ice_prev);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "H_flux_t")
      {
        // heat flux at the top [W/m2]
        try
        {
          res_series[plt](at) = plotter.calc_heat_flux_top(th_change_top, at>0);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "H_flux_b")
      {
        // heat flux at the bot [W/m2]
        try
        {
          res_series[plt](at) = plotter.calc_heat_flux_bot(th_change_bot, at>0);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "qv_flux_t")
      {
        // moisture flux at the top [kg/kg * m/s]
        try
        {
          res_series[plt](at) = plotter.calc_moist_flux_top(rv_change_top, at>0);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }
      else if (plt == "qv_flux_b")
      {
        // moisture flux at the bot [kg/kg * m/s]
        try
        {
          res_series[plt](at) = plotter.calc_moist_flux_bot(rv_change_bot, at>0);
        }
        catch(...) {if(at==first_timestep) data_found[plt]=0;}
      }

      else assert(false);
    } // var loop
  } // ------- end of time loop ------

  // processing done after reading whole time series
  for (auto &plt : plots.series)
  {
    std::cerr << "post processing: " << plt << std::endl;
    // if no data was found, skip to the next var, dont save the data=0 as it was confusing
    if(!data_found[plt])
      continue;

    bool plot_std_dev = 0;

    if (plt == "ract_com")
    {
      res_series[plt] /= 1000.;
      res_pos *= 60.;
    }
    else if (plt == "th_com")
    {
      res_series[plt] /= 1000.;
      res_pos *= 60.;
    }
    else if (plt == "ract_avg")
    {
      plot_std_dev = true;
      res_pos *= 60.;
    }
    else if (plt == "ract_std_dev")
    {
      res_pos *= 60.;
    }
    else if (plt == "cloud_avg_act_conc")
    {
      plot_std_dev = true;
      res_pos *= 60.;
    }
    else if (plt == "cloud_std_dev_act_conc")
    {
      res_pos *= 60.;
    }
    else if (plt == "cloud_avg_supersat")
    {
      plot_std_dev = true;
      res_pos *= 60.;
    }
    else if (plt == "cloud_std_dev_supersat")
    {
      res_pos *= 60.;
    }
    else if (plt == "cloud_avg_act_rad")
    {
      res_pos *= 60.;
    }
    else if (plt == "cloud_avg_std_dev_act_rad")
    {
      res_pos *= 60.;
    }
    else if (plt == "sd_conc")
    {
      plot_std_dev = true;
      res_pos *= 60.;
    }
    else if (plt == "cl_sd_conc")
    {
      plot_std_dev = true;
      res_pos *= 60.;
    }
    else if (plt == "cl_sd_conc_act")
    {
      plot_std_dev = true;
      res_pos *= 60.;
    }
    else if (plt == "tot_water")
    {
      res_pos *= 60.;
      res_series[plt] *= 1e3;
    }
    else if (plt == "com_vel")
    {
      res_pos *= 60.;
    }
    else if (plt == "com_supersat")
    {
      res_pos *= 60.;
      res_series[plt] *= 100.; // to get %
    }
    else if (plt == "com_mom0")
    {
      res_pos *= 60.;
      res_series[plt] /= 1e6;
    }
    else if (plt == "com_mom1")
    {
      res_pos *= 60.;
      res_series[plt] *= 1e6;
    }
    else if (plt == "com_mom2")
    {
      res_pos *= 60.;
      res_series[plt] *= 1e6;
    }
    else if (plt == "com_sd_conc")
    {
      res_pos *= 60.;
    }
    else if (plt == "RH_max")
    {
      res_pos *= 60.;
    }
    else if (plt == "cloud_top_height")
    {
      res_pos *= 60.;
    }
    else if (plt == "tot_tke" || plt == "tot_tke" || plt == "sgs_tke" || plt == "uw_resolved_tke")
    {
      res_pos *= 60.;
    }
    else if (plt == "lwp")
    {
      res_series[plt] *= (n["z"] - 1) * n["dz"]; // top and bottom cells are smaller
    }
    else if (plt == "rwp")
    {
      res_series[plt] *= (n["z"] - 1) * n["dz"]; // top and bottom cells are smaller
    }
    else if (plt == "cwp")
    {
      res_series[plt] *= (n["z"] - 1) * n["dz"]; // top and bottom cells are smaller
    }
    else if (plt == "er")
    {
      // central difference, in cm
      Range nofirstlast = Range(1, last_timestep-1);
      auto res_series_tmp = res_series[plt].copy();
      res_series[plt](nofirstlast) = where(res_series_tmp(nofirstlast+1) > 0., (res_series_tmp(nofirstlast+1) - res_series_tmp(nofirstlast-1)) * n["dz"] * 1e2 / (2 * n["dt"] * n["outfreq"])  + D * (res_series_tmp(nofirstlast) - 0.5) * n["dz"] * 1e2, 0.);

      // larger stencil
//      Range notwo = Range(2, last_timestep-2);
   //   res_series[plt](notwo) = where(res_series_tmp(notwo+1) > 0., ( 2. / 3. * (res_series_tmp(notwo+1) - res_series_tmp(notwo-1)) + 1. / 12. * (res_series_tmp(notwo+2) - res_series_tmp(notwo-2)) ) * n["dz"] * 1e2 / (n["dt"] * n["outfreq"])  + D * (res_series_tmp(notwo) - 0.5) * n["dz"] * 1e2, 0.);

      //res_series[plt](0) = 0.;
      res_series[plt](0) = (res_series_tmp(1) - res_series_tmp(0)) * n["dz"] * 1e2 / (n["dt"] * n["outfreq"])  + D * (res_series_tmp(0) - 0.5) * n["dz"] * 1e2;
      res_series[plt](last_timestep) = (res_series_tmp(last_timestep) - res_series_tmp(last_timestep-1)) * n["dz"] * 1e2 / (n["dt"] * n["outfreq"])  + D * (res_series_tmp(last_timestep) - 0.5) * n["dz"] * 1e2;
    }
    else if (plt == "Qv")
    {
      res_pos *= 3600.;
    }
    else if (plt == "LWC")
    {
      res_pos *= 3600.;
      res_series[plt] *= 1e3; // g/kg
    }
    else if (plt == "N_drop_removal")
    {
      res_pos *= 3600.;
      res_series[plt] *= 60; // per minute
    }
    else if (plt == "N_ice_removal")
    {
      res_pos *= 3600.;
      res_series[plt] *= 60; // per minute
    }
    else if (plt == "LWC_gm-3")
    {
      res_pos *= 3600.;
      res_series[plt] *= 1e3; // g/m^3
    }
    else if (plt == "IWC_gm-3")
    {
      res_pos *= 3600.;
      res_series[plt] *= 1e3; // g/m^3
    }
    else if (plt == "T")
    {
      res_pos *= 3600.;
    }
    else if (plt == "S_drop")
    {
      res_pos *= 3600.;
    }
    else if (plt == "RH")
    {
      res_pos *= 3600.;
    }
    else if (plt == "N_drop")
    {
      res_pos *= 3600.;
      res_series[plt] /= 1e6; // 1/m^3 -> 1/cm^3
    }
    else if (plt == "N_ice")
    {
      res_pos *= 3600.;
      res_series[plt] /= 1e6; // 1/m^3 -> 1/cm^3
    }
    else if (plt == "N_aerosol")
    {
      res_pos *= 3600.;
      res_series[plt] /= 1e6; // 1/m^3 -> 1/cm^3
    }
    else if (plt == "r_dmean")
    {
      res_pos *= 3600.;
      res_series[plt] *= 1e6; // m -> um
    }
    else if (plt == "r_imean")
    {
      res_pos *= 3600.;
      res_series[plt] *= 1e6; // m -> um
    }
    else if (plt == "r_mean2")
    {
      res_pos *= 3600.;
      res_series[plt] *= 1e6; // m -> um
    }
    else if (plt == "disp_r")
    {
      res_pos *= 3600.;
    }
    else if (plt == "Sigma2_Qv")
    {
      res_pos *= 3600.;
    }
    else if (plt == "Sigma2_S_drop")
    {
      res_pos *= 3600.;
    }

    // set labels for the gnuplot plot
    gnuplot_series_set_labels(gp, plt);

    gp << "plot '-' with l";
    if(plot_std_dev)
      gp << ", '-' w l, '-' w l";
    else if(plt == "cl_acnv25_dycoms" || plt == "cl_acnv25_rico" || plt == "cl_accr25_dycoms" || plt == "cl_accr25_rico")
      gp << ", '-' w l";
    gp << " \n";

    std::cout << plt << " " << res_pos << res_series[plt] << res_series_std_dev[plt] << std::endl;
    gp.send1d(boost::make_tuple(res_pos, res_series[plt]));
    oprof_file << plt << endl ;
    oprof_file << res_series[plt] ;
    if(plot_std_dev)
    {
      oprof_file << res_series_std_dev[plt] ;
      res_series[plt] = res_series[plt] + res_series_std_dev[plt];
      gp.send1d(boost::make_tuple(res_pos, res_series[plt]));
      res_series[plt] = res_series[plt] - 2*res_series_std_dev[plt];
      gp.send1d(boost::make_tuple(res_pos, res_series[plt]));
    }
    if(plt == "cl_acnv25_dycoms" || plt == "cl_acnv25_rico" || plt == "cl_accr25_dycoms" || plt == "cl_accr25_rico")
    {
      // acnv/accr rate averaged since the start of the simulation
      int nt = last_timestep - first_timestep + 1;
      Array<double, 1> res_series_acc_sum(nt);
      res_series_acc_sum = 0;
      for(int t=1; t<nt; ++t)
        res_series_acc_sum(t) = (res_series_acc_sum(t-1) + res_series[plt](t));
      for(int t=1; t<nt; ++t)
        res_series_acc_sum(t) /= double(t);

      gp.send1d(boost::make_tuple(res_pos, res_series_acc_sum));
      oprof_file << "acc_" << plt << endl ;
      oprof_file << res_series_acc_sum ;
    }
   // plot(gp, res);
  } // var loop
}
