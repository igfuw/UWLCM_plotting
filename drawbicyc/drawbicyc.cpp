#include "include/plot_series.hpp"
#include "include/plot_prof.hpp"
#include "include/plot_fields.hpp"
#include "include/plot_qv_qc_2_6_10_min.hpp"
#include "include/plot_lgrngn_spec.hpp"

int main(int argc, char** argv)
{
  // make args global
  ac=argc;
  av=argv;

  // general opts
  opts_main.add_options()
    ("profs", po::value<bool>()->default_value(true), "plot profiles?")
    ("series", po::value<bool>()->default_value(true) , "plot series?")
    ("fields", po::value<bool>()->default_value(false) , "plot fields?")
    ("spectra", po::value<bool>()->default_value(false) , "plot spectra?")
    ("qv_qc_2_6_10_min", po::value<bool>()->default_value(false) , "plot comparison of qv and qc fields at 2, 6 and 10 min?")
    ("dir", po::value<std::string>()->required() , "directory containing out_lgrngn")
    ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn")
    ("type", po::value<std::string>()->required(), "one of: dycoms, moist_thermal, rico, cumulus_congestus, gccn_ccn_conc")//, base_prflux_vs_clhght")
  ;

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).options(opts_main).allow_unregistered().run(), vm); //     ignores unknown

  // checking if all required options present
  po::notify(vm);

  // handling the "micro" option
  std::string micro = vm["micro"].as<std::string>();

  // handling the "type" option
  std::string type = vm["type"].as<std::string>();
  if(type != "dycoms" && type != "moist_thermal" && type != "rico" && type != "pi_chamber" && type != "cumulus_congestus" && type != "pi_chamber_icmw" && type != "gccn_ccn_conc")
    throw std::runtime_error("Unrecognized 'type' option, only dycoms, rico, moist_thermal, pi_chamber, pi_chamber_icmw, cumulus_congestus, gccn_ccn_conc  available now");//, base_prflux_vs_clhght available now");

  // should profiles be normalized by inversion height
  const bool normalize_prof = type == "dycoms";

  // parse dir name
  std::string
    dir = vm["dir"].as<std::string>(),
    h5  = dir;// + "out_" + micro;

  // reading required plot types
  bool flag_series = vm["series"].as<bool>(),
       flag_profiles = vm["profs"].as<bool>(),
       flag_fields = vm["fields"].as<bool>(),
       flag_lgrngn_spec = vm["spectra"].as<bool>(),
       flag_qv_qc_2_6_10_min = vm["qv_qc_2_6_10_min"].as<bool>();

  // detecting input data dimensionality
  H5::H5File h5f(h5 + "/const.h5", H5F_ACC_RDONLY);
  H5::DataSet h5d = h5f.openDataSet("G");
  H5::DataSpace h5s = h5d.getSpace();
  int NDims = h5s.getSimpleExtentNdims();

  // selecting type of cloud mask
  mask_type_t mask_type = mask_type_t::unset;
  if(type == "dycoms")
  {
    std::cout << "Using Dycoms_rf02 cloud mask." << std::endl;
    mask_type = mask_type_t::Dycoms_rf02;
  }
  else
  {
    std::cout << "Using Rico11 cloud mask." << std::endl;
    mask_type = mask_type_t::Rico11;
  }
  
  // detecting if subgrid model was on
//  std::cerr << "checking of sgs model was used... ";
  bool sgs = h5f.nameExists("sgs");// true;
  /*
  try 
  {
    auto h5g = h5f.openGroup("sgs");
  }
  catch (...)
  {
    sgs = false;
  }
  */
 // std::cerr << sgs << std::endl;

  Plots plots(type, sgs);

  if(NDims == 2)
  {

    if(flag_series)   plot_series(PlotterMask<2>(h5, micro, mask_type_t::Rico11), plots, type);
    if(flag_profiles) plot_profiles(PlotterMask<2>(h5, micro, mask_type_t::Rico11), plots, type, normalize_prof);
//    if(flag_fields)   plot_fields(PlotterMask<2>(h5, micro), plots, type);
//    if(flag_qv_qc_2_6_10_min)   plot_qv_qc_2_6_10_min(PlotterMask<2>(h5, micro));
  }
  else if(NDims == 3)
  {
    if(flag_series)   plot_series(PlotterMask<3>(h5, micro, mask_type), plots, type);
    if(flag_profiles) plot_profiles(PlotterMask<3>(h5, micro, mask_type), plots, type, normalize_prof);
//    if(flag_fields)   plot_fields(PlotterMask<3>(h5, micro), plots, type);
//    if(flag_qv_qc_2_6_10_min)   plot_qv_qc_2_6_10_min(PlotterMask<2>(h5, micro));
//    if(flag_lgrngn_spec) {
//      plot_lgrngn_spec_positions(PlotterMask<3>(h5, "lgrngn"));
//      plot_lgrngn_spec(PlotterMask<3>(h5, "lgrngn"));
//    }
  }
  else
    assert(false && "need 2d or 3d input data");

return 0;
} // main
