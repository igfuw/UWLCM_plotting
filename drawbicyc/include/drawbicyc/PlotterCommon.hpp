#pragma once

#include <blitz/array.h>
#include <blitz/array/stencil-et.h>
#include <H5Cpp.h>
#include <map>

class PlotterCommon
{
  public:
  using arr_prof_t = blitz::Array<float,1>;

  const string file;
  std::map<std::string, double> map;
  std::map<std::string, arr_prof_t> map_prof;
  blitz::Array<float, 1> timesteps;
  double CellVol, DomainSurf, DomainVol;

  protected:
  H5::H5File h5f;
  H5::DataSet h5d;
  H5::Group h5g;
  H5::DataSpace h5s;

  void h5load(
    const string &file,
    const string &dataset,
    bool srfc = false
  )
  {
    notice_macro("about to close current file")
    h5f.close();

    notice_macro("about to open file: " << file)
    h5f.openFile(file, H5F_ACC_RDONLY);

    notice_macro("about to read dataset: " << dataset)
    h5d = h5f.openDataSet(dataset);
    h5s = h5d.getSpace();
  }

  template <class gp_t>
  void plot(gp_t &gp)
  {
    //gp << "set cbtics format \"%.2tE%+03T\"\n";
    gp << "set cbtics font \", 8\"\n";
  //  gp << "set rmargin 2cm\n";
  }

  public:

  float puddle_liq_vol(int at)
  {
    notice_macro("about to close current file")
    h5f.close();
  
    string timestep_file = file + "/timestep" + zeropad(at, 10) + ".h5";
    notice_macro("about to open file: " << timestep_file)
    h5f.openFile(timestep_file, H5F_ACC_RDONLY);
  
    notice_macro("about to read group: puddle")
    h5g = h5f.openGroup("puddle");

    float ret;
    auto attr = h5g.openAttribute("liquid_volume");
    attr.read(attr.getDataType(), &ret);
    return ret;
  }

  float puddle_prtcl_no(int at)
  {
    notice_macro("about to close current file")
    h5f.close();
  
    string timestep_file = file + "/timestep" + zeropad(at, 10) + ".h5";
    notice_macro("about to open file: " << timestep_file)
    h5f.openFile(timestep_file, H5F_ACC_RDONLY);
  
    notice_macro("about to read group: puddle")
    h5g = h5f.openGroup("puddle");

    float ret;
    auto attr = h5g.openAttribute("particle_number");
    attr.read(attr.getDataType(), &ret);
    return ret;
  }

  float acnv25(int at)
  {
    notice_macro("about to close current file")
    h5f.close();
  
    string timestep_file = file + "/timestep" + zeropad(at, 10) + ".h5";
    notice_macro("about to open file: " << timestep_file)
    h5f.openFile(timestep_file, H5F_ACC_RDONLY);
  
    notice_macro("about to read group: /")
    h5g = h5f.openGroup("/");

    float ret;
    auto attr = h5g.openAttribute("acnv25");
    attr.read(attr.getDataType(), &ret);
    return ret;
  }

  //ctor
  PlotterCommon(const string &file):
    file(file)
  {
    // init dt and outfreq
    {
      notice_macro("about to open file: " << file << "/const.h5")
      h5f.openFile(file + "/const.h5", H5F_ACC_RDONLY);

      notice_macro("about to read group: advection")
      h5g = h5f.openGroup("advection");

      float dt;
      {
        auto attr = h5g.openAttribute("dt");
        attr.read(attr.getDataType(), &dt);
      }
      map["dt"] = dt;

      // read number of timesteps
      hsize_t n;
      h5load(file + "/const.h5", "T");
      h5s.getSimpleExtentDims(&n, NULL);
      this->map["t"] = n;
      // read timesteps
      timesteps.resize(n);
      h5d.read(timesteps.data(), H5::PredType::NATIVE_FLOAT);

      // read environmental pressure profile
      h5load(file + "/const.h5", "p_e");
      h5s.getSimpleExtentDims(&n, NULL);
      map_prof.emplace("p_e", arr_prof_t(n));
      h5d.read(map_prof["p_e"].data(), H5::PredType::NATIVE_FLOAT);

      // read SGS mixing length profile
      h5load(file + "/const.h5", "mix_len");
      h5s.getSimpleExtentDims(&n, NULL);
      map_prof.emplace("mix_len", arr_prof_t(n));
      h5d.read(map_prof["mix_len"].data(), H5::PredType::NATIVE_FLOAT);

      // read output frequency
      float outfreq;
      {
        auto root_group = h5f.openGroup("/");
        auto attr = root_group.openAttribute("user_params outfreq");
        attr.read(attr.getDataType(), &outfreq);
      }
      map["outfreq"] = outfreq;
    }
  }
};

