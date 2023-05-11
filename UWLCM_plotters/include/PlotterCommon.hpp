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
  double CellVol, DomainSurf, DomainVol, CellVol_ref;

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
    if(h5f.getFileName() != file)
    {
      notice_macro("about to close current file: " << h5f.getFileName())
      h5f.close();

      notice_macro("about to open file: " << file)
      h5f.openFile(file, H5F_ACC_RDONLY);
    }

    notice_macro("about to read dataset: " << dataset)
    h5d = h5f.openDataSet(dataset);
    h5s = h5d.getSpace();
  }

  float h5load_attr(const string &file, const string &attr_name, const string &group_name)
  {
    if(h5f.getFileName() != file)
    {
      notice_macro("about to close current file: " << h5f.getFileName())
      h5f.close();
  
      notice_macro("about to open file: " << file)
      h5f.openFile(file, H5F_ACC_RDONLY);
    }
  
    notice_macro(std::string("about to read group: " + group_name))
    h5g = h5f.openGroup(group_name);

    float ret;
    notice_macro(std::string("about to open attribute: " + attr_name))
    auto attr = h5g.openAttribute(attr_name);
    notice_macro(std::string("about to read attribute value"))
    attr.read(attr.getDataType(), &ret);
    notice_macro(std::string("attribute value read: ") + std::to_string(ret))
    return ret;
  }

  template <class gp_t>
  void plot(gp_t &gp)
  {
    //gp << "set cbtics format \"%.2tE%+03T\"\n";
    gp << "set cbtics font \", 8\"\n";
  //  gp << "set rmargin 2cm\n";
  }

  public:

  float h5load_attr_timestep(int at, const std::string attr_name, const std::string group_name = "/")
  {
    string timestep_file = file + "/timestep" + zeropad(at, 10) + ".h5";
    return h5load_attr(timestep_file, attr_name, group_name);
  }

  //ctor
  PlotterCommon(const string &file):
    file(file)
  {
    // init h5f
    notice_macro("about to open file: " << file << "/const.h5")
    h5f.openFile(file + "/const.h5", H5F_ACC_RDONLY);

    // init dt and outfreq
    {
      map["dt"] = h5load_attr(file + "/const.h5", "dt", "advection");
      map["outfreq"] = h5load_attr(file + "/const.h5", "outfreq", "user_params");
      map["MPI_compiler"] = h5load_attr(file + "/const.h5", "MPI compiler (true/false)", "MPI details");

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

      // read dry air density profile
      h5load(file + "/const.h5", "rhod");
      h5s.getSimpleExtentDims(&n, NULL);
      map_prof.emplace("rhod", arr_prof_t(n));
      h5d.read(map_prof["rhod"].data(), H5::PredType::NATIVE_FLOAT);

      // read dry air density profile on refined grid
      try
      {
        h5load(file + "/const.h5", "refined rhod");
        h5s.getSimpleExtentDims(&n, NULL);
        map_prof.emplace("refined rhod", arr_prof_t(n));
        h5d.read(map_prof["refined rhod"].data(), H5::PredType::NATIVE_FLOAT);
      }
      catch(...) // for pre-refinement simulations, use rhod as refined rhod
      {
        map_prof.emplace("refined rhod", map_prof["rhod"].copy());
	std::cerr << "refined rhod as copy of rhod: " << map_prof["refined rhod"];
      }
    }
  }
};

