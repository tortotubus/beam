#include "elff/io/CXX/vtkHDFPolyData.hpp"
#include "elff/io/CXX/vtkPolyData.hpp"

#include "elff/io/C/vtkHDFPolyData.hpp"

namespace ELFF {
namespace IO {
namespace CXX {

vtkHDFPolyData::vtkHDFPolyData(std::string filename)
  : vtk_polydata(C::vtk_polydata_init(0, 0, 0, 0, 0))
  , filename(filename)
{
}

vtkHDFPolyData::vtkHDFPolyData(std::string filename,
                               C::vtkPolyData c_vtk_polydata)
  : vtk_polydata(c_vtk_polydata)
  , filename(filename)
{
}

vtkHDFPolyData::vtkHDFPolyData(std::string filename,
                               CXX::vtkPolyData& cxx_vtk_polydata)
  : vtk_polydata(cxx_vtk_polydata.to_c_struct())
  , filename(filename)
{
}

vtkHDFPolyData::~vtkHDFPolyData()
{
  C::vtk_polydata_free(&vtk_polydata);
}

void
vtkHDFPolyData::write_new_static(bool overwrite)
{
  C::vtkHDFPolyData vtk_hdf_pd =
    C::vtk_HDF_polydata_init_static(filename.c_str(), overwrite, &vtk_polydata);

  C::vtk_HDF_polydata_close(&vtk_hdf_pd);
}

void
vtkHDFPolyData::write_new_static(bool overwrite, int compression_level)
{
  (void)compression_level;
  write_new_static(overwrite);
}

void
vtkHDFPolyData::write_new_transient(bool overwrite, float time)
{
  C::vtkHDFPolyData vtk_hdf_pd =
    C::vtk_HDF_polydata_init_transient(filename.c_str(), overwrite, &vtk_polydata, time);

  C::vtk_HDF_polydata_close(&vtk_hdf_pd);
}

void
vtkHDFPolyData::write_new_transient(bool overwrite,
                                    float time,
                                    int compression_level)
{
  (void)compression_level;
  write_new_transient(overwrite, time);
}


void
vtkHDFPolyData::append_transient(float time)
{
  C::vtkHDFPolyData vtk_hdf_pd =
    C::vtk_HDF_polydata_append_transient(filename.c_str(), &vtk_polydata, time);

  C::vtk_HDF_polydata_close(&vtk_hdf_pd);
}

void
vtkHDFPolyData::append_transient(float time, int compression_level)
{
  (void)compression_level;
  append_transient(time);
}

}
}
}
