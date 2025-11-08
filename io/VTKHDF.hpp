#pragma once

#include "config/config.hpp"
#include "general/error.hpp"

#include <H5Cpp.h>
#include <string>
#include <cstring>

namespace beam {
typedef enum
{
  IMAGEDATA,
  POLYDATA,
  UNSTRUCTUREDGRID,
  HYPERTREEGRID,
  OVERLAPPINGAMR,
  PARTITIONEDDATASETCOLLECTION,
  MULTIBLOCKDATASET
} VTKHDFType;

typedef struct
{
  int major;
  int minor;
} VTKHDFVersion;

class VTKHDF
{
protected:
  std::string filename;
  VTKHDFType type;
  VTKHDFVersion version;

  H5::H5File file;
  H5::Group grp_vtkhdf;


  std::string type_to_string(VTKHDFType type) {
    switch (type) {
      case POLYDATA: return "PolyData";
      case IMAGEDATA: return "ImageData";
      case UNSTRUCTUREDGRID: return "UnstructuredGrid";
      case OVERLAPPINGAMR: return "OverlappingAMR";
      case PARTITIONEDDATASETCOLLECTION: return "PartitionedDataSetCollection";
      case MULTIBLOCKDATASET: return "MultiBlockDataSet";
      default: BEAM_ABORT("VTKHDF Type not known.\n");
    }
  }

  static H5::Group ensure_group(H5::H5File &file, std::string path) {
    if (H5Lexists(file.getId(), path.c_str(), H5P_DEFAULT) > 0) 
      return file.openGroup(path);

    return file.createGroup(path);
  }

  static void ensure_attr_fixed_ascii_scalar(H5::Group &group, std::string name, std::string value) {
    if (H5Aexists(group.getId(), name.c_str()) > 0) { 
      // Attribute exists
    } else { 
      // Attribute does not exist
      H5::StrType datatype (H5::PredType::C_S1, value.size());
      datatype.setStrpad(H5T_STR_NULLPAD);
      datatype.setCset(H5T_CSET_ASCII);
      
      H5::DataSpace space(H5S_SCALAR);

      auto attribute = group.createAttribute(name, datatype, space);
      attribute.write(datatype, value);
      return;
    }
  }

  void ensure_common_schema() {
    // Create/overwrite file
    try { 
      file = H5::H5File(filename, H5F_ACC_TRUNC); 
    } catch (const H5::FileIException&) { 
      file = H5::H5File(filename, H5F_ACC_EXCL);
    }

    // Group: /VTKHDF
    grp_vtkhdf = ensure_group(file, "/VTKHDF");

    // Attribute: /VTKHDF/Type
    ensure_attr_fixed_ascii_scalar(grp_vtkhdf, "Type", type_to_string(type));

    // Attribute: /VTKHDF/Version
    long long attr_version_value[2] = {version.major, version.minor};
    hsize_t attr_version_dims[1] = {2};
    H5::DataSpace attr_version_space(1, attr_version_dims);
    H5::Attribute attr_version = grp_vtkhdf.createAttribute("Version", H5::PredType::STD_I64LE, attr_version_space);
    attr_version.write(H5::PredType::NATIVE_LLONG, attr_version_value);
  }

public:
  VTKHDF(std::string filename, VTKHDFType type)
    : filename(filename)
    , type(type)
    , version({ .major = 2, .minor = 0 }) {
      ensure_common_schema();
    };
  VTKHDF(std::string filename, VTKHDFType type, VTKHDFVersion version)
    : filename(filename)
    , type(type)
    , version(version) {
      ensure_common_schema();
    };
};
}