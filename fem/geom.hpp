// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#pragma once

// #include "../config/config.hpp"
// #include "../general/hash.hpp"
// #include "../linalg/densemat.hpp"
// #include "intrules.hpp"

// #include <memory>
// #include <unordered_map>

#include "../config/config.hpp"

namespace beam {
namespace fem {
/** Types of domains for integration rules and reference finite elements:
    Geometry::POINT    - a point
    Geometry::SEGMENT  - the interval [0,1]
    Geometry::TRIANGLE - triangle with vertices (0,0), (1,0), (0,1)
    Geometry::SQUARE   - the unit square (0,1)x(0,1)
    Geometry::TETRAHEDRON - w/ vert. (0,0,0),(1,0,0),(0,1,0),(0,0,1)
    Geometry::CUBE - the unit cube
    Geometry::PRISM - w/ vert. (0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,0,1),(0,1,1)
    Geometry::PYRAMID - w/ vert. (0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1)
*/
class Geometry {
public:
  enum Type {
    INVALID = -1,
    POINT = 0,
    SEGMENT,
    TRIANGLE,
    // SQUARE,
    // TETRAHEDRON,
    // CUBE,
    // PRISM,
    // PYRAMID,
    NUM_GEOMETRIES
  };

  static const int NumGeom = NUM_GEOMETRIES;
  static const int MaxDim = 3;
  static const int NumBdrArray[NumGeom];
  static const char *Name[NumGeom];
  static const real_t Volume[NumGeom];
  static const int Dimension[NumGeom];
  static const int DimStart[MaxDim + 2]; // including MaxDim+1
  static const int NumVerts[NumGeom];
  static const int NumEdges[NumGeom];
  static const int NumFaces[NumGeom];

  // Structure that holds constants describing the Geometries.
  template <Type Geom> struct Constants;

private:
  // IntegrationRule *GeomVert[NumGeom];
  // IntegrationPoint GeomCenter[NumGeom];
  // DenseMatrix *GeomToPerfGeomJac[NumGeom];
  // DenseMatrix *PerfGeomToGeomJac[NumGeom];

public:
  Geometry();
  ~Geometry();

  // const IntegrationRule *GetVertices(int GeomType) const;
  // const IntegrationPoint &GetCenter(int GeomType) const { return GeomCenter[GeomType]; }
  // static void GetRandomPoint(int GeomType, IntegrationPoint &ip);
  // static bool CheckPoint(int GeomType, const IntegrationPoint &ip);
  // static bool CheckPoint(int GeomType, const IntegrationPoint &ip, real_t eps);
  // static bool ProjectPoint(int GeomType, const IntegrationPoint &beg, IntegrationPoint &end);
  // static bool ProjectPoint(int GeomType, IntegrationPoint &ip);
  // const DenseMatrix &GetGeomToPerfGeomJac(int GeomType) const { return *GeomToPerfGeomJac[GeomType]; }
  // const DenseMatrix *GetPerfGeomToGeomJac(int GeomType) const { return PerfGeomToGeomJac[GeomType]; }
  // void GetPerfPointMat(int GeomType, DenseMatrix &pm) const;
  // void JacToPerfJac(int GeomType, const DenseMatrix &J, DenseMatrix &PJ) const;
  // static bool IsTensorProduct(Type geom) { return geom == SEGMENT || geom == SQUARE || geom == CUBE; }

  // static Type TensorProductGeometry(int dim) {
  //   switch (dim) {
  //   case 0:
  //     return POINT;
  //   case 1:
  //     return SEGMENT;
  //   case 2:
  //     return SQUARE;
  //   case 3:
  //     return CUBE;
  //   default:
  //     MFEM_ABORT("Invalid dimension.");
  //     return INVALID;
  //   }
  // }

  // static int GetInverseOrientation(Type geom_type, int orientation);
  // int NumBdr(int GeomType) const { return NumBdrArray[GeomType]; }
};

template <> struct Geometry::Constants<Geometry::POINT> {
  static const int Dimension = 0;
  static const int NumVert = 1;

  static const int NumOrient = 1;
  static const int Orient[NumOrient][NumVert];
  static const int InvOrient[NumOrient];
};

template <> struct Geometry::Constants<Geometry::SEGMENT> {
  static const int Dimension = 1;
  static const int NumVert = 2;
  static const int NumEdges = 1;
  static const int Edges[NumEdges][2];

  static const int NumOrient = 2;
  static const int Orient[NumOrient][NumVert];
  static const int InvOrient[NumOrient];
};

template <> struct Geometry::Constants<Geometry::TRIANGLE> {
  static const int Dimension = 2;
  static const int NumVert = 3;
  static const int NumEdges = 3;
  static const int Edges[NumEdges][2];
  // Upper-triangular part of the local vertex-to-vertex graph.
  struct VertToVert {
    static const int I[NumVert];
    static const int J[NumEdges][2]; // {end,edge_idx}
  };
  static const int NumFaces = 1;
  static const int FaceVert[NumFaces][NumVert];
  // For a given base tuple v={v0,v1,v2}, the orientation of a permutation
  // u={u0,u1,u2} of v, is an index 'j' such that u[i]=v[Orient[j][i]].
  // The static method Mesh::GetTriOrientation, computes the index 'j' of the
  // permutation that maps the second argument 'test' to the first argument
  // 'base': test[Orient[j][i]]=base[i].
  static const int NumOrient = 6;
  static const int Orient[NumOrient][NumVert];
  // The inverse of orientation 'j' is InvOrient[j].
  static const int InvOrient[NumOrient];
};

// template <> struct Geometry::Constants<Geometry::SQUARE> {
//   static const int Dimension = 2;
//   static const int NumVert = 4;
//   static const int NumEdges = 4;
//   static const int Edges[NumEdges][2];
//   // Upper-triangular part of the local vertex-to-vertex graph.
//   struct VertToVert {
//     static const int I[NumVert];
//     static const int J[NumEdges][2]; // {end,edge_idx}
//   };
//   static const int NumFaces = 1;
//   static const int FaceVert[NumFaces][NumVert];
//   static const int NumOrient = 8;
//   static const int Orient[NumOrient][NumVert];
//   static const int InvOrient[NumOrient];
// };

// template <> struct
//    Geometry::Constants<Geometry::TETRAHEDRON>
// {
//    static const int Dimension = 3;
//    static const int NumVert = 4;
//    static const int NumEdges = 6;
//    static const int Edges[NumEdges][2];
//    static const int NumFaces = 4;
//    static const int FaceTypes[NumFaces];
//    static const int MaxFaceVert = 3;
//    static const int FaceVert[NumFaces][MaxFaceVert];
//    // Upper-triangular part of the local vertex-to-vertex graph.
//    struct VertToVert
//    {
//       static const int I[NumVert];
//       static const int J[NumEdges][2]; // {end,edge_idx}
//    };
//    static const int NumOrient = 24;
//    static const int Orient[NumOrient][NumVert];
//    static const int InvOrient[NumOrient];
// };

// template <> struct Geometry::Constants<Geometry::CUBE> {
//   static const int Dimension = 3;
//   static const int NumVert = 8;
//   static const int NumEdges = 12;
//   static const int Edges[NumEdges][2];
//   static const int NumFaces = 6;
//   static const int FaceTypes[NumFaces];
//   static const int MaxFaceVert = 4;
//   static const int FaceVert[NumFaces][MaxFaceVert];
//   // Upper-triangular part of the local vertex-to-vertex graph.
//   struct VertToVert {
//     static const int I[NumVert];
//     static const int J[NumEdges][2]; // {end,edge_idx}
//   };
// };

// template <> struct
//    Geometry::Constants<Geometry::PRISM>
// {
//    static const int Dimension = 3;
//    static const int NumVert = 6;
//    static const int NumEdges = 9;
//    static const int Edges[NumEdges][2];
//    static const int NumFaces = 5;
//    static const int FaceTypes[NumFaces];
//    static const int MaxFaceVert = 4;
//    static const int FaceVert[NumFaces][MaxFaceVert];
//    // Upper-triangular part of the local vertex-to-vertex graph.
//    struct VertToVert
//    {
//       static const int I[NumVert];
//       static const int J[NumEdges][2]; // {end,edge_idx}
//    };
// };

// template <> struct
//    Geometry::Constants<Geometry::PYRAMID>
// {
//    static const int Dimension = 3;
//    static const int NumVert = 5;
//    static const int NumEdges = 8;
//    static const int Edges[NumEdges][2];
//    static const int NumFaces = 5;
//    static const int FaceTypes[NumFaces];
//    static const int MaxFaceVert = 4;
//    static const int FaceVert[NumFaces][MaxFaceVert];
//    // Upper-triangular part of the local vertex-to-vertex graph.
//    struct VertToVert
//    {
//       static const int I[NumVert];
//       static const int J[NumEdges][2]; // {end,edge_idx}
//    };
// };
}
} // namespace beam
