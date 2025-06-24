#!/usr/bin/env python3
import numpy as np
import h5py as h5

from vtkmodules.vtkCommonCore import vtkDoubleArray, vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData, vtkCellArray
from vtkmodules.vtkFiltersGeometry import vtkGeometryFilter
from vtkmodules.vtkFiltersCore import vtkAppendFilter
from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
import vtkmodules.numpy_interface.dataset_adapter as dsa
from vtkmodules.util.numpy_support import numpy_to_vtk as npvtk, vtk_to_numpy as vtknp

# cell‐type names in the VTKHDF layout
connectivities = ['Vertices','Lines','Polygons','Strips']

# ─── Set up the empty HDF structure ────────────────────────────────────────
def generate_geometry_structure(root):
    root.attrs['Version'] = (2,0)
    ascii_type = b'PolyData'
    root.attrs.create('Type', ascii_type, dtype=h5.string_dtype('ascii', len(ascii_type)))

    root.create_dataset('NumberOfPoints', (0,),  maxshape=(None,),   dtype='i8')
    root.create_dataset('Points',         (0,3), maxshape=(None,3), dtype='f4')

    for name in connectivities:
        grp = root.create_group(name)
        grp.create_dataset('NumberOfConnectivityIds', (0,), maxshape=(None,), dtype='i8')
        grp.create_dataset('NumberOfCells',           (0,), maxshape=(None,), dtype='i8')
        grp.create_dataset('Offsets',                 (0,), maxshape=(None,), dtype='i8')
        grp.create_dataset('Connectivity',            (0,), maxshape=(None,), dtype='i8')

    pd = root.create_group('PointData')
    pd.create_dataset('Warping', (0,3), maxshape=(None,3), dtype='f4')
    pd.create_dataset('Normals', (0,3), maxshape=(None,3), dtype='f4')

    cd = root.create_group('CellData')
    cd.create_dataset('Materials', (0,), maxshape=(None,), dtype='i8')

def generate_step_structure(root):
    steps = root.create_group('Steps')
    steps.attrs['NSteps'] = 0
    steps.create_dataset('Values', (0,), maxshape=(None,), dtype='f4')

    for name in ['PartOffsets','NumberOfParts','PointOffsets']:
        steps.create_dataset(name, (0,), maxshape=(None,), dtype='i8')
    for name in ['CellOffsets','ConnectivityIdOffsets']:
        steps.create_dataset(name, (0,4), maxshape=(None,4), dtype='i8')

    pdo = steps.create_group('PointDataOffsets')
    pdo.create_dataset('Warping', (0,), maxshape=(None,), dtype='i8')
    pdo.create_dataset('Normals', (0,), maxshape=(None,), dtype='i8')

    cdo = steps.create_group('CellDataOffsets')
    cdo.create_dataset('Materials', (0,), maxshape=(None,), dtype='i8')

def append_dataset(dset, array):
    orig = dset.shape[0]
    dset.resize(orig + array.shape[0], axis=0)
    dset[orig:] = array

def append_data(root, poly_data, newStep):
    # 1) Increment time‐step
    steps = root['Steps']
    steps.attrs['NSteps'] += 1
    append_dataset(steps['Values'], np.array([newStep],dtype='f4'))

    # 2) Append point data
    npts = poly_data.GetNumberOfPoints()
    pts  = vtknp(poly_data.GetPoints().GetData())
    append_dataset(root['NumberOfPoints'], np.array([npts],dtype='i8'))
    append_dataset(root['Points'], pts)

    # 3) Append each cell‐type (only Vertices is nonempty here)
    cell_arrays = {
        'Vertices': poly_data.GetVerts(),
        'Lines':    poly_data.GetLines(),
        'Polygons': poly_data.GetPolys(),
        'Strips':   poly_data.GetStrips()
    }
    for name, ca in cell_arrays.items():
        conn = vtknp(ca.GetConnectivityArray())
        offs = vtknp(ca.GetOffsetsArray())
        append_dataset(root[name]['NumberOfConnectivityIds'], np.array([conn.size],dtype='i8'))
        append_dataset(root[name]['Connectivity'],            conn)
        append_dataset(root[name]['NumberOfCells'],           np.array([offs.size-1],dtype='i8'))
        append_dataset(root[name]['Offsets'],                 offs)

    # 4) Append point‐data arrays (Warping, Normals)
    warping = vtknp(poly_data.GetPointData().GetArray('Warping'))
    normals = vtknp(poly_data.GetPointData().GetArray('Normals'))
    append_dataset(root['PointData']['Warping'], warping)
    append_dataset(root['PointData']['Normals'], normals)

    # 5) Append cell‐data (Materials)
    mats = vtknp(poly_data.GetCellData().GetArray('Materials'))
    append_dataset(root['CellData']['Materials'], mats)

# ─── Build a vtkPolyData with exactly three vertices ───────────────────────
def generate_three_vertices():
    # a) Points
    pts = vtkPoints()
    pts.InsertNextPoint(0.0, 0.0, 0.0)
    pts.InsertNextPoint(1.0, 0.0, 0.0)
    pts.InsertNextPoint(0.0, 1.0, 0.0)

    # b) Vertices cell‐array
    verts = vtkCellArray()
    for i in range(3):
        verts.InsertNextCell(1)
        verts.InsertCellPoint(i)

    # c) Assemble
    pd = vtkPolyData()
    pd.SetPoints(pts)
    pd.SetVerts(verts)

    # d) Dummy point‐data arrays
    warping = npvtk(np.zeros((3,3),dtype='f4'))
    warping.SetName('Warping')
    pd.GetPointData().AddArray(warping)

    normals = vtkDoubleArray()
    normals.SetNumberOfComponents(3)
    normals.SetName('Normals')
    for _ in range(3):
        normals.InsertNextTuple((0,0,1))
    pd.GetPointData().AddArray(normals)

    # e) Dummy cell‐data array (one material per vertex)
    mats = npvtk(np.zeros((3,),dtype='i8'))
    mats.SetName('Materials')
    pd.GetCellData().AddArray(mats)

    return pd

# ─── Main: write both .vtkhdf and a .vtp “piece” ─────────────────────────
if __name__ == '__main__':
    # 1) Create HDF file
    with h5.File('three_vertices.vtkhdf','w') as f:
        root = f.create_group('VTKHDF')
        generate_geometry_structure(root)
        generate_step_structure(root)

        three_pd = generate_three_vertices()
        append_data(root, three_pd, newStep=0.0)

    print("Wrote three_vertices.vtkhdf")

    # 2) Write matching XML piece for quick visualization
    #    (this writes only the vertices, no connectivity beyond that)
    writer = vtkXMLPolyDataWriter()
    writer.SetFileName('three_vertices_piece.vtp')
    writer.SetInputData(generate_three_vertices())
    writer.Write()
    print("Wrote three_vertices_piece.vtp")
