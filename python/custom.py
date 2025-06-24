import h5py
import numpy as np

# File name to write
filename = "three_vertices_custom.vtkhdf"

# Data for the three vertices example
n_points = 3
points = np.array([
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0]
], dtype=np.float32)

# PointData arrays
normals = np.tile([0.0, 0.0, 1.0], (n_points, 1)).astype(np.float32)
warping = np.zeros((n_points, 3), dtype=np.float32)

# CellData: 3 vertex‐cells → one material ID per cell
materials = np.zeros((n_points,), dtype=np.int64)

# Vertex cell arrays
connectivity = np.arange(n_points, dtype=np.int64)          # [0,1,2]
n_cells = np.array([n_points], dtype=np.int64)              # [3]
n_conn_ids = np.array([connectivity.size], dtype=np.int64)  # [3]
# Offsets now length = n_cells + 1 → [0,1,2,3]
offsets = np.arange(n_points + 1, dtype=np.int64)

# Helper to create an empty cell‐type group
def create_empty_cell_group(parent, name):
    grp = parent.create_group(name)
    grp.create_dataset("Connectivity", data=np.zeros((0,), dtype=np.int64), maxshape=(None,))
    grp.create_dataset("NumberOfCells", data=np.array([0], dtype=np.int64), maxshape=(None,))
    grp.create_dataset("NumberOfConnectivityIds", data=np.array([0], dtype=np.int64), maxshape=(None,))
    grp.create_dataset("Offsets", data=np.zeros((1,), dtype=np.int64), maxshape=(None,))
    return grp

# Write the HDF5 file
with h5py.File(filename, "w") as f:
    vtk = f.create_group("VTKHDF")

    # Attributes
    vtk.attrs.create("Type", b"PolyData",
                     dtype=h5py.string_dtype("ascii", 8))
    vtk.attrs.create("Version", np.array([2, 0], dtype=np.int64))

    # Point counts and coordinates
    vtk.create_dataset("NumberOfPoints", data=np.array([n_points], dtype=np.int64),
                       maxshape=(None,), dtype=np.int64)
    vtk.create_dataset("Points", data=points,
                       maxshape=(None, 3), dtype=np.float32)

    # PointData group
    # pd = vtk.create_group("PointData")
    # pd.create_dataset("Normals", data=normals,
    #                   maxshape=(None, 3), dtype=np.float32)
    # pd.create_dataset("Warping", data=warping,
    #                   maxshape=(None, 3), dtype=np.float32)

    # CellData group
    # cd = vtk.create_group("CellData")
    # cd.create_dataset("Materials", data=materials,
    #                   maxshape=(None,), dtype=np.int64)

    # Empty cell‐type groups: Lines, Polygons, Strips
    for name in ["Lines", "Polygons", "Strips"]:
        create_empty_cell_group(vtk, name)

    # Vertices group with real data
    vg = vtk.create_group("Vertices")
    vg.create_dataset("Connectivity", data=connectivity,
                      maxshape=(None,), dtype=np.int64)
    vg.create_dataset("NumberOfCells", data=n_cells,
                      maxshape=(None,), dtype=np.int64)
    vg.create_dataset("NumberOfConnectivityIds", data=n_conn_ids,
                      maxshape=(None,), dtype=np.int64)
    vg.create_dataset("Offsets", data=offsets,
                      maxshape=(None,), dtype=np.int64)

print(f"Wrote {filename}")
