#include "foreach_cell_bfs.h"
#include "vtkType.h"

/**
 * @brief This struct and the associated members collect descriptive information from basilisk to describe the tree as a
 * HyperTreeGrid.
 */
typedef struct {
    size_t max_vertices;    /**< Number of cells in the tree level with the greatest number of cells */
    int64_t depth_per_tree; /**< Maximum depth of our tree */

    Bit_t *descriptors;      /**< Breadth-first search of our tree. See @ref vtk_hdf_get_descriptors */
    size_t descriptors_size; /**< Number of bytes in @ref descriptors */
    int64_t n_descriptors;   /**< Number of bits in @ref descriptors */

    int64_t number_of_cells;                 /**< Number of cells in the tree */
    int64_t *number_of_cells_per_tree_depth; /**< Number of the cells (per-level) in the tree */
    int64_t number_of_trees;                 /**< Number of trees; in basilisk this is always 1. */
    int64_t tree_ids;                        /**< Ids of each tree; we only have one tree with id 0. */

    bool has_mask;    /**< If we are using a mask or not */
    Bit_t *mask;      /**< Binary mask to hide cells in the tree */
    size_t mask_size; /**< Number of bytes in @ref mask */
    size_t n_masks;   /**< Number of bits in @ref mask */

    double *x;  /**< Bounds of the tree on the x-axis. */
    double *y;  /**< Bounds of the tree on the y-axis. */
    double *z;  /**< Bounds of the tree on the z-axis. */
    size_t n_x; /**< Number of bounds of the tree on the x-axis. See @ref x */
    size_t n_y; /**< Number of bounds of the tree on the y-axis. See @ref y */
    size_t n_z; /**< Number of bounds of the tree on the z-axis. See @ref z */

} vtkHDFHyperTreeGridData;

/**
 * @brief Get the X/Y/Z corners of the tree.
 * 
 * @memberof vtkHDFHyperTreeGridData
 */
void vtk_hdf_get_coordinates(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {

#if dimension >= 1
    vtk_hdf_hypertreegrid_data->n_x = 2;
    vtk_hdf_hypertreegrid_data->x = malloc(2 * sizeof(double));
    vtk_hdf_hypertreegrid_data->x[0] = X0;
    vtk_hdf_hypertreegrid_data->x[1] = X0 + L0;
#else
    vtk_hdf_hypertreegrid_data->n_x = 1;
    vtk_hdf_hypertreegrid_data->x = malloc(1 * sizeof(double));
    vtk_hdf_hypertreegrid_data->x[0] = 0;
#endif

#if dimension >= 2
    vtk_hdf_hypertreegrid_data->n_y = 2;
    vtk_hdf_hypertreegrid_data->y = malloc(2 * sizeof(double));
    vtk_hdf_hypertreegrid_data->y[0] = Y0;
    vtk_hdf_hypertreegrid_data->y[1] = Y0 + L0;
#else
    vtk_hdf_hypertreegrid_data->n_y = 1;
    vtk_hdf_hypertreegrid_data->y = malloc(1 * sizeof(double));
    vtk_hdf_hypertreegrid_data->y[0] = 0;
#endif

#if dimension >= 3
    vtk_hdf_hypertreegrid_data->n_z = 2;
    vtk_hdf_hypertreegrid_data->z = malloc(2 * sizeof(double));
    vtk_hdf_hypertreegrid_data->z[0] = Z0;
    vtk_hdf_hypertreegrid_data->z[1] = Z0 + L0;
#else
    vtk_hdf_hypertreegrid_data->n_z = 1;
    vtk_hdf_hypertreegrid_data->z = malloc(1 * sizeof(double));
    vtk_hdf_hypertreegrid_data->z[0] = 0;
#endif
}

/**
 * @brief Count the number of cells on each level of the (binary/quad/oct)tree. Also compute the maximum
 *
 * @memberof vtkHDFHyperTreeGridData
 */
void vtk_hdf_get_number_of_vertices_per_depth(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {

    int depth = grid->maxdepth;

    // int64_t depth_per_tree_old = vtk_hdf_hypertreegrid_data->depth_per_tree;
    vtk_hdf_hypertreegrid_data->depth_per_tree = (int64_t)(depth + 1);

    // Check if we have already counted this previously
    if (vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth != NULL) {
        // Free and reallocate
        free(vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth);
        vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth =
            malloc(vtk_hdf_hypertreegrid_data->depth_per_tree * sizeof(int64_t));

    } else {
        // If we havent, allocate the array
        vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth =
            malloc(vtk_hdf_hypertreegrid_data->depth_per_tree * sizeof(int64_t));
    }

    // Set the count to zero on each level
    memset(vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth, 0,
           vtk_hdf_hypertreegrid_data->depth_per_tree * sizeof(int64_t));

    // Count the cells on each level
    foreach_cell_BFS() { vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth[level]++; }

    // Compute the largest cell count on any level and total number of cells
    vtk_hdf_hypertreegrid_data->max_vertices = 0;
    vtk_hdf_hypertreegrid_data->number_of_cells = 0;

    for (size_t i = 0; i < vtk_hdf_hypertreegrid_data->depth_per_tree; i++) {
        int64_t nv = vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth[i];
        if (nv > vtk_hdf_hypertreegrid_data->max_vertices)
            vtk_hdf_hypertreegrid_data->max_vertices = nv;
        vtk_hdf_hypertreegrid_data->number_of_cells += nv;
    }
}

/**
 * @brief Create a mask for the grid to hide remote leaf cells.
 *
 * @memberof vtkHDFHyperTreeGridData
 */
void vtk_hdf_get_local_mask(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {

    size_t number_of_cells = vtk_hdf_hypertreegrid_data->number_of_cells;
    // size_t max_depth = vtk_hdf_hypertreegrid_data->depth_per_tree;

    size_t n_masks = number_of_cells;
    vtk_hdf_hypertreegrid_data->n_masks = n_masks;
    size_t mask_size = (size_t)(((n_masks + 7) / 8));
    vtk_hdf_hypertreegrid_data->mask_size = mask_size;

    vtk_hdf_hypertreegrid_data->mask = malloc(mask_size * sizeof(Bit_t));
    memset(vtk_hdf_hypertreegrid_data->mask, 0, mask_size);

    uint8_t bit_count = 0;
    size_t byte_count = 0;

    foreach_cell_BFS() {
        if (!(is_local(cell)) && is_leaf(cell)) {
            bool val = true;
            vtk_hdf_hypertreegrid_data->mask[byte_count] |= (uint8_t)(val << (7 - bit_count));
        }

        bit_count++;

        if (bit_count == 8) {
            bit_count = 0;
            byte_count++;
        }
    }

    vtk_hdf_hypertreegrid_data->has_mask = true;
}

/**
 * @brief Traverse the (binary-/quad-/oct-)tree using breadth-first search and record if the node is refined or not.
 *
 * This function uses the custom @ref foreach_cell_BFS() iterator macro to perform breadth first search. If the nth
 * visited node is not a leaf (i.e. it is refined), then we set the corresponding bit in our byte array to 1 (true).
 * Otherwise we set it to 0 (false).
 *
 * @memberof vtkHDFHyperTreeGridData
 */
void vtk_hdf_get_descriptors(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {

    // Get the total number of cells
    size_t total_vertices = vtk_hdf_hypertreegrid_data->number_of_cells;

    // Get the maximum depth of the tree
    size_t max_depth = vtk_hdf_hypertreegrid_data->depth_per_tree;

    // We decide that the total number of descriptors (bits) we need is the number of cells in all but the deepest
    // level: Since this is the deepest level then it follows that any nodes there must be leaves, and so any leaves
    // there are implied by the 1 from their parents
    size_t n_descriptors = total_vertices - vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth[max_depth - 1];

    // Store the result in our object; we need this later
    vtk_hdf_hypertreegrid_data->n_descriptors = (int64_t)n_descriptors;

    // Calculate the number of bytes needed to pack n number of descriptors
    size_t descriptors_size = (size_t)(((n_descriptors + 7) / 8));

    // Store this result
    vtk_hdf_hypertreegrid_data->descriptors_size = descriptors_size;

    // Heap allocate the byte array
    vtk_hdf_hypertreegrid_data->descriptors = malloc(descriptors_size * sizeof(Bit_t));
    memset(vtk_hdf_hypertreegrid_data->descriptors, 0, descriptors_size);

    // Keep a count of the bit + byte where we are in the byte array
    uint8_t bit_count = 0;
    size_t byte_count = 0;

    // Breadth-first search
    foreach_cell_BFS() {

        // If the node is refined (not a leaf), set the correct bit of the byte to 1
        if (!is_leaf(cell))
            vtk_hdf_hypertreegrid_data->descriptors[byte_count] |= (uint8_t)(1 << (7 - bit_count));

        // Increment bit
        bit_count++;

        // Increment byte
        if (bit_count == 8) {
            bit_count = 0;
            byte_count++;
        }
    }
}

/**
 * @brief Desctructor for the @ref vtkHDFHyperTreeGridData struct.
 *
 * @memberof vtkHDFHyperTreeGridData
 */
void vtk_hdf_hypertreegrid_data_free(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {
    if (vtk_hdf_hypertreegrid_data) {
        free(vtk_hdf_hypertreegrid_data->x);
        free(vtk_hdf_hypertreegrid_data->y);
        free(vtk_hdf_hypertreegrid_data->z);
        free(vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth);
        free(vtk_hdf_hypertreegrid_data->descriptors);
        free(vtk_hdf_hypertreegrid_data->mask);
    }
}

/**
 * @brief Constructor for the @ref vtkHDFHyperTreeGridData struct. 
 *
 * @memberof vtkHDFHyperTreeGridData
 */
vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data_init(void) {
    vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data = calloc(1, sizeof(vtkHDFHyperTreeGridData));

    if (!vtk_hdf_hypertreegrid_data) {
        // handle malloc error
    };

    if (!grid) {
        fprintf(stderr, "Grid has not yet been initialized: Please initialize the "
                        "grid first.\n");
        exit(1);
    }

    vtk_hdf_hypertreegrid_data->number_of_trees = 1;
    vtk_hdf_hypertreegrid_data->tree_ids = 0;

    vtk_hdf_get_coordinates(vtk_hdf_hypertreegrid_data);
    vtk_hdf_get_number_of_vertices_per_depth(vtk_hdf_hypertreegrid_data);
    vtk_hdf_get_descriptors(vtk_hdf_hypertreegrid_data);

    #if MPI_SINGLE_FILE
        vtk_hdf_get_local_mask(vtk_hdf_hypertreegrid_data);
    #endif

    return vtk_hdf_hypertreegrid_data;
}
