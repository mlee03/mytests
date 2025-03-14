/***********************************************************************
 *                   GNU Lesser General Public License
 *
 * This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
 *
 * FRE-NCtools is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * FRE-NCtools is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FRE-NCTools.  If not, see
 * <http://www.gnu.org/licenses/>.
 **********************************************************************/
#ifndef CREATE_XGRID_GPU_H_
#define CREATE_XGRID_GPU_H_

#include "globals_gpu.h"

int get_upbound_nxcells_2dx2d_gpu_(int *nlon_input_cells_in,  int *nlat_input_cells_in,
                                   int *nlon_output_cells_in,  int *nlat_output_cells_in,
                                   double *input_grid_lon,  double *input_grid_lat,
                                   double *output_grid_lon,  double *output_grid_lat,
                                   double *output_cell_lon_min, double *output_cell_lon_max,
                                   double *output_cell_lat_min, double *output_cell_lat_max,
                                   double *output_cell_lon_cent,double *skip_input_cells,
                                   int *approx_nxcells_per_ij1, int *ij2_start, int *ij2_end);

int get_upbound_nxcells_2dx2d_gpu(int *nlon_input_cells_in,  int *nlat_input_cells_in,
                                  int *nlon_output_cells_in,  int *nlat_output_cells_in,
                                  double *input_grid_lon,  double *input_grid_lat,
                                  double *output_grid_lon,  double *output_grid_lat,
                                  double *output_cell_lon_min, double *output_cell_lon_max,
                                  double *output_cell_lat_min, double *output_cell_lat_max,
                                  double *output_cell_lon_cent,double *skip_input_cells,
                                  int *approx_nxcells_per_ij1, int *ij2_start, int *ij2_end);

int create_xgrid_2dx2d_order1_gpu_(int *nlon_input_cells_in, int *nlat_input_cells_in,
                                   int *nlon_output_cells_in, int *nlat_output_cells_in,
                                   double *input_grid_lon, double *input_grid_lat,
                                   double *output_grid_lon, double *output_grid_lat, int nxcells_ish,
                                   double *output_cells_lon_min, double *output_cells_lon_max,
                                   double *output_cells_lat_min, double *output_cells_lat_max,
                                   double *output_cells_lon_cent, double *output_cells_area, int *output_cells_nvertices,
                                   double *output_cells_lon_vertices, double *output_cells_lat_vertices,
                                   double *mask_input_grid, int *approx_nxcells_ij1, int *ij2_start, int *ij2_end,
                                   int *input_parent_index, int *output_parent_index, double *xarea);

int create_xgrid_2dx2d_order1_gpu(int *nlon_input_cells_in, int *nlat_input_cells_in,
                                  int *nlon_output_cells_in, int *nlat_output_cells_in,
                                  double *input_grid_lon, double *input_grid_lat,
                                  double *output_grid_lon, double *output_grid_lat, int nxcells_ish,
                                  double *output_cells_lon_min, double *output_cells_lon_max,
                                  double *output_cells_lat_min, double *output_cells_lat_max,
                                  double *output_cells_lon_cent, double *output_cells_area, int *output_cells_nvertices,
                                  double *output_cells_lon_vertices, double *output_cells_lat_vertices,
                                  double *mask_input_grid, int *approx_nxcells_ij1, int *ij2_start, int *ij2_end,
                                  int *input_parent_index, int *output_parent_index, double *xarea);

#endif
