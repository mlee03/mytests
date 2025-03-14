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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <openacc.h>
#include "general_utils_gpu.h"
#include "create_xgrid_gpu.h"
#include "create_xgrid_utils_gpu.h"
#include "globals_gpu.h"

/*******************************************************************************
void get_upbound_nxcells_2dx2d_gpu
This function computes the upperbound to nxgrid.  This upper bound will be used
to malloc arrays used in create_xgrid
*******************************************************************************/
int get_upbound_nxcells_2dx2d_gpu_(int *nlon_input_cells_in, int *nlat_input_cells_in,
                                   int *nlon_output_cells_in, int *nlat_output_cells_in,
                                   double *input_grid_lon, double *input_grid_lat,
                                   double *output_grid_lon, double *output_grid_lat,
                                   double *output_cell_lon_min, double *output_cell_lon_max,
                                   double *output_cell_lat_min, double *output_cell_lat_max,
                                   double *output_cell_lon_cent, double *skip_input_cells,                 
                                   int *approx_xcells_per_ij1, int *ij2_start, int *ij2_end)
{
  return get_upbound_nxcells_2dx2d_gpu(nlon_input_cells_in, nlat_input_cells_in, nlon_output_cells_in, nlat_output_cells_in,
                                       input_grid_lon, input_grid_lat, output_grid_lon, output_grid_lat,
                                       output_cell_lon_min, output_cell_lon_max, output_cell_lat_min, output_cell_lat_max,
                                       output_cell_lon_cent, skip_input_cells, approx_xcells_per_ij1, ij2_start, ij2_end);
}


int get_upbound_nxcells_2dx2d_gpu(int *nlon_input_cells_in, int *nlat_input_cells_in,
                                  int *nlon_output_cells_in, int *nlat_output_cells_in,
                                  double *input_grid_lon, double *input_grid_lat,
                                  double *output_grid_lon, double *output_grid_lat,
                                  double *output_cell_lon_min, double *output_cell_lon_max,
                                  double *output_cell_lat_min, double *output_cell_lat_max,
                                  double *output_cell_lon_cent, double *skip_input_cells,                 
                                  int *approx_xcells_per_ij1, int *ij2_start, int *ij2_end)
{
  
  int nlon_input_cells = *nlon_input_cells_in;
  int nlat_input_cells = *nlat_input_cells_in;
  int nlon_output_cells = *nlon_output_cells_in;
  int nlat_output_cells = *nlat_output_cells_in;
  
  int input_grid_ncells  = (nlon_input_cells)*(nlat_input_cells);
  int output_grid_ncells = (nlon_output_cells)*(nlat_output_cells);
  int input_grid_npts    = (nlon_input_cells+1)*(nlat_input_cells+1);
  int output_grid_npts   = (nlon_output_cells+1)*(nlat_output_cells+1);
  
  int upbound_nxcells=0;
  
#pragma acc data present(output_grid_lon[:output_grid_npts],          \
                         output_grid_lat[:output_grid_npts],          \
                         input_grid_lon[:input_grid_npts],            \
                         input_grid_lat[:input_grid_npts],            \
                         output_cell_lon_min[:output_grid_ncells],     \
                         output_cell_lat_min[:output_grid_ncells],     \
                         output_cell_lon_max[:output_grid_ncells],     \
                         output_cell_lat_max[:output_grid_ncells],     \
                         output_cell_lon_cent[:output_grid_ncells],    \
                         approx_xcells_per_ij1[:input_grid_ncells],    \
                         ij2_start[:input_grid_ncells],                \
                         ij2_end[:input_grid_ncells],                  \
                         skip_input_cells[:input_grid_ncells])
#pragma acc data copyin(input_grid_ncells, output_grid_ncells)
#pragma acc data copy(upbound_nxcells)
#pragma acc parallel loop independent reduction(+:upbound_nxcells)
  for( int ij1=0 ; ij1<input_grid_ncells ; ij1++) {
    if( skip_input_cells[ij1] > MASK_THRESH ) {
      
      int i_approx_xcells_per_ij1=0;
      int ij2_min=output_grid_ncells, ij2_max=0;
      double input_cell_lon_vertices[MV], input_cell_lat_vertices[MV];

      get_cell_vertices_gpu(ij1, nlon_input_cells, input_grid_lon, input_grid_lat,
                            input_cell_lon_vertices, input_cell_lat_vertices);

      double input_cell_lat_min = minval_double_gpu(4, input_cell_lat_vertices);
      double input_cell_lat_max = maxval_double_gpu(4, input_cell_lat_vertices);
      int nvertices = fix_lon_gpu(input_cell_lon_vertices, input_cell_lat_vertices, 4, M_PI);
      double input_cell_lon_min = minval_double_gpu(nvertices, input_cell_lon_vertices);
      double input_cell_lon_max = maxval_double_gpu(nvertices, input_cell_lon_vertices);
      double input_cell_lon_cent = avgval_double_gpu(nvertices, input_cell_lon_vertices);

      approx_xcells_per_ij1[ij1]=0;

#pragma acc loop independent reduction(+:upbound_nxcells) reduction(+:i_approx_xcells_per_ij1) \
                             reduction(min:ij2_min) reduction(max:ij2_max)

      for(int ij2=0; ij2<output_grid_ncells; ij2++) {

        double dlon_cent, output_cell_lon_min_i, output_cell_lon_max_i;
        double rotate=0.0;
        
        if(output_cell_lat_min[ij2] >= input_cell_lat_max) continue;
        if(output_cell_lat_max[ij2] <= input_cell_lat_min) continue;
        
        dlon_cent = output_cell_lon_cent[ij2] - input_cell_lon_cent;
        if(dlon_cent < -M_PI) rotate = +TPI;
        if(dlon_cent > M_PI)  rotate = -TPI;

        // adjust according to input_grid_lon_cent
        output_cell_lon_min_i = output_cell_lon_min[ij2] + rotate;
        output_cell_lon_max_i = output_cell_lon_max[ij2] + rotate;
        
        //output_cell_lon should in the same range as input_cell_lon after lon_fix,
        // so no need to consider cyclic condition
        if(output_cell_lon_min_i >= input_cell_lon_max ) continue;
        if(output_cell_lon_max_i <= input_cell_lon_min ) continue;
        
        //Note, the check for AREA_RATIO_THRESH has been removed
        //Thus, the computed value of upbound_nxcells will be equal to or greater than nxgrid
        i_approx_xcells_per_ij1++;
        upbound_nxcells++;
        ij2_min = min(ij2_min, ij2);
        ij2_max = max(ij2_max, ij2);
        
      } //ij2
      approx_xcells_per_ij1[ij1] = i_approx_xcells_per_ij1;
      ij2_start[ij1] = ij2_min ;
      ij2_end[ij1]   = ij2_max;

    } //mask
  } //ij1

  return upbound_nxcells;

}

int create_xgrid_2dx2d_order1_gpu_(int *nlon_input_cells_in, int *nlat_input_cells_in,
                                   int *nlon_output_cells_in, int *nlat_output_cells_in,
                                   double *input_grid_lon, double *input_grid_lat,
                                   double *output_grid_lon, double *output_grid_lat, int nxcells_ish,
                                   double *output_cells_lon_min, double *output_cells_lon_max,
                                   double *output_cells_lat_min, double *output_cells_lat_max,
                                   double *output_cells_lon_cent, double *output_cells_area, int *output_cells_nvertices,
                                  double *output_cells_lon_vertices, double *output_cells_lat_vertices, 
                                   double *mask_input_grid, int *approx_nxcells_ij1, int *ij2_start, int *ij2_end,
                                   int *input_parent_index, int *output_parent_index, double *xarea)
{  
  int nxcells = create_xgrid_2dx2d_order1_gpu(nlon_input_cells_in, nlat_input_cells_in, nlon_output_cells_in,
                                              nlat_output_cells_in, input_grid_lon, input_grid_lat, output_grid_lon,
                                              output_grid_lat, nxcells_ish, output_cells_lon_min, output_cells_lon_max,
                                              output_cells_lat_min, output_cells_lat_max, output_cells_lon_cent,
                                              output_cells_area, output_cells_nvertices, output_cells_lon_vertices,
                                              output_cells_lat_vertices, mask_input_grid, approx_nxcells_ij1, ij2_start,
                                              ij2_end, input_parent_index, output_parent_index, xarea);
  return nxcells;
}

int create_xgrid_2dx2d_order1_gpu(int *nlon_input_cells_in, int *nlat_input_cells_in,
                                  int *nlon_output_cells_in, int *nlat_output_cells_in,
                                  double *input_grid_lon, double *input_grid_lat,
                                  double *output_grid_lon, double *output_grid_lat, int nxcells_ish,
                                  double *output_cells_lon_min, double *output_cells_lon_max,
                                  double *output_cells_lat_min, double *output_cells_lat_max,
                                  double *output_cells_lon_cent, double *output_cells_area, int *output_cells_nvertices,
                                  double *output_cells_lon_vertices, double *output_cells_lat_vertices, 
                                  double *mask_input_grid, int *approx_nxcells_ij1, int *ij2_start, int *ij2_end,
                                  int *input_parent_index, int *output_parent_index, double *xarea)
{

  if(nxcells_ish<1) return 0;

  int nxcells=0;

  int nlon_input_cells = *nlon_input_cells_in;
  int nlat_input_cells = *nlat_input_cells_in;
  int nlon_output_cells = *nlon_output_cells_in;
  int nlat_output_cells = *nlat_output_cells_in;
  
  int input_grid_ncells  = nlon_input_cells*nlat_input_cells;
  int output_grid_ncells = nlon_output_cells*nlat_output_cells;
  int input_grid_npts    = (nlon_input_cells+1)*(nlat_input_cells+1);
  int output_grid_npts   = (nlon_output_cells+1)*(nlat_output_cells+1);

  int *input_parent_index_i  = NULL ; input_parent_index_i  = (int *)malloc(nxcells_ish*sizeof(int));
  int *output_parent_index_i = NULL ; output_parent_index_i = (int *)malloc(nxcells_ish*sizeof(int));
  double *xarea_i            = NULL ; xarea_i =  (double *)malloc(nxcells_ish*sizeof(double));
  int *nxcells_ij1           = NULL ; nxcells_ij1 = (int *)malloc(input_grid_ncells*sizeof(int));

  printf("herehere %d %d\n", input_grid_ncells, output_grid_ncells);
  exit(0);
  
#pragma acc enter data create(input_parent_index_i[:nxcells_ish], output_parent_index_i[:nxcells_ish], \
                              xarea_i[:nxcells_ish], nxcells_ij1[:input_grid_ncells])
  
#pragma acc data present(output_grid_lon[:output_grid_npts],         \
                         output_grid_lat[:output_grid_npts],         \
                         input_grid_lon[:input_grid_npts],           \
                         input_grid_lat[:input_grid_npts],           \
                         output_cells_lon_min[:output_grid_ncells],  \
                         output_cells_lon_max[:output_grid_ncells],  \
                         output_cells_lat_min[:output_grid_ncells],  \
                         output_cells_lat_max[:output_grid_ncells],  \
                         output_cells_lon_cent[:output_grid_ncells], \
                         output_cells_nvertices[:output_grid_ncells],\
                         output_cells_lon_vertices[MAX_V*output_grid_ncells],\
                         output_cells_lat_vertices[MAX_V*output_grid_ncells],\
                         output_cells_area[:output_grid_ncells], \
                         approx_nxcells_ij1[:input_grid_ncells], \
                         ij2_start[:input_grid_ncells],          \
                         ij2_end[:input_grid_ncells],            \
                         mask_input_grid[:input_grid_ncells],    \
                         nxcells_ij1[:input_grid_ncells],        \
                         input_parent_index_i[:nxcells_ish],     \
                         output_parent_index_i[:nxcells_ish],    \
                         xarea_i[:nxcells_ish]) copyin(input_grid_ncells, output_grid_ncells) copy(nxcells)
#pragma acc parallel loop reduction(+:nxcells)
  for(int ij1=0; ij1<input_grid_ncells; ij1++) {
    if(mask_input_grid[ij1] > MASK_THRESH)  {

      double input_cell_lon_vertices[MV], input_cell_lat_vertices[MV];
      int cume_nxcells=0, ixcell=0;

      get_cell_vertices_gpu(ij1, nlon_input_cells, input_grid_lon, input_grid_lat,
                            input_cell_lon_vertices, input_cell_lat_vertices);
      double input_cell_lat_min = minval_double_gpu(4, input_cell_lat_vertices);
      double input_cell_lat_max = maxval_double_gpu(4, input_cell_lat_vertices);
      int nvertices1 = fix_lon_gpu(input_cell_lon_vertices, input_cell_lat_vertices, 4, M_PI);
      double input_cell_lon_min = minval_double_gpu(nvertices1, input_cell_lon_vertices);
      double input_cell_lon_max = maxval_double_gpu(nvertices1, input_cell_lon_vertices);
      double input_cell_lon_cent = avgval_double_gpu(nvertices1, input_cell_lon_vertices);
      double input_cell_area = poly_area_gpu(input_cell_lon_vertices, input_cell_lat_vertices, nvertices1);

#pragma acc loop seq
      for(int i=1; i<=ij1 ; i++) cume_nxcells += approx_nxcells_ij1[i-1];
      nxcells_ij1[ij1]=0;

#pragma acc loop seq reduction(+:ixcell)
      for(int ij2=ij2_start[ij1]; ij2<=ij2_end[ij1]; ij2++) {

        int nvertices2, xvertices=1;
        double dlon_cent, output_cell_lon_min_i, output_cell_lon_max_i, output_cell_area_i;
        double output_cell_lon_vertices_i[MAX_V], output_cell_lat_vertices_i[MAX_V];
        double xcell_lon_vertices[MV], xcell_lat_vertices[MV];

        double rotate=0.0;

        if(output_cells_lat_min[ij2] >= input_cell_lat_max) continue;
        if(output_cells_lat_max[ij2] <= input_cell_lat_min) continue;

        /* adjust according to input_grid_lon_cent*/
        output_cell_lon_min_i = output_cells_lon_min[ij2];
        output_cell_lon_max_i = output_cells_lon_max[ij2];
        nvertices2 = output_cells_nvertices[ij2];
        output_cell_area_i = output_cells_area[ij2];

        dlon_cent = output_cells_lon_cent[ij2] - input_cell_lon_cent;

        if(dlon_cent < -M_PI) rotate = TPI;
        if(dlon_cent > M_PI)  rotate = -TPI;

        output_cell_lon_min_i += rotate;
        output_cell_lon_max_i += rotate;
        for (int l=0; l<nvertices2; l++) {
          output_cell_lon_vertices_i[l] = output_cells_lon_vertices[ij2*MAX_V+l] + rotate;
          output_cell_lat_vertices_i[l] = output_cells_lat_vertices[ij2*MAX_V+l];
        }

        //output_cell_lon should in the same range as input_cell_lon after lon_fix,
        // so no need to consider cyclic condition
        if(output_cell_lon_min_i >= input_cell_lon_max ) continue;
        if(output_cell_lon_max_i <= input_cell_lon_min ) continue;
        if ( (xvertices = clip_2dx2d_gpu( input_cell_lon_vertices, input_cell_lat_vertices, nvertices1,
                                          output_cell_lon_vertices_i, output_cell_lat_vertices_i, nvertices2,
                                          xcell_lon_vertices, xcell_lat_vertices)) > 0 ){
          double xcell_area = poly_area_gpu(xcell_lon_vertices, xcell_lat_vertices, xvertices);
          if( xcell_area/min(input_cell_area, output_cell_area_i) > AREA_RATIO_THRESH ) {
            xarea_i[cume_nxcells+ixcell] = xcell_area;
            input_parent_index_i[cume_nxcells+ixcell]  = ij1;
            output_parent_index_i[cume_nxcells+ixcell] = ij2;
            ixcell++;
          }
        }
      }
      nxcells+=ixcell;
      nxcells_ij1[ij1]=ixcell;
    }
  }

  free(input_parent_index_i) ; input_parent_index_i = NULL;
  free(output_parent_index_i); output_parent_index_i = NULL;
  free(nxcells_ij1)          ; nxcells_ij1 = NULL;
  free(xarea_i)              ; xarea_i = NULL;

  return nxcells;

};/* get_xgrid_2Dx2D_order1 */

