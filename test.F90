program main

  use openacc
  use iso_c_binding
  implicit none

  integer, parameter :: nlon=10, nlat=10
  integer :: n=nlon*nlat
  real(c_double), allocatable :: lon(:), lat(:)
  real(c_double), allocatable :: mask(:)
  real(c_double), allocatable :: approx_nxcells(:), ij2_start(:), ij2_end(:)
  
  real(c_double), allocatable :: lon_min(:), lon_max(:), lat_min(:), lat_max(:)  
  real(c_double), allocatable :: lon_cent(:), area(:), lon_vertices(:), lat_vertices(:)
  integer, allocatable :: nvertices(:)

  integer, allocatable :: input_parent_index(:)
  integer, allocatable :: output_parent_index(:)
  real(c_double), allocatable :: xarea(:)
  
  integer i, j, ij

  integer :: nxcells_ish, nxcells
  integer :: get_upbound_nxcells_2dx2d_gpu
  integer :: create_xgrid_2dx2d_order1_gpu
  
  allocate(lon((nlon+1)*(nlat+1)))
  allocate(lat((nlon+1)*(nlat+1)))
  allocate(lon_min(n))
  allocate(lon_max(n))
  allocate(lat_min(n))
  allocate(lat_max(n))
  allocate(lon_cent(n))
  allocate(area(n))
  allocate(nvertices(n))
  allocate(lon_vertices(8*n))
  allocate(lat_vertices(8*n))

  allocate(mask(n))
  allocate(approx_nxcells(n))
  allocate(ij2_start(n))
  allocate(ij2_end(n))
  
  ij=1
  do i=1, nlat+1
     do j=1, nlon+1
        lon(ij) = real(i*3.14*.001,c_double)
        lat(ij) = real(i*3.14*.001,c_double)
        ij = ij+1
     end do
  end do

  mask = 1.0

!$acc enter data copyin(lon(:(nlon+1)*(nlat+1)), lat(:(nlon+1)*(nlat+1)), mask(:n))
!$acc enter data create(approx_nxcells(:n), ij2_start(:n), ij2_end(:n))  
  call get_grid_cell_struct_gpu(nlon, nlat, lon, lat, lon_min, lon_max, lat_min, lat_max, &
       lon_cent, area, nvertices, lon_vertices, lat_vertices)
  
  nxcells_ish = get_upbound_nxcells_2dx2d_gpu(nlon, nlat, nlon, nlat, lon, lat, lon, lat, &
       lon_min, lon_max, lat_min, lat_max, lon_cent, mask, approx_nxcells, ij2_start, ij2_end)

  allocate(input_parent_index(nxcells_ish))
  allocate(output_parent_index(nxcells_ish))
  allocate(xarea(nxcells_ish))

  write(*,*) 'here', nxcells_ish
  
!$acc enter data create(input_parent_index(:nxcells_ish), output_parent_index(:nxcells_ish), xarea(:nxcells_ish))
  
  nxcells = create_xgrid_2dx2d_order1_gpu(nlon, nlat, nlon, nlat, lon, lat, lon, lat, nxcells_ish, &
       lon_min, lon_max, lat_min, lat_max, lon_cent, area, nvertices, lon_vertices, lat_vertices, &
       mask, approx_nxcells, ij2_start, ij2_end, input_parent_index, output_parent_index, xarea)  
  
end program main
