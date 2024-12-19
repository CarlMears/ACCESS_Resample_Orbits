   
    module equirectangular_map_constants
        use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64

        integer(int32),parameter                :: CELL_CENTERED = 1
        integer(int32),parameter                :: EDGE_CENTERED = 2

        integer(int32),parameter                :: SPHERICAL     = 1
        integer(int32),parameter                :: RECTANGULAR   = 2
        integer(int32),parameter                :: CYLINDRICAL   = 3 ! not implemented yet

        real(real32),parameter                  :: SMALL_NUMBER  = 1.0e-3



    ! explanation of the grid type flag:
    !
    !  There are basically two types of evenly spaced rectangular gridded map --
    !
    !        1. Those for which the value represents the value in a cell - I'll call this cell_centered
    !              for these maps, the lower left value represents the value in the cell defined
    !             by long = (0.0:delta_x), lat = (-90,-90+delta_y), and is associated with at point at
    !              (delta_x/2.0, -90+delta_x/2.0)
    !        2. Those for which the value represents the value on the vertex between cells - I'll call this edge_centered
    !              for these maps, the lower left value represents the value at the point 
    !              by long = 0.0, lat = (-90).  A peculiar feature of this sort of grid is
    !              that the row at -90.0 must be identical, since it all refers to the same point!
    !        For now, all routines in the module assume that maps obey spherical boundary conditions
    !        For now, all maps are not redundant at the vertical edges
    !        
    !        This means that for both types of grids, 360.0/(delta_x) longitude points are needed, and the only different is
    !        whether the values are offset by 0.5*delta_x
    !
    !        For cell_centered grids, 180/(delta_y) latitude points are needed, and for 
    !        Edge centered grids 180/(delta_y) + 1 latitude points are needed
    !
    !        note that mixed grid types are allowed, with the lat and long grids being different types

    end module equirectangular_map_constants