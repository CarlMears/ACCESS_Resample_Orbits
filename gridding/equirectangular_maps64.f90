module equirectangular_maps64
    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use NAN_support, only: nan_f32, nan_f64
    
    use equirectangular_map_constants

    implicit none

    type map_str64
        real(real32)                            :: min_x         ! minimum x value
        real(real32)                            :: max_x         ! maximum x value
        real(real32)                            :: delta_x       ! spacing points in x direction
        integer(int32)                          :: num_x         ! number of x points
        integer(int32)                          :: type_x        ! type of grid in x direction (see explantion below)
        real(real32)                            :: min_y         ! minimum y value
        real(real32)                            :: max_y         ! maximum x value
        real(real32)                            :: delta_y         ! spacing points in x direction
        integer(int32)                          :: num_y            ! number of x points
        integer(int32)                          :: type_y         ! type of grid in x direction (see explantion below)
        integer(int32)                          :: type_bc         ! type of boundary conditions
        real(real64),dimension(:,:),allocatable :: dat             ! map data
        real(real32),dimension(:,:),allocatable :: lons
        real(real32),dimension(:,:),allocatable :: lats
    end type map_str64

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
        

contains
    
    subroutine set_up_map64(min_x,delta_x,num_x,type_x,min_y,delta_y,num_y,type_y,type_bc,mp,error)

        use NAN_support

        implicit none

        real(real32),intent(IN)       :: min_x
        real(real32),intent(IN)       :: delta_x
        integer(int32),intent(IN)     :: num_x
        integer(int32),intent(IN)     :: type_x
        real(real32),intent(IN)       :: min_y
        real(real32),intent(IN)       :: delta_y
        integer(int32),intent(IN)     :: num_y 
        integer(int32),intent(IN)     :: type_y
        integer(int32),intent(IN)     :: type_bc
        type(map_str64),intent(INOUT)    :: mp
        integer(int32),intent(OUT)    :: error
        
        integer(int32)                :: ilon,ilat

        if (ieee_is_finite(min_x)) then
            mp%min_x = min_x
        else
            mp%min_x = nan_f32
            error = 1
            return
        endif

        if ((ieee_is_finite(delta_x)) .and.(delta_x .gt. 0.0)) then
            mp%delta_x = delta_x
        else
            mp%delta_x = nan_f32
            error = 1
            return
        endif

        if (num_x .gt. 0) then
            mp%num_x = num_x
        else
            mp%num_x = 0
            error = 1
            return
        endif

        if (type_x .eq. CELL_CENTERED) then
            mp%type_x = type_x
        else
            if (type_x .eq. EDGE_CENTERED) then
                mp%type_x = type_x
            else
                mp%type_x = 0
                error = 1
                return
            endif
        endif


        if (ieee_is_finite(min_y)) then
            mp%min_y = min_y
        else
            mp%min_y = nan_f32
            error = 1
            return
        endif

        if ((ieee_is_finite(delta_y)) .and.(delta_y .gt. 0.0)) then
            mp%delta_y = delta_y
        else
            mp%delta_y = nan_f32
            error = 1
            return
        endif

        if (num_y .gt. 0) then
            mp%num_y = num_y
        else
            mp%num_y = 0
            error = 1
            return
        endif

        if (type_y .eq. CELL_CENTERED) then
            mp%type_y = type_y
        else
            if (type_y .eq. EDGE_CENTERED) then
                mp%type_y = type_y
            else
                mp%type_y = 0
                error = 1
                return
            endif
        endif

        select case(type_bc)    
            case(SPHERICAL)
                mp%type_bc = SPHERICAL
            case(RECTANGULAR)
                mp%type_bc = RECTANGULAR
            case default
                mp%type_bc = 0
                error = 1
                return
        end select

        
        if (mp%type_bc .eq. SPHERICAL) then 

            ! make sure things make sense for spherical boundaries

            if (modulo(mp%num_x,2) .eq. 1) then
                print *,'Error, must have an even number of x value for spherical boundaries'
                Error = -4
                return
            endif

            if (mp%min_x .ne. 0.0) then
                print *,'Error, must have x_min = 0 for spherical boundaries'
                Error = -1
                return
            endif

            if (abs(mp%delta_x*mp%num_x - 360.0) .gt. SMALL_NUMBER) then
                print *,'Error, delta_x and num_x doesnt make sense for spherical boundaries'
                Error = -8
                return
            endif

            if (abs(mp%min_y+90.0) .gt. SMALL_NUMBER) then
                print *,'Error, must have y_min = -90.0 for spherical boundaries '
                Error = -1
                return
            endif
            select case(type_y)
                case (CELL_CENTERED)
                    if (abs(mp%delta_y*mp%num_y-180.0) .gt. SMALL_NUMBER) then
                        print *,'Error, delta_y and num_y make no sense for spherical boundaries'
                        error = -8
                        return
                    endif
                case (EDGE_CENTERED)
                    if (abs(mp%delta_y*(mp%num_y-1)-180.0) .gt. SMALL_NUMBER) then
                        print *,'Error, delta_y and num_y make no sense for spherical boundaries'
                        error = -8
                        return
                    endif
            end select
            

        endif


        select case(mp%type_bc)
            case(SPHERICAL)
                mp%max_x = 360.0 
                mp%max_y = 90.0
            case(RECTANGULAR)
                select case(mp%type_x)
                    case(EDGE_CENTERED)
                        mp%max_x = min_x+delta_x*(num_x -1)
                    case(CELL_CENTERED)
                        mp%max_x = min_x+delta_x*num_x
                    case default
                        mp%max_x = nan_f32
                end select

                select case(mp%type_y)
                    case(EDGE_CENTERED)
                        mp%max_y = min_y +delta_y*(num_y-1)
                    case(CELL_CENTERED)
                        mp%max_y = min_y + delta_y*num_y
                    case default
                        mp%max_y = nan_f32
                end select
            case default
                mp%max_x = nan_f32
                mp%max_y = nan_f32
        end select

        ! allocate space for the map

        allocate(mp%dat(mp%num_x,mp%num_y),STAT = error)
        allocate(mp%lons(mp%num_x,mp%num_y),STAT = error)
        allocate(mp%lats(mp%num_x,mp%num_y),STAT = error)
 
        ! fill the lat and lon arrays

        do ilon = 1,mp%num_x
            mp%lons(ilon,:) = mp%min_x + (ilon-1)*mp%delta_x
        enddo

        do ilat = 1,mp%num_y
            mp%lats(:,ilat) = mp%min_y + (ilat-1)*mp%delta_y
        enddo

        return

    end subroutine set_up_map64

    subroutine zero_map64(mp,error)

        implicit none

        type(map_str64),intent(INOUT)        :: mp
        integer(int32),intent(INOUT)            :: error

        mp%dat = 0.0D0
        error = ior(error,0)

    end subroutine zero_map64

    subroutine deallocate_map64(mp,error)
        
        implicit none

        type(map_str64),intent(INOUT)        :: mp
        integer(int32),intent(OUT)         :: error

        mp%min_x    = nan_f32
        mp%max_x    = nan_f32
        mp%delta_x  = nan_f32
        mp%num_x    = 0
        mp%type_x   = 0
        mp%min_y    = nan_f32
        mp%max_y    = nan_f32
        mp%delta_y    = nan_f32
        mp%num_y       = 0
        mp%type_y   = 0
        mp%type_bc  = 0

        deallocate(mp%dat,STAT = error)
        deallocate(mp%lons,STAT = error)
        deallocate(mp%lats,STAT = error)

        return

    end subroutine deallocate_map64

    subroutine get_map_value64(mp,xx,yy,value,error)

        ! this subroutine returns the interpolated value at a point x,y
        ! note that the interpolation is somewhat different depending on 
        ! what type of grid is used for the map structure

        implicit none

        type(map_str64),intent(IN)        :: mp       ! map structure to be interpolated
        real(real32),intent(IN)                :: xx        ! longitude
        real(real32),intent(IN)                :: yy        ! latitude
        real(real64),intent(OUT)                :: value    ! the interpolated value
        integer(int32),intent(OUT)            :: error

        integer(int32)                        :: ix_left,ix_right        ! map indices of points to the right (left) of desired point
        integer(int32)                        :: iy_lower,iy_upper    ! map indices of points above (below) the desired point
        real(real32)                            :: dx,dy                ! distance of the desired point to the
                                                                ! right (above) the lower left grid point

        real(real32)                            :: wt_left,wt_right
        real(real32)                            :: wt_lower,wt_upper
        real(real64)                            :: tot,tot_wt,wt
        real(real32)                            :: x,y,v

        error = 0

        ! modulo x and y to be within bounds

        if ((.not. ieee_is_finite(xx)) .or. (.not. ieee_is_finite(yy)) .or. &
            (mp%num_x .eq. 0) .or. (mp%num_y .eq. 0)) then
            value = nan_f64
            return
        endif

        x = modulo(xx-mp%min_x,mp%max_x - mp%min_x) + mp%min_x
        y = modulo(yy-mp%min_y,mp%max_y - mp%min_y) + mp%min_y

        ! find x indices and distance
        select case (mp%type_x)

            case (CELL_CENTERED) 

                ix_left        = floor((x-mp%min_x)/mp%delta_x - 0.5) + 1
                ix_right    = ix_left + 1
                dx            = x - mp%min_x - ((ix_left - 0.5) * mp%delta_x)

            case (EDGE_CENTERED) 
                ix_left  = floor((x-mp%min_x)/mp%delta_x) + 1
                ix_right = ix_left + 1
                dx       = x - mp%min_x - ((ix_left - 1) * mp%delta_x)         
            case default
                error = 1
                value = nan_f32
                return
        end select

        ! find y indices and distance

        select case(mp%type_y)

            case (CELL_CENTERED)
                
                iy_lower    = floor((y-mp%min_y)/mp%delta_y - 0.5) + 1
                iy_upper    = iy_lower + 1
                dy            = y - mp%min_y - ((iy_lower - 0.5) * mp%delta_y)
            case (EDGE_CENTERED)
                iy_lower    = floor((y-mp%min_y)/mp%delta_y) + 1
                iy_upper    = iy_lower + 1
                dy            = y - mp%min_y - ((iy_lower - 1) * mp%delta_y)
            case default
                error = 1
                return
        end select

        ! now perform the interpolation

        wt_right = dx/mp%delta_x
        wt_left  = 1.0 - wt_right

        wt_upper = dy/mp%delta_y
        wt_lower = 1.0 - wt_upper

        tot = 0.0 
        tot_wt = 0.0

        v = map_value64(mp,ix_left,iy_lower,error) 
        if (ieee_is_finite(v)) then
            wt = wt_left * wt_lower
            tot = tot +( v * wt)
            tot_wt = tot_wt+wt
        endif 

        v = map_value64(mp,ix_left,iy_upper,error) 
        if (ieee_is_finite(v)) then
            wt = wt_left * wt_upper
            tot = tot + (v * wt)
            tot_wt = tot_wt+wt
        endif

        v = map_value64(mp,ix_right,iy_upper,error) 
        if (ieee_is_finite(v)) then
            wt = wt_right*wt_upper
            tot = tot + (v * wt)
            tot_wt = tot_wt + wt
        endif
        
        v = map_value64(mp,ix_right,iy_lower,error) 
        if (ieee_is_finite(v)) then
            wt = wt_right*wt_lower 
            tot = tot + (v * wt)
            tot_wt = tot_wt+wt
        endif 

        value = tot/tot_wt

        return

    end subroutine get_map_value64

    subroutine get_map_value_no_interp64(mp,xx,yy,value,error)

        ! this subroutine returns the non interpolated value closest to a point x,y

        implicit none

        type(map_str64),intent(IN)        :: mp       ! map structure to be interpolated
        real(real32),intent(IN)                :: xx        ! longitude
        real(real32),intent(IN)                :: yy        ! latitude
        real(real64),intent(OUT)                :: value    ! the interpolated value
        integer(int32),intent(OUT)            :: error

        integer(int32)                        :: ix                    ! map index of point closest to the desired point
        integer(int32)                        :: iy                    ! map index of point closest to the desired point

        real(real32)                            :: x,y
        error = 0


        if ((.not. ieee_is_finite(xx)) .or. (.not. ieee_is_finite(yy)) .or. &
            (mp%num_x .eq. 0) .or. (mp%num_y .eq. 0)) then
            value = nan_f32
            return
        endif

        ! modulo x and y to be within bounds

        x = modulo(xx-mp%min_x,mp%max_x - mp%min_x) + mp%min_x
        y = modulo(yy-mp%min_y,mp%max_y - mp%min_y) + mp%min_y

        ! find x indices and distance
        select case (mp%type_x)

            case (CELL_CENTERED) 
                ix        = floor((x-mp%min_x)/mp%delta_x) + 1
            case (EDGE_CENTERED) 
                ix  = floor((x-mp%min_x)/mp%delta_x + 0.5) + 1
            case default
                error = 1
                value = nan_f32
                return
        end select

        ! find y indices and distance

        select case(mp%type_y)

            case (CELL_CENTERED)
                iy    = floor((y-mp%min_y)/mp%delta_y) + 1
            case (EDGE_CENTERED)
                iy    = floor((y-mp%min_y)/mp%delta_y + 0.5) + 1
            case default
                error = 1
                return
        end select

        ! now get the data

        value = map_value64(mp,ix,iy,error)
         
        return

    end subroutine get_map_value_no_interp64


    subroutine wrap_indices64(mp,ix,iy,ix_wrapped,iy_wrapped,error)

        implicit none

        type(map_str64),intent(IN)        :: mp
        integer(int32),intent(IN)            :: ix
        integer(int32),intent(IN)            :: iy
        integer(int32),intent(OUT)            :: ix_wrapped
        integer(int32),intent(OUT)            :: iy_wrapped
        integer(int32),intent(OUT)            :: error

        logical                            :: reversed

        select case (mp%type_bc)

            case(SPHERICAL)

                ! y values are reflected across poles

                iy_wrapped = modulo(((iy-1)+4*mp%num_y),(2*mp%num_y))+1
                reversed = .false.
                if (iy_wrapped .gt. mp%num_y) then  ! fold   
                      iy_wrapped = (2*mp%num_y)+1 - iy_wrapped
                      reversed = .true.
                endif

                if (reversed) then
                    ix_wrapped = ix + mp%num_x/2    ! reflect across pole --  modulo will be taken care of in next part
                else
                    ix_wrapped = ix
                endif


                ! x values are wrapped 

                ix_wrapped = modulo(ix_wrapped-1,mp%num_x) + 1

            case (RECTANGULAR)

                if ((ix .lt. 1) .or. (ix .gt. mp%num_x) .or. &
                    (iy .lt. 1) .or. (iy .gt. mp%num_y)) then
                   ix_wrapped = 0
                   iy_wrapped = 0     
                   error = 1
                else
                    ix_wrapped = ix
                    iy_wrapped = iy
                endif
                return

            case default
                print *,'Error: Map Boundary Conditions Undefined'
                error = -1
                ix_wrapped = 0
                iy_wrapped = 0     
        end select

    end subroutine wrap_indices64


    real(real64) function map_value64(mp,ix,iy,error)

        ! this subroutine returns the value of the map at xi,yi even
        ! if xi and yi are out of range by applying the appropriate
        ! boundary conditions

        implicit none

        type(map_str64),intent(IN):: mp
        integer(int32),intent(IN)    :: ix
        integer(int32),intent(IN)    :: iy
        integer(int32),intent(OUT)    :: error

        integer(int32)                :: ix_wrapped
        integer(int32)                :: iy_wrapped

        error = 0

        call wrap_indices64(mp,ix,iy,ix_wrapped,iy_wrapped,error) 

        if (error .eq. 0) then
            map_value64 = mp%dat(ix_wrapped,iy_wrapped)
        else
            map_value64 = nan_f64
        endif


        return

    end function map_value64

    subroutine get_grid_location64(mp,ix,iy,x,y,error)

        ! this subroutine returns the location in real map space of
        ! the grid point ix,iy

        implicit none

        type(map_str64),intent(IN)            :: mp        ! map structure
        integer(int32),intent(IN)                :: ix        ! x grid position
        integer(int32),intent(IN)                :: iy        ! y grid position
        real(real32),intent(OUT)                    :: x        ! x position in real space
        real(real32),intent(OUT)                    :: y        ! y position in real space
        integer(int32),intent(OUT)                :: error    ! error flag

        error = 0

        if ((ix .lt. 1) .or. (ix .gt. mp%num_x) .or. &
            (iy .lt.    1) .or. (iy .gt. mp%num_y)) then
           error = 1
           x = nan_f32
           y = nan_f32
           return
        endif

        select case (mp%type_x)

            case(CELL_CENTERED)
                x = (ix - 0.5)*mp%delta_x+mp%min_x
            case(EDGE_CENTERED)
                x = (ix - 1)*mp%delta_x+mp%min_x
            case default
                error = 2
                x = nan_f32
        end select

        select case (mp%type_y)
            
            case(CELL_CENTERED)
                y = (iy - 0.5)*mp%delta_y+mp%min_y
            case(EDGE_CENTERED)
                y = (iy - 1.0)*mp%delta_y+mp%min_y
            case default
                error = 4
                y = nan_f32
        end select

        return

    end subroutine get_grid_location64

    subroutine get_vertex_location64(mp,ix,iy,x,y,error)

        ! this subroutine returns the location in real map space of
        ! the vertex point ix,iy

        implicit none

        type(map_str64),intent(IN)            :: mp        ! map structure
        integer(int32),intent(IN)                :: ix        ! x grid position
        integer(int32),intent(IN)                :: iy        ! y grid position
        real(real32),intent(OUT)                    :: x        ! x position in real space
        real(real32),intent(OUT)                    :: y        ! y position in real space
        integer(int32),intent(OUT)                :: error    ! error flag

        error = 0

        if ((ix .lt. 0) .or. (ix .gt. mp%num_x) .or. &
            (iy .lt.    0) .or. (iy .gt. mp%num_y)) then
           error = 1
           x = nan_f32
           y = nan_f32
           return
        endif

        select case (mp%type_x)

            case(CELL_CENTERED)
                x = ix*mp%delta_x+mp%min_x
            case(EDGE_CENTERED)
                x = (ix - 0.5)*mp%delta_x+mp%min_x
            case default
                error = 2
                x = nan_f32
        end select

        select case (mp%type_y)
            
            case(CELL_CENTERED)
                y = iy*mp%delta_y+mp%min_y
            case(EDGE_CENTERED)
                if (iy .gt. 0) then
                    y = (iy - 0.5)*mp%delta_y+mp%min_y
                else
                    y = nan_f32
                    error = 1
                    return
                endif
            case default
                error = 4
                y = nan_f32
        end select

        return

    end subroutine get_vertex_location64



    subroutine get_grid_indices64(mp,x,y,ix,iy,error)

        ! returns the indices of the grid point in whose cell the
        ! point x,y is located

        implicit none

        type(map_str64),intent(IN)            :: mp   ! map structure
        real(real32),intent(IN)                    :: x
        real(real32),intent(IN)                    :: y
        integer(int32),intent(OUT)                :: ix
        integer(int32),intent(OUT)                :: iy
        integer(int32),intent(OUT)                :: error

        ! check ranges

        if ((x .lt. mp%min_x) .or. (x .gt. mp%max_x) .or. &
            (y .lt. mp%min_y) .or. (y .gt. mp%max_y)) then
            error = 1
            ix = 0
            iy = 0
            return
        endif

        ! find x index

        select case (mp%type_x)

            case(CELL_CENTERED)
                ix = floor((x-mp%min_x)/mp%delta_x)+1
                if (ix .eq. mp%num_x+1) ix = mp%num_x  ! this is really only for x = max_x
            case(EDGE_CENTERED)
                ix = floor((x-mp%min_x)/mp%delta_x+0.5)+1
                if (ix .eq. mp%num_x+1) then
                    select case (mp%type_bc)
                        case (SPHERICAL)
                            ix = 1
                        case (RECTANGULAR)
                            ix = mp%num_x
                            ! I don't think this should evere happen
                        case default
                            ix = 0
                            error = 1
                    end select
                endif
            case default
                error = 2
                ix = 0
        end select

        select case (mp%type_y)

            case(CELL_CENTERED)
                iy = floor((y-mp%min_y)/mp%delta_y)+1
                if (iy .eq. mp%num_y+1) iy = mp%num_y  ! this is really only for x = max_x
            case(EDGE_CENTERED)
                iy = floor((y-mp%min_y)/mp%delta_y+0.5)+1
                if (iy .eq. mp%num_y+1) then
                    iy = mp%num_y
                endif
            case default
                error = 2
                iy = 0
        end select

        return

    end subroutine get_grid_indices64

    subroutine get_vertex_indices64(mp,x,y,ix,iy,error)

        ! returns the indices of the vertex point nearest x,y.  The vertex point is
        ! where the edges of the cells cross

        implicit none

        type(map_str64),intent(IN)            :: mp   ! map structure
        real(real32),intent(IN)                    :: x
        real(real32),intent(IN)                    :: y
        integer(int32),intent(OUT)                :: ix
        integer(int32),intent(OUT)                :: iy
        integer(int32),intent(OUT)                :: error

        ! check ranges

        if ((x .lt. mp%min_x) .or. (x .gt. mp%max_x) .or. &
            (y .lt. mp%min_y) .or. (y .gt. mp%max_y)) then
            error = 1
            ix = 0
            iy = 0
            return
        endif

        ! find x vertex index

        select case (mp%type_x)

            case(CELL_CENTERED)
                ix = floor((x-mp%min_x)/mp%delta_x-0.5)+1
            case(EDGE_CENTERED)
                ix = floor((x-mp%min_x)/mp%delta_x)+1
            case default
                error = 2
                ix = 0
        end select

        select case (mp%type_y)

            case(CELL_CENTERED)
                iy = floor((y-mp%min_y)/mp%delta_y-0.5)+1
            case(EDGE_CENTERED)
                iy = floor((y-mp%min_y)/mp%delta_y)+1
            case default
                error = 2
                iy = 0
        end select

        return

    end subroutine get_vertex_indices64

    subroutine calc_relative_position64(mp,x,y,ix,iy,x_rel,y_rel,error)

        ! calculates the relative positon of the point x,y relative to the
        ! vertex ix,iy.  Since x, y is constrained to be within the
        ! bounds of the map, the modulo arithmatic is relatively simple,
        ! as long as x,y id within half a delta of the vertex.  This is
        ! assumed to be the case here

        implicit none

        type(map_str64),intent(IN)            :: mp   ! map structure
        real(real32),intent(IN)                    :: x
        real(real32),intent(IN)                    :: y
        integer(int32),intent(IN)                :: ix
        integer(int32),intent(IN)                :: iy
        real(real32),intent(OUT)                    :: x_rel
        real(real32),intent(OUT)                    :: y_rel
        integer(int32),intent(OUT)                :: error

        real(real32)                                :: xv,yv

        ! check ranges

        if ((x .lt. mp%min_x) .or. (x .gt. mp%max_x) .or. &
            (y .lt. mp%min_y) .or. (y .gt. mp%max_y)) then
            error = 1
            x_rel = nan_f32 
            y_rel = nan_f32
            return
        endif

        ! find location of vertex in real space

        call get_vertex_location64(mp,ix,iy,xv,yv,error)

        select case(mp%type_bc)

            case(SPHERICAL)
                x_rel = modulo(x-xv+180,360.0) - 180.0
                y_rel = y - yv

            case(RECTANGULAR)
                x_rel = x - xv
                y_rel = y - yv

            case default
                x_rel = nan_f32
                y_rel = nan_f32
        end select

        if ((x_rel .gt. mp%delta_x/2.0) .or. (abs(y_rel) .gt. mp%delta_y/2.0)) then
            print *,'Warning rel_x, or rel_y in calc_rel_position too big'
            print *,'Rel X = ',x_rel,' Rel Y = ',y_rel
        endif

        return

    end subroutine calc_relative_position64



    subroutine downsample_map64(oldmap,num_x,type_x,num_y,type_y,newmap,error)

        ! this subroutine downsamples a map object to a lower resolution

        type(map_str64),intent(IN)        ::    oldmap    ! original map to be downsampled
        integer(int32),intent(IN)            ::    num_x    ! number of x grid points in new map
        integer(int32),intent(IN)            ::  type_x    ! x grid type for new map
        integer(int32),intent(IN)            ::    num_y    ! number of y grid points in new map
        integer(int32),intent(IN)            ::  type_y    ! y grid type for new map
        type(map_str64),intent(OUT)        ::  newmap    ! new map structure, resampled to new
                                                    ! lower resolution
        integer(int32),intent(OUT)            ::  error   ! error flag

        real(real32)                            ::  new_dx
        real(real32)                            ::  new_dy
        real(real32)                            ::  x
        real(real32)                            ::  y
        real(real32)                            ::  x_rel
        real(real32)                            ::  y_rel
        integer(int32)                        ::  ix
        integer(int32)                        ::  iy
        integer(int32)                        ::    ix_new
        integer(int32)                        ::    iy_new
        integer(int32)                        ::    ix_bx
        integer(int32)                        ::    iy_bx
        integer(int32)                        ::  ix_wrapped
        integer(int32)                        ::  iy_wrapped
        real(real32),dimension(2,2)            ::    box_wts
        real(real32)                            ::  v

        real(real64),dimension(:,:),pointer    ::  map_sum
        real(real64),dimension(:,:),pointer    ::    map_sum_wts
    
        error = 0

        ! check new nums to make sure we're downsampling

        if ((num_x .lt. 1) .or.    (num_y .lt. 1) .or.  &
            (num_y .gt. oldmap%num_y) .or. &
            (num_x .gt. oldmap%num_x)) then
                error = 1
                return
        endif

        ! set up the new map structure

        select case(oldmap%type_bc)
            case(SPHERICAL)    
                select case(type_x)
                    case(EDGE_CENTERED)
                        new_dx = oldmap%delta_x/(real(num_x)/real(oldmap%num_x))
                    case(CELL_CENTERED)
                        new_dx = oldmap%delta_x/(real(num_x)/real(oldmap%num_x))
                end select

                select case(type_y)
                    case(EDGE_CENTERED)
                        select case(oldmap%type_y)
                            case(EDGE_CENTERED)
                                new_dy = oldmap%delta_y/(real(num_y-1)/real(oldmap%num_y-1))
                            case(CELL_CENTERED)
                                new_dy = oldmap%delta_y/(real(num_y-1)/real(oldmap%num_y))
                        end select
                    case(CELL_CENTERED)
                        select case(oldmap%type_y)
                            case(EDGE_CENTERED)
                                new_dy = oldmap%delta_y/(real(num_y)/real(oldmap%num_y-1))
                            case(CELL_CENTERED)
                                new_dy = oldmap%delta_y/(real(num_y)/real(oldmap%num_y))
                        end select
                end select
            case(RECTANGULAR)
                select case(type_x)
                    case(EDGE_CENTERED)
                        select case(oldmap%type_x)
                            case(EDGE_CENTERED)
                                new_dy = oldmap%delta_x/(real(num_x-1)/real(oldmap%num_x-1))
                            case(CELL_CENTERED)
                                new_dy = oldmap%delta_x/(real(num_x-1)/real(oldmap%num_x))
                        end select
                    case(CELL_CENTERED)
                        select case(oldmap%type_y)
                            case(EDGE_CENTERED)
                                new_dy = oldmap%delta_x/(real(num_x)/real(oldmap%num_x-1))
                            case(CELL_CENTERED)
                                new_dy = oldmap%delta_x/(real(num_x)/real(oldmap%num_x))
                        end select
                end select

                select case(type_y)
                    case(EDGE_CENTERED)
                        select case(oldmap%type_y)
                            case(EDGE_CENTERED)
                                new_dy = oldmap%delta_y/(real(num_y-1)/real(oldmap%num_y-1))
                            case(CELL_CENTERED)
                                new_dy = oldmap%delta_y/(real(num_y-1)/real(oldmap%num_y))
                        end select
                    case(CELL_CENTERED)
                        select case(oldmap%type_y)
                            case(EDGE_CENTERED)
                                new_dy = oldmap%delta_y/(real(num_y)/real(oldmap%num_y-1))
                            case(CELL_CENTERED)
                                new_dy = oldmap%delta_y/(real(num_y)/real(oldmap%num_y))
                        end select
                end select
        end select

        call set_up_map64(oldmap%min_x,new_dx,num_x,type_x, &
                        oldmap%min_y,new_dy,num_y,type_y, &
                        oldmap%type_bc, newmap,error)
        if (error .ne.  0) return

        ! allocate two temporary matrices 

        allocate(map_sum(num_x,num_y),STAT = error)
        if (error .ne. 0) return
        allocate(map_sum_wts(num_x,num_y),STAT = error)
        if (error .ne. 0) return


        map_sum     = 0.0D0
        map_sum_wts = 0.0D0



        ! perform downsampling

        ! loop over cells in old map, and place properly wieghted amounts into cell in new grid

        do ix = 1,oldmap%num_x
            do iy = 1,oldmap%num_y

                ! find position of cell in old map

                call get_grid_location64(oldmap,ix,iy,x,y,error)
                
                ! find which vertex in the new map the center of the old cell falls nearest

                call get_vertex_indices64(newmap,x,y,ix_new,iy_new,error)

                ! find position relative to vertex

                call calc_relative_position64(newmap,x,y,ix_new,iy_new,x_rel,y_rel,error)

                ! calculate area overlap of the old cell with the four quadrants centered at the
                ! vertex.

                call calc_box_weights64(x_rel,y_rel,oldmap%delta_x,oldmap%delta_y,box_wts)

            !    if (box_wts(2,2) .gt. 0.0) then
            !        print *,'fadsfadsfads'
            !    endif

                ! loop over the 4 cell adjacent to the vertex.  Ignore NAN values in old map.

                do ix_bx = 0,1
                    do iy_bx = 0,1
                        v = map_value64(oldmap,ix,iy,error)
                        if (ieee_is_finite(v)) then
                            
                            call wrap_indices64(newmap,ix_new+ix_bx,iy_new+iy_bx,ix_wrapped,iy_wrapped,error)

!                            write(6,333)x,y,ix_wrapped,iy_wrapped,v,box_wts(ix_bx+1,iy_bx+1)
333                            format(1x,2f9.2,i4,i4,2f9.2)

                            map_sum(ix_wrapped,iy_wrapped) =                        &
                                        map_sum(ix_wrapped,iy_wrapped) +            &
                                        (v * box_wts(ix_bx+1,iy_bx+1))

                             map_sum_wts(ix_wrapped,iy_wrapped) =                    &
                                        map_sum_wts(ix_wrapped,iy_wrapped) +        &
                                        box_wts(ix_bx+1,iy_bx+1)

                        
                        endif
                    enddo
                enddo

            enddo
        enddo

        ! loop over cell in new map, computing average value in each cell
        ! if a given cell has no wieght, set it to NAN


        do ix_new = 1,num_x
            do iy_new = 1,num_y

                if (map_sum_wts(ix_new,iy_new) .gt. 0.0) then
                    newmap%dat(ix_new,iy_new) =  map_sum(ix_new,iy_new)/map_sum_wts(ix_new,iy_new)
                else
                    newmap%dat(ix_new,iy_new) = nan_f32  
                endif
            enddo
        enddo

        ! release memory used for the temporary arrays
        

        deallocate(map_sum,STAT = error)
        deallocate(map_sum_wts,STAT = error)

        print *,error


        return

    end subroutine downsample_map64
    
    
    subroutine calc_box_weights64(x,y,len_x,len_y,box_wts)

        implicit none

        real(real32),intent(IN)                    :: x
        real(real32),intent(IN)                    :: y
        real(real32),intent(IN)                    :: len_x
        real(real32),intent(IN)                    :: len_y
        real(real32),dimension(2,2),intent(OUT)    :: box_wts

        real(real32)                                :: llx
        real(real32)                                :: lly
        real(real32)                                :: urx
        real(real32)                                :: ury

        real(real32)                                :: tot_wt



        llx = x - (len_x/2.0)  ! position of lower left corner (x coordinate)
        lly = y - (len_y/2.0)  ! position of upper right corner (y coordinate)

        urx = x + (len_x/2.0)  ! position of upper right corner (x coordinate)
        ury = y + (len_y/2.0)  ! position of upper right corner (y coordinate)

        ! calculate area of each rectangle -- trim rectangle at axes if partly outside of quandrant

        box_wts(1,1) = (min(llx,0.0) - min(urx,0.0)) * (min(lly,0.0) - min(ury,0.0))       !  lower left quadrant;  remember < is the minimum operator
        box_wts(1,2) = (min(llx,0.0) - min(urx,0.0)) * (max(lly,0.0) - max(ury,0.0))    !  upper left
        box_wts(2,2) = (max(llx,0.0) - max(urx,0.0)) * (max(lly,0.0) - max(ury,0.0))    !  upper right
        box_wts(2,1) = (max(llx,0.0) - max(urx,0.0)) * (min(lly,0.0) - min(ury,0.0))    !  lower right

        ! normalize to total weight = 1.0

        tot_wt = box_wts(1,1)+box_wts(1,2)+box_wts(2,2)+box_wts(2,1)
        box_wts = box_wts/tot_wt

        return
    end subroutine calc_box_weights64

    subroutine write_map_netcdf64(mp,nc_filename,error,title)
      
        use io_nc, only: handle_nc_err, minmax_iso8601
        use netcdf

        implicit none

        type(map_str64),intent(IN)            :: mp
        character(*),intent(IN)             :: nc_filename
        integer(int32),intent(OUT)          :: error
        character(*),intent(in),optional    :: title

        integer(int32)  :: ncid,varid,dim_lon,dim_lat

        error = 0
        call handle_nc_err(nf90_create(nc_filename, ior(NF90_CLOBBER, NF90_NETCDF4), ncid))
        
        ! Define global attributes
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-1.8,ACDD-1.3"))
        if (present(title)) call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "title", "AMSR2"))
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "institution", "REMSS"))
        !call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "date_created", timestamp))
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_name", "Remote Sensing Systems"))
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_email", "support@remss.com"))
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_url", "http://www.remss.com"))
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lat_min", mp%min_y))
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lat_max", mp%max_y))
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lon_min", mp%min_x))
        call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lon_max", mp%max_x))

        call handle_nc_err(nf90_def_dim(ncid, "longitude", mp.num_x, dim_lon))
        call handle_nc_err(nf90_def_dim(ncid, "latitude", mp.num_y, dim_lat))

        call handle_nc_err(nf90_def_var(ncid, "longitude", NF90_FLOAT,dim_lon, varid))
        call handle_nc_err(nf90_put_var(ncid, varid, mp%lons(:,1)))

        call handle_nc_err(nf90_def_var(ncid, "latitude", NF90_FLOAT, dim_lat, varid))
        call handle_nc_err(nf90_put_var(ncid, varid, mp%lats(1,:)))

        call handle_nc_err(nf90_def_var(ncid, "Data", NF90_DOUBLE, [dim_lon,dim_lat], varid, &
                            deflate_level=2, shuffle=.true.))
        call handle_nc_err(nf90_put_var(ncid, varid, mp%dat))

        call handle_nc_err(nf90_inq_varid(ncid, "latitude", varid))
        call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "latitude"))
        call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_north"))

        call handle_nc_err(nf90_inq_varid(ncid, "longitude", varid))
        call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "longitude"))
        call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_east"))

        call handle_nc_err(nf90_close(ncid))

    end subroutine write_map_netcdf64

    subroutine write_map64(mp,unit_num,error)

        implicit none

        type(map_str64),intent(IN)            :: mp
        integer(int32),intent(IN)          :: unit_num
        integer(int32),intent(OUT)         :: error

        error = 0



        write(unit_num) mp%min_x
        write(unit_num) mp%max_x
        write(unit_num) mp%delta_x
        write(unit_num) mp%num_x
        write(unit_num) mp%min_y
        write(unit_num) mp%max_y
        write(unit_num) mp%delta_y
        write(unit_num) mp%num_y
        write(unit_num) mp%dat

        return
    end subroutine write_map64

    subroutine read_map64(unit_num,mp,error)

        implicit none

        type(map_str64),intent(OUT)            :: mp
        integer(int32),intent(IN)                :: unit_num
        integer(int32),intent(OUT)                :: error

        real(real32)                                :: min_x
        real(real32)                                :: max_x
        real(real32)                                :: delta_x
        integer(int32)                            :: num_x
        real(real32)                                :: min_y
        real(real32)                                :: max_y
        real(real32)                                :: delta_y
        integer(int32)                            :: num_y
        integer(int32)                            :: type_x
        integer(int32)                            :: type_y

        ! try to figure out axis type based on limits

        read(unit_num) min_x
        read(unit_num) max_x
        read(unit_num) delta_x
        read(unit_num) num_x
        read(unit_num) min_y
        read(unit_num) max_y
        read(unit_num) delta_y
        read(unit_num) num_y


        if (abs(max_x - (min_x+delta_x*num_x)) .lt. SMALL_NUMBER) then
            type_x = CELL_CENTERED
        else
            if (abs(max_x - (min_x+delta_x*(num_x-1))) .lt. SMALL_NUMBER) then
                type_x = EDGE_CENTERED
            else
                print *,'Error -- Unable to determine x axis type'
                type_x = 0
                error = 1
                return
            endif
        endif

        if (abs(max_y - (min_y+delta_y*num_y)) .lt. SMALL_NUMBER) then
            type_y = CELL_CENTERED
        else
            if (abs(max_y - (min_y+delta_y*(num_y-1))) .lt. SMALL_NUMBER) then
                type_y = EDGE_CENTERED
                type_x = EDGE_CENTERED
            else
                print *,'Error -- Unable to determine y axis type'
                type_y = 0
                error = 1
                return
            endif
        endif

        error = 0



        call set_up_map64(min_x,delta_x,num_x,type_x,min_y,delta_y,num_y,type_y,SPHERICAL,mp,error)

        read(unit_num) mp%dat

        return
    end subroutine read_map64

end module equirectangular_maps64







