module CF_metadata

    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite

    use io_nc, only: handle_nc_err, minmax_iso8601
    use netcdf
    implicit none

    private
    public CF_1_8_attr,init_CF_1_8_attr,put_CF_1_8_attr
    
    type CF_1_8_attr
        real(real32),dimension(2)                :: actual_range
        real(real32)                             :: add_offset
        character(len=:),allocatable             :: ancillary_variables
        character(len=:),allocatable             :: axis
        character(len=:),allocatable             :: bounds
        character(len=:),allocatable             :: calendar
        character(len=:),allocatable             :: cell_measures
        character(len=:),allocatable             :: cell_methods
        character(len=:),allocatable             :: cf_role
        character(len=:),allocatable             :: climatology
        character(len=:),allocatable             :: comment
        character(len=:),allocatable             :: compress
        character(len=:),allocatable             :: computed_standard_name
        character(len=:),allocatable             :: Conventions
        character(len=:),allocatable             :: coordinates
        character(len=:),allocatable             :: external_variables
        real(real32)                             :: FillValue
        character(len=:),allocatable             :: feature_type
        integer(int32)                           :: flag_mask
        character(len=:),allocatable             :: flag_meanings
        integer(int32)                           :: flag_values
        character(len=:),allocatable             :: formula_terms
        character(len=:),allocatable             :: geometry
        character(len=:),allocatable             :: geometry_type
        character(len=:),allocatable             :: grid_mapping
        character(len=:),allocatable             :: history
        character(len=:),allocatable             :: instance_dimension
        character(len=:),allocatable             :: institution
        character(len=:),allocatable             :: interior_ring
        integer(int32)                           :: leap_month
        integer(int32)                           :: leap_year
        character(len=:),allocatable             :: long_name
        real(real32)                             :: missing_value
        integer(int32),dimension(12)             :: month_lengths
        character(len=:),allocatable             :: node_coordiantes
        character(len=:),allocatable             :: node_count
        character(len=:),allocatable             :: nodes
        character(len=:),allocatable             :: part_node_count
        character(len=:),allocatable             :: positive
        character(len=:),allocatable             :: references
        character(len=:),allocatable             :: sample_dimension
        real(real32)                             :: scale_factor
        character(len=:),allocatable             :: source
        real(real32)                             :: standard_error_multiplier
        character(len=:),allocatable             :: standard_name
        character(len=:),allocatable             :: title
        character(len=:),allocatable             :: units
        real(real32)                             :: valid_max
        real(real32)                             :: valid_min
        real(real32),dimension(2)                :: valid_range
    end type CF_1_8_attr


contains

    subroutine init_CF_1_8_attr(a)
        
        type(CF_1_8_attr),intent(inout)   :: a 

        real(real32) :: nan_f32
        !real(real64) :: nan_f64

        nan_f32 = ieee_value(0.0_real32, ieee_quiet_nan)    
        !nan_f64 = ieee_value(0.0_real64, ieee_quiet_nan)    

        a%actual_range = nan_f32
        a%add_offset = nan_f32
        a%ancillary_variables = ''
        a%axis = ''
        a%bounds = ''
        a%calendar = ''
        a%cell_measures = ''
        a%cell_methods = ''
        a%cf_role = ''
        a%climatology = ''
        a%comment = ''
        a%compress = ''
        a%computed_standard_name = ''
        a%Conventions = ''
        a%coordinates = ''
        a%external_variables = ''
        a%FillValue = nan_f32
        a%feature_type = ''
        a%flag_mask = 0
        a%flag_meanings = ''
        a%flag_values = 0
        a%formula_terms = ''
        a%geometry = ''
        a%geometry_type = ''
        a%grid_mapping = ''
        a%history = ''
        a%instance_dimension = ''
        a%institution = ''
        a%interior_ring = ''
        a%leap_month = 0
        a%leap_year = 0
        a%long_name = ''
        a%missing_value = nan_f32
        a%month_lengths = 0
        a%node_coordiantes = ''
        a%node_count = ''
        a%nodes = ''
        a%part_node_count = ''
        a%positive = ''
        a%references = ''
        a%sample_dimension = ''
        a%scale_factor = nan_f32
        a%source = ''
        a%standard_error_multiplier = nan_f32
        a%standard_name = ''
        a%title = ''
        a%units = ''
        a%valid_max = nan_f32
        a%valid_min = nan_f32
        a%valid_range = nan_f32

    end subroutine init_CF_1_8_attr

    subroutine put_CF_1_8_attr(a,ncid,varid)
        type(CF_1_8_attr),intent(in)   :: a 
        integer :: ncid
        integer,optional :: varid

        if (.not. present(varid))then
            varid = NF90_GLOBAL
        endif

        if ((ieee_is_finite(a%actual_range(1))) .and.(ieee_is_finite(a%actual_range(2)))) then
            print *,'actual range not working yet...'
        endif

        if (ieee_is_finite(a%add_offset)) then
            call handle_nc_err(nf90_put_att(ncid, varid, "add_offset", a%add_offset))
        endif

        if (len(a%ancillary_variables) > 0) then
            call handle_nc_err(nf90_put_att(ncid, varid, "ancillary_variables", a%ancillary_variables))
        endif

        ! a%axis = ''
        ! a%bounds = ''
        ! a%calendar = ''
        ! a%cell_measures = ''
        ! a%cell_methods = ''
        ! a%cf_role = ''
        ! a%climatology = ''
        ! a%comment = ''
        ! a%compress = ''
        ! a%computed_standard_name = ''
        ! a%Conventions = ''
        ! a%coordinates = ''
        ! a%external_variables = ''
        if (ieee_is_finite(a%FillValue)) then
            print *, ncid
            print *, varid
            print *, a%FillValue

            call handle_nc_err(nf90_put_att(ncid, varid, "FillValue", a%FillValue))
        endif
        ! a%feature_type = ''
        ! a%flag_mask = 0
        ! a%flag_meanings = ''
        ! a%flag_values = 0
        ! a%formula_terms = ''
        ! a%geometry = ''
        ! a%geometry_type = ''
        ! a%grid_mapping = ''
        ! a%history = ''
        ! a%instance_dimension = ''
        ! a%institution = ''
        ! a%interior_ring = ''
        ! a%leap_month = 0
        ! a%leap_year = 0
        if (len(a%long_name) > 0) then
            call handle_nc_err(nf90_put_att(ncid, varid, "long_name", a%long_name))
        endif
        if (ieee_is_finite(a%missing_value)) then
            call handle_nc_err(nf90_put_att(ncid, varid, "missing_value", a%missing_value))
        endif


    ! a%month_lengths = 0
    ! a%node_coordiantes = ''
    ! a%node_count = ''
    ! a%nodes = ''
    ! a%part_node_count = ''
    ! a%positive = ''
    ! a%references = ''
    ! a%sample_dimension = ''
        if (ieee_is_finite(a%scale_factor)) then
            call handle_nc_err(nf90_put_att(ncid, varid, "scale_factor", a%scale_factor))
        endif
        if (len(a%standard_name) > 0) then
            call handle_nc_err(nf90_put_att(ncid, varid,'standard_name',a%standard_name))
        endif
        if (len(a%title) > 0) then
            call handle_nc_err(nf90_put_att(ncid, varid,'title',a%title))
        endif
        if (len(a%units) > 0) then
            call handle_nc_err(nf90_put_att(ncid, varid,'units',a%units))
        endif
        if (ieee_is_finite(a%valid_max)) then
            call handle_nc_err(nf90_put_att(ncid, varid, 'valid_max', a%valid_max))
        endif
        if (ieee_is_finite(a%valid_min)) then
            call handle_nc_err(nf90_put_att(ncid, varid, 'valid_min', a%valid_min))
        endif
        
    end subroutine put_CF_1_8_attr
end module CF_metadata
    

  