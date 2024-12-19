module access_py_interfaces
    use, intrinsic :: iso_c_binding, only: c_int, c_bool
    use python_types_rewrite, only: iop_type

    implicit none
    private

    type, bind(C) :: py_process_resample_orbit
        integer(c_int) :: iorbit
        integer(c_int) :: footprint_size
        integer(c_int), dimension(20) :: channel_list
        logical(c_bool) :: write_swath_tbs
        logical(c_bool) :: write_swath_geoloc
        logical(c_bool) :: overwrite_existing_files
        type(iop_type) :: pathname_tb_output
    end type py_process_resample_orbit

    public :: py_process_resample_orbit
end module access_py_interfaces