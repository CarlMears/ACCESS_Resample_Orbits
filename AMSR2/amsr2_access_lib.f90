module amsr2_lib
    use, intrinsic :: iso_c_binding, only: c_int, c_bool
    use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    implicit none
    private
contains
    subroutine process_and_resample_amsr2_l2b_orbit_PY_HOOK(py_args) bind(C, name="process_and_resample_amsr2_l2b_orbit_PY_HOOK")
        !DEC$ ATTRIBUTES DLLEXPORT :: process_and_resample_amsr2_l2b_orbit_PY_HOOK
        use, intrinsic :: iso_c_binding, only: c_f_pointer
        use python_types_rewrite, only: python_type
        use access_py_interfaces, only: py_process_resample_orbit
        use amsr2_access_entry, only: process_and_resample_amsr2_l2b_orbit

        type(py_process_resample_orbit), intent(in) :: py_args

        ! write(*,*) "Would call with filename output: ", python_type%string(py_args%pathname_tb_output)
        ! return
        call process_and_resample_amsr2_l2b_orbit( &
            py_args%iorbit, py_args%channel_list, &
            py_args%write_swath_tbs, py_args%write_swath_geoloc, py_args%overwrite_existing_files, &
            python_type%string(py_args%pathname_tb_output) &
        )
    end subroutine process_and_resample_amsr2_l2b_orbit_PY_HOOK
end module amsr2_lib