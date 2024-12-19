!------------------------------------------------------------------------------
!>  @file logger.f90
!>  @author Michael Densberger and Carl Mears
!>  @brief A collection of functions to log messages to the screen and files within Fortran.
!>  @details 
!!  
!------------------------------------------------------------------------------

module logger
use, intrinsic :: iso_fortran_env, only: int32, int64
!use printer, only: printf
implicit none
    private
    
    character(len=*), parameter :: module_name='logger' !> Module name.
    logical :: initialized=.false. !> Init status.
    
    integer(int64), parameter  :: timestamp_length=22 !> Length of timestamp string used when logging.
    character(len=timestamp_length) :: timestamp 
    
    integer(int64) :: log_unit !> Log unit to be used for log file. Automatically set via new newunit during open call.
    logical :: log_to_screen, log_to_file !> Logical values set by the user.
    
    integer(int32), parameter :: debug_level=0 !> Debug messages.
    integer(int32), parameter :: info_level=1 !> Messages a user might want to see.
    integer(int32), parameter :: error_level=2 !> Errors that are recoverable.
    integer(int32), parameter :: crit_level=3 !> Errors that are not recoverable.
    character(len=*), parameter :: debug_str = "dbg"
    character(len=*), parameter :: info_str = "inf"
    character(len=*), parameter :: error_str = "err"
    character(len=*), parameter :: crit_str = "crt"
    character(len=*), dimension(4), parameter :: error_strs = [debug_str, info_str, error_str, crit_str]
    integer(int32) :: screen_filter_level !> Suppress log messages to the screen below this level.
    integer(int32) :: file_filter_level !> Suppress log messages to the file below this level.
    character(len=4096) :: log_msg_string   !> An empty string that users may modify for formatted log messages
    !> Ex: 
    !> integer :: example_int 
    !> write(log_msg_string,"(a,i0.3)") "My example int is: ", example_int
    !> call log_debug(log_msg_string)
    public debug_level, &
    info_level, &
    error_level, &
    crit_level, &
    log_msg_string, &
    initialize_logger, close_log, &
    set_filter, set_screen_filter, set_file_filter, &
    log_line, log_debug, log_info, log_error, log_crit
    
contains
    !---------------------------------------------------------------------------  
    !>  @brief Inititalize the logger.
    !>  @param[in], character(len=*), program_name: The name of the program and log file.
    !>  @param[in], logical, use_log_file: True to log to a file as well.
    !>  @param[in], character(len=*), program_name: The path to the log file.
    !>  @param[in], integer(int32), optional, log_level: The level to log at.
    !---------------------------------------------------------------------------  
    subroutine initialize_logger(program_name, use_log_file, path, log_level, clobber)
        character(len=*), intent(in) :: program_name
        logical, intent(in), optional :: use_log_file
        character(len=*), intent(in), optional :: path
        integer(int32), intent(in), optional :: log_level
        logical,intent(in), optional  :: clobber
        integer(int64) :: stat
        
        
        character(len=4096) :: full_path 
        logical :: file_exists  
        
        
        stat=0
        screen_filter_level=0 ! Default log all messages.
        file_filter_level=0
        log_to_screen=.true.
        
        if (present(log_level)) then
            screen_filter_level=log_level ! Use provided value if preset.
            file_filter_level=log_level
        end if
        
        if (present(use_log_file) .and. use_log_file .eqv. .true.) then
            log_to_file=.true.
        else
            log_to_file=.false.
        end if   

        if (log_to_file) then
            if (present(path)) then
                full_path=trim(path)
            else
                ! TODO: get cwd
                full_path=trim(program_name) // '.log'
            end if
            
            inquire(file=trim(full_path), exist=file_exists)
            
            if (file_exists .eqv. .true.) then
                if (present(clobber) .and. clobber .eqv. .true.) then
                    open(newunit=log_unit, file=full_path, status='replace', iostat=stat, action='write') 
                else
                    open(newunit=log_unit, file=full_path, status='old', iostat=stat, position='append', action='write') 
                endif
            else 
                open(newunit=log_unit, file=full_path, status='new', iostat=stat, action='write') 
            end if 
        end if
        
        initialized=.true.
    end subroutine initialize_logger
    
    !---------------------------------------------------------------------------  
    !>  @brief Set the internal timestamp string to the current time.
    !---------------------------------------------------------------------------  
    subroutine update_prefix()
        integer(int64), dimension(8) :: values
        
        call date_and_time(VALUES=values) ! Update with current time        
        write(timestamp,1000) values(1), '-', values(2), '-', values(3), ' ', values(5), ':', values(6), ':', values(7), ' - '
        1000 format (i4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a) 
    end subroutine update_prefix
    
    !---------------------------------------------------------------------------  
    !>  @brief Translate the current integer log level to a string log level.
    !---------------------------------------------------------------------------  
    subroutine get_str_level(level, str)
        integer(int32), optional, intent(in) :: level
        character(len=:), allocatable, intent(out) :: str
        
        allocate(character(len=len('@')+len(error_strs(level+1))+len('-')) :: str)
        str = '@' // error_strs(level+1) // ':'
    end subroutine get_str_level
    
    !---------------------------------------------------------------------------  
    !>  @brief Close the log.
    !---------------------------------------------------------------------------  
    subroutine close_log()
        close(log_unit)
    end subroutine
    
    !---------------------------------------------------------------------------  
    !>  @brief Set both the screen and the file filter at the same time.
    !---------------------------------------------------------------------------  
    subroutine set_filter(level)
        integer(int32) :: level
        
        screen_filter_level=level
        file_filter_level=level
    end subroutine

    !---------------------------------------------------------------------------  
    !>  @brief Set the screen filter.
    !---------------------------------------------------------------------------  
    subroutine set_screen_filter(level)
        integer(int32) :: level
        
        screen_filter_level=level
    end subroutine
    
    !---------------------------------------------------------------------------  
    !>  @brief Set the file filter.
    !---------------------------------------------------------------------------  
    subroutine set_file_filter(level)
        integer(int32) :: level
        
        file_filter_level=level
    end subroutine
    
    !---------------------------------------------------------------------------  
    !>  @brief Log a string at level: debug
    !---------------------------------------------------------------------------  
    subroutine log_debug(line)
        character(len=*), intent(in) :: line
        
        call log_line(line, debug_level)
    end subroutine log_debug
    
    !---------------------------------------------------------------------------  
    !>  @brief Log a string at level: info
    !---------------------------------------------------------------------------  
    subroutine log_info(line)
        character(len=*), intent(in) :: line
        
        call log_line(line, info_level)
    end subroutine log_info
    
    !---------------------------------------------------------------------------  
    !>  @brief Log a string at level: error
    !---------------------------------------------------------------------------  
    subroutine log_error(line)
        character(len=*), intent(in) :: line
        
        call log_line(line, error_level)
    end subroutine log_error
    
    !---------------------------------------------------------------------------  
    !>  @brief Log a string at level: critical
    !---------------------------------------------------------------------------  
    subroutine log_crit(line)
        character(len=*), intent(in) :: line
        
        call log_line(line, crit_level)
    end subroutine log_crit
    
    !---------------------------------------------------------------------------  
    !>  @brief Log a string.
    !---------------------------------------------------------------------------  
    subroutine log_line(line, level)
        character(len=*), intent(in) :: line
        integer(int32), optional, intent(in) :: level
        
        integer(int32) :: level_
        character(len=:), allocatable :: level_str
        character(len=:), allocatable :: line_to_log
        
        if (present(level)) then
            level_ = level
        else
            level_ = debug_level
        end if
        
        if (line .eq. '\n') then
            allocate(character(len=0) :: line_to_log)
            line_to_log=''
            timestamp='' 
        else 
            if (present(level)) then
                call get_str_level(level, level_str)
                allocate(character(len=len(level_str // line)) :: line_to_log)
                line_to_log=level_str // line
                call update_prefix()
            else 
                allocate(character(len=len(line)) :: line_to_log)
                line_to_log=line
                call update_prefix()
            end if
        end if
        
        if (initialized .eqv. .true.) then
            ! Write to console  if enabled
            !if ((log_to_screen .eqv. .true.) .and. (level_ .ge. screen_filter_level)) call printf(timestamp // trim(line_to_log))
            if (log_to_screen .eqv. .true. .and. level_ .ge. screen_filter_level) write(*,'(1x,a)') timestamp // trim(line_to_log)
            
            ! Write to log if enabled
            if ((log_to_file .eqv. .true.) .and. (level_ .ge. file_filter_level)) write(log_unit,'(1x,a)') timestamp // trim(line_to_log)
        else             
            ! Write to console only, log file not set.
            !call printf(timestamp // trim(line_to_log))
            write(*,'(1x,a)') timestamp // trim(line_to_log)
        end if
        
        deallocate(line_to_log)
    end subroutine
end module logger