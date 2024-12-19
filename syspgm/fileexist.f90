!September 15, 2006, 05:06:37 PM    version changed on Sept 20, 2006.  Routine now only waits for 100 sec rather than 1000


module FileExist 

    public file_exist,file_exists

contains
    subroutine file_exist(filename, lexist)
    implicit none
    character(*) filename
    logical(4) lexist
    integer(4) ierr1,ierr2,igo,itry

    inquire(file=filename, exist=lexist) !i trust inquire when it says file exists
    if(lexist) return

!     check to make sure does not exist

!    do itry=1,100
    do itry=1,10

        open( 901,file=filename,status='new',iostat=ierr1)    !return 10 if file exist, returns 29 or 30 if cannot find path
        close(901,status='delete',iostat=ierr2)

        if(ierr2.ne.0) then
            write(*,*) filename
            stop 'file_exist unable to delete test file, pgm stopped'
        endif

        if(ierr1.eq.0) return  !able to write a file, file must not exist
        if(ierr1.eq.9) return  !do not have premission to write a file, assume file does not exist

        if(ierr1.eq.10) then      !file must exist
        lexist=.true.
        return
        endif

        call sleepqq(10000) !wait 10 sec

    enddo
       
    write(*,1001) ierr1,filename
    1001 format('path cannot be found in file_exist ',i5,1x,a80)
    write(*,*) 'enter 0 to say file does not exist, enter 1 to stop program'
    read(*,*) igo
    if(igo.eq.1) stop 'pgm stopped by user'
    lexist=.false.

    return
    end subroutine file_exist 

    logical function file_exists(filename)  
        character(*) filename

        logical  :: lexist
        call file_exist(filename,lexist)
        file_exists=lexist
        return
    end function file_exists
end module FileExist

