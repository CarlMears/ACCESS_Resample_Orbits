!   7/13/2021 convert to module

module open_file_routines

    public open_direct_w_error,open_direct,open_binary,open_binary_buffered
    public open_text,open_text_append,openbig_try,openbig

contains
!   3/7/06 version change on 3/28/07.  open_binary_buffered routine added
    subroutine open_direct_w_error(iunit,afile,astatus,aaction,ashare,irecl,iseconds,error)
        implicit none

        character*(*) afile,astatus,aaction,ashare
        integer(4) iunit,irecl,iseconds,itry,ierror,igo,kstart,error

        kstart=1; itry=0
        error = 0

        10 continue

        open(iunit,file=afile,status=astatus,recl=irecl,form='unformatted',access='direct',action=aaction,share=ashare,iostat=ierror)

        if(ierror.ne.0) then
        if(kstart.eq.1) then
            kstart=0
            write(*,1001) ierror,iseconds,afile
        endif

        itry=itry+1
        if(itry-1.eq.iseconds) then
            write(*,1001) ierror,iseconds,afile
            1001 format('open_direct having problems: ',i7,i10,a100)
            write(*,*) 'returning error'
            error = -1
            return
        endif

        call sleepqq(1000) !check every 1 seconds
        goto 10
        endif

        return
    end subroutine open_direct_w_error
    
    subroutine open_direct(iunit,afile,astatus,aaction,ashare,irecl,iseconds)
        implicit none

        character*(*) afile,astatus,aaction,ashare
        integer(4) iunit,irecl,iseconds,itry,ierror,igo,kstart

        kstart=1; itry=0

        10 continue

        open(iunit,file=afile,status=astatus,recl=irecl,form='unformatted',access='direct',action=aaction,share=ashare,iostat=ierror)

        if(ierror.ne.0) then

        if(kstart.eq.1) then
            kstart=0
            write(*,1001) ierror,iseconds,afile
        endif

        itry=itry+1
        if(itry-1.eq.iseconds) then
            write(*,1001) ierror,iseconds,afile
            1001   format('open_direct having problems: ',i7,i10,a100)
            write(*,*) 'enter 1 to keep trying, enter 2 to stop program'
            read(*,*) igo
            if(igo.eq.2) stop 'program stopped by user'
            itry=0
            goto 10
        endif

        call sleepqq(1000) !check every 1 seconds
        goto 10
        endif
        return
    end subroutine open_direct


    subroutine open_binary(iunit,afile,astatus,aaction,ashare,iseconds)
        implicit none

        character*(*) afile,astatus,aaction,ashare
        integer(4) iunit,iseconds,itry,ierror,igo,kstart

        kstart=1; itry=0

        10 continue

        open(iunit,file=afile,status=astatus,form='binary',access='sequential',action=aaction,share=ashare,iostat=ierror)

        if(ierror.ne.0) then

        if(kstart.eq.1) then
        kstart=0
        write(*,1001) ierror,iseconds,afile
        endif

        itry=itry+1
        if(itry-1.eq.iseconds) then
            write(*,1001) ierror,iseconds,afile
            1001 format('open_binary having problems: ',i7,i10,a100)
            write(*,*) 'enter 1 to keep trying, enter 2 to stop program'
            read(*,*) igo
            if(igo.eq.2) stop 'program stopped by user'
            itry=0
            goto 10
        endif

        call sleepqq(1000) !check every 1 seconds
        goto 10
        endif

        return
    end subroutine open_binary


    subroutine open_binary_buffered(iunit,afile,astatus,aaction,ashare,iseconds)
        implicit none

        character*(*) afile,astatus,aaction,ashare
        integer(4) iunit,iseconds,itry,ierror,igo,kstart

        kstart=1; itry=0

        10 continue

        open(iunit,file=afile,status=astatus,form='binary',access='sequential',action=aaction,share=ashare,iostat=ierror, &
            blocksize=32000,buffercount=64,buffered='yes')

        if(ierror.ne.0) then

        if(kstart.eq.1) then
        kstart=0
        write(*,1001) ierror,iseconds,afile
        endif

        itry=itry+1
        if(itry-1.eq.iseconds) then
            write(*,1001) ierror,iseconds,afile
            1001 format('open_binary having problems: ',i7,i10,a100)
            write(*,*) 'enter 1 to keep trying, enter 2 to stop program'
            read(*,*) igo
            if(igo.eq.2) stop 'program stopped by user'
            itry=0
            goto 10
        endif
        call sleepqq(1000) !check every 1 seconds
        goto 10
        endif

        return
    end subroutine open_binary_buffered

    subroutine open_text(iunit,afile,astatus,aaction,ashare,iseconds)
        implicit none

        character*(*) afile,astatus,aaction,ashare
        integer(4) iunit,iseconds,itry,ierror,igo,kstart

        kstart=1; itry=0

        10 continue

        open(iunit,file=afile,status=astatus,form='formatted',access='sequential',action=aaction,share=ashare,iostat=ierror)

        if(ierror.ne.0) then

        if(kstart.eq.1) then
        kstart=0
        write(*,1001) ierror,iseconds,afile
        endif

        itry=itry+1
        if(itry-1.eq.iseconds) then
            write(*,1001) ierror,iseconds,afile
            1001 format('open_text having problems: ',i7,i10,a100)
            write(*,*) 'enter 1 to keep trying, enter 2 to stop program'
            read(*,*) igo
            if(igo.eq.2) stop 'program stopped by user'
            itry=0
            goto 10
        endif

        call sleepqq(1000) !check every 1 seconds
        goto 10
        endif

        return
    end subroutine open_text


    subroutine open_text_append(iunit,afile,astatus,aaction,ashare,iseconds)
        implicit none

        character*(*) afile,astatus,aaction,ashare
        integer(4) iunit,iseconds,itry,ierror,igo,kstart

        kstart=1; itry=0

        10 continue

        open(iunit,file=afile,status=astatus,form='formatted',access='sequential',action=aaction,share=ashare,position='append', &
        iostat=ierror)

        if(ierror.ne.0) then
        if(kstart.eq.1) then
        kstart=0
        write(*,1001) ierror,iseconds,afile
        endif

        itry=itry+1
        if(itry-1.eq.iseconds) then
            write(*,1001) ierror,iseconds,afile
            1001 format('open_text_append having problems: ',i7,i10,a100)
            write(*,*) 'enter 1 to keep trying, enter 2 to stop program'
            read(*,*) igo
            if(igo.eq.2) stop 'program stopped by user'
            itry=0
            goto 10
        endif

        call sleepqq(1000) !check every 1 seconds
        goto 10
        endif

        return
    end subroutine open_text_append

    subroutine openbig_try(iunit,filename,s)
        use msflib                                        
        implicit none
    
        character*(*) filename,s
        integer(4) iunit,itry,ierror,igo
    
        itry=0
    
        10 continue
    
        open(iunit,file=filename,status=s,blocksize=32000,access='sequential',form='binary',share='denywr',iostat=ierror)
    
        if(ierror.ne.0) then
            itry=itry+1
            if(itry.eq.10) then
                itry=0
                write(*,*) filename
                write(*,*) 'openbig_try having problems',ierror
                write(*,*) 'enter 1 to keep trying, enter 2 to stop program'
                read(*,*) igo
                if(igo.eq.2) stop 'program stopped by user'
                goto 10
            endif
            call sleepqq(60000) !check every 1 minutes
            goto 10
        endif
    
        return
    end subroutine openbig_try

    subroutine openbig(iunit,filename,s)

        integer(4)   :: iunit
        character*(*) filename,s
        
        integer(4)   :: ierror
        integer(4)   :: igo
  
        5 continue
  
        open(iunit,file=filename,status=s,blocksize=32000,access='sequential',form='binary', &
            share='denynone',iostat=ierror,err=10)
        return
  
        10 continue
        write(*,1001) filename,iunit,s,ierror
        1001 format('error in openbig ',a50,1x,i3,1x,a10,i7)
        write(*,*) 'enter 0 to try again, enter 1 to stop'
        read(*,*) igo
        if(igo.eq.1) stop 'program stopped by user'
        goto 5
    end subroutine openbig
end module open_file_routines

