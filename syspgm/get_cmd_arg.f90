subroutine get_cmd_arg(numarg,acmd)
    use dflib 
	implicit none

	character(256),intent(out) :: acmd(100)
	integer(4),intent(out)     ::  numarg
	
	integer(2) iarg,istatus

    numarg=nargs( )
    numarg=numarg-1 !narg also include the exe aname, so numarg is one minus
	if(numarg.le.100)then
          do iarg=1,numarg
            call getarg(iarg, acmd(iarg),istatus)
	      if(istatus.le.0) then
	          write(*,1001) numarg,istatus
                stop
            endif
            acmd(iarg)=trim(acmd(iarg))  !probably do not need this
	    enddo
	    return
	else
        write(*,1001) numarg,istatus
 1001   format('error in get_cmd_arg, enter anything, program will stop ',2i10)
        stop
    endif
end


	

