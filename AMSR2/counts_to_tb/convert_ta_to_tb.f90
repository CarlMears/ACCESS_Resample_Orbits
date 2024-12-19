!     this routine converts ta to tb and stores results in the ta array 
        subroutine convert_ta_to_tb
            use l2_module                                             
            implicit none
            
            integer(4) icel,iscan
            integer(4) ifreq,ich1,ich2
            real(4) tax(2),tbx(2)
    
            do iscan=1,numscan
                do icel=1,maxcel

                    !     native ta with no resampling applied
                    do ifreq=1,6
                        ich1=2*(ifreq-1) + 1
                        ich2=ich1+1
                        tax=ta(icel,iscan,ich1:ich2)
                        call tatb(2,ifreq, tax,tbx)  !2 means go from ta to tb
                        ta(icel,iscan,ich1:ich2)=tbx
                    enddo  !ifreq
    
                    !     ta maps to 7 ghz footprint
                    do ifreq=1,6
                        ich1=2*(ifreq-1) + 1
                        ich2=ich1+1
                        tax=tar(icel,iscan,ich1:ich2)
                        call tatb(2,ifreq, tax,tbx)  !2 means go from ta to tb
                        tar(icel,iscan,ich1:ich2)=tbx
                    enddo  !ifreq
    
                    !     ta maps to 11 ghz footprint
                    do ifreq=3,6
                        ich1=2*(ifreq-1) + 9
                        ich2=ich1+1
                        tax=tar(icel,iscan,ich1:ich2)
                        call tatb(2,ifreq, tax,tbx)  !2 means go from ta to tb
                        tar(icel,iscan,ich1:ich2)=tbx
                    enddo  !ifreq

                    !     ta maps to 19 ghz footprint, no remapping of 19 ghz
                    do ifreq=5,6
                        ich1=2*(ifreq-1) + 13
                        ich2=ich1+1
                        tax=tar(icel,iscan,ich1:ich2)
                        call tatb(2,ifreq, tax,tbx)  !2 means go from ta to tb
                        tar(icel,iscan,ich1:ich2)=tbx
                    enddo  !ifreq
    
                enddo !icel
            enddo !iscan

            return
        end subroutine convert_ta_to_tb
    
        subroutine convert_ta_to_tb_native_only
            use l2_module                                             
            implicit none
            
            integer(4) icel,iscan
            integer(4) ifreq,ich1,ich2
            real(4) tax(2),tbx(2)
            
            do iscan=1,numscan
                do icel=1,maxcel
                    !native ta with no resampling applied
                    do ifreq=1,6  !used to be 6
                        ich1=2*(ifreq-1) + 1
                        ich2=ich1+1
                        tax=ta(icel,iscan,ich1:ich2)
                        call tatb(2,ifreq, tax,tbx)  !2 means go from ta to tb
                        ta(icel,iscan,ich1:ich2)=tbx
                    enddo  !ifreq
                enddo !icel
            enddo !iscan

            do iscan=1,numscan*2
                do icel=1,maxcel_89
                    !native ta with no resampling applied
                    if (modulo(iscan,2)==0) then
                        ifreq=7
                    else
                        ifreq=8
                    endif
                    tax=ta_89(icel,iscan,:)
                    call tatb(2,ifreq, tax,tbx)  !2 means go from ta to tb
                    ta_89(icel,iscan,:)=tbx
                enddo !icel
            enddo !iscan

            return
        end subroutine convert_ta_to_tb_native_only
