    subroutine fd_tbsss_clm(itime,xlat,xlon, tbsss) 
        use l2_module
        implicit none

        integer(4) itime,i1,i2,j1,j2,k1,k2,ich
        real(4) xlat,xlon,tbsss(10),brief,a1,a2,b1,b2,c1,c2                          
                                                                       
        !      do time,lat,lon interpolation                                    
                                                                       
        !      2629800 is the average num sec in a month, 86400*365.25/12   

        brief=(itime-1314900)/2629800.d0                                  
        if(brief.lt.0) then                                               
            i1=12                                                             
            i2=1                                                              
            a1=-brief                                                         
            a2=1.-a1                                                          
        else                                                              
            i1=1+brief                                                        
            i2=i1+1                                                           
            if(i2.eq.13) i2=1                                                 
            a1=i1-brief                                                       
            a2=1.-a1                                                          
        endif  
      
        brief=xlat+89.5
        j1=1+brief                                                        
        j2=j1+1                                                           
        b1=j1-brief                                                       
        b2=1-b1
        if(j1.eq.  0) j1=  1
        if(j2.eq.181) j2=180

        brief=xlon-0.5 
        k1=1+brief                                                       
        k2=k1+1                                                           
        c1=k1-brief                                                       
        c2=1-c1                                                          
        if(k1.eq.  0) k1=360
        if(k2.eq.361) k2=  1 

 
        do ich=1,10                                                                  
            tbsss(ich)= &            
                a1*b1*(c1*real(isal_corr(ich,k1,j1,i1))+c2*real(isal_corr(ich,k2,j1,i1)))+  &               
                a1*b2*(c1*real(isal_corr(ich,k1,j2,i1))+c2*real(isal_corr(ich,k2,j2,i1)))+  &               
                a2*b1*(c1*real(isal_corr(ich,k1,j1,i2))+c2*real(isal_corr(ich,k2,j1,i2)))+  &               
                a2*b2*(c1*real(isal_corr(ich,k1,j2,i2))+c2*real(isal_corr(ich,k2,j2,i2))) 
        enddo
        tbsss=0.01*tbsss-0.5
                                                                 
        return                                                            
    end subroutine fd_tbsss_clm                                                              

