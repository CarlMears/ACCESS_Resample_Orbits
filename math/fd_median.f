!     k is desired position, to find median k=(n+1)/2
!     n is number of values in arr to be considered
!     arr in integer(4) arrry that contains numbers

      function fd_median(k,n,arr_in)   ! Numerical Recipees routine called select
	implicit none
	integer(4) k,n
	integer(4) arr_in(*),arr(n),fd_median,a,temp
      integer(4) i,ir,j,l,mid

	arr=arr_in(1:n)

      l=1
      ir=n
    1 if(ir-l.le.1)then
        if(ir-l.eq.1)then
          if(arr(ir).lt.arr(l))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
          endif
        endif
        fd_median=arr(k)
        return
      else
        mid=(l+ir)/2
        temp=arr(mid)
        arr(mid)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)

          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
    3   continue
          i=i+1
        if(arr(i).lt.a)goto 3
    4   continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
    5   arr(l+1)=arr(j)
        arr(j)=a
        if(j.ge.k)ir=j-1
        if(j.le.k)l=i

      endif
      goto 1
      end function fd_median



