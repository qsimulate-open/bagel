c
c Author: Toru Shiozaki
c Date  : September 2008
c
c Assuming 16kB L1 cache or more
c The arrays will occupy 12.8kB
c
#include <config.h>
      subroutine mytranspose_inplace(h,m,n)
      implicit none
      integer m,n
      double precision h(*)
#ifdef HAVE_MKL_H
#include <mkl.fi>
      call mkl_dimatcopy("c","t",m,n,1.0d0,h,m,n)
#else
      call abort
#endif
      return
      end subroutine
c
      subroutine mytranspose(h,m,n,vec)
      implicit none
      integer i,j,m,n,mresidual,nresidual,mlim,nlim
      double precision h(m,n),vec(n,m)
#ifdef HAVE_MKL_H
#include <mkl.fi>
      call mkl_domatcopy("c","t",m,n,1.0d0,h,m,vec,n)
#else
      if((m.eq.1).or.(n.eq.1)) then
c      quick return if possible
       call dcopy(m*n,h,1,vec,1)
       return
      endif
      mresidual=mod(m,10)
      nresidual=mod(n,10)
      mlim=m-mresidual-1
      nlim=n-nresidual-1
      do i=1,nlim,10
       do j=1,mlim,10
        vec(i  ,j  )=h(j  ,i  )
        vec(i+1,j  )=h(j  ,i+1)
        vec(i+2,j  )=h(j  ,i+2)
        vec(i+3,j  )=h(j  ,i+3)
        vec(i+4,j  )=h(j  ,i+4)
        vec(i+5,j  )=h(j  ,i+5)
        vec(i+6,j  )=h(j  ,i+6)
        vec(i+7,j  )=h(j  ,i+7)
        vec(i+8,j  )=h(j  ,i+8)
        vec(i+9,j  )=h(j  ,i+9)
        vec(i  ,j+1)=h(j+1,i  )
        vec(i+1,j+1)=h(j+1,i+1)
        vec(i+2,j+1)=h(j+1,i+2)
        vec(i+3,j+1)=h(j+1,i+3)
        vec(i+4,j+1)=h(j+1,i+4)
        vec(i+5,j+1)=h(j+1,i+5)
        vec(i+6,j+1)=h(j+1,i+6)
        vec(i+7,j+1)=h(j+1,i+7)
        vec(i+8,j+1)=h(j+1,i+8)
        vec(i+9,j+1)=h(j+1,i+9)
        vec(i  ,j+2)=h(j+2,i  )
        vec(i+1,j+2)=h(j+2,i+1)
        vec(i+2,j+2)=h(j+2,i+2)
        vec(i+3,j+2)=h(j+2,i+3)
        vec(i+4,j+2)=h(j+2,i+4)
        vec(i+5,j+2)=h(j+2,i+5)
        vec(i+6,j+2)=h(j+2,i+6)
        vec(i+7,j+2)=h(j+2,i+7)
        vec(i+8,j+2)=h(j+2,i+8)
        vec(i+9,j+2)=h(j+2,i+9)
        vec(i  ,j+3)=h(j+3,i  )
        vec(i+1,j+3)=h(j+3,i+1)
        vec(i+2,j+3)=h(j+3,i+2)
        vec(i+3,j+3)=h(j+3,i+3)
        vec(i+4,j+3)=h(j+3,i+4)
        vec(i+5,j+3)=h(j+3,i+5)
        vec(i+6,j+3)=h(j+3,i+6)
        vec(i+7,j+3)=h(j+3,i+7)
        vec(i+8,j+3)=h(j+3,i+8)
        vec(i+9,j+3)=h(j+3,i+9)
        vec(i  ,j+4)=h(j+4,i  )
        vec(i+1,j+4)=h(j+4,i+1)
        vec(i+2,j+4)=h(j+4,i+2)
        vec(i+3,j+4)=h(j+4,i+3)
        vec(i+4,j+4)=h(j+4,i+4)
        vec(i+5,j+4)=h(j+4,i+5)
        vec(i+6,j+4)=h(j+4,i+6)
        vec(i+7,j+4)=h(j+4,i+7)
        vec(i+8,j+4)=h(j+4,i+8)
        vec(i+9,j+4)=h(j+4,i+9)
        vec(i  ,j+5)=h(j+5,i  )
        vec(i+1,j+5)=h(j+5,i+1)
        vec(i+2,j+5)=h(j+5,i+2)
        vec(i+3,j+5)=h(j+5,i+3)
        vec(i+4,j+5)=h(j+5,i+4)
        vec(i+5,j+5)=h(j+5,i+5)
        vec(i+6,j+5)=h(j+5,i+6)
        vec(i+7,j+5)=h(j+5,i+7)
        vec(i+8,j+5)=h(j+5,i+8)
        vec(i+9,j+5)=h(j+5,i+9)
        vec(i  ,j+6)=h(j+6,i  )
        vec(i+1,j+6)=h(j+6,i+1)
        vec(i+2,j+6)=h(j+6,i+2)
        vec(i+3,j+6)=h(j+6,i+3)
        vec(i+4,j+6)=h(j+6,i+4)
        vec(i+5,j+6)=h(j+6,i+5)
        vec(i+6,j+6)=h(j+6,i+6)
        vec(i+7,j+6)=h(j+6,i+7)
        vec(i+8,j+6)=h(j+6,i+8)
        vec(i+9,j+6)=h(j+6,i+9)
        vec(i  ,j+7)=h(j+7,i  )
        vec(i+1,j+7)=h(j+7,i+1)
        vec(i+2,j+7)=h(j+7,i+2)
        vec(i+3,j+7)=h(j+7,i+3)
        vec(i+4,j+7)=h(j+7,i+4)
        vec(i+5,j+7)=h(j+7,i+5)
        vec(i+6,j+7)=h(j+7,i+6)
        vec(i+7,j+7)=h(j+7,i+7)
        vec(i+8,j+7)=h(j+7,i+8)
        vec(i+9,j+7)=h(j+7,i+9)
        vec(i  ,j+8)=h(j+8,i  )
        vec(i+1,j+8)=h(j+8,i+1)
        vec(i+2,j+8)=h(j+8,i+2)
        vec(i+3,j+8)=h(j+8,i+3)
        vec(i+4,j+8)=h(j+8,i+4)
        vec(i+5,j+8)=h(j+8,i+5)
        vec(i+6,j+8)=h(j+8,i+6)
        vec(i+7,j+8)=h(j+8,i+7)
        vec(i+8,j+8)=h(j+8,i+8)
        vec(i+9,j+8)=h(j+8,i+9)
        vec(i  ,j+9)=h(j+9,i  )
        vec(i+1,j+9)=h(j+9,i+1)
        vec(i+2,j+9)=h(j+9,i+2)
        vec(i+3,j+9)=h(j+9,i+3)
        vec(i+4,j+9)=h(j+9,i+4)
        vec(i+5,j+9)=h(j+9,i+5)
        vec(i+6,j+9)=h(j+9,i+6)
        vec(i+7,j+9)=h(j+9,i+7)
        vec(i+8,j+9)=h(j+9,i+8)
        vec(i+9,j+9)=h(j+9,i+9)
       enddo
       do j=mlim+2,m
        vec(i  ,j  )=h(j  ,i  )
        vec(i+1,j  )=h(j  ,i+1)
        vec(i+2,j  )=h(j  ,i+2)
        vec(i+3,j  )=h(j  ,i+3)
        vec(i+4,j  )=h(j  ,i+4)
        vec(i+5,j  )=h(j  ,i+5)
        vec(i+6,j  )=h(j  ,i+6)
        vec(i+7,j  )=h(j  ,i+7)
        vec(i+8,j  )=h(j  ,i+8)
        vec(i+9,j  )=h(j  ,i+9)
       enddo
      enddo
      do i=nlim+2,n
       do j=1,mlim,10
        vec(i  ,j  )=h(j  ,i  )
        vec(i  ,j+1)=h(j+1,i  )
        vec(i  ,j+2)=h(j+2,i  )
        vec(i  ,j+3)=h(j+3,i  )
        vec(i  ,j+4)=h(j+4,i  )
        vec(i  ,j+5)=h(j+5,i  )
        vec(i  ,j+6)=h(j+6,i  )
        vec(i  ,j+7)=h(j+7,i  )
        vec(i  ,j+8)=h(j+8,i  )
        vec(i  ,j+9)=h(j+9,i  )
       enddo
       do j=mlim+2,m
        vec(i  ,j  )=h(j  ,i  )
       enddo
      enddo
#endif
      return
      end
c
