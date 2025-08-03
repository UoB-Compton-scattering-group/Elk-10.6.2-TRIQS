
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine lopzflmn(lmax,n,ld,zflm,zlflm1,zlflm2,zlflm3)
implicit none
! arguments
integer, intent(in) :: lmax,n,ld
complex(8), intent(in) :: zflm(ld,n)
complex(8), intent(out) :: zlflm1(ld,n),zlflm2(ld,n),zlflm3(ld,n)
! local variables
integer l,m,lm,i
real(8) t1
complex(8), parameter :: zi=(0.d0,1.d0),zmi=(0.d0,-1.d0)
complex(8) z1
lm=0
do l=0,lmax
  do m=-l,l
    lm=lm+1
    if (m == -l) then
      zlflm1(lm,1:n)=0.d0
      zlflm2(lm,1:n)=0.d0
    end if
    if (m < l) then
      t1=0.5d0*sqrt(dble((l-m)*(l+m+1)))
      do i=1,n
        z1=t1*zflm(lm,i)
        zlflm1(lm+1,i)=z1
        zlflm2(lm+1,i)=zmi*z1
      end do
    end if
    if (m > -l) then
      t1=0.5d0*sqrt(dble((l+m)*(l-m+1)))
      do i=1,n
        z1=t1*zflm(lm,i)
        zlflm1(lm-1,i)=zlflm1(lm-1,i)+z1
        zlflm2(lm-1,i)=zlflm2(lm-1,i)+zi*z1
      end do
    end if
    if (m /= 0) then
      zlflm3(lm,1:n)=dble(m)*zflm(lm,1:n)
    else
      zlflm3(lm,1:n)=0.d0
    end if
  end do
end do
end subroutine

