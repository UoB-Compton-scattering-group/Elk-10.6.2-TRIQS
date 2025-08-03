
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine pades(ns,r,ni,zi,ui,no,zo,uo)
implicit none
! arguments
integer, intent(in) :: ns
real(8), intent(in) :: r
integer, intent(in) :: ni
complex(8), intent(in) :: zi(ni),ui(ni)
integer, intent(in) :: no
complex(8), intent(in) :: zo(no)
complex(8), intent(out) :: uo(no)
! local variables
integer i
real(8), parameter :: pi=3.1415926535897932385d0
real(8) t1,t2
complex(8) z1
! automatic arrays
complex(8) u1(ni),u2(no)
if (ns < 1) then
  write(*,*)
  write(*,'("Error(pades): ns < 1 : ",I8)') ns
  write(*,*)
  stop
end if
if (ns == 1) then
  call pade(ni,zi,ui,no,zo,uo)
  return
end if
uo(1:no)=0.d0
do i=1,ns
  t1=dble(i-1)/dble(ns)
  t2=6.d0*pi*t1
  z1=r*t1*cmplx(cos(t2),sin(t2),8)
  u1(1:ni)=ui(1:ni)+z1
  call pade(ni,zi,u1,no,zo,u2)
  uo(1:no)=uo(1:no)+u2(1:no)-z1
end do
t1=1.d0/dble(ns)
uo(1:no)=t1*uo(1:no)
end subroutine

