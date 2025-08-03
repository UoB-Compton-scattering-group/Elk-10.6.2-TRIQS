
! Copyright (C) 2025 Wenhan Chen, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine dkevulr(evecu,dk)
use modmain
use modulr
implicit none
! arguments
complex(8), intent(in) :: evecu(nstulr,nstulr)
real(8), intent(out) :: dk(nstulr)
! local variables
integer ikpa,ist,jst,i
real(8) v(3),sm
do jst=1,nstulr
  v(:)=0.d0
  do ikpa=1,nkpa
    sm=0.d0
    do ist=1,nstsv
      i=(ikpa-1)*nstsv+ist
      sm=sm+evecu(i,jst)%re**2+evecu(i,jst)%im**2
    end do
    v(1:3)=v(1:3)+sm*vqc(1:3,ikpa)
  end do
  dk(jst)=sqrt(v(1)**2+v(2)**2+v(3)**2)
end do
end subroutine

