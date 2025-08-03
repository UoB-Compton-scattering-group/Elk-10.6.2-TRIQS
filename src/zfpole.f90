
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure complex(8) function zfpole(c,z)
use modgw
implicit none
! arguments
complex(8), intent(in) :: c(*),z
! local variables
integer i,j
complex(8) z1
zfpole=c(1)
i=2
do j=1,npole
  z1=c(i)+z
  if (abs(z1%re)+abs(z1%im) > 1.d-8) then
    zfpole=zfpole+c(i+1)/z1
  end if
  i=i+2
end do
end function

