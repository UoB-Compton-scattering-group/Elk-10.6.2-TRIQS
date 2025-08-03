
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeengyu(fnum)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
write(fnum,*)
write(fnum,'("Energies per unit cell:")')
write(fnum,'(" Fermi",T30,": ",G22.12)') efermi
write(fnum,'(" sum of eigenvalues",T30,": ",G22.12)') evalsum
end subroutine

