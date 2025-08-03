
! Copyright (C) 2022 J. K. Dewhurst, S. Sharma and E. Harris-Lee.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modtdhfc

! energy cut-off for TDHFC states around the Fermi energy
real(8) ecutthc
! total number of TDHFC states
integer nthc
! number of states for each k-point
integer, allocatable :: nthck(:)
! index to used states for each k-point
integer, allocatable :: idxthc(:,:)

end module

