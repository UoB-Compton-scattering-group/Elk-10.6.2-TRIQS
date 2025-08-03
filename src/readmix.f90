
! Copyright (C) 2025 Eddie Harris-Lee, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readmix(iscl0,nwork,work)
use modmain
implicit none
! arguments
integer, intent(out) :: iscl0
integer, intent(in) :: nwork
real(8), intent(out) :: work(nwork)
! local variables
integer ios,nwork_
open(80,file='MIXWORK'//trim(filext),form='UNFORMATTED',action='READ', &
 status='OLD',iostat=ios)
if (ios /= 0) then
  write(*,*)
  write(*,'("Error(readmix): error opening ",A)') 'MIXWORK'//trim(filext)
  write(*,*)
  stop
end if
read(80) iscl0
read(80) nwork_
if (nwork /= nwork_) then
  write(*,*)
  write(*,'("Error(readmix): differing nwork")')
  write(*,'(" current     : ",I12)') nwork
  write(*,'(" MIXWORK.OUT : ",I12)') nwork_
  write(*,*)
  stop
end if
read(80) work
close(80)
end subroutine

