
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writefsm(fnum)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer is,ia,ias
write(fnum,*)
if (any(abs(fsmtype) == [1,3,4,6])) then
  write(fnum,'("FSM global effective field",T30,": ",3G18.10)') bfsmc(1:ndmag)
end if
if (any(abs(fsmtype) == [2,3,5,6])) then
  write(fnum,'("FSM local muffin-tin effective fields :")')
  do is=1,nspecies
    write(fnum,'(" species : ",I4," (",A,")")') is,trim(spsymb(is))
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(fnum,'("  atom ",I4,T30,": ",3G18.10)') ia,bfsmcmt(1:ndmag,ias)
    end do
  end do
end if
end subroutine

