
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine acgwse(ist,jst,sem,wr,ser)
use modmain
use modgw
implicit none
! arguments
integer ist,jst
complex(8), intent(in) :: sem(nstsv,nstsv,0:nwfm)
real(8), intent(in) :: wr(nwplot)
complex(8), intent(out) :: ser(nstsv,nstsv,nwplot)
! automatic arrays
complex(8) zm(0:nwfm),zwr(nwplot),zr(nwplot)
zm(0:nwfm)=sem(ist,jst,0:nwfm)
zwr(1:nwplot)=wr(1:nwplot)
select case(actype)
case(1)
! fit a multipole model
  call acpole(zm,zwr,zr)
case(10)
! stabilised Pade approximant
  call pades(nspade,swidth,nwfm+1,wfm,zm,nwplot,zwr,zr)
case default
  write(*,*)
  write(*,'("Error(acgwse): actype not defined : ",I8)') actype
  write(*,*)
  stop
end select
ser(ist,jst,1:nwplot)=zr(1:nwplot)
end subroutine

