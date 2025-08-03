
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine momentu
use modmain
use modulr
use modtest
implicit none
! local variables
integer idm,is,ias
integer nrc,nrci,ifq,ir
real(8) t1
! automatic arrays
real(8) rfft(nqpt)
complex(8) zfft(nfqrz)
! external functions
complex(8), external :: zfmtint
if (.not.spinpol) return
! calculate muffin-tin moments
mommttot(:)=0.d0
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    do ifq=1,nfqrz
      zfft(ifq)=zfmtint(nrc,nrci,wr2cmt(:,is),magqmt(:,ias,idm,ifq))
    end do
    mommt(idm,ias)=dble(zfft(1))
    mommttot(idm)=mommttot(idm)+mommt(idm,ias)
    call rzfftifc(3,ngridq,1,rfft,zfft)
    mommtru(idm,ias,1:nqpt)=rfft(1:nqpt)
  end do
end do
! find the interstitial and total moments
t1=omega/dble(ngtc)
do idm=1,ndmag
  do ifq=1,nfqrz
    zfft(ifq)=sum(magqir(1:ngtc,idm,ifq)*cfrc(1:ngtc))
  end do
  momir(idm)=t1*dble(zfft(1))
  momtot(idm)=mommttot(idm)+momir(idm)
  call rzfftifc(3,ngridq,1,rfft,zfft)
  do ir=1,nqpt
    momirru(idm,ir)=t1*rfft(ir)
    momtotru(idm,ir)=sum(mommtru(idm,1:natmtot,ir))+momirru(idm,ir)
  end do
end do
! total moment magnitude
if (ncmag) then
  momtotm=sqrt(momtot(1)**2+momtot(2)**2+momtot(3)**2)
else
  momtotm=abs(momtot(1))
end if
! write the muffin-tin moments to test file
call writetest(770,'ULR muffin-tin moments',nv=ndmag*natmtot*nqpt,tol=2.d-2, &
 rva=mommtru)
end subroutine

