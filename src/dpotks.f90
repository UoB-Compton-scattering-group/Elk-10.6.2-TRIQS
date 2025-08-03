
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dpotks
use modmain
use modphonon
use modomp
implicit none
! local variables
integer is,ias,np,nthd
! convert density derivative to spherical coordinates
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP SCHEDULE(DYNAMIC) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call zbshtip(nrmt(is),nrmti(is),drhomt(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! compute the exchange-correlation potential derivative
call dpotxc
! convert density derivative to spherical harmonics
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP SCHEDULE(DYNAMIC) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call zfshtip(nrmt(is),nrmti(is),drhomt(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! generate the Coulomb potential derivative
call dpotcoul
! add to the Kohn-Sham potential derivative
do ias=1,natmtot
  is=idxis(ias)
  np=npmt(is)
  dvsmt(1:np,ias)=dvsmt(1:np,ias)+dvclmt(1:np,ias)
end do
dvsir(1:ngtot)=dvsir(1:ngtot)+dvclir(1:ngtot)
! remove the gradient part of the potential derivative for displaced muffin-tin
np=npmt(isph)
dvsmt(1:np,iasph)=dvsmt(1:np,iasph)+gvsmt(1:np)
end subroutine

