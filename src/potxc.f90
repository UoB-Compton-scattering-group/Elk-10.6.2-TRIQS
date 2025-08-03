
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potxc
! !INTERFACE:
subroutine potxc(tsh,xctype_,rhomt_,rhoir_,magmt_,magir_,taumt_,tauir_,exmt_, &
 exir_,ecmt_,ecir_,vxcmt_,vxcir_,bxcmt_,bxcir_,wxcmt_,wxcir_)
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Computes the exchange-correlation potential and energy density. In the
!   muffin-tin, the density is transformed from spherical harmonic coefficients
!   $\rho_{lm}$ to spherical coordinates $(\theta,\phi)$ with a backward
!   spherical harmonic transformation (SHT). Once calculated, the
!   exchange-correlation potential and energy density are transformed with a
!   forward SHT.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh
integer, intent(in) :: xctype_(3)
real(8), intent(in) :: rhomt_(npmtmax,natmtot),rhoir_(ngtot)
real(8), intent(in) :: magmt_(npmtmax,natmtot,ndmag),magir_(ngtot,ndmag)
real(8), intent(in) :: taumt_(npmtmax,natmtot,nspinor),tauir_(ngtot,nspinor)
real(8), intent(out) :: exmt_(npmtmax,natmtot),exir_(ngtot)
real(8), intent(out) :: ecmt_(npmtmax,natmtot),ecir_(ngtot)
real(8), intent(out) :: vxcmt_(npmtmax,natmtot),vxcir_(ngtot)
real(8), intent(out) :: bxcmt_(npmtmax,natmtot,ndmag),bxcir_(ngtot,ndmag)
real(8), intent(out) :: wxcmt_(npmtmax,natmtot),wxcir_(ngtot)
! local variables
integer ias,nthd
call holdthd(natmtot+1,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
! muffin-tin exchange-correlation potential, field and energy density
!$OMP DO SCHEDULE(DYNAMIC)
do ias=1,natmtot
  call potxcmt(tsh,ias,xctype_,rhomt_,magmt_,taumt_,exmt_,ecmt_,vxcmt_,bxcmt_, &
   wxcmt_)
end do
!$OMP END DO NOWAIT
! interstitial exchange-correlation potential, field and energy density
!$OMP SINGLE
call potxcir(xctype_,rhoir_,magir_,tauir_,exir_,ecir_,vxcir_,bxcir_,wxcir_)
!$OMP END SINGLE
! symmetrise the muffin-tin exchange-correlation potentials and magnetic fields
if (tsh) then
!$OMP SECTIONS
!$OMP SECTION
  call symrfmt(nrmt,nrmti,npmt,npmtmax,vxcmt_)
!$OMP SECTION
  if (spinpol) call symrvfmt(.true.,ncmag,nrmt,nrmti,npmt,npmtmax,bxcmt_)
!$OMP END SECTIONS
end if
!$OMP END PARALLEL
call freethd(nthd)
end subroutine
!EOC

