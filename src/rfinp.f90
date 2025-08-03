
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rfinp
! !INTERFACE:
real(8) function rfinp(rfmt1,rfir1,rfmt2,rfir2)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   rfmt1 : first function in real spherical harmonics for all muffin-tins
!           (in,real(npmtmax,natmtot))
!   rfir1 : first real interstitial function in real-space (in,real(ngtot))
!   rfmt2 : second function in real spherical harmonics for all muffin-tins
!           (in,real(npmtmax,natmtot))
!   rfir2 : second real interstitial function in real-space (in,real(ngtot))
! !DESCRIPTION:
!   Calculates the inner product of two real functions over the entire unit
!   cell. The input muffin-tin functions should have angular momentum cut-off
!   {\tt lmaxo}. In the interstitial region, the integrand is multiplied with
!   the characteristic function, $\tilde{\Theta}({\bf r})$, to remove the
!   contribution from the muffin-tin. See routines {\tt rfmtinp} and
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created July 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rfmt1(npmtmax,natmtot),rfir1(ngtot)
real(8), intent(in) :: rfmt2(npmtmax,natmtot),rfir2(ngtot)
! local variables
integer is,ias,nthd
real(8) t1
! external functions
real(8), external :: rfmtinp
rfinp=0.d0
call holdthd(natmtot+1,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is,t1) REDUCTION(+:rfinp) &
!$OMP NUM_THREADS(nthd)
!$OMP DO SCHEDULE(DYNAMIC)
do ias=1,natmtot
  is=idxis(ias)
! muffin-tin contribution
  rfinp=rfinp+rfmtinp(nrmt(is),nrmti(is),wr2mt(:,is),rfmt1(:,ias),rfmt2(:,ias))
end do
!$OMP END DO NOWAIT
! interstitial contribution
!$OMP SINGLE
t1=sum(rfir1(1:ngtot)*rfir2(1:ngtot)*cfunir(1:ngtot))
rfinp=rfinp+t1*(omega/ngtot)
!$OMP END SINGLE
!$OMP END PARALLEL
call freethd(nthd)
end function
!EOC

