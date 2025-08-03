
! Copyright (C) 2002-2011 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine wavefcr(tsh,lrstp,is,ia,ist,m,ld,wfcr)
use modmain
implicit none
! arguments
logical, intent(in) :: tsh
integer, intent(in) :: lrstp,is,ia,ist
! pass in m-1/2
integer, intent(in) :: m
integer, intent(in) :: ld
complex(4), intent(out) :: wfcr(ld,2)
! local variables
integer ias,nr,nri,ir
integer k,l,lm1,lm2
integer i,i1,i2
real(8) c1,c2,t0,t1,t2
l=lsp(ist,is)
k=ksp(ist,is)
if (l > lmaxo) then
  wfcr(:,:)=0.e0
  return
end if
ias=idxas(ia,is)
! calculate the Clebsch-Gordon coefficients
t1=sqrt(dble(l+m+1)/dble(2*l+1))
t2=sqrt(dble(l-m)/dble(2*l+1))
if (k == l+1) then
  c1=t1
  c2=t2
else
  c1=t2
  c2=-t1
end if
if (abs(m) <= l) then
  lm1=l*(l+1)+m+1
else
  lm1=0
end if
if (abs(m+1) <= l) then
  lm2=l*(l+1)+m+2
else
  lm2=0
end if
nr=nrmt(is)
nri=nrmti(is)
! zero the wavefunction
if (lrstp == 1) then
  wfcr(1:npmt(is),1:2)=0.e0
else
  wfcr(1:npcmt(is),1:2)=0.e0
end if
!----------------------------------!
!     inner part of muffin-tin     !
!----------------------------------!
if (l > lmaxi) goto 10
if (tsh) then
  i1=lm1
  i2=lm2
else
  i=0
end if
do ir=1,nri,lrstp
! major component of radial wavefunction
  t0=rwfcr(ir,1,ist,ias)*rlmt(ir,-1,is)
  if (tsh) then
    if (lm1 > 0) wfcr(i1,1)=t0*c1
    if (lm2 > 0) wfcr(i2,2)=t0*c2
    i1=i1+lmmaxi
    i2=i2+lmmaxi
  else
    t1=t0*c1
    t2=t0*c2
    if (lm1 > 0) wfcr(i+1:i+lmmaxi,1)=t1*zbshti(1:lmmaxi,lm1)
    if (lm2 > 0) wfcr(i+1:i+lmmaxi,2)=t2*zbshti(1:lmmaxi,lm2)
    i=i+lmmaxi
  end if
end do
!----------------------------------!
!     outer part of muffin-tin     !
!----------------------------------!
10 continue
if (lrstp == 1) then
  i=lmmaxi*nrmti(is)
else
  i=lmmaxi*nrcmti(is)
end if
if (tsh) then
  i1=i+lm1
  i2=i+lm2
end if
do ir=nri+lrstp,nr,lrstp
  t0=rwfcr(ir,1,ist,ias)*rlmt(ir,-1,is)
  if (tsh) then
    if (lm1 > 0) wfcr(i1,1)=t0*c1
    if (lm2 > 0) wfcr(i2,2)=t0*c2
    i1=i1+lmmaxo
    i2=i2+lmmaxo
  else
    t1=t0*c1
    t2=t0*c2
    if (lm1 > 0) wfcr(i+1:i+lmmaxo,1)=t1*zbshto(1:lmmaxo,lm1)
    if (lm2 > 0) wfcr(i+1:i+lmmaxo,2)=t2*zbshto(1:lmmaxo,lm2)
    i=i+lmmaxo
  end if
end do
end subroutine

