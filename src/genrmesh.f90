
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genrmesh
! !INTERFACE:
subroutine genrmesh
! !USES:
use modmain
use modvars
! !DESCRIPTION:
!   Generates the coarse and fine radial meshes for each atomic species in the
!   crystal. Also determines which points are in the inner part of the
!   muffin-tin using the value of {\tt fracinr}.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,nr,nrc,ir,irc,l
real(8) t1
! estimate the number of radial mesh points to infinity
nrspmax=1
do is=1,nspecies
! logarithmic mesh
  t1=log(rmaxsp(is)/rminsp(is))/log(rmt(is)/rminsp(is))
  t1=dble(nrmt(is)-1)*t1
  nrsp(is)=nint(t1)+1
  nrspmax=max(nrspmax,nrsp(is))
end do
! compute and store (R_mt)ˡ
if (allocated(rmtl)) deallocate(rmtl)
allocate(rmtl(0:lmaxo+3,nspecies))
do is=1,nspecies
  do l=0,lmaxo+3
    rmtl(l,is)=rmt(is)**l
  end do
end do
! generate the radial meshes
if (allocated(rsp)) deallocate(rsp)
allocate(rsp(nrspmax,nspecies))
if (allocated(rlmt)) deallocate(rlmt)
allocate(rlmt(nrmtmax,-lmaxo-1:lmaxo+2,nspecies))
if (allocated(wr2mt)) deallocate(wr2mt)
allocate(wr2mt(nrmtmax,nspecies))
if (allocated(wprmt)) deallocate(wprmt)
allocate(wprmt(4,nrmtmax,nspecies))
if (allocated(wcrmt)) deallocate(wcrmt)
allocate(wcrmt(12,nrmtmax,nspecies))
do is=1,nspecies
! generate logarithmic radial mesh
  t1=log(rmt(is)/rminsp(is))/dble(nrmt(is)-1)
  do ir=1,nrsp(is)
    rsp(ir,is)=rminsp(is)*exp(dble(ir-1)*t1)
  end do
! calculate rˡ on the fine radial mesh
  nr=nrmt(is)
  do l=-lmaxo-1,lmaxo+2
    rlmt(1:nr,l,is)=rsp(1:nr,is)**l
  end do
! determine the weights for spline integration on the fine radial mesh
  call wsplint(nr,rsp(:,is),wr2mt(:,is))
! multiply by r²
  wr2mt(1:nr,is)=wr2mt(1:nr,is)*rlmt(1:nr,2,is)
! determine the weights for partial integration on fine radial mesh
  call wsplintp(nr,rsp(:,is),wprmt(:,:,is))
! determine the weights for the spline coefficients
  call wspline(nr,rsp(:,is),wcrmt(:,:,is))
end do
! determine the fraction of the muffin-tin radius which defines the inner part
if (fracinr < 0.d0) fracinr=sqrt(dble(lmmaxi)/dble(lmmaxo))
! set up the coarse radial meshes and find the inner part of the muffin-tin
! where rho is calculated with lmaxi
if (allocated(rcmt)) deallocate(rcmt)
allocate(rcmt(nrcmtmax,nspecies))
if (allocated(rlcmt)) deallocate(rlcmt)
allocate(rlcmt(nrcmtmax,-lmaxo-1:lmaxo+2,nspecies))
if (allocated(wr2cmt)) deallocate(wr2cmt)
allocate(wr2cmt(nrcmtmax,nspecies))
if (allocated(wprcmt)) deallocate(wprcmt)
allocate(wprcmt(4,nrcmtmax,nspecies))
if (allocated(wcrcmt)) deallocate(wcrcmt)
allocate(wcrcmt(12,nrcmtmax,nspecies))
do is=1,nspecies
  t1=fracinr*rmt(is)
  nrmti(is)=1
  nrcmti(is)=1
  do ir=1,nrmt(is),lradstp
    irc=(ir-1)/lradstp+1
    rcmt(irc,is)=rsp(ir,is)
    if (rsp(ir,is) < t1) then
      nrmti(is)=ir
      nrcmti(is)=irc
    end if
  end do
! store rˡ on the coarse radial mesh
  do l=-lmaxo-1,lmaxo+2
    do ir=1,nrmt(is),lradstp
      irc=(ir-1)/lradstp+1
      rlcmt(irc,l,is)=rlmt(ir,l,is)
    end do
  end do
! determine the weights for spline integration on the coarse radial mesh
  nrc=nrcmt(is)
  call wsplint(nrc,rcmt(:,is),wr2cmt(:,is))
! multiply by r²
  wr2cmt(1:nrc,is)=wr2cmt(1:nrc,is)*rlcmt(1:nrc,2,is)
! determine the weights for partial integration on coarse radial mesh
  call wsplintp(nrc,rcmt(:,is),wprcmt(:,:,is))
! determine the weights for the spline coefficients
  call wspline(nrc,rcmt(:,is),wcrcmt(:,:,is))
end do
! write to VARIABLES.OUT
if (wrtvars) then
  call writevars('nrsp',nv=nspecies,iva=nrsp)
  call writevars('nrmt',nv=nspecies,iva=nrmt)
  call writevars('nrmti',nv=nspecies,iva=nrmti)
  call writevars('lradstp',iv=lradstp)
  call writevars('nrcmt',nv=nspecies,iva=nrcmt)
  call writevars('nrcmti',nv=nspecies,iva=nrcmti)
  do is=1,nspecies
    call writevars('rsp',nv=nrmt(is),rva=rsp(:,is))
  end do
end if
end subroutine
!EOC

