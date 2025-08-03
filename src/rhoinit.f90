
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhoinit
! !INTERFACE:
subroutine rhoinit
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Initialises the crystal charge density. Inside the muffin-tins it is set to
!   the spherical atomic density. In the interstitial region it is taken to be
!   constant such that the total charge is correct. Requires that the atomic
!   densities have already been calculated.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer lmax,is,ia,ias,nthd
integer nr,nri,nro,nrs,iro,ir
integer nrc,nrci,irco,irc
integer l,lm,i0,i1,ig,ifg
real(8) x,sm,t1,t2
complex(8) z1,z2
! automatic arrays
real(8) ffg(ngvc),wr(nrspmax),jl(0:lmaxi,nrcmtmax)
complex(8) zfmt(npcmtmax),zfft(ngtc)
lmax=min(lmaxi,1)
! compute the superposition of all the atomic density tails
zfft(1:ngtc)=0.d0
call holdthd(nspecies,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ffg,wr,nr,nrs,nro,ig) &
!$OMP PRIVATE(t1,t2,sm,ir,x,ia,ias,ifg) &
!$OMP REDUCTION(+:zfft) &
!$OMP SCHEDULE(DYNAMIC) &
!$OMP NUM_THREADS(nthd)
do is=1,nspecies
  nr=nrmt(is)
  nrs=nrsp(is)
  nro=nrs-nr+1
! determine the weights for the radial integral
  call wsplint(nro,rsp(nr,is),wr(nr))
  do ig=1,ngvc
    t1=gc(ig)
! spherical bessel function j_0(x) times the atomic density tail
    if (t1 > epslat) then
      t2=1.d0/t1
      sm=0.d0
      do ir=nr,nrs
        x=t1*rsp(ir,is)
        sm=sm+t2*sin(x)*rhosp(ir,is)*rsp(ir,is)*wr(ir)
      end do
    else
      sm=sum(rhosp(nr:nrs,is)*(rsp(nr:nrs,is)**2)*wr(nr:nrs))
    end if
! apply low-pass filter
    t1=sm*exp(-4.d0*(gc(ig)/gmaxvr)**2)
    ffg(ig)=(fourpi/omega)*t1
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ig=1,ngvc
      ifg=igfc(ig)
      zfft(ifg)=zfft(ifg)+ffg(ig)*conjg(sfacg(ig,ias))
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! compute the tails in each muffin-tin
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(jl,zfmt,is,nrc,nrci) &
!$OMP PRIVATE(irco,ig,ifg,irc,x) &
!$OMP PRIVATE(z1,z2,lm,l,i0,i1) &
!$OMP SCHEDULE(DYNAMIC) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  irco=nrci+1
  zfmt(1:npcmt(is))=0.d0
  do ig=1,ngvc
    ifg=igfc(ig)
    do irc=1,nrc
      x=gc(ig)*rcmt(irc,is)
      call sbessel(lmax,x,jl(:,irc))
    end do
    z1=zfft(ifg)*sfacg(ig,ias)
    do l=0,lmax
      do lm=l**2+1,(l+1)**2
        z2=z1*conjg(ylmg(lm,ig))
        i1=lmmaxi*(nrci-1)+lm
        zfmt(lm:i1:lmmaxi)=zfmt(lm:i1:lmmaxi)+jl(l,1:nrci)*z2
        i0=i1+lmmaxi
        i1=lmmaxo*(nrc-irco)+i0
        zfmt(i0:i1:lmmaxo)=zfmt(i0:i1:lmmaxo)+jl(l,irco:nrc)*z2
      end do
    end do
  end do
  call ztorfmt(nrc,nrci,zfmt,rhomt(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! convert the density from a coarse to a fine radial mesh
call rfmtctof(rhomt)
! add the atomic charge density and the excess charge in each muffin-tin
t1=chgexs/omega
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  iro=nri+1
  i1=lmmaxi*(nri-1)+1
  rhomt(1:i1:lmmaxi,ias)=rhomt(1:i1:lmmaxi,ias)+(t1+rhosp(1:nri,is))*y00i
  i0=i1+lmmaxi
  i1=lmmaxo*(nr-iro)+i0
  rhomt(i0:i1:lmmaxo,ias)=rhomt(i0:i1:lmmaxo,ias)+(t1+rhosp(iro:nr,is))*y00i
end do
! deallocate the global array rhosp as it is not used again
deallocate(rhosp)
! interstitial density determined from the atomic tails and excess charge
call zfftifc(3,ngdgc,1,zfft)
do ir=1,ngtc
  rhoir(ir)=dble(zfft(ir))+t1
! make sure that the density is always positive
  if (rhoir(ir) < 1.d-10) rhoir(ir)=1.d-10
end do
! convert the interstitial density from coarse to fine grid
call rfirctof(rhoir,rhoir)
end subroutine
!EOC

