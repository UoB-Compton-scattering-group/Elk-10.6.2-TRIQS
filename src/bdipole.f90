
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bdipole
! !INTERFACE:
subroutine bdipole
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the magnetic dipole field arising from the spin and orbital
!   current. The total current density is
!   $$ {\bf j}({\bf r}) = {\rm Im}\sum_{i{\bf k}}^{\rm occ}
!    \varphi_{i{\bf k}}^{\dag}({\bf r})\nabla\varphi_{i{\bf k}}({\bf r})
!    -\frac{1}{c}{\bf A}\,\rho({\bf r})
!    +\frac{g_s}{4}\nabla\times{\bf m}({\bf r}), $$
!   where $g_s$ is the electron spin $g$-factor. The vector potential arising
!   from ${\bf j}({\bf r})$ is calculated by
!   $$ {\bf A}({\bf r})=\frac{1}{c}\int d^3r'\,\frac{{\bf j}({\bf r}')}
!    {|{\bf r}-{\bf r}'|}, $$
!   using the Poisson equation solver {\tt zpotcoul}. Finally, the magnetic
!   field is determined from ${\bf B}({\bf r})=\nabla\times{\bf A}({\bf r})$.
!   This field is included as a Zeeman term in the second-variational
!   Hamiltonian:
!   $$ \hat{H}\rightarrow\hat{H}+\frac{g_s}{4c}{\bf B}\cdot\boldsymbol\sigma. $$
!
! !REVISION HISTORY:
!   Created April 2018 (T. Mueller)
!EOP
!BOC
implicit none
! local variables
integer idm,is,ias
integer nrc,nrci,np,npc
real(8) cb,t1,t2
! automatic arrays
real(8) rfmt(npcmtmax)
! allocatable arrays
real(8), allocatable :: rvfmt(:,:,:),rvfir(:,:)
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
! external functions
real(8), external :: rfmtint
! prefactor for the spin dipole magnetic field
cb=gfacte/(4.d0*solsc)
! allocate local arrays
allocate(rvfmt(npmtmax,natmtot,3),rvfir(ngtot,3))
allocate(zrhomt(npmtmax,natmtot),zrhoir(ngtot))
allocate(zvclmt(npmtmax,natmtot),zvclir(ngtot))
! compute the curl of the magnetisation density, i.e. the magnetisation current
call curlrvf(magmt,magir,rvfmt,rvfir)
! multiply by prefactor
do idm=1,3
  do ias=1,natmtot
    is=idxis(ias)
    np=npmt(is)
    rvfmt(1:np,ias,idm)=cb*rvfmt(1:np,ias,idm)
  end do
end do
rvfir(1:ngtot,1:3)=cb*rvfir(1:ngtot,1:3)
! add the gauge-invariant current density if required
if (tjr.and.(iscl >= 1)) then
  call genjr
  t1=1.d0/solsc
  do idm=1,3
    do ias=1,natmtot
      is=idxis(ias)
      np=npmt(is)
      rvfmt(1:np,ias,idm)=rvfmt(1:np,ias,idm)+t1*jrmt(1:np,ias,idm)
    end do
  end do
  rvfir(1:ngtot,1:3)=rvfir(1:ngtot,1:3)+t1*jrir(1:ngtot,1:3)
end if
do idm=1,3
! transform to complex spherical harmonics
  do ias=1,natmtot
    is=idxis(ias)
    call rtozfmt(nrmt(is),nrmti(is),rvfmt(:,ias,idm),zrhomt(:,ias))
  end do
! solve Poisson's equation in the muffin-tin to find the A-field
  call genzvclmt(nrmt,nrmti,nrmtmax,rlmt,wprmt,npmtmax,zrhomt,zvclmt)
  zrhoir(1:ngtot)=rvfir(1:ngtot,idm)
! solve in the entire unit cell
  call zpotcoul(0,nrmt,nrmti,npmt,nrmtmax,rlmt,ngridg,igfft,ngvec,gc,gclg, &
   ngvec,jlgrmt,ylmg,sfacg,zrhoir,npmtmax,zvclmt,zvclir)
! convert muffin-tin A-field to real spherical harmonics
  do ias=1,natmtot
    is=idxis(ias)
    call ztorfmt(nrmt(is),nrmti(is),zvclmt(:,ias),rvfmt(:,ias,idm))
  end do
! store the real part of the interstitial A-field
  rvfir(1:ngtot,idm)=dble(zvclir(1:ngtot))
end do
! compute the curl of A to obtain the dipole B-field
call curlrvf(rvfmt,rvfir,bdmt,bdir)
! scale dipole B-field if required (by scaling the prefactor)
cb=cb*bdipscf
! add to the Kohn-Sham field
do idm=1,3
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
! convert to coarse radial mesh
    call rfmtftoc(nrc,nrci,bdmt(:,ias,idm),rfmt)
! convert to spherical coordinates
    call rbshtip(nrc,nrci,rfmt)
    bsmt(1:npc,ias,idm)=bsmt(1:npc,ias,idm)+cb*rfmt(1:npc)
! store the average dipole field in each muffin-tin
    t1=rfmtint(nrmt(is),nrmti(is),wr2mt(:,is),bdmt(:,ias,idm))
    t2=(4.d0/3.d0)*pi*rmt(is)**3
    bdmta(idm,ias)=t1/t2
  end do
end do
do idm=1,3
  bsir(1:ngtot,idm)=bsir(1:ngtot,idm)+cb*bdir(1:ngtot,idm)
end do
deallocate(rvfmt,rvfir)
deallocate(zrhomt,zrhoir,zvclmt,zvclir)
end subroutine
!EOC

