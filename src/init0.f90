
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: init0
! !INTERFACE:
subroutine init0
! !USES:
use modmain
use modxcifc
use moddftu
use modtddft
use modphonon
use modulr
use modgw
use modtest
use modvars
use modmpi
use modomp
! !DESCRIPTION:
!   Performs basic consistency checks as well as allocating and initialising
!   global variables not dependent on the $k$-point set.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
logical hybrid_
integer is,ia,ias,npc
integer idu,ist,nr,l,n,i
integer xcspin_
real(8) hybridc_,t1
real(8) ts0,ts1

!-------------------------------!
!     zero timing variables     !
!-------------------------------!
timeinit=0.d0
timemat=0.d0
timefv=0.d0
timesv=0.d0
timerho=0.d0
timepot=0.d0
timefor=0.d0
call timesec(ts0)

!------------------------------------!
!     angular momentum variables     !
!------------------------------------!
if (lmaxo > lmaxapw) then
  write(*,*)
  write(*,'("Error(init0): lmaxo > lmaxapw : ",2I8)') lmaxo,lmaxapw
  write(*,*)
  stop
end if
lmaxi=min(lmaxi,lmaxo)
lmmaxapw=(lmaxapw+1)**2
lmmaxi=(lmaxi+1)**2
lmmaxo=(lmaxo+1)**2
! check DOS lmax is within range
lmaxdos=min(lmaxdos,lmaxo)
! write to VARIABLES.OUT
if (wrtvars) then
  call writevars('lmaxapw',iv=lmaxapw)
  call writevars('lmaxi',iv=lmaxi)
  call writevars('lmaxo',iv=lmaxo)
end if

!------------------------------------!
!     index to atoms and species     !
!------------------------------------!
natmmax=0
ias=0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=ias+1
    idxas(ia,is)=ias
    idxis(ias)=is
    idxia(ias)=ia
  end do
! maximum number of atoms over all species
  natmmax=max(natmmax,natoms(is))
end do
! total number of atoms
natmtot=ias
! number of phonon branches
nbph=3*natmtot
! write to VARIABLES.OUT
if (wrtvars) then
  call writevars('nspecies',iv=nspecies)
  call writevars('natoms',nv=nspecies,iva=natoms)
  call writevars('spsymb',nv=nspecies,sva=spsymb)
  call writevars('spname',nv=nspecies,sva=spname)
  call writevars('spzn',nv=nspecies,rva=spzn)
end if

!------------------------!
!     spin variables     !
!------------------------!
if (spinsprl) then
  spinpol=.true.
  spinorb=.false.
  if (any(task == [5,51,52,53,61,62,63,700,701])) then
    write(*,*)
    write(*,'("Error(init0): spin-spirals do not work with task ",I4)') task
    write(*,*)
    stop
  end if
  if (xctype(1) < 0) then
    write(*,*)
    write(*,'("Error(init0): spin-spirals do not work with the OEP method")')
    write(*,*)
    stop
  end if
  if (tefvit) then
    write(*,*)
    write(*,'("Error(init0): spin-spirals do not work with iterative &
     &diagonalisation")')
    write(*,*)
    stop
  end if
end if
! de-phasing required only for spin-spirals
if (.not.spinsprl) ssdph=.false.
! spin-orbit coupling, B-field-orbit coupling, fixed spin moment, spin spirals
! or spin-polarised cores requires a spin-polarised calculation
if (spinorb.or.bforb.or.(fsmtype /= 0).or.spinsprl.or.spincore) spinpol=.true.
! number of spinor components and maximum allowed occupancy
if (spinpol) then
  nspinor=2
  occmax=1.d0
else
  nspinor=1
  occmax=2.d0
end if
! number of spin-dependent first-variational functions per state and map from
! second- to first-variational spin index
if (spinsprl) then
  nspnfv=2
  jspnfv(1)=1
  jspnfv(2)=2
else
  nspnfv=1
  jspnfv(1)=1
  jspnfv(2)=1
end if
! no calculation of second-variational eigenvectors by default
tevecsv=.false.
! spin-polarised calculations require second-variational eigenvectors
if (spinpol) tevecsv=.true.
! Hartree-Fock/RDMFT/TDDFT/GW/ULR requires second-variational eigenvectors
if (any(task == [5,10,170,300,460,461,462,463,600,601,620,700,701,720]).or. &
 ksgwrho) then
  tevecsv=.true.
end if
! get exchange-correlation functional data
call getxcdata(xctype,xcdescr,xcspin,xcgrad,hybrid0,hybridc)
if ((spinpol).and.(xcspin == 0)) then
  write(*,*)
  write(*,'("Error(init0): requested spin-polarised run with &
   &spin-unpolarised")')
  write(*,'(" exchange-correlation functional")')
  write(*,*)
  stop
end if
! set flag for hybrid functional
if (task == 5) then
  hybrid=hybrid0
else
  hybrid=.false.
end if
! check for collinearity in the z-direction and set the dimension of the
! magnetisation and exchange-correlation vector fields
if (spinpol) then
  ndmag=1
  if ((abs(bfieldc0(1)) > epslat).or.(abs(bfieldc0(2)) > epslat)) ndmag=3
  do is=1,nspecies
    do ia=1,natoms(is)
      if ((abs(bfcmt0(1,ia,is)) > epslat).or. &
          (abs(bfcmt0(2,ia,is)) > epslat)) ndmag=3
    end do
  end do
! spin-orbit coupling is non-collinear in general
  if (spinorb) ndmag=3
! source-free fields and spin-spirals must be non-collinear
  if (nosource.or.spinsprl) then
    ndmag=3
    cmagz=.false.
  end if
! force collinear magnetism along the z-axis if required
  if (cmagz) ndmag=1
else
  ndmag=0
end if
! set the non-collinear flag
ncmag=(ndmag == 3)
! check for meta-GGA with non-collinearity
if (any(xcgrad == [3,4,5,6]).and.ncmag) then
  write(*,*)
  write(*,'("Error(init0): meta-GGA is not valid for non-collinear magnetism")')
  write(*,*)
  stop
end if
if (tbdip.and.(.not.ncmag)) then
  write(*,*)
  write(*,'("Error(init0): non-collinear magnetism required for inclusion of &
   &the dipole field")')
  write(*,*)
  stop
end if
! spin-polarised cores
if (.not.spinpol) spincore=.false.
if (fsmtype /= 0) then
! set fixed spin moment effective field to zero
  bfsmc(:)=0.d0
! set muffin-tin FSM fields to zero
  if (allocated(bfsmcmt)) deallocate(bfsmcmt)
  allocate(bfsmcmt(3,natmtot))
  bfsmcmt(:,:)=0.d0
  if (mp_mpi.and.(mixtype /= 1)) then
    write(*,'("Info(init0): mixtype changed to 1 for FSM calculation")')
  end if
  mixtype=1
end if
! number of independent spin components of the f_xc spin tensor
if (spinpol) then
  if (ncmag) then
    nscfxc=10
  else
    nscfxc=3
  end if
else
  nscfxc=1
end if
! set the magnetic fields to the initial values
bfieldc(:)=bfieldc0(:)
bfcmt(:,:,:)=bfcmt0(:,:,:)
! if reducebf < 1 then reduce the external magnetic fields immediately for
! non-self-consistent calculations or resumptions
if (reducebf < 1.d0-1.d-4) then
  if (all(task /= [0,1,2,3,5,28,200,201,350,351])) then
    bfieldc(:)=0.d0
    bfcmt(:,:,:)=0.d0
  end if
end if
if (tmwrite.or.(ftmtype /= 0).or.(task == 400)) then
  if (.not.spinorb) then
    write(*,*)
    write(*,'("Error(init0): tensor moments require spin-orbit coupling &
     &enabled")')
    write(*,'(" set spinorb=.true.")')
    write(*,*)
    stop
  end if
end if
! generate the fixed tensor moment density matrices if required
if (ftmtype /= 0) call gendmftm
! write to VARIABLES.OUT
if (wrtvars) then
  call writevars('nspinor',iv=nspinor)
  call writevars('ndmag',iv=ndmag)
  call writevars('sxcscf',rv=sxcscf)
  call writevars('bfieldc0',nv=3,rva=bfieldc0)
end if

!----------------------------------!
!     crystal structure set up     !
!----------------------------------!
! use reference lattice vectors if required
if (tavref) then
  avec0(:,:)=avec(:,:)
  avec(:,:)=avecref(:,:)
end if
! generate the reciprocal lattice vectors, unit cell and Brillouin zone volumes
call reciplat(avec,bvec,omega,omegabz)
! inverse of the lattice vector matrix
call r3minv(avec,ainv)
! inverse of the reciprocal vector matrix
call r3minv(bvec,binv)
! Cartesian coordinates of the spin-spiral vector
call r3mv(bvec,vqlss,vqcss)
do is=1,nspecies
  do ia=1,natoms(is)
! map atomic lattice coordinates to [0,1)
    call r3frac(epslat,atposl(:,ia,is))
! determine atomic Cartesian coordinates
    call r3mv(avec,atposl(:,ia,is),atposc(:,ia,is))
  end do
end do
! check for overlapping muffin-tins and adjust radii if required
call checkmt
! compute the total muffin-tin volume (M. Meinert)
omegamt=0.d0
do is=1,nspecies
  omegamt=omegamt+dble(natoms(is))*(fourpi/3.d0)*rmt(is)**3
end do
! input q-vector in Cartesian coordinates
call r3mv(bvec,vecql,vecqc)

!-------------------------------!
!     vector fields E and A     !
!-------------------------------!
! static electric field
tefield=.false.
if (any(abs(efieldc(:)) > epslat)) then
! no shift of the atomic positions
  tshift=.false.
! electric field vector in lattice coordinates
  call r3mv(ainv,efieldc,efieldl)
! potential at maximum distance
  vmaxefc=dmaxefc*norm2(efieldc(1:3))
! allocate array for average electric field in each muffin-tin
  if (allocated(efcmt)) deallocate(efcmt)
  allocate(efcmt(3,natmtot))
! set the E-field flag
  tefield=.true.
end if
! static vector potential
tafield=.false.
if (any(abs(afieldc(:)) > epslat)) then
  tafield=.true.
! A-field in lattice coordinates
  call r3mv(ainv,afieldc,afieldl)
! vector potential added in second-variational step
  tevecsv=.true.
end if
! static spin-dependent vector potential
tafsp=.false.
if (spinpol.and.(any(abs(afspc(:,:)) > epslat))) then
  tafsp=.true.
  tevecsv=.true.
end if
! time-dependent vector potential
tafieldt=.false.
if (any(task == [460,461,462,463,480,481,485])) then
! read time-dependent A-field from file
  call readafieldt
  tafieldt=.true.
! zero the induced A-field and its time derivative
  afindt(:,:)=0.d0
end if
! write to VARIABLES.OUT
if (wrtvars) then
  call writevars('efieldc',nv=3,rva=efieldc)
  call writevars('afieldc',nv=3,rva=afieldc)
  call writevars('afspc',nv=9,rva=afspc)
end if

!---------------------------------!
!     crystal symmetry set up     !
!---------------------------------!
call symmetry

!-----------------------!
!     radial meshes     !
!-----------------------!
nrmtmax=1
nrcmtmax=1
do is=1,nspecies
! make the muffin-tin mesh commensurate with lradstp
  nrmt(is)=nrmt(is)-mod(nrmt(is)-1,lradstp)
  nrmtmax=max(nrmtmax,nrmt(is))
! number of coarse radial mesh points
  nrcmt(is)=(nrmt(is)-1)/lradstp+1
  nrcmtmax=max(nrcmtmax,nrcmt(is))
end do
! set up atomic and muffin-tin radial meshes
call genrmesh
! number of points in packed muffin-tins
npmtmax=1
npcmtmax=1
do is=1,nspecies
  npmti(is)=lmmaxi*nrmti(is)
  npmt(is)=npmti(is)+lmmaxo*(nrmt(is)-nrmti(is))
  npmtmax=max(npmtmax,npmt(is))
  npcmti(is)=lmmaxi*nrcmti(is)
  npcmt(is)=npcmti(is)+lmmaxo*(nrcmt(is)-nrcmti(is))
  npcmtmax=max(npcmtmax,npcmt(is))
end do
! index to first muffin-tin point in packed array over all atoms
if (allocated(ipcmt)) deallocate(ipcmt)
allocate(ipcmt(natmtot))
ipcmt(1)=1
npc=npcmt(1)
npcmttot=npc
do ias=2,natmtot
  is=idxis(ias)
  ipcmt(ias)=ipcmt(ias)+npc
  npc=npcmt(is)
  npcmttot=npcmttot+npc
end do

!--------------------------------------!
!     charges and number of states     !
!--------------------------------------!
chgzn=0.d0
chgcrtot=0.d0
chgval=0.d0
nstspmax=0
nstcr=0
do is=1,nspecies
! nuclear charge
  chgzn=chgzn+spzn(is)*natoms(is)
! find the maximum number of atomic states
  nstspmax=max(nstspmax,nstsp(is))
! compute the electronic charge for each species, as well as the total core and
! valence charge
  spze(is)=0.d0
  chgcr(is)=0.d0
  do ist=1,nstsp(is)
    spze(is)=spze(is)+occsp(ist,is)
    if (spcore(ist,is)) then
      chgcr(is)=chgcr(is)+occsp(ist,is)
      nstcr=nstcr+2*ksp(ist,is)*natoms(is)
    else
      chgval=chgval+occsp(ist,is)*natoms(is)
    end if
  end do
  chgcrtot=chgcrtot+chgcr(is)*natoms(is)
end do
! add excess charge
chgval=chgval+chgexs
! total charge
chgtot=chgcrtot+chgval
if (chgtot < 1.d-8) then
  write(*,*)
  write(*,'("Error(init0): zero total charge")')
  write(*,*)
  stop
end if
! subtract small charge to ensure consistent Fermi energy for insulators
chgval=chgval-epsocc
! effective Wigner radius
rwigner=(3.d0/(fourpi*(chgtot/omega)))**(1.d0/3.d0)
! write to VARIABLES.OUT
if (wrtvars) then
  call writevars('spze',nv=nspecies,rva=spze)
  call writevars('chgcr',nv=nspecies,rva=chgcr)
  call writevars('chgexs',rv=chgexs)
  call writevars('chgval',rv=chgtot)
end if

!-------------------------!
!     G-vector arrays     !
!-------------------------!
! determine gkmax from rgkmax
if (nspecies == 0) isgkmax=-2
select case(isgkmax)
case(:-4)
! use largest muffin-tin radius
  gkmax=rgkmax/maxval(rmt(1:nspecies))
case(-3)
! use smallest muffin-tin radius
  gkmax=rgkmax/minval(rmt(1:nspecies))
case(-2)
! use the fixed value of 2.0
  gkmax=rgkmax/2.d0
case(-1,0)
! use average muffin-tin radius
  t1=sum(natoms(1:nspecies)*rmt(1:nspecies))/dble(natmtot)
  gkmax=rgkmax/t1
case(1:)
! use user-specified muffin-tin radius
  if (isgkmax <= nspecies) then
    gkmax=rgkmax/rmt(isgkmax)
  else
    write(*,*)
    write(*,'("Error(init0): isgkmax > nspecies : ",2I8)') isgkmax,nspecies
    write(*,*)
    stop
  end if
end select
! generate the G-vectors
call gengvec
! apply strain to A-, B- and G-vectors if required
call strainabg
! write number of G-vectors to test file
call writetest(900,'number of G-vectors',iv=ngvec)
! Poisson solver pseudocharge density constant
if (nspecies > 0) then
  t1=0.25d0*gmaxvr*maxval(rmt(1:nspecies))
else
  t1=0.25d0*gmaxvr*2.d0
end if
npsd=max(nint(t1),1)
lnpsd=lmaxo+npsd+1
! generate the Coulomb Green's function in G-space = 4π/G²
call gengclg
! compute the spherical Bessel functions j_l(|G|R_mt)
if (allocated(jlgrmt)) deallocate(jlgrmt)
allocate(jlgrmt(0:lnpsd,ngvec,nspecies))
call genjlgprmt(lnpsd,ngvec,gc,ngvec,jlgrmt)
! generate the spherical harmonics of the G-vectors
call genylmg
! allocate structure factor array for G-vectors
if (allocated(sfacg)) deallocate(sfacg)
allocate(sfacg(ngvec,natmtot))
! generate structure factors for G-vectors
call gensfacgp(ngvec,vgc,ngvec,sfacg)
! generate the smooth step function form factors
if (allocated(ffacg)) deallocate(ffacg)
allocate(ffacg(ngtot,nspecies))
call genffacgp(ngtot,gc,ngtot,ffacg)
! generate the smooth characteristic function
call gencfun
! G-vector variables for coarse grid with |G| < 2 gkmax
call gengvc
! generate the characteristic function on the coarse grid
call gencfrc
! write to VARIABLES.OUT
if (wrtvars) then
  call writevars('avec',nv=9,rva=avec)
  call writevars('bvec',nv=9,rva=bvec)
  call writevars('omega',rv=omega)
  do is=1,nspecies
    call writevars('atposl',n1=is,nv=3*natoms(is),rva=atposl(:,:,is))
  end do
  do is=1,nspecies
    call writevars('atposc',n1=is,nv=3*natoms(is),rva=atposc(:,:,is))
  end do
  call writevars('vqlss',nv=3,rva=vqlss)
  call writevars('vqcss',nv=3,rva=vqcss)
  call writevars('gmaxvr',rv=gmaxvr)
  call writevars('ngridg',nv=3,iva=ngridg)
  call writevars('intgv',nv=6,iva=intgv)
  call writevars('ngvec',iv=ngvec)
  call writevars('ivg',nv=3*ngtot,iva=ivg)
  call writevars('igfft',nv=ngtot,iva=igfft)
end if

!-------------------------!
!     atoms and cores     !
!-------------------------!
! determine the nuclear Coulomb potential
if (allocated(vcln)) deallocate(vcln)
allocate(vcln(nrspmax,nspecies))
do is=1,nspecies
  nr=nrsp(is)
  call potnucl(ptnucl,nr,rsp(:,is),spzn(is),vcln(:,is))
  vcln(1:nr,is)=vcln(1:nr,is)*y00i
end do
! solve the Kohn-Sham-Dirac equations for all atoms
call allatoms
! allocate core state occupancy and eigenvalue arrays and set to default
if (allocated(occcr)) deallocate(occcr)
allocate(occcr(nstspmax,natmtot))
if (allocated(evalcr)) deallocate(evalcr)
allocate(evalcr(nstspmax,natmtot))
do ias=1,natmtot
  is=idxis(ias)
  do ist=1,nstsp(is)
    occcr(ist,ias)=occsp(ist,is)
    evalcr(ist,ias)=evalsp(ist,is)
  end do
end do
! allocate core state radial wavefunction array
if (allocated(rwfcr)) deallocate(rwfcr)
allocate(rwfcr(nrspmax,2,nstspmax,natmtot))
! number of core spin channels
if (spincore) then
  nspncr=2
else
  nspncr=1
end if
! allocate core state charge density array
if (allocated(rhocr)) deallocate(rhocr)
allocate(rhocr(nrmtmax,natmtot,nspncr))

!-------------------------------------------------------------!
!     charge density, potentials and exchange-correlation     !
!-------------------------------------------------------------!
! combined target array for density and magnetisation
if (allocated(rhmg)) deallocate(rhmg)
n=npmtmax*natmtot+ngtot
if (spinpol) n=n+(npmtmax*natmtot+ngtot)*ndmag
allocate(rhmg(n))
! associate pointer arrays with target
rhomt(1:npmtmax,1:natmtot) => rhmg(1:)
i=npmtmax*natmtot+1
rhoir(1:ngtot) => rhmg(i:)
if (spinpol) then
  i=i+ngtot
  magmt(1:npmtmax,1:natmtot,1:ndmag) => rhmg(i:)
  i=i+npmtmax*natmtot*ndmag
  magir(1:ngtot,1:ndmag) => rhmg(i:)
end if
if (any(task == [371,372,373]).or.tafield.or.tdjr1d.or.tdjr2d.or.tdjr3d) then
  tjr=.true.
end if
! allocate current density arrays
if (allocated(jrmt)) deallocate(jrmt)
if (allocated(jrir)) deallocate(jrir)
if (tjr) then
  allocate(jrmt(npmtmax,natmtot,3),jrir(ngtot,3))
end if
! Coulomb potential
if (allocated(vclmt)) deallocate(vclmt)
allocate(vclmt(npmtmax,natmtot))
if (allocated(vclir)) deallocate(vclir)
allocate(vclir(ngtot))
! exchange energy density
if (allocated(exmt)) deallocate(exmt)
allocate(exmt(npmtmax,natmtot))
if (allocated(exir)) deallocate(exir)
allocate(exir(ngtot))
! correlation energy density
if (allocated(ecmt)) deallocate(ecmt)
allocate(ecmt(npmtmax,natmtot))
if (allocated(ecir)) deallocate(ecir)
allocate(ecir(ngtot))
! exchange-correlation potential
if (allocated(vxcmt)) deallocate(vxcmt)
allocate(vxcmt(npmtmax,natmtot))
if (allocated(vxcir)) deallocate(vxcir)
allocate(vxcir(ngtot))
! exchange-correlation and dipole magnetic fields
if (allocated(bxcmt)) deallocate(bxcmt)
if (allocated(bxcir)) deallocate(bxcir)
if (allocated(bdmt)) deallocate(bdmt)
if (allocated(bdir)) deallocate(bdir)
if (allocated(bdmta)) deallocate(bdmta)
if (spinpol) then
  allocate(bxcmt(npmtmax,natmtot,ndmag),bxcir(ngtot,ndmag))
  if (tbdip) then
    allocate(bdmt(npmtmax,natmtot,ndmag),bdir(ngtot,ndmag))
    allocate(bdmta(ndmag,natmtot))
    bdmta(1:ndmag,1:natmtot)=0.d0
  end if
end if
! combined target array for Kohn-Sham potential and magnetic field
if (allocated(vsbs)) deallocate(vsbs)
n=npmtmax*natmtot+ngtc
if (spinpol) n=n+(npcmtmax*natmtot+ngtc)*ndmag
allocate(vsbs(n))
! associate pointer arrays with target
vsmt(1:npmtmax,1:natmtot) => vsbs(1:)
i=npmtmax*natmtot+1
vsirc(1:ngtc) => vsbs(i:)
if (spinpol) then
  i=i+ngtc
  bsmt(1:npcmtmax,1:natmtot,1:ndmag) => vsbs(i:)
  i=i+npcmtmax*natmtot*ndmag
  bsirc(1:ngtc,1:ndmag) => vsbs(i:)
! allocate the Kohn-Sham magnetic field on the intersitial grid
  if (allocated(bsir)) deallocate(bsir)
  allocate(bsir(ngtot,ndmag))
end if
! interstitial Kohn-Sham potential
if (allocated(vsir)) deallocate(vsir)
allocate(vsir(ngtot))
! interstitial Kohn-Sham potential in G-space
if (allocated(vsig)) deallocate(vsig)
allocate(vsig(ngvc))
! kinetic energy density and meta-GGA exchange-correlation potential
if (allocated(taumt)) deallocate(taumt)
if (allocated(tauir)) deallocate(tauir)
if (allocated(taucr)) deallocate(taucr)
if (allocated(wxcmt)) deallocate(wxcmt)
if (allocated(wxcir)) deallocate(wxcir)
if (any(xcgrad == [3,4,5,6])) then
  allocate(taumt(npmtmax,natmtot,nspinor),tauir(ngtot,nspinor))
  allocate(taucr(npmtmax,natmtot,nspinor))
  allocate(wxcmt(npmtmax,natmtot),wxcir(ngtot))
! approximate kinetic energy density functional used to compute the functional
! derivative δτ(r')/δρ(r) for meta-GGA
  call getxcdata(ktype,kdescr,xcspin_,kgrad,hybrid_,hybridc_)
end if
! spin-orbit coupling radial function
if (allocated(socfr)) deallocate(socfr)
if (spinorb) then
  allocate(socfr(nrcmtmax,natmtot))
end if
! allocate muffin-tin charge and moment arrays
if (allocated(chgcrlk)) deallocate(chgcrlk)
allocate(chgcrlk(natmtot))
if (allocated(chgmt)) deallocate(chgmt)
allocate(chgmt(natmtot))
if (allocated(mommt)) deallocate(mommt)
allocate(mommt(3,natmtot))
! check if scaled spin exchange-correlation should be used
tssxc=(abs(sxcscf-1.d0) > 1.d-6)
! spin-spiral phase factors
if (ssdph) then
  if (allocated(zqss)) deallocate(zqss)
  allocate(zqss(natmtot))
  do ias=1,natmtot
    is=idxis(ias)
    ia=idxia(ias)
    t1=-0.5d0*dot_product(vqcss(1:3),atposc(1:3,ia,is))
    zqss(ias)=cmplx(cos(t1),sin(t1),8)
  end do
end if
! mixing vector: either density/magnetisation or potential/magnetic field
if (mixrho) then
  vmixer => rhmg
else
  vmixer => vsbs
end if
! zero the mixing vector
vmixer(:)=0.d0

!-------------------------!
!     force variables     !
!-------------------------!
if (tforce) then
  if (allocated(forcehf)) deallocate(forcehf)
  allocate(forcehf(3,natmtot))
  if (allocated(forcetot)) deallocate(forcetot)
  allocate(forcetot(3,natmtot))
end if

!-------------------------------------------------!
!     DFT+U and fixed tensor moment variables     !
!-------------------------------------------------!
if ((dftu /= 0).or.(ftmtype /= 0)) then
! density matrix elements in each muffin-tin
  if (allocated(dmatmt)) deallocate(dmatmt)
  allocate(dmatmt(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
! potential matrix elements in each muffin-tin
  if (allocated(vmatmt)) deallocate(vmatmt)
  allocate(vmatmt(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
! zero the potential matrix
  vmatmt(:,:,:,:,:)=0.d0
! matrix elements in spherical coordinates for TDDFT+U
  if (any(task == [460,461,462,463,478])) then
    if (allocated(vmatmti)) deallocate(vmatmti)
    allocate(vmatmti(lmmaxi,lmmaxi,nspinor,nspinor,natmtot))
    if (allocated(vmatmto)) deallocate(vmatmto)
    allocate(vmatmto(lmmaxo,lmmaxo,nspinor,nspinor,natmtot))
  end if
! require the potential matrix elements be calculated
  tvmatmt=.true.
! flags for non-zero muffin-tin potential matrices
  if (allocated(tvmmt)) deallocate(tvmmt)
  allocate(tvmmt(0:lmaxdm,natmtot))
  tvmmt(:,:)=.false.
! require second-variational eigenvectors
  tevecsv=.true.
end if
if (dftu /= 0) then
  if (any(task == [5,300,600,601,610,620,640])) then
    write(*,*)
    write(*,'("Error(init0): DFT+U does not work with task ",I8)') task
    write(*,*)
    stop
  end if
! DFT+U energy for each atom
  if (allocated(engyadu)) deallocate(engyadu)
  allocate(engyadu(natmmax,ndftu))
! flag the muffin-tin potential matrices which are non-zero
  do idu=1,ndftu
    is=isldu(1,idu)
    if (is > nspecies) then
      write(*,*)
      write(*,'("Error(init0): invalid species number : ",I8)') is
      write(*,*)
      stop
    end if
    l=isldu(2,idu)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      tvmmt(l,ias)=.true.
    end do
  end do
! zero the initial values of screening length
  lamdu0(:)=0.d0
! write to VARIABLES.OUT
  if (wrtvars) then
    call writevars('udufix',nv=ndftu,rva=udufix)
  end if
end if
if (ftmtype /= 0) then
! allocate and zero the fixed tensor moment potential array
  if (allocated(vmftm)) deallocate(vmftm)
  allocate(vmftm(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
  vmftm(:,:,:,:,:)=0.d0
! flag the muffin-tin potential matrices which are non-zero
  do i=1,ntmfix
    is=itmfix(1,i)
    ia=itmfix(2,i)
    ias=idxas(ia,is)
    l=itmfix(3,i)
    tvmmt(l,ias)=.true.
  end do
end if

!-----------------------!
!     miscellaneous     !
!-----------------------!
! determine nuclear radii and volumes
call nuclei
! determine the nuclear-nuclear energy
call energynn
! get smearing function description
call getsdata(stype,sdescr)
! get mixing type description
call getmixdata(mixtype,mixdescr)
! generate the spherical harmonic transform (SHT) matrices
call genshtmat
! find the maximum size of the spherical Bessel function array over all species
call findnjcmax
! allocate 1D plotting arrays
if (allocated(dvp1d)) deallocate(dvp1d)
allocate(dvp1d(nvp1d))
if (allocated(vplp1d)) deallocate(vplp1d)
allocate(vplp1d(3,npp1d))
if (allocated(dpp1d)) deallocate(dpp1d)
allocate(dpp1d(npp1d))
! initial self-consistent loop number
iscl=1
tlast=.false.
! set the Fermi energy to zero
efermi=0.d0
! set the temperature from the smearing width
tempk=swidth/kboltz

call timesec(ts1)
timeinit=timeinit+ts1-ts0

end subroutine
!EOC

