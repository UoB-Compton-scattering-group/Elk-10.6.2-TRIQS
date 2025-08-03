
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init2
use modmain
use modrdm
use modphonon
use modgw
use modvars
use modmpi
implicit none
! local variables
logical lsym(48)
integer isym,iv(3)
real(8) boxl(3,0:3)
real(8) ts0,ts1

call timesec(ts0)

!---------------------!
!     q-point set     !
!---------------------!
! check if the system is an isolated molecule
if (molecule) ngridq(:)=1
! store the point group symmetries for reducing the q-point set
if (reduceq == 0) then
  nsymqpt=1
  symqpt(:,:,1)=symlat(:,:,1)
else
  lsym(:)=.false.
  do isym=1,nsymcrys
    lsym(lsplsymc(isym))=.true.
  end do
  nsymqpt=0
  do isym=1,nsymlat
    if (lsym(isym)) then
      nsymqpt=nsymqpt+1
      symqpt(:,:,nsymqpt)=symlat(:,:,isym)
    end if
  end do
end if
if (any(task == [105,180,185,320,330,331])) then
! equal k- and q-point grids for nesting function, BSE and linear-reposnse TDDFT
  ngridq(:)=ngridk(:)
else if (any(task == [5,300,600,601,620]).or.(xctype(1) < 0).or.ksgwrho) then
! allow the q-point grid to be smaller than the k-point grid for OEP,
! Hartree-Fock, RDMFT and GW
  if (any(ngridq(:) < 1)) ngridq(:)=ngridk(:)
else
  ngridq(:)=abs(ngridq(:))
end if
! check that the q-point and k-point grids are commensurate for some tasks
if (any(task == [5,205,240,241,300,600,601,620]).or.(xctype(1) < 0).or. &
 ksgwrho) then
  iv(:)=mod(ngridk(:),ngridq(:))
  if (any(iv(:) /= 0)) then
    write(*,*)
    write(*,'("Error(init2): k-point grid incommensurate with q-point grid")')
    write(*,'(" ngridk : ",3I6)') ngridk
    write(*,'(" ngridq : ",3I6)') ngridq
    write(*,*)
    stop
  end if
end if
! allocate the q-point arrays
if (allocated(ivqiq)) deallocate(ivqiq)
allocate(ivqiq(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
if (allocated(ivqiqnr)) deallocate(ivqiqnr)
allocate(ivqiqnr(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
nqptnr=ngridq(1)*ngridq(2)*ngridq(3)
if (allocated(ivq)) deallocate(ivq)
allocate(ivq(3,nqptnr))
if (allocated(vql)) deallocate(vql)
allocate(vql(3,nqptnr))
if (allocated(vqc)) deallocate(vqc)
allocate(vqc(3,nqptnr))
if (allocated(wqpt)) deallocate(wqpt)
allocate(wqpt(nqptnr))
! set up the q-point box (offset should always be zero)
boxl(:,:)=0.d0
boxl(1,1)=1.d0; boxl(2,2)=1.d0; boxl(3,3)=1.d0
! generate the q-point set
! (note that the vectors vql and vqc are in the first Brillouin zone)
call genppts(.true.,nsymqpt,symqpt,ngridq,nqptnr,epslat,bvec,boxl,nqpt,ivqiq, &
 ivqiqnr,ivq,vql,vqc,wqpt,wqptnr)
! write the q-points to QPOINTS.OUT
if (mp_mpi) call writeqpts
! write to VARIABLES.OUT
if (wrtvars) then
  call writevars('nsymqpt',iv=nsymqpt)
  call writevars('symqpt',nv=9*nsymqpt,iva=symqpt)
  call writevars('ngridq',nv=3,iva=ngridq)
  call writevars('nqpt',iv=nqpt)
  call writevars('ivqiq',nv=nqptnr,iva=ivqiq)
  call writevars('ivq',nv=3*nqptnr,iva=ivq)
  call writevars('vql',nv=3*nqptnr,rva=vql)
  call writevars('wqpt',nv=nqpt,rva=wqpt)
end if

!--------------------------------------------------------!
!     OEP, Hartree-Fock, RDMFT, BSE and GW variables     !
!--------------------------------------------------------!
if (any(task == [5,180,185,300,320,330,331,600,601,620]).or. &
 (xctype(1) < 0).or.ksgwrho) then
! determine the regularised Coulomb Green's function for small q
  call gengclq
! output the Coulomb Green's function to GCLQ.OUT
  if (mp_mpi) call writegclq
! initialise OEP variables
  if (xctype(1) < 0) call initoep
end if
if (task == 300) then
  if (allocated(vclmat)) deallocate(vclmat)
  allocate(vclmat(nstsv,nstsv,nkpt))
  if (allocated(dkdc)) deallocate(dkdc)
  allocate(dkdc(nstsv,nstsv,nkpt))
end if

call timesec(ts1)
timeinit=timeinit+ts1-ts0

end subroutine

