
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ephcouple
use modmain
use modphonon
use modmpi
use modomp
implicit none
! local variables
integer iq,ik,jk,jkq
integer ist,jst,isym,ip
integer is,ia,ias,js,jas
integer nr,nri,np,i,n,nthd
real(8) vl(3),de,x
real(8) t1,t2,t3,t4
! allocatable arrays
real(8), allocatable :: gq(:,:),y(:)
complex(8), allocatable :: ev(:,:),a(:,:)
complex(8), allocatable :: dvmt(:,:,:),dvir(:,:)
complex(8), allocatable :: zfmt(:),gzfmt(:,:,:)
complex(8), allocatable :: ephmat(:,:,:)
! external functions
real(8), external :: sdelta
! increase the angular momentum cut-off on the inner part of the muffin-tin
lmaxi0=lmaxi
lmaxi=max(lmaxi,4)
! initialise universal variables
call init0
call init1
call init2
call initph
! allocate global arrays
if (allocated(dvsbs)) deallocate(dvsbs)
n=npmtmax*natmtot
allocate(dvsbs(n))
dvsmt(1:npmtmax,1:natmtot) => dvsbs(1:)
if (allocated(dvsir)) deallocate(dvsir)
allocate(dvsir(ngtot))
! allocate local arrays
allocate(gq(nbph,nqpt),y(nbph))
allocate(ev(nbph,nbph),a(nbph,nbph))
allocate(dvmt(npcmtmax,natmtot,nbph),dvir(ngtc,nbph))
allocate(zfmt(npmtmax),gzfmt(npmtmax,3,natmtot))
! read in the density and potentials from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! set the speed of light >> 1 (non-relativistic approximation)
solsc=sol*100.d0
! new file extension for eigenvector files with c >> 1
filext='_EPH.OUT'
! generate the first- and second-variational eigenvectors and eigenvalues
call linengy
call genapwlofr
call gensocfr
call genevfsv
! precise determination of the Fermi energy
swidth0=swidth
swidth=1.d-5
call occupy
swidth=swidth0
! restore the speed of light
solsc=sol
! compute the gradients of the Kohn-Sham potential for the rigid-ion term
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  call rtozfmt(nr,nri,vsmt(:,ias),zfmt)
  call gradzfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),zfmt,npmtmax,gzfmt(:,:,ias))
end do
! energy window for calculating the electron-phonon matrix elements
if (task == 240) then
  de=4.d0*swidth
else
  de=1.d6
end if
! loop over phonon q-points
do iq=1,nqpt
  if (mp_mpi) write(*,'("Info(ephcouple): ",I6," of ",I6," q-points")') iq,nqpt
! diagonalise the dynamical matrix
  call dynev(dynq(:,:,iq),wphq(:,iq),ev)
! generate the matrix for converting between Cartesian and phonon coordinates
  call genmcph(wphq(:,iq),ev,a)
  i=0
  do is=1,nspecies
    nr=nrmt(is)
    nri=nrmti(is)
    np=npmt(is)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ip=1,3
        i=i+1
! read in the Cartesian change in Kohn-Sham potential
        call readdvs(iq,is,ia,ip,dvsmt,dvsir)
! add the rigid-ion term
        dvsmt(1:np,ias)=dvsmt(1:np,ias)-gzfmt(1:np,ip,ias)
        do jas=1,natmtot
          js=idxis(jas)
! convert to coarse radial mesh
          call zfmtftoc(nrcmt(js),nrcmti(js),dvsmt(:,jas),dvmt(:,jas,i))
! apply the radial integration weights
          call zfmtwr(nrcmt(js),nrcmti(js),wr2cmt(:,js),dvmt(:,jas,i))
        end do
! multiply the interstitial potential with the characteristic function and
! convert to coarse grid
        call zfirftoc(dvsir,dvir(:,i))
      end do
    end do
  end do
  y(1:nbph)=0.d0
  call holdthd(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ephmat,jk,vl,isym,jkq) &
!$OMP PRIVATE(t1,t2,t3,t4,ist,jst,x,i) &
!$OMP REDUCTION(+:y) &
!$OMP NUM_THREADS(nthd)
  allocate(ephmat(nstsv,nstsv,nbph))
!$OMP DO SCHEDULE(DYNAMIC)
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi) /= lp_mpi) cycle
! equivalent reduced k-point
    jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! compute the electron-phonon coupling matrix elements
    call genephmat(iq,ik,de,a,dvmt,dvir,ephmat)
! write the matrix elements to file if required
    if (task == 241) call putephmat(iq,ik,ephmat)
! k+q-vector in lattice coordinates
    vl(1:3)=vkl(1:3,ik)+vql(1:3,iq)
! index to k+q-vector
    call findkpt(vl,isym,jkq)
    t1=pi*wkptnr*occmax
! loop over second-variational states
    do ist=1,nstsv
      x=(evalsv(ist,jkq)-efermi)/swidth
      t2=t1*sdelta(stype,x)/swidth
      do jst=1,nstsv
! loop over phonon branches
        do i=1,nbph
          x=(wphq(i,iq)+evalsv(jst,jk)-evalsv(ist,jkq))/swidth
          t3=t2*sdelta(stype,x)/swidth
          t4=dble(ephmat(ist,jst,i))**2+aimag(ephmat(ist,jst,i))**2
          y(i)=y(i)+wphq(i,iq)*t3*t4
        end do
      end do
    end do
! end loop over k-points
  end do
!$OMP END DO
  deallocate(ephmat)
!$OMP END PARALLEL
  call freethd(nthd)
! store in phonon linewidths array
  gq(1:nbph,iq)=y(1:nbph)
! end loop over phonon q-points
end do
! add gq from each MPI process
if (np_mpi > 1) then
  n=nbph*nqpt
  call mpi_allreduce(mpi_in_place,gq,n,mpi_double_precision,mpi_sum,mpicom, &
   ierror)
end if
! restore the default file extension
filext='.OUT'
if (mp_mpi) then
! write the phonon linewidths to file
  call writegamma(gq)
! write electron-phonon coupling constants to file
  call writelambda(gq)
  if (task == 241) then
    write(*,*)
    write(*,'("Info(ephcouple):")')
    write(*,'(" wrote electron-phonon matrix elements to EPHMAT.OUT")')
  end if
end if
! deallocate global arrays
deallocate(dvsbs)
! deallocate local arrays
deallocate(gq,y,ev,a,dvmt,dvir,zfmt,gzfmt)
! restore original input parameters
lmaxi=lmaxi0
end subroutine

