
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine epsinv
use modmain
use modmpi
use modomp
implicit none
! local variables
integer iq,ik,ig,iw
integer n,nthd
! automatic arrays
integer(omp_lock_kind) lock(nwrf)
real(8) vgqc(3,ngrf),gqc(ngrf),gclgq(ngrf)
! allocatable arrays
real(8), allocatable :: jlgqr(:,:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:),epsi(:,:,:)
complex(4), allocatable :: vchi0(:,:,:)
! allocate local arrays
allocate(jlgqr(njcmax,nspecies,ngrf))
allocate(ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot))
allocate(epsi(ngrf,ngrf,nwrf),vchi0(ngrf,ngrf,nwrf))
! initialise the OpenMP locks
do iw=1,nwrf
  call omp_init_lock(lock(iw))
end do
if (mp_mpi) write(*,*)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! loop over q-points
do iq=1,nqpt
  if (mp_mpi) write(*,'("Info(epsinv): ",I6," of ",I6," q-points")') iq,nqpt
! generate the G+q-vectors and related functions
  call gengqf(ngrf,vqc(:,iq),vgqc,gqc,jlgqr,ylmgq,sfacgq)
! generate the regularised Coulomb Green's function in G+q-space
  call gengclgq(.true.,iq,ngrf,gqc,gclgq)
! use the symmetric form
  gclgq(1:ngrf)=sqrt(gclgq(1:ngrf))
! zero the response function
  vchi0(1:ngrf,1:ngrf,1:nwrf)=0.e0
  call holdthd(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP SCHEDULE(DYNAMIC) &
!$OMP NUM_THREADS(nthd)
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi) /= lp_mpi) cycle
! compute v¹⸍² χ₀ v¹⸍²
    call genvchi0(.false.,ik,lock,vql(:,iq),gclgq,jlgqr,ylmgq,sfacgq,ngrf,vchi0)
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
! add vchi0 from each process and redistribute
  if (np_mpi > 1) then
    n=ngrf*ngrf*nwrf
    call mpi_allreduce(mpi_in_place,vchi0,n,mpi_complex,mpi_sum,mpicom,ierror)
  end if
! negate and add δ(G,G')
  epsi(1:ngrf,1:ngrf,1:nwrf)=-vchi0(1:ngrf,1:ngrf,1:nwrf)
  do ig=1,ngrf
    epsi(ig,ig,1:nwrf)=epsi(ig,ig,1:nwrf)+1.d0
  end do
!-------------------------------------!
!     invert epsilon over G-space     !
!-------------------------------------!
  call holdthd(nwrf,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP SCHEDULE(DYNAMIC) &
!$OMP NUM_THREADS(nthd)
  do iw=1,nwrf
    call zminv(ngrf,epsi(:,:,iw))
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
! write inverse RPA epsilon to EPSINV.OUT
  if (mp_mpi) call putepsinv(iq,epsi)
! end loop over q-points
end do
! destroy the OpenMP locks
do iw=1,nwrf
  call omp_destroy_lock(lock(iw))
end do
deallocate(jlgqr,ylmgq,sfacgq,epsi,vchi0)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

