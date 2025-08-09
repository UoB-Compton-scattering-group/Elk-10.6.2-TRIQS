
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dmftdmat
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik,lp,nthd
! initialise universal variables
call init0
call init1
call init3
! read Fermi energy from file
call readfermi
! get the eigenvalues from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
end do
!get DFT+DMFT density matrix
call occupy
call getdmatdmft
! compute the DFT+DMFT density matrices and write the natural orbitals and occupation
! numbers to EVECSV.OUT and OCCSV.OUT, respectively
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP SCHEDULE(DYNAMIC) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi) /= lp_mpi) cycle
  call dmftdmatk(ik)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! broadcast occupation number array to every MPI process
if (np_mpi > 1) then
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(occsv(:,ik),nstsv,mpi_double_precision,lp,mpicom,ierror)
  end do
end if
if (mp_mpi) then
! write the occupation numbers to file
  do ik=1,nkpt
    call putoccsv(filext,ik,occsv(:,ik))
  end do
  write(*,*)
  write(*,'("Info(dmftdmat):")')
  write(*,'(" DFT+DMFT density matrices determined for each k-point")')
  write(*,*)
  write(*,'(" Natural orbitals and occupation numbers written to")')
  write(*,'(" EVECSV.OUT and OCCSV.OUT, respectively")')
end if
end subroutine

