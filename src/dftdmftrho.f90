
! Copyright (C) 2018 A. Davydov, P. Elliott, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine dftdmftrho
use modmain
use modgw
use modmpi
use modomp
implicit none
! local variables
logical exist
integer ik,nthd
integer nmix,nwork,lp
real(8) dv
! allocatable arrays
real(8), allocatable :: work(:)
! initialise universal variables
call init0
call init1
call init2
call init3
! read density and potentials from file
call readstate
call genvsig
! read the Fermi energy
call readfermi
! re-calculate the ground-state
call gencore
call linengy
call genapwlofr
call gensocfr
call genevfsv
call occupy

if (mp_mpi) then
! check if interface has been generated from a potential previous loop 
  inquire(file="DMFT_INFO.OUT", exist=exist)
! append DMFT_INFO.OUT
  if (exist) then
      open(60,file='DMFT_INFO.OUT',form='FORMATTED',position='append')
! create new DMFT_INFO.OUT file
  else
    open(60,file='DMFT_INFO.OUT',form='FORMATTED')
  end if
! open FERMIDOS.OUT
  open(62,file='FERMIDOS.OUT',form='FORMATTED')
! open MOMENT.OUT if required
  if (spinpol) open(63,file='MOMENT.OUT',form='FORMATTED')
! open GAP.OUT
  open(64,file='GAP.OUT',form='FORMATTED')
! open RMSDVS.OUT
  open(65,file='RMSDVS.OUT',form='FORMATTED')
! open MOMENTM.OUT
  if (spinpol) open(68,file='MOMENTM.OUT',form='FORMATTED')
  call writebox(60,"Calculating DFT+DMFT density")
end if
! size of mixing vector
nmix=size(vmixer)
! determine the size of the mixer work array
nwork=-1
call mixerifc(mixtype,nmix,vmixer,dv,nwork,vmixer)
allocate(work(nwork))
! initialise the mixer
iscl=0
call mixerifc(mixtype,nmix,vmixer,dv,nwork,work)
iscl=1
if (mp_mpi) then
  write(60,*)
  write(60,'("Previous Kohn-Sham Fermi energy : ",G18.10)') efermi
end if

!read in dmft density matrix
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
  call gwdmatk(ik)
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
! determine the DFT+DMFT density and magnetisation
call rhomag
! compute the Kohn-Sham potential and magnetic field before potential mixing
if (.not.mixrho) call potks(.true.)
! mix the old density/magnetisation or potential/field with the new
call mixerifc(mixtype,nmix,vmixer,dv,nwork,work)
! compute the Kohn-Sham potential and magnetic field after density mixing
if (mixrho) call potks(.true.)
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! compute the energy components
call energy
! calculate the Kohn-Sham energy eigenvectors, eigenvalues
! and occupations from the updated DFT+DMFT density
call genapwlofr
call gensocfr
call genevfsv
call occupy

if (mp_mpi) then
! write the Kohn-Sham occupation numbers to file
  do ik=1,nkpt
    call putoccsv(filext,ik,occsv(:,ik))
  end do
  call writeeval
  call writefermi
! write STATE.OUT file
  call writestate
! adnj edit - output DMFT energy components. Note that the total
!             energy will not be used as a convergence criteria
  !if(task.ne.630) 
  call writeengy(60)
! end edit
  write(60,*)
  write(60,'("Density of states at Fermi energy : ",G18.10)') fermidos
  write(60,'(" (states/Hartree/unit cell)")')
  write(60,*)
  write(60,'("Estimated indirect band gap : ",G18.10)') bandgap(1)
  write(60,'(" from k-point ",I6," to k-point ",I6)') ikgap(1),ikgap(2)
  write(60,'("Estimated direct band gap   : ",G18.10)') bandgap(2)
  write(60,'(" at k-point ",I6)') ikgap(3)
! output charges and moments
  call writechg(60)
  if (spinpol) call writemom(60)
  flush(60)
! write DOS at Fermi energy to FERMIDOS.OUT
  write(62,'(G18.10)') fermidos
  flush(62)
  if (spinpol) then
! write total moment to MOMENT.OUT
    write(63,'(3G18.10)') momtot(1:ndmag)
    flush(63)
! write total moment magnitude to MOMENTM.OUT
    write(68,'(G18.10)') momtotm
    flush(68)
  end if
! write estimated Kohn-Sham indirect band gap
  write(64,'(G22.12)') bandgap(1)
  flush(64)
end if
! check for convergence
if (mp_mpi) then
  write(60,*)
  write(60,'("RMS change in Kohn-Sham potential (target) : ",G18.10," (",&
    &G18.10,")")') dv,epspot
  flush(60)
  write(65,'(G18.10)') dv
  flush(65)
end if
if (dv.lt.epspot) then
  if (mp_mpi) then
    write(60,*)
    write(60,'("Convergence targets achieved")')
  end if
end if
! reset the OpenMP thread variables
call omp_reset

! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
if (mp_mpi) then
  call writebox(60,"Calculated DFT+DMFT density")
! close the DMFT_INFO.OUT file
  close(60)
! close the FERMIDOS.OUT file
  close(62)
! close the MOMENT.OUT and MOMENTM.OUT files
  if (spinpol) then
    close(63); close(68)
  end if
! close the GAP.OUT file
  close(64)
! close the RMSDVS.OUT file
  close(65)
end if
end subroutine

