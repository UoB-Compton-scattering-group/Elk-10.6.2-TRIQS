
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnzh(n,ld,a,w)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: n,ld
complex(8), intent(inout) :: a(ld,n)
real(8), intent(out) :: w(n)
! local variables
integer nb,lwork,lrwork
integer nts,info,nthd
! automatic arrays
integer iwork(3+5*n)
! allocatable arrays
real(8), allocatable :: rwork(:)
complex(8), allocatable :: work(:)
! external functions
integer, external :: ilaenv
! find the optimal blocksize for allocating the work array
nb=ilaenv(1,'ZHETRD','U',n,-1,-1,-1)
lwork=max(2*n+n**2,(nb+1)*n)
lrwork=1+5*n+2*n**2
allocate(work(lwork),rwork(lrwork))
! enable MKL parallelism
call holdthd(maxthdmkl,nthd)
nts=mkl_set_num_threads_local(nthd)
call zheevd('V','U',n,a,ld,w,work,lwork,rwork,lrwork,iwork,3+5*n,info)
nts=mkl_set_num_threads_local(0)
call freethd(nthd)
if (info /= 0) then
  write(*,*)
  write(*,'("Error(eveqnzh): diagonalisation failed")')
  write(*,'(" ZHEEVD returned INFO = ",I8)') info
  write(*,*)
  stop
end if
deallocate(rwork,work)
end subroutine

