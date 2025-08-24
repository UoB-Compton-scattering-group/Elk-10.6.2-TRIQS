! Copyright (C) 2025 A. D. N. James. 
! This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine dmftdmatk(ik)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: ik
! local variables
integer ist,jst,iw,i,nthd
real(8) e,t1
complex(8) z1
! automatic arrays
complex(8) gs(nstsv),g(nstsv,nstsv),ge(4,nstsv,nstsv)
complex(8) evecsv(nstsv,nstsv),d(nstsv,nstsv),a(nstsv,nstsv)
! allocatable arrays
complex(8), allocatable :: se(:,:,:)
! external functions
complex(8), external :: gwtails
!adnj edit - sets the density matrix from dmft output
d(:,:)=dmatkdmft(:,:,ik)
! diagonalise the density matrix for the natural orbitals and occupation numbers
call eveqnzh(nstsv,nstsv,d,occsv(:,ik))
occsv(1:nstsv,ik)=occsv(1:nstsv,ik)*occmax
! get the second-variational eigenvectors from file
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! apply unitary transformation to the third-variational states so that they
! refer to the first-variational basis
call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,d,nstsv,zzero,a,nstsv)
! write the density matrix to file as second-variational eigenvectors
call putevecsv(filext,ik,a)
end subroutine

