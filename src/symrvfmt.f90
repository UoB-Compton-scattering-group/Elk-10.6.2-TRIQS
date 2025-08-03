
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symrvfmt(tspin,tnc,nrmt_,nrmti_,npmt_,ld,rvfmt)
use modmain
implicit none
! arguments
logical, intent(in) :: tspin,tnc
integer, intent(in) :: nrmt_(nspecies),nrmti_(nspecies),npmt_(nspecies),ld
real(8), intent(inout) :: rvfmt(ld,natmtot,*)
! local variables
integer nd,isym,lspl,lspn
integer is,ia,ja,ias,jas
integer nr,nri,np,i
real(8) sc(3,3),t0,t1
real(8) x1,x2,x3,y1,y2,y3
! allocatable arrays
real(8), allocatable :: rvfmt1(:,:,:),rvfmt2(:,:)
! dimension of the vector field
if (tnc) then
  nd=3
else
  nd=1
end if
allocate(rvfmt1(ld,natmmax,nd),rvfmt2(ld,nd))
t0=1.d0/dble(nsymcrys)
do is=1,nspecies
  nr=nrmt_(is)
  nri=nrmti_(is)
  np=npmt_(is)
! make copy of vector field for all atoms of current species
  do i=1,nd
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      rvfmt1(1:np,ia,i)=rvfmt(1:np,ias,i)
    end do
  end do
  do ia=1,natoms(is)
! only symmetrise first equivalent atom (rotate into others)
    if (.not.tfeqat(ia,is)) cycle
    ias=idxas(ia,is)
    rvfmt(1:np,ias,1:nd)=0.d0
! begin loop over crystal symmetries
    do isym=1,nsymcrys
! equivalent atom
      ja=ieqatom(ia,is,isym)
! parallel transport of vector field
      lspl=lsplsymc(isym)
      do i=1,nd
        call rotrfmt(symlatc(:,:,lspl),nr,nri,rvfmt1(:,ja,i),rvfmt2(:,i))
      end do
      if (tspin) then
! global spin proper rotation matrix in Cartesian coordinates
        lspn=lspnsymc(isym)
        sc(1:3,1:3)=symlatd(lspn)*symlatc(1:3,1:3,lspn)
      else
! set spin rotation equal to spatial rotation
        lspn=lspl
        sc(1:3,1:3)=symlatc(1:3,1:3,lspl)
      end if
! global spin rotation of vector field
      if (tnc) then
! non-collinear case
        do i=1,np
          x1=rvfmt2(i,1); x2=rvfmt2(i,2); x3=rvfmt2(i,3)
          y1=sc(1,1)*x1+sc(1,2)*x2+sc(1,3)*x3
          y2=sc(2,1)*x1+sc(2,2)*x2+sc(2,3)*x3
          y3=sc(3,1)*x1+sc(3,2)*x2+sc(3,3)*x3
          rvfmt(i,ias,1)=rvfmt(i,ias,1)+y1
          rvfmt(i,ias,2)=rvfmt(i,ias,2)+y2
          rvfmt(i,ias,3)=rvfmt(i,ias,3)+y3
        end do
      else
! collinear case
        t1=sc(3,3)
        rvfmt(1:np,ias,1)=rvfmt(1:np,ias,1)+t1*rvfmt2(1:np,1)
      end if
! end loop over crystal symmetries
    end do
! normalise
    do i=1,nd
      rvfmt(1:np,ias,i)=t0*rvfmt(1:np,ias,i)
    end do
! rotate into equivalent atoms
    do ja=1,natoms(is)
      if (eqatoms(ia,ja,is).and.(ia /= ja)) then
        isym=findloc(ieqatom(ia,is,1:nsymcrys),ja,1)
        jas=idxas(ja,is)
! parallel transport of vector field (using operation inverse)
        lspl=isymlat(lsplsymc(isym))
        do i=1,nd
          call rotrfmt(symlatc(:,:,lspl),nr,nri,rvfmt(:,ias,i),rvfmt(:,jas,i))
        end do
        if (tspin) then
! inverse of global proper rotation matrix in Cartesian coordinates
          lspn=isymlat(lspnsymc(isym))
          sc(1:3,1:3)=symlatd(lspn)*symlatc(1:3,1:3,lspn)
        else
! set spin rotation equal to spatial rotation
          lspn=lspl
          sc(1:3,1:3)=symlatc(1:3,1:3,lspl)
        end if
! global spin rotation of vector field
        if (tnc) then
! non-collinear case
          do i=1,np
            x1=rvfmt(i,jas,1); x2=rvfmt(i,jas,2); x3=rvfmt(i,jas,3)
            y1=sc(1,1)*x1+sc(1,2)*x2+sc(1,3)*x3
            y2=sc(2,1)*x1+sc(2,2)*x2+sc(2,3)*x3
            y3=sc(3,1)*x1+sc(3,2)*x2+sc(3,3)*x3
            rvfmt(i,jas,1)=y1; rvfmt(i,jas,2)=y2; rvfmt(i,jas,3)=y3
          end do
        else
! collinear case
          t1=sc(3,3)
          rvfmt(1:np,jas,1)=t1*rvfmt(1:np,jas,1)
        end if
      end if
    end do
! end loop over atoms and species
  end do
end do
deallocate(rvfmt1,rvfmt2)
end subroutine

