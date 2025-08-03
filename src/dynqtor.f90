
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynqtor(dq,dr)
use modmain
use modphonon
implicit none
! arguments
complex(8), intent(in) :: dq(nbph,nbph,nqpt)
real(8), intent(out) :: dr(nbph,nbph,nqptnr)
! local variables
integer intr(2,3),iq,n
integer isym,lspl,ir
integer i1,i2,i3,j1,j2,j3
real(8) v1(3),v2(3),s(3,3),t1
complex(8) z1
! automatic arrays
complex(8) dqs(nbph,nbph)
intr(2,:)=ngridq(:)/2
intr(1,:)=intr(2,:)-ngridq(:)+1
dr(:,:,:)=0.d0
! loop over q-vectors
do j1=0,ngridq(1)-1
  v1(1)=dble(j1)/dble(ngridq(1))
  do j2=0,ngridq(2)-1
    v1(2)=dble(j2)/dble(ngridq(2))
    do j3=0,ngridq(3)-1
      v1(3)=dble(j3)/dble(ngridq(3))
      iq=ivqiq(j1,j2,j3)
! map v1 to the first Brillouin zone
      call vecfbz(epslat,bvec,v1)
! rotate and add the dynamical matrix of the reduced q-point with all symmetries
      n=0
      dqs(:,:)=0.d0
      do isym=1,nsymcrys
        lspl=lsplsymc(isym)
        s(:,:)=dble(symlat(:,:,lspl))
        call r3mtv(s,vql(:,iq),v2)
        call vecfbz(epslat,bvec,v2)
        t1=abs(v1(1)-v2(1))+abs(v1(2)-v2(2))+abs(v1(3)-v2(3))
        if (t1 < epslat) then
          call dynsymapp(isym,vql(:,iq),dq(:,:,iq),dqs)
          n=n+1
        end if
      end do
      if (n == 0) then
        write(*,*)
        write(*,'("Error(dynqtor): vector ",3G18.10)') v1
        write(*,'(" cannot be mapped to reduced q-point set")')
        write(*,*)
        stop
      end if
      t1=1.d0/dble(n)
      dqs(:,:)=t1*dqs(:,:)
! subtract the non-analytic term if required
      if (tphnat) call dynqnat(-1,v1,dqs)
! loop over R-vectors
      ir=0
      do i3=intr(1,3),intr(2,3)
        do i2=intr(1,2),intr(2,2)
          do i1=intr(1,1),intr(2,1)
            ir=ir+1
            t1=twopi*(v1(1)*dble(i1)+v1(2)*dble(i2)+v1(3)*dble(i3))
            z1=cmplx(cos(t1),sin(t1),8)
            dr(:,:,ir)=dr(:,:,ir)+dble(z1*dqs(:,:))
          end do
        end do
      end do
    end do
  end do
end do
t1=1.d0/dble(nqptnr)
dr(:,:,:)=t1*dr(:,:,:)
end subroutine

