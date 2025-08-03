
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bfieldfsm
! !INTERFACE:
subroutine bfieldfsm
! !USES:
use modmain
! !DESCRIPTION:
!   Updates the effective magnetic field, ${\bf B}_{\rm FSM}$, required for
!   fixing the spin moment to a given value, ${\bf M}_{\rm FSM}$. This is done
!   by adding a vector to the field which is proportional to the difference
!   between the moment calculated in the $i$th self-consistent loop and the
!   required moment:
!   $$ {\bf B}_{\rm FSM}^{i+1}={\bf B}_{\rm FSM}^i+\lambda\left({\bf M}^i
!    -{\bf M}_{\rm FSM}\right), $$
!   where $\lambda$ is a scaling factor.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias
real(8) v1(3),v2(3),t1,t2
if ((.not.spinpol).or.(fsmtype == 0)) return
! fixed spin direction not valid for collinear magnetism
if ((.not.ncmag).and.(fsmtype < 0)) return
! determine the global effective field
if ((abs(fsmtype) == 1).or.(abs(fsmtype) == 3)) then
  if (ncmag) then
    v1(:)=momtot(:)
  else
    v1(:)=0.d0
    v1(3)=momtot(1)
  end if
  v2(:)=v1(:)-momfix(:)
  if (ncmag) then
    bfsmc(:)=bfsmc(:)+taufsm*v2(:)
  else
    bfsmc(1)=bfsmc(1)+taufsm*v2(3)
  end if
! make sure that the constraining field is perpendicular to the fixed moment
! for fixed direction calculations (Y. Kvashnin and LN)
  if (fsmtype < 0) call r3vo(momfix,bfsmc)
end if
! determine the muffin-tin fields for fixed local moments
if ((abs(fsmtype) == 2).or.(abs(fsmtype) == 3)) then
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! if any component is >= 1000 then do not fix the moment
      if (any(abs(mommtfix(1:3,ia,is)) >= 1000.d0)) cycle
      if (ncmag) then
        v1(:)=mommt(:,ias)
      else
        v1(:)=0.d0
        v1(3)=mommt(1,ias)
      end if
      v2(:)=v1(:)-mommtfix(:,ia,is)
      if (ncmag) then
        bfsmcmt(:,ias)=bfsmcmt(:,ias)+taufsm*v2(:)
      else
        bfsmcmt(1,ias)=bfsmcmt(1,ias)+taufsm*v2(3)
      end if
! fixed spin direction
      if (fsmtype < 0) call r3vo(mommtfix(:,ia,is),bfsmcmt(:,ias))
    end do
  end do
end if
! global fixed spin magnitude
if ((fsmtype == 4).or.(fsmtype == 6)) then
  t1=taufsm*(momtotm-momfixm)
  bfsmc(1:ndmag)=bfsmc(1:ndmag)+t1*momtot(1:ndmag)
end if
! local fixed spin magnitude
if ((fsmtype == 5).or.(fsmtype == 6)) then
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      t1=mommtfixm(ia,is)
      if (t1 < 0.d0) cycle
      t2=norm2(mommt(1:ndmag,ias))
      t1=taufsm*(t2-t1)
      bfsmcmt(1:ndmag,ias)=bfsmcmt(1:ndmag,ias)+t1*mommt(1:ndmag,ias)
    end do
  end do
end if
end subroutine
!EOC

