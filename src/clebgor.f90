
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: clebgor
! !INTERFACE:
real(8) function clebgor(j1,j2,j3,m1,m2,m3)
! !INPUT/OUTPUT PARAMETERS:
!   j1, j2, j3 : angular momentum quantum numbers (in,integer)
!   m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the Clebsch-Gordon coefficients using the Wigner $3j$-symbols
!   $$ C(J_1 J_2 J_3 | m_1 m_2 m_3)=(-1)^{J_1-J_2+m_3}\sqrt{2J_3+1}
!    \begin{pmatrix} J_1 & J_2 & J_3 \\ m_1 & m_2 & -m_3 \end{pmatrix}. $$
!   Suitable for $J_i\le 50$. See {\tt wigner3j}.
!
! !REVISION HISTORY:
!   Created September 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: j1,j2,j3
integer, intent(in) :: m1,m2,m3
! external functions
real(8), external :: wigner3j
if ((j1 < 0).or.(j2 < 0).or.(j3 < 0).or.(abs(m1) > j1).or.(abs(m2) > j2) &
 .or.(abs(m3) > j3)) then
  write(*,*)
  write(*,'("Error(clebgor): non-physical arguments :")')
  write(*,'("j1 = ",I8," j2 = ",I8," j3 = ",I8)') j1,j2,j3
  write(*,'("m1 = ",I8," m2 = ",I8," m3 = ",I8)') m1,m2,m3
  write(*,*)
  stop
end if
if ((j1 == 0).and.(j2 == 0).and.(j3 == 0)) then
  clebgor=1.d0
  return
end if
if ((j1 > 50).or.(j2 > 50).or.(j3 > 50)) then
  write(*,*)
  write(*,'("Error(clebgor): angular momenta out of range : ",3I8)') j1,j2,j3
  write(*,*)
  stop
end if
if ((m1+m2 /= m3).or.(j1+j2 < j3).or.(j2+j3 < j1).or.(j1+j3 < j2)) then
  clebgor=0.d0
  return
end if
clebgor=sqrt(dble(2*j3+1))*wigner3j(j1,j2,j3,m1,m2,-m3)
if (mod(j1-j2+m3,2) /= 0) clebgor=-clebgor
end function
!EOC

