
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wfmtsv_sp(tsh,is,ias,nst,idx,ngp,apwalm,evecfv,evecsv,ld,wfmt)
use modmain
use modomp
implicit none
! arguments
logical, intent(in) :: tsh
integer, intent(in) :: is,ias,nst,idx(*),ngp(nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
integer, intent(in) :: ld
complex(4), intent(out) :: wfmt(ld,nspinor,nst)
! local variables
logical tasv
integer io,ilo,ispn,jspn
integer nrc,nrci,nrco,iro
integer l,lm,npc,npci
integer n,p,i,j,k,nthd
complex(8) zq(2),z1
! automatic arrays
complex(8) x(nstfv,nspnfv),y(nlmwf(is),nspinor,nst)
! external functions
complex(8), external :: zdotu
iro=nrmti(is)+lradstp
nrc=nrcmt(is)
nrci=nrcmti(is)
nrco=nrc-nrci
npc=npcmt(is)
npci=npcmti(is)
! de-phasing factor for spin-spirals
if (ssdph) then
  zq(1)=zqss(ias)
  zq(2)=conjg(zq(1))
end if
! check if all the second-variational wavefunctions should be calculated
tasv=(idx(1) == 0)
call holdthd(nst,nthd)
!-----------------------!
!     APW functions     !
!-----------------------!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(p,l,lm,io,ispn,jspn) &
!$OMP PRIVATE(n,i,j,k,z1,ilo) &
!$OMP NUM_THREADS(nthd)
p=0
do l=0,lmaxo
  do lm=l**2+1,(l+1)**2
    do io=1,apword(l,is)
      p=p+1
      if (tevecsv) then
        do jspn=1,nspnfv
          n=ngp(jspn)
!$OMP DO SCHEDULE(DYNAMIC)
          do j=1,nstfv
            x(j,jspn)=zdotu(n,evecfv(:,j,jspn),1,apwalm(:,io,lm,ias,jspn),1)
          end do
!$OMP END DO
        end do
! loop only over required states
!$OMP DO SCHEDULE(DYNAMIC)
        do j=1,nst
! index to state in evecsv
          if (tasv) then; k=j; else; k=idx(j); end if
          y(p,1,j)=zdotu(nstfv,evecsv(1,k),1,x,1)
          if (spinpol) then
            jspn=jspnfv(2)
            y(p,2,j)=zdotu(nstfv,evecsv(nstfv+1,k),1,x(1,jspn),1)
          end if
        end do
!$OMP END DO
      else
!$OMP DO SCHEDULE(DYNAMIC)
        do j=1,nst
          if (tasv) then; k=j; else; k=idx(j); end if
          y(p,1,j)=zdotu(ngp(1),evecfv(:,k,1),1,apwalm(:,io,lm,ias,1),1)
        end do
!$OMP END DO
      end if
    end do
  end do
end do
!$OMP DO SCHEDULE(DYNAMIC)
do j=1,nst
  wfmt(1:npc,1:nspinor,j)=0.e0
  do ispn=1,nspinor
    p=0
    do l=0,lmaxo
      do lm=l**2+1,(l+1)**2
        i=npci+lm
        do io=1,apword(l,is)
          p=p+1
          z1=y(p,ispn,j)
          if (abs(z1%re)+abs(z1%im) < 1.d-10) cycle
          if (ssdph) z1=z1*zq(ispn)
          if (l <= lmaxi) then
            call cfzrf(nrci,z1,apwfr(1,1,io,l,ias),lmmaxi,wfmt(lm,ispn,j))
          end if
          call cfzrf(nrco,z1,apwfr(iro,1,io,l,ias),lmmaxo,wfmt(i,ispn,j))
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
p=0
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do lm=l**2+1,(l+1)**2
    p=p+1
    i=idxlo(lm,ilo,ias)
    if (tevecsv) then
      do jspn=1,nspnfv
        n=ngp(jspn)
        x(1:nstfv,jspn)=evecfv(n+i,1:nstfv,jspn)
      end do
!$OMP DO SCHEDULE(DYNAMIC)
      do j=1,nst
        if (tasv) then; k=j; else; k=idx(j); end if
        y(p,1,j)=zdotu(nstfv,evecsv(1,k),1,x,1)
        if (spinpol) then
          jspn=jspnfv(2)
          y(p,2,j)=zdotu(nstfv,evecsv(nstfv+1,k),1,x(1,jspn),1)
        end if
      end do
!$OMP END DO
    else
      do j=1,nst
        if (tasv) then; k=j; else; k=idx(j); end if
        y(p,1,j)=evecfv(ngp(1)+i,k,1)
      end do
    end if
  end do
end do
!$OMP DO SCHEDULE(DYNAMIC)
do j=1,nst
  do ispn=1,nspinor
    p=0
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      do lm=l**2+1,(l+1)**2
        p=p+1
        i=npci+lm
        z1=y(p,ispn,j)
        if (abs(z1%re)+abs(z1%im) < 1.d-10) cycle
        if (ssdph) z1=z1*zq(ispn)
        if (l <= lmaxi) then
          call cfzrf(nrci,z1,lofr(1,1,ilo,ias),lmmaxi,wfmt(lm,ispn,j))
        end if
        call cfzrf(nrco,z1,lofr(iro,1,ilo,ias),lmmaxo,wfmt(i,ispn,j))
      end do
    end do
! convert to spherical coordinates if required
    if (.not.tsh) call cbshtip(nrc,nrci,wfmt(:,ispn,j))
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
return

contains

pure subroutine cfzrf(n,z,rf,ld,cf)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: z
real(8), intent(in) :: rf(lradstp,n)
integer, intent(in) :: ld
complex(4), intent(inout) :: cf(ld,n)
cf(1,1:n)=cf(1,1:n)+z*rf(1,1:n)
end subroutine

end subroutine

