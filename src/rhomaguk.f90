
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhomaguk(ik0,lock,evecu)
use modmain
use modulr
use modomp
implicit none
! arguments
integer, intent(in) :: ik0
integer(omp_lock_kind), intent(inout) :: lock(nqpt)
complex(8), intent(in) :: evecu(nstulr,nstulr)
! local variables
integer ik,ikpa,jkpa
integer nst,ist,jst,i,j
integer ngk0,is,ias
integer npc,ir,nthd
real(8) ts0,ts1
real(4) wo
! automatic arrays
integer idx(nstsv)
complex(8) zfft(nqpt)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: evectv(:,:,:),evecsvt(:,:)
complex(4), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
call timesec(ts0)
! central k-point
ik=(ik0-1)*nkpa+1
! number of G+k-vectors for central k-point
ngk0=ngk(1,ik)
! get the eigenvectors from file
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
call match(ngk0,vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),apwalm)
allocate(evectv(nstsv,nstsv,nqpt))
call holdthd(nqpt,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfft,evecsvt,wfmt,wfir) &
!$OMP PRIVATE(ikpa,jkpa,ist,jst,i,j) &
!$OMP PRIVATE(ir,wo,ias,is,npc) &
!$OMP NUM_THREADS(nthd)
allocate(evecsvt(nstsv,nstsv))
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfir(ngtc,nspinor,nstsv))
! loop over long-range states in subsets of size nstsv
do jkpa=1,nkpa
! number of and index to occupied states in subset
!$OMP SINGLE
  nst=0
  do jst=1,nstsv
    j=(jkpa-1)*nstsv+jst
    if (abs(occulr(j,ik0)) < epsocc) cycle
    nst=nst+1
    idx(nst)=jst
  end do
!$OMP END SINGLE
  if (nst == 0) cycle
!$OMP DO SCHEDULE(DYNAMIC)
  do jst=1,nst
    j=(jkpa-1)*nstsv+idx(jst)
    do ist=1,nstsv
      zfft(1:nqpt)=0.d0
      do ikpa=1,nkpa
        i=(ikpa-1)*nstsv+ist
! store the long-range state in FFT Q-space
        zfft(iqfft(ikpa))=evecu(i,j)
      end do
! Fourier transform to R-space
      call zfftifc(3,ngridq,1,zfft)
      evectv(ist,jst,1:nqpt)=zfft(1:nqpt)
    end do
  end do
!$OMP END DO
! parallel loop over R-points
!$OMP DO SCHEDULE(DYNAMIC)
  do ir=1,nqpt
! convert third-variational states to second-variational states
    call zgemm('N','N',nstsv,nst,nstsv,zone,evecsv,nstsv,evectv(:,:,ir), &
     nstsv,zzero,evecsvt,nstsv)
! generate the wavefunctions in single-precision
    call genwfsv_sp(.false.,.false.,nst,[0],ngdgc,igfc,ngk0,igkig(:,1,ik), &
     apwalm,evecfv,evecsvt,wfmt,ngtc,wfir)
! loop over second-variational states
    do jst=1,nst
      j=(jkpa-1)*nstsv+idx(jst)
      wo=occulr(j,ik0)*wkpt(ik)
! add to the density and magnetisation
      call omp_set_lock(lock(ir))
! muffin-tin part
      do ias=1,natmtot
        is=idxis(ias)
        npc=npcmt(is)
        if (spinpol) then
          if (ncmag) then
            call rmk1(npc,wo,wfmt(:,ias,1,jst),wfmt(:,ias,2,jst), &
             rhormt(:,ias,ir),magrmt(:,ias,1,ir),magrmt(:,ias,2,ir), &
             magrmt(:,ias,3,ir))
          else
            call rmk2(npc,wo,wfmt(:,ias,1,jst),wfmt(:,ias,2,jst), &
             rhormt(:,ias,ir),magrmt(:,ias,1,ir))
          end if
        else
          call rmk3(npc,wo,wfmt(:,ias,1,jst),rhormt(:,ias,ir))
        end if
      end do
! interstitial part
      if (spinpol) then
        if (ncmag) then
          call rmk1(ngtc,wo,wfir(:,1,jst),wfir(:,2,jst),rhorir(:,ir), &
           magrir(:,1,ir),magrir(:,2,ir),magrir(:,3,ir))
        else
          call rmk2(ngtc,wo,wfir(:,1,jst),wfir(:,2,jst),rhorir(:,ir), &
           magrir(:,1,ir))
        end if
      else
        call rmk3(ngtc,wo,wfir(:,1,jst),rhorir(:,ir))
      end if
      call omp_unset_lock(lock(ir))
    end do
! end parallel loop over R-points
  end do
!$OMP END DO
end do
deallocate(evecsvt,wfmt,wfir)
!$OMP END PARALLEL
call freethd(nthd)
deallocate(apwalm,evecfv,evecsv,evectv)
call timesec(ts1)
!$OMP ATOMIC
timerho=timerho+ts1-ts0
return

contains

pure subroutine rmk1(n,wo,wf1,wf2,rho,mag1,mag2,mag3)
implicit none
! arguments
integer, intent(in) :: n
real(4), intent(in) :: wo
complex(4), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag1(n),mag2(n),mag3(n)
! local variables
integer i
real(4) wo2,t1,t2
real(4) a1,b1,a2,b2
wo2=2.e0*wo
!$OMP SIMD PRIVATE(a1,b1,a2,b2,t1,t2) SIMDLEN(8)
do i=1,n
  a1=real(wf1(i)); b1=aimag(wf1(i))
  a2=real(wf2(i)); b2=aimag(wf2(i))
  t1=a1**2+b1**2; t2=a2**2+b2**2
  mag1(i)=mag1(i)+wo2*(a1*a2+b1*b2)
  mag2(i)=mag2(i)+wo2*(a1*b2-b1*a2)
  mag3(i)=mag3(i)+wo*(t1-t2)
  rho(i)=rho(i)+wo*(t1+t2)
end do
end subroutine

pure subroutine rmk2(n,wo,wf1,wf2,rho,mag)
implicit none
! arguments
integer, intent(in) :: n
real(4), intent(in) :: wo
complex(4), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag(n)
! local variables
integer i
real(4) t1,t2
!$OMP SIMD PRIVATE(t1,t2) SIMDLEN(8)
do i=1,n
  t1=real(wf1(i))**2+aimag(wf1(i))**2
  t2=real(wf2(i))**2+aimag(wf2(i))**2
  mag(i)=mag(i)+wo*(t1-t2)
  rho(i)=rho(i)+wo*(t1+t2)
end do
end subroutine

pure subroutine rmk3(n,wo,wf,rho)
implicit none
! arguments
integer, intent(in) :: n
real(4), intent(in) :: wo
complex(4), intent(in) :: wf(n)
real(8), intent(inout) :: rho(n)
rho(1:n)=rho(1:n)+wo*(real(wf(1:n))**2+aimag(wf(1:n))**2)
end subroutine

end subroutine

