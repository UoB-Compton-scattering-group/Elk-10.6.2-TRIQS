
! Copyright (C) 2002-2009 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modmain

!----------------------------!
!     lattice parameters     !
!----------------------------!
! lattice vectors stored column-wise
real(8) avec(3,3),avec0(3,3),davec(3,3)
! inverse of lattice vector matrix
real(8) ainv(3,3)
! reciprocal lattice vectors
real(8) bvec(3,3),bvec0(3,3)
! inverse of reciprocal lattice vector matrix
real(8) binv(3,3),binv0(3,3)
! unit cell volume
real(8) omega,omega0
! Brillouin zone volume
real(8) omegabz
! any vector with length less than epslat is considered zero
real(8) epslat

!--------------------------!
!     atomic variables     !
!--------------------------!
! maximum allowed species
integer, parameter :: maxspecies=8
! maximum allowed atoms per species
integer, parameter :: maxatoms=200
! number of species
integer nspecies
! number of atoms for each species
integer natoms(maxspecies),natoms0(maxspecies)
! maximum number of atoms over all the species
integer natmmax
! total number of atoms
integer natmtot,natmtot0
! index to atoms and species
integer idxas(maxatoms,maxspecies)
! inverse atoms and species indices
integer idxis(maxatoms*maxspecies),idxis0(maxatoms*maxspecies)
integer idxia(maxatoms*maxspecies)
! molecule is .true. is the system is an isolated molecule
logical molecule
! primcell is .true. if primitive unit cell is to be found automatically
logical primcell,primcell0
! atomic positions in lattice coordinates
real(8) atposl(3,maxatoms,maxspecies),atposl0(3,maxatoms,maxspecies)
real(8) datposl(3,maxatoms,maxspecies)
! atomic positions in Cartesian coordinates
real(8) atposc(3,maxatoms,maxspecies),atposc0(3,maxatoms,maxspecies)
! magnitude of random displacements added to the atomic positions
real(8) rndatposc
! tatdisp is .true. if small amplitude atomic displacements are to be included
! when calculating the Coulomb potential
logical :: tatdisp=.false.
! trdatdv is .true. if the atomic displacements and velocities are to be read
! from file
logical trdatdv
! atomic displacements and velocities in Cartesian coordinates
real(8) atdvc(3,0:1,maxatoms,maxspecies)
! atomic damping force coefficient
real(8) atdfc

!----------------------------------!
!     atomic species variables     !
!----------------------------------!
! species files path
character(256) sppath
! species filenames
character(256) spfname(maxspecies)
! species name
character(64) spname(maxspecies)
! species symbol
character(64) spsymb(maxspecies)
! species nuclear charge
real(8) spzn(maxspecies)
! ptnucl is .true. if the nuclei are to be treated as point charges, if .false.
! the nuclei have a finite spherical distribution
logical ptnucl
! nuclear radius
real(8) rnucl(maxspecies)
! nuclear volume
real(8) volnucl(maxspecies)
! number of radial mesh points to nuclear radius
integer nrnucl(maxspecies)
! Thomson radius
real(8) rtmsn(maxspecies)
! Thomson volume
real(8) voltmsn(maxspecies)
! number of radial mesh points to Thomson radius
integer nrtmsn(maxspecies)
! nuclear Coulomb potential
real(8), allocatable :: vcln(:,:)
! species electronic charge
real(8) spze(maxspecies)
! species mass
real(8) spmass(maxspecies)
! smallest radial point for each species
real(8) rminsp(maxspecies)
! effective infinity for species
real(8) rmaxsp(maxspecies)
! number of radial points to effective infinity for each species
integer nrsp(maxspecies)
! maximum nrsp over all the species
integer nrspmax
! maximum allowed states for each species
integer, parameter :: maxstsp=40
! number of states for each species
integer nstsp(maxspecies)
! maximum nstsp over all the species
integer nstspmax
! core-valence cut-off energy for species file generation
real(8) ecvcut
! semi-core-valence cut-off energy for species file generation
real(8) esccut
! state principle quantum number for each species
integer nsp(maxstsp,maxspecies)
! state l value for each species
integer lsp(maxstsp,maxspecies)
! state k value for each species
integer ksp(maxstsp,maxspecies)
! spcore is .true. if species state is core
logical spcore(maxstsp,maxspecies)
! total number of core states
integer nstcr
! state eigenvalue for each species
real(8) evalsp(maxstsp,maxspecies)
! state occupancy for each species
real(8) occsp(maxstsp,maxspecies)
! species radial mesh to effective infinity
real(8), allocatable :: rsp(:,:)
! species charge density
real(8), allocatable :: rhosp(:,:)
! species self-consistent potential
real(8), allocatable :: vrsp(:,:)
! exchange-correlation type for atomic species (the converged ground-state of
! the crystal does not depend on this choice)
integer xctsp(3)

!---------------------------------------------------------------!
!     muffin-tin radial mesh and angular momentum variables     !
!---------------------------------------------------------------!
! scale factor for number of muffin-tin points
real(8) nrmtscf,dnrmtscf
! number of muffin-tin radial points for each species
integer nrmt(maxspecies)
! maximum nrmt over all the species
integer nrmtmax
! muffin-tin radius scale factor
real(8) rmtscf
! order of averaging applied to the muffin-tin radii
integer mrmtav
! optional default muffin-tin radius for all atoms
real(8) rmtall
! minimum allowed distance between muffin-tin surfaces
real(8) rmtdelta
! muffin-tin radii
real(8) rmt(maxspecies),rmt0(maxspecies)
! trmt0 is .true. if the original muffin-tin radii rmt0 are to be retained
! between tasks
logical trmt0
! (R_mt)ˡ for l up to lmaxo+3
real(8), allocatable :: rmtl(:,:)
! total muffin-tin volume
real(8) omegamt
! radial step length for coarse mesh
integer lradstp
! number of coarse radial mesh points
integer nrcmt(maxspecies)
! maximum nrcmt over all the species
integer nrcmtmax
! coarse muffin-tin radial mesh
real(8), allocatable :: rcmt(:,:)
! r^l on fine radial mesh
real(8), allocatable :: rlmt(:,:,:)
! r^l on coarse radial mesh
real(8), allocatable :: rlcmt(:,:,:)
! weights for spline integration on fine radial mesh multiplied by r²
real(8), allocatable :: wr2mt(:,:)
! weights for spline partial integration on fine radial mesh
real(8), allocatable :: wprmt(:,:,:)
! weights for spline coefficients on fine radial mesh
real(8), allocatable :: wcrmt(:,:,:)
! weights for spline integration on coarse radial mesh multiplied by r²
real(8), allocatable :: wr2cmt(:,:)
! weights for spline partial integration on coarse radial mesh
real(8), allocatable :: wprcmt(:,:,:)
! weights for spline coefficients on coarse radial mesh
real(8), allocatable :: wcrcmt(:,:,:)
! maximum allowable angular momentum for augmented plane waves
integer, parameter :: maxlapw=30
! maximum angular momentum for augmented plane waves
integer lmaxapw,dlmaxapw
! (lmaxapw+1)²
integer lmmaxapw
! maximum angular momentum on the outer part of the muffin-tin
integer lmaxo,dlmaxo
! (lmaxo+1)²
integer lmmaxo
! maximum angular momentum on the inner part of the muffin-tin
integer lmaxi,lmaxi0
! (lmaxi+1)²
integer lmmaxi
! fraction of muffin-tin radius which constitutes the inner part
real(8) fracinr
! number of fine/coarse radial points on the inner part of the muffin-tin
integer nrmti(maxspecies),nrcmti(maxspecies)
! number of fine/coarse points in packed muffin-tins
integer npmti(maxspecies),npmt(maxspecies)
integer npcmti(maxspecies),npcmt(maxspecies)
! maximum number of points over all packed muffin-tins
integer npmtmax,npcmtmax
! total number of muffin-tin points for all atoms
integer npcmttot
! index to first muffin-tin point in packed array for all atoms
integer, allocatable :: ipcmt(:)
! smoothing order used when calculating gradients in the muffin-tin
integer msmgmt

!--------------------------------!
!     spin related variables     !
!--------------------------------!
! spinpol is .true. for spin-polarised calculations
logical spinpol,spinpol0
! spinorb is .true. for spin-orbit coupling
logical spinorb,spinorb0
! scale factor of spin-orbit coupling term in Hamiltonian
real(8) socscf
! bforb is .true. for external B-field-orbit coupling
logical bforb
! bfdmag is .true. for external B-field diamagnetic coupling
logical bfdmag
! dimension of magnetisation and magnetic vector fields (1 or 3)
integer ndmag
! ncmag is .true. if the magnetisation is non-collinear, i.e. when ndmag = 3
logical ncmag
! if cmagz is .true. then collinear magnetism along the z-axis is enforced
logical cmagz,cmagz0
! fixed spin moment type
!  0      : none
!  1 (-1) : total moment (direction)
!  2 (-2) : individual muffin-tin moments (direction)
!  3 (-3) : total and muffin-tin moments (direction)
!  4      : total moment magnitude
!  5      : individual muffin-tin moment magnitudes
!  6      : total and muffin-tin moment magnitudes
integer fsmtype,fsmtype0
! fixed total spin magnetic moment
real(8) momfix(3),momfix0(3),dmomfix(3)
! fixed total spin magnetic moment magnitude
real(8) momfixm
! fixed spin moment global effective field in Cartesian coordinates
real(8) bfsmc(3)
! muffin-tin fixed spin moments
real(8) mommtfix(3,maxatoms,maxspecies),mommtfix0(3,maxatoms,maxspecies)
! muffin-tin fixed spin moment magnitudes
real(8) mommtfixm(maxatoms,maxspecies)
! muffin-tin fixed spin moment effective fields in Cartesian coordinates
real(8), allocatable :: bfsmcmt(:,:)
! fixed spin moment field step size
real(8) taufsm
! second-variational spinor dimension (1 or 2)
integer nspinor
! global external magnetic field in Cartesian coordinates
real(8) bfieldc(3)
! initial field
real(8) bfieldc0(3),bfieldc00(3),dbfieldc0(3)
! external magnetic field in each muffin-tin in Cartesian coordinates
real(8) bfcmt(3,maxatoms,maxspecies)
! initial field as read in from input file
real(8) bfcmt0(3,maxatoms,maxspecies),bfcmt00(3,maxatoms,maxspecies)
! magnitude of random vectors added to muffin-tin fields
real(8) rndbfcmt
! external magnetic fields are multiplied by reducebf after each s.c. loop
real(8) reducebf,reducebf0
! small change in magnetic field used for calculating the magnetoelectric tensor
real(8) deltabf
! spinsprl is .true. if a spin-spiral is to be calculated
logical spinsprl,spinsprl0
! ssdph is .true. if the muffin-tin spin-spiral magnetisation is de-phased
logical ssdph
! spin-spiral phase factor for each atom
complex(8), allocatable :: zqss(:)
! number of spin-dependent first-variational functions per state
integer nspnfv
! map from second- to first-variational spin index
integer jspnfv(2)
! spin-spiral q-vector in lattice coordinates
real(8) vqlss(3),dvqlss(3)
! spin-spiral q-vector in Cartesian coordinates
real(8) vqcss(3)
! current q-point in spin-spiral supercell calculation
integer iqss
! number of primitive unit cells in spin-spiral supercell
integer nscss
! number of fixed spin direction points on the sphere for finding the magnetic
! anisotropy energy (MAE)
integer npmae0,npmae
! (theta,phi) coordinates for each MAE direction
real(8), allocatable :: tpmae(:,:)

!----------------------------------------------------!
!     static electric field and vector potential     !
!----------------------------------------------------!
! tefield is .true. if a polarising constant electric field is applied
logical tefield
! electric field vector in Cartesian coordinates
real(8) efieldc(3)
! electric field vector in lattice coordinates
real(8) efieldl(3)
! average electric field in Cartesian coordinates in each muffin-tin
real(8), allocatable :: efcmt(:,:)
! maximum distance over which the electric field is applied
real(8) dmaxefc
! maximum allowed absolute value of the potential generated by efieldc
real(8) vmaxefc
! tafield is .true. if a constant vector potential is applied
logical tafield
! vector potential A-field in Cartesian coordinates which couples to the
! paramagnetic current
real(8) afieldc(3),afieldc0(3),dafieldc(3)
! A-field in lattice coordinates
real(8) afieldl(3)
! tafsp is .true. if a constant spin-dependent vector potential is applied
logical tafsp
! spin-dependent vector potential (3 x 3 tensor) in Cartesian coordinates
real(8) afspc(3,3),dafspc(3,3)

!----------------------------!
!     symmetry variables     !
!----------------------------!
! type of symmetry allowed for the crystal
!  0 : only the identity element is used
!  1 : full symmetry group is used
!  2 : only symmorphic symmetries are allowed
integer symtype
! number of Bravais lattice point group symmetries
integer nsymlat
! Bravais lattice point group symmetries
integer symlat(3,3,48)
! determinants of lattice symmetry matrices (1 or -1)
integer symlatd(48)
! index to inverses of the lattice symmetries
integer isymlat(48)
! lattice point group symmetries in Cartesian coordinates
real(8) symlatc(3,3,48)
! tshift is .true. if atomic basis is allowed to be shifted
logical tshift,tshift0
! tsyminv is .true. if the crystal has inversion symmetry
logical tsyminv
! maximum of symmetries allowed
integer, parameter :: maxsymcrys=192
! number of crystal symmetries
integer nsymcrys
! crystal symmetry translation vector in lattice and Cartesian coordinates
real(8) vtlsymc(3,maxsymcrys),vtcsymc(3,maxsymcrys)
! tv0symc is .true. if the translation vector is zero
logical tv0symc(maxsymcrys)
! spatial rotation element in lattice point group for each crystal symmetry
integer lsplsymc(maxsymcrys)
! global spin rotation element in lattice point group for each crystal symmetry
integer lspnsymc(maxsymcrys)
! equivalent atom index for each crystal symmetry
integer, allocatable :: ieqatom(:,:,:)
! eqatoms(ia,ja,is) is .true. if atoms ia and ja are equivalent
logical, allocatable :: eqatoms(:,:,:)
! tfeqat is .true. if this is the first atom in a subset of equivalent atoms
logical, allocatable :: tfeqat(:,:)
! number of site symmetries
integer, allocatable :: nsymsite(:)
! site symmetry spatial rotation element in lattice point group
integer, allocatable :: lsplsyms(:,:)
! site symmetry global spin rotation element in lattice point group
integer, allocatable :: lspnsyms(:,:)

!----------------------------!
!     G-vector variables     !
!----------------------------!
! G-vector cut-off for interstitial potential and density
real(8) gmaxvr,dgmaxvr
! G-vector grid sizes
integer ngridg(3),ngridg0(3)
! G-vector grid sizes for coarse grid with |G| < 2 gkmax
integer ngdgc(3)
! total number of G-vectors
integer ngtot,ngtot0
! total number of G-vectors for coarse grid
integer ngtc
! integer grid intervals for each direction
integer intgv(2,3)
! number of G-vectors with |G| < gmaxvr
integer ngvec
! number of G-vectors for coarse grid with |G| < 2 gkmax
integer ngvc
! G-vector integer coordinates (i1,i2,i3)
integer, allocatable :: ivg(:,:),ivg0(:,:)
! map from (i1,i2,i3) to G-vector index
integer, allocatable :: ivgig(:,:,:)
! number of prime factors for the G-vector FFT
integer npfftg
! map from G-vector index to FFT array
integer, allocatable :: igfft(:),igfft0(:)
! number of prime factors for the coarse G-vector FFT
integer npfftgc
! map from G-vector index to FFT array for coarse grid
integer, allocatable :: igfc(:)
! number of complex FFT elements for real-complex transforms
integer nfgrz
! number of elements on the coarse grid
integer nfgrzc
! map from real-complex FFT index to G-point index
integer, allocatable :: igrzf(:)
! map on the coarse grid
integer, allocatable :: igrzfc(:)
! G-vectors in Cartesian coordinates
real(8), allocatable :: vgc(:,:)
! length of G-vectors
real(8), allocatable :: gc(:)
! Coulomb Green's function in G-space = 4π / G²
real(8), allocatable :: gclg(:)
! spherical Bessel functions jₗ(|G|Rₘₜ)
real(8), allocatable :: jlgrmt(:,:,:)
! spherical harmonics of the G-vectors
complex(8), allocatable :: ylmg(:,:)
! structure factors for the G-vectors
complex(8), allocatable :: sfacg(:,:)
! smooth step function form factors for all species and G-vectors
real(8), allocatable :: ffacg(:,:)
! characteristic function in G-space: 0 inside the muffin-tins and 1 outside
complex(8), allocatable :: cfunig(:)
! characteristic function in real-space: 0 inside the muffin-tins and 1 outside
real(8), allocatable :: cfunir(:)
! characteristic function in real-space for coarse grid
real(8), allocatable :: cfrc(:)

!---------------------------!
!     k-point variables     !
!---------------------------!
! autokpt is .true. if the k-point set is determined automatically
logical autokpt,autokpt0
! radius of sphere used to determine k-point density when autokpt is .true.
real(8) radkpt
! k-point grid sizes
integer ngridk(3),ngridk0(3),dngridk(3)
! k-point offset
real(8) vkloff(3),vkloff0(3)
! type of reduction to perform on k-point set
!  0 : no reduction
!  1 : reduce with full crystal symmetry group
!  2 : reduce with symmorphic symmetries only
integer reducek,reducek0
! number of point group symmetries used for k-point reduction
integer nsymkpt
! point group symmetry matrices used for k-point reduction
integer symkpt(3,3,48)
! total number of reduced k-points
integer nkpt
! total number of non-reduced k-points
integer nkptnr
! locations of k-points on integer grid
integer, allocatable :: ivk(:,:)
! map from integer grid to reduced k-point index
integer, allocatable :: ivkik(:,:,:)
! map from integer grid to non-reduced k-point index
integer, allocatable :: ivkiknr(:,:,:)
! k-points in lattice coordinates
real(8), allocatable :: vkl(:,:)
! k-points in Cartesian coordinates
real(8), allocatable :: vkc(:,:)
! reduced k-point weights
real(8), allocatable :: wkpt(:)
! weight of each non-reduced k-point
real(8) wkptnr
! k-point at which to determine effective mass tensor
real(8) vklem(3)
! displacement size for computing the effective mass tensor
real(8) deltaem
! number of displacements in each direction
integer ndspem
! number of k-points subdivision used for calculating the polarisation phase
integer nkspolar

!------------------------------!
!     G+k-vector variables     !
!------------------------------!
! species for which the muffin-tin radius will be used for calculating gkmax
integer isgkmax
! smallest muffin-tin radius times gkmax
real(8) rgkmax,drgkmax
! maximum |G+k| cut-off for APW functions
real(8) gkmax
! number of G+k-vectors for augmented plane waves
integer, allocatable :: ngk(:,:)
! maximum number of G+k-vectors over all k-points
integer ngkmax
! index from G+k-vectors to G-vectors
integer, allocatable :: igkig(:,:,:)
! G+k-vectors in lattice coordinates
real(8), allocatable :: vgkl(:,:,:,:)
! G+k-vectors in Cartesian coordinates
real(8), allocatable :: vgkc(:,:,:,:)
! length of G+k-vectors
real(8), allocatable :: gkc(:,:,:)
! structure factors for the G+k-vectors
complex(8), allocatable :: sfacgk(:,:,:,:)

!---------------------------!
!     q-point variables     !
!---------------------------!
! q-point grid sizes
integer ngridq(3)
! integer grid intervals for the q-points
integer intq(2,3)
! type of reduction to perform on q-point set (see reducek)
integer reduceq
! number of point group symmetries used for q-point reduction
integer nsymqpt
! point group symmetry matrices used for q-point reduction
integer symqpt(3,3,48)
! total number of reduced q-points
integer nqpt
! total number of non-reduced q-points
integer nqptnr
! locations of q-points on integer grid
integer, allocatable :: ivq(:,:)
! map from integer grid to reduced index
integer, allocatable :: ivqiq(:,:,:)
! map from integer grid to non-reduced index
integer, allocatable :: ivqiqnr(:,:,:)
! number of prime factors for the q-vector FFT
integer npfftq
! map from q-vector index to complex-complex FFT array
integer, allocatable :: iqfft(:)
! number of complex FFT elements for real-complex transforms
integer nfqrz
! map from q-point index to real-complex FFT index
integer, allocatable :: ifqrz(:)
! map from real-complex FFT index to q-point index
integer, allocatable :: iqrzf(:)
! q-points in lattice coordinates
real(8), allocatable :: vql(:,:)
! q-points in Cartesian coordinates
real(8), allocatable :: vqc(:,:)
! q-point weights
real(8), allocatable :: wqpt(:)
! weight for each non-reduced q-point
real(8) wqptnr
! regularised Coulomb Green's function in q-space
real(8), allocatable :: gclq(:)
! if t0gclq0 is .true. then the Coulomb Green's function at q = 0 is set to zero
logical t0gclq0

!-----------------------------------------------------!
!     spherical harmonic transform (SHT) matrices     !
!-----------------------------------------------------!
! trotsht is .true. if the spherical cover used for the SHT is to be rotated
logical :: trotsht=.false.
! spherical cover rotation matrix
real(8) rotsht(3,3)
! real backward SHT matrix for lmaxi
real(8), allocatable :: rbshti(:,:)
! real forward SHT matrix for lmaxi
real(8), allocatable :: rfshti(:,:)
! real backward SHT matrix for lmaxo
real(8), allocatable :: rbshto(:,:)
! real forward SHT matrix for lmaxo
real(8), allocatable :: rfshto(:,:)
! complex backward SHT matrix for lmaxi
complex(8), allocatable :: zbshti(:,:)
! complex forward SHT matrix for lmaxi
complex(8), allocatable :: zfshti(:,:)
! complex backward SHT matrix for lmaxo
complex(8), allocatable :: zbshto(:,:)
! complex forward SHT matrix for lmaxo
complex(8), allocatable :: zfshto(:,:)
! single-precision copies of the complex SHT matrices
complex(4), allocatable :: cbshti(:,:),cfshti(:,:)
complex(4), allocatable :: cbshto(:,:),cfshto(:,:)

!---------------------------------------------------------------!
!     density, potential and exchange-correlation variables     !
!---------------------------------------------------------------!
! exchange-correlation functional type
integer xctype(3)
! exchange-correlation functional description
character(264) xcdescr
! exchange-correlation functional spin requirement
integer xcspin
! exchange-correlation functional density gradient requirement
!  0 : no gradients
!  1 : gradients required for GGA with no post-processing: |∇ρ|, ∇²ρ,
!      (∇ρ).(∇|∇ρ|)
!  2 : gradients required for GGA with post-processing: |∇ρ|²
!  3 : as 2 but with the laplacian, ∇²ρ
!  4 : as 2 but with the kinetic energy density, τ
!  5 : as 4 but with the laplacian, ∇²ρ
!  6 : as 4 but for potential-only meta-GGA functionals
integer xcgrad
! small constant used to stabilise non-collinear GGA
real(8) dncgga
! kinetic energy density functional type
integer ktype(3)
! kinetic energy density functional description
character(264) kdescr
! kinetic energy density gradient requirement (see xcgrad)
integer kgrad
! combined target array for rhomt, rhoir, magmt and magir
real(8), allocatable, target :: rhmg(:)
! muffin-tin and interstitial charge density
real(8), pointer, contiguous :: rhomt(:,:),rhoir(:)
! muffin-tin and interstitial magnetisation vector field
real(8), pointer, contiguous :: magmt(:,:,:),magir(:,:)
! trhonorm is .true. if the density is to be normalised after every iteration
logical trhonorm
! tjr is .true. if the current density j(r) is to be calculated
logical tjr,tjr0
! muffin-tin and interstitial gauge-invariant current density vector field
real(8), allocatable :: jrmt(:,:,:),jrir(:,:)
! muffin-tin and interstitial Coulomb potential
real(8), allocatable :: vclmt(:,:),vclir(:)
! Poisson solver pseudocharge density constant
integer npsd
! lmaxo+npsd+1
integer lnpsd
! muffin-tin and interstitial exchange energy density
real(8), allocatable :: exmt(:,:),exir(:)
! muffin-tin and interstitial correlation energy density
real(8), allocatable :: ecmt(:,:),ecir(:)
! muffin-tin and interstitial exchange-correlation potential
real(8), allocatable :: vxcmt(:,:),vxcir(:)
! muffin-tin and interstitial exchange-correlation magnetic field
real(8), allocatable :: bxcmt(:,:,:),bxcir(:,:)
! muffin-tin and interstitial magnetic dipole field
real(8), allocatable :: bdmt(:,:,:),bdir(:,:)
! average dipole field in each muffin-tin
real(8), allocatable :: bdmta(:,:)
! tbdip is .true. if the spin and current dipole fields are to be added to the
! Kohn-Sham magnetic field
logical tbdip
! dipole magnetic field scaling factor (default 1)
real(8) bdipscf
! combined target array for vsmt, vsir, bsmt and bsir
real(8), allocatable, target :: vsbs(:)
! muffin-tin Kohn-Sham effective potential
real(8), pointer, contiguous :: vsmt(:,:)
! interstitial Kohn-Sham effective potential
real(8), allocatable :: vsir(:)
! vsir multiplied by the characteristic function and stored on a coarse grid
real(8), pointer, contiguous :: vsirc(:)
! muffin-tin Kohn-Sham effective magnetic field in spherical coordinates and on
! a coarse radial mesh
real(8), pointer, contiguous :: bsmt(:,:,:)
! interstitial Kohn-Sham effective magnetic field
real(8), allocatable :: bsir(:,:)
! bsir multiplied by the characteristic function and stored on a coarse grid
real(8), pointer, contiguous :: bsirc(:,:)
! G-space interstitial Kohn-Sham effective potential
complex(8), allocatable :: vsig(:)
! nosource is .true. if the field is to be made source-free
logical nosource
! tssxc is .true. if scaled spin exchange-correlation is to be used
logical tssxc
! spin exchange-correlation scaling factor
real(8) sxcscf,dsxcscf
! spin-orbit coupling radial function
real(8), allocatable :: socfr(:,:)
! kinetic energy density
real(8), allocatable :: taumt(:,:,:),tauir(:,:)
! core kinetic energy density
real(8), allocatable :: taucr(:,:,:)
! meta-GGA exchange-correlation potential
real(8), allocatable :: wxcmt(:,:),wxcir(:)
! Tran-Blaha '09 constant c [Phys. Rev. Lett. 102, 226401 (2009)]
real(8) c_tb09
! tc_tb09 is .true. if the Tran-Blaha constant has been read in
logical tc_tb09
! if trdstate is .true. the density and potential can be read from STATE.OUT
logical :: trdstate=.false.
! temperature in degrees Kelvin
real(8) tempk
! if mixrho is .true. then the (density, magnetisation) is mixed, otherwise the
! (potential, magnetic field)
logical mixrho
! mixing vector: either (density, magnetisation) or (potential, magnetic field)
real(8), pointer, contiguous :: vmixer(:)

!--------------------------!
!     mixing variables     !
!--------------------------!
! type of mixing to use for the potential
integer mixtype
! mixing type description
character(64) mixdescr
! if mixsave is .true. then the mixer work array is saved after each self-
! consistent loop and will be read in at the beginning of a restart
logical mixsave
! adaptive mixing parameters (formerly beta0 and betamax)
real(8) amixpm(2)
! subspace dimension for Broyden mixing
integer mixsdb
! Broyden mixing parameters alpha and w0
real(8) broydpm(2)

!----------------------------------------------!
!     charge, moment and current variables     !
!----------------------------------------------!
! tolerance for error in total charge
real(8) epschg
! total nuclear charge
real(8) chgzn
! core charges
real(8) chgcr(maxspecies)
! total core charge
real(8) chgcrtot
! core leakage charge
real(8), allocatable :: chgcrlk(:)
! total valence charge
real(8) chgval
! excess charge
real(8) chgexs,dchgexs
! total charge
real(8) chgtot
! calculated total charge
real(8) chgcalc
! interstitial region charge
real(8) chgir
! muffin-tin charges
real(8), allocatable :: chgmt(:)
! total muffin-tin charge
real(8) chgmttot
! effective Wigner radius
real(8) rwigner
! total moment
real(8) momtot(3)
! total moment magnitude
real(8) momtotm
! interstitial region moment
real(8) momir(3)
! muffin-tin moments
real(8), allocatable :: mommt(:,:)
! total muffin-tin moment
real(8) mommttot(3)
! total gauge-invariant current and its magnitude
real(8) jtot(3),jtotm

!-----------------------------------------!
!     APW and local-orbital variables     !
!-----------------------------------------!
! maximum allowable APW order
integer, parameter :: maxapword=3
! polynomial order used for APW radial derivatives
integer, parameter :: npapw=8
! APW order
integer apword(0:maxlapw,maxspecies)
! maximum of apword over all angular momenta and species
integer apwordmax
! total number of APW coefficients (l, m and order) for each species
integer lmoapw(maxspecies)
! energy step size used for APW numerical derivatives
real(8) deapw
! APW initial linearisation energies
real(8) apwe0(maxapword,0:maxlapw,maxspecies)
! APW linearisation energies
real(8), allocatable :: apwe(:,:,:)
! APW derivative order
integer apwdm(maxapword,0:maxlapw,maxspecies)
! apwve is .true. if the linearisation energies are allowed to vary
logical apwve(maxapword,0:maxlapw,maxspecies)
! APW radial functions
real(8), allocatable :: apwfr(:,:,:,:,:)
! derivate of radial functions at the muffin-tin surface multiplied by R_MT²/2
real(8), allocatable :: apwdfr(:,:,:)
! maximum number of local-orbitals
integer, parameter :: maxlorb=200
! maximum allowable local-orbital order
integer, parameter :: maxlorbord=4
! polynomial order used for local-orbital radial derivatives
integer, parameter :: nplorb=8
! number of local-orbitals
integer nlorb(maxspecies)
! maximum nlorb over all species
integer nlomax
! total number of local-orbitals
integer nlotot
! local-orbital order
integer lorbord(maxlorb,maxspecies)
! maximum lorbord over all species
integer lorbordmax
! local-orbital angular momentum
integer lorbl(maxlorb,maxspecies)
! maximum lorbl over all species
integer lolmax
! (lolmax+1)²
integer lolmmax
! energy step size used for local-orbital numerical derivatives
real(8) delorb
! local-orbital initial energies
real(8) lorbe0(maxlorbord,maxlorb,maxspecies)
! index which arranges the local-orbitals in ascending order of energy
integer idxelo(maxlorb,maxspecies)
! local-orbital energies
real(8), allocatable :: lorbe(:,:,:)
! local-orbital derivative order
integer lorbdm(maxlorbord,maxlorb,maxspecies)
! lorbve is .true. if the linearisation energies are allowed to vary
logical lorbve(maxlorbord,maxlorb,maxspecies)
! local-orbital radial functions
real(8), allocatable :: lofr(:,:,:,:)
! band energy search tolerance
real(8) epsband
! maximum allowed change in energy during band energy search; enforced only if
! default energy is less than zero
real(8) demaxbnd
! minimum default linearisation energy over all APWs and local-orbitals
real(8) e0min
! if autolinengy is .true. then the fixed linearisation energies are set to the
! Fermi energy minus dlefe
logical autolinengy
! difference between linearisation and Fermi energies when autolinengy is .true.
real(8) dlefe
! lorbcnd is .true. if conduction state local-orbitals should be added
logical lorbcnd
! conduction state local-orbital order
integer lorbordc
! excess order of the APW and local-orbital functions
integer nxoapwlo
! excess local orbitals
integer nxlo
! number of (l,m) components used for generating the muffin-tin wavefunctions
integer nlmwf(maxspecies)

!-------------------------------------------!
!     overlap and Hamiltonian variables     !
!-------------------------------------------!
! overlap and Hamiltonian matrices sizes at each k-point
integer, allocatable :: nmat(:,:)
! maximum nmat over all k-points
integer nmatmax
! index to the position of the local-orbitals in the H and O matrices
integer, allocatable :: idxlo(:,:,:)
! APW-local-orbital overlap integrals
real(8), allocatable :: oalo(:,:,:)
! local-orbital-local-orbital overlap integrals
real(8), allocatable :: ololo(:,:,:)
! APW-APW Hamiltonian integrals
real(8), allocatable :: haa(:,:,:,:,:,:)
! local-orbital-APW Hamiltonian integrals
real(8), allocatable :: hloa(:,:,:,:,:)
! local-orbital-local-orbital Hamiltonian integrals
real(8), allocatable :: hlolo(:,:,:,:)
! complex Gaunt coefficient array
complex(8), allocatable :: gntyry(:,:,:)
! tefvr is .true. if the first-variational eigenvalue equation is to be solved
! as a real symmetric problem
logical tefvr
! tefvit is .true. if the first-variational eigenvalue equation is to be solved
! iteratively
logical tefvit
! minimum and maximum allowed number of eigenvalue equation iterations
integer minitefv,maxitefv
! eigenvalue mixing parameter for iterative solver
real(8) befvit
! iterative solver convergence tolerance
real(8) epsefvit

!--------------------------------------------!
!     eigenvalue and occupancy variables     !
!--------------------------------------------!
! number of empty states per atom and spin
real(8) nempty0,dnempty0
! number of empty states
integer nempty
! number of first-variational states
integer nstfv
! number of second-variational states
integer nstsv
! smearing type
integer stype
! smearing function description
character(64) sdescr
! smearing width
real(8) swidth,swidth0
! autoswidth is .true. if the smearing width is to be determined automatically
logical autoswidth
! effective mass used in smearing width formula
real(8) mstar
! maximum allowed occupancy (1 or 2)
real(8) occmax
! convergence tolerance for occupation numbers
real(8) epsocc
! second-variational occupation numbers
real(8), allocatable :: occsv(:,:)
! Fermi energy for second-variational states
real(8) efermi
! tscissor is .true. if the scissor correction is non-zero
logical tscissor
! scissor correction applied to eigenvalues and momentum matrix elements
real(8) scissor
! density of states at the Fermi energy
real(8) fermidos
! estimated indirect and direct band gaps
real(8) bandgap(2)
! k-points of indirect and direct gaps
integer ikgap(3)
! error tolerance for the first-variational eigenvalues
real(8) evaltol
! second-variational eigenvalues
real(8), allocatable :: evalsv(:,:)
! tevecsv is .true. if second-variational eigenvectors are calculated
logical tevecsv
! maximum number of k-point and states indices in user-defined list
integer, parameter :: maxkst=20
! number of k-point and states indices in user-defined list
integer nkstlist
! user-defined list of k-point and state indices
integer kstlist(2,maxkst)

!------------------------------!
!     core state variables     !
!------------------------------!
! occupation numbers for core states
real(8), allocatable :: occcr(:,:)
! eigenvalues for core states
real(8), allocatable :: evalcr(:,:)
! radial wavefunctions for core states
real(8), allocatable :: rwfcr(:,:,:,:)
! radial charge density for core states
real(8), allocatable :: rhocr(:,:,:)
! spincore is .true. if the core is to be treated as spin-polarised
logical spincore
! number of core spin-channels
integer nspncr

!--------------------------!
!     energy variables     !
!--------------------------!
! core, valence and total occupied eigenvalue sum
real(8) evalsum
! electron kinetic energy
real(8) engykn
! core electron kinetic energy
real(8) engykncr
! nuclear-nuclear energy
real(8) engynn
! electron-nuclear energy
real(8) engyen
! Hartree energy
real(8) engyhar
! Coulomb energy (E_nn + E_en + E_H)
real(8) engycl
! electronic Coulomb potential energy
real(8) engyvcl
! Madelung term
real(8) engymad
! exchange-correlation potential energy
real(8) engyvxc
! exchange-correlation effective field energy
real(8) engybxc
! energy of external global magnetic field
real(8) engybext
! exchange energy
real(8) engyx
! correlation energy
real(8) engyc
! electronic entropy
real(8) entrpy
! entropic contribution to free energy
real(8) engyts
! total energy
real(8) engytot

!--------------------------------------------!
!     force, stress and strain variables     !
!--------------------------------------------!
! tforce is .true. if force should be calculated
logical tforce,tforce0
! Hellmann-Feynman force on each atom
real(8), allocatable :: forcehf(:,:)
! total force on each atom
real(8), allocatable :: forcetot(:,:)
! previous total force on each atom
real(8), allocatable :: forcetotp(:,:)
! maximum force magnitude over all atoms
real(8) forcemax
! maximum allowed force magnitude; if this force is reached for any atom then
! all forces are rescaled so that the maximum force magnitude is this value
real(8) maxforce
! tfav0 is .true. if the average force should be zero in order to prevent
! translation of the atomic basis
logical tfav0,tfav00
! atomic position optimisation type
!  0 : no optimisation
!  1 : unconstrained optimisation
integer atpopt
! maximum number of atomic position optimisation steps
integer maxatpstp
! default step size parameter for atomic position optimisation
real(8) tau0atp
! step size parameters for each atom
real(8), allocatable :: tauatp(:)
! number of strain tensors
integer nstrain
! current strain tensor
integer :: istrain=0
! strain tensors
real(8) strain(3,3,9)
! small displacement parameter multiplied by the strain tensor for computing the
! stress tensor; also used for calculating the piezoelectric tensor
real(8) deltast
! symmetry reduced stress tensor components
real(8) stress(9)
! previous stress tensor
real(8) stressp(9)
! stress tensor component magnitude maximum
real(8) stressmax
! reference lattice vectors for generating the G-vectors and derived quantities
real(8) avecref(3,3)
! tavref is .true. if avecref is non-zero
logical tavref
! lattice vector optimisation type
!  0 : no optimisation
!  1 : unconstrained optimisation
!  2 : iso-volumetric optimisation
integer latvopt
! maximum number of lattice vector optimisation steps
integer maxlatvstp
! default step size parameter for lattice vector optimisation
real(8) tau0latv
! step size for each stress tensor component acting on the lattice vectors
real(8) taulatv(9)

!--------------------------------------------------------!
!     self-consistent loop and convergence variables     !
!--------------------------------------------------------!
! maximum number of self-consistent loops
integer maxscl,maxscl0
! current self-consistent loop number
integer iscl
! tlast is .true. if the calculation is on the last self-consistent loop
logical tlast
! tstop is .true. if the STOP file exists
logical tstop
! trestart is .true. if the code should be completely restarted
logical trestart
! number of self-consistent loops after which STATE.OUT is written
integer nwrite
! Kohn-Sham potential convergence tolerance
real(8) epspot
! energy convergence tolerance
real(8) epsengy
! force convergence tolerance
real(8) epsforce
! stress tensor convergence tolerance
real(8) epsstress

!----------------------------------------------------------!
!     density of states, optics and response variables     !
!----------------------------------------------------------!
! number of energy intervals in the DOS/optics function plot
integer nwplot
! fine k-point grid size for integration of functions in the Brillouin zone
integer ngrkf
! smoothing level for DOS/optics function plot
integer nswplot
! energy interval for DOS/optics function plot
real(8) wplot(2)
! maximum angular momentum for the partial DOS plot
integer lmaxdos
! dosocc is .true. if the DOS is to be weighted by the occupancy
logical dosocc
! tpdos is .true. if the partial DOS should be calculated
logical tpdos
! dosmsum is .true. if the partial DOS is to be summed over m
logical dosmsum
! dosssum is .true. if the partial DOS is to be summed over spin
logical dosssum
! number of optical matrix components required
integer noptcomp
! required optical matrix components
integer optcomp(3,27)
! intraband is .true. if the intraband term is to be added to the optical matrix
logical intraband
! lmirep is .true. if the (l,m) band characters should correspond to the
! irreducible representations of the site symmetries
logical lmirep
! spin-quantisation axis in Cartesian coordinates used when plotting the
! spin-resolved DOS (z-axis by default)
real(8) sqados(3)
! q-vector in lattice and Cartesian coordinates for calculating the matrix
! elements ⟨i,k+q| exp(iq⋅r) |j,k⟩
real(8) vecql(3),vecqc(3)
! maximum initial-state energy allowed in ELNES transitions
real(8) emaxelnes
! structure factor energy window
real(8) wsfac(2)

!-------------------------------------!
!     1D/2D/3D plotting variables     !
!-------------------------------------!
! number of vertices in 1D plot
integer nvp1d
! total number of points in 1D plot
integer npp1d
! starting point for 1D plot
integer ip01d
! vertices in lattice coordinates for 1D plot
real(8), allocatable :: vvlp1d(:,:)
! distance to vertices in 1D plot
real(8), allocatable :: dvp1d(:)
! plot vectors in lattice coordinates for 1D plot
real(8), allocatable :: vplp1d(:,:)
! distance to points in 1D plot
real(8), allocatable :: dpp1d(:)
! corner vectors of 2D plot in lattice coordinates
real(8) vclp2d(3,0:2)
! grid sizes of 2D plot
integer np2d(2)
! corner vectors of 3D plot in lattice coordinates
real(8) vclp3d(3,0:3)
! grid sizes of 3D plot
integer np3d(3)

!-------------------------------------------------------------!
!     OEP, Hartree-Fock and Kohn-Sham inversion variables     !
!-------------------------------------------------------------!
! maximum number of core states over all species
integer ncrmax
! maximum number of OEP iterations
integer maxitoep
! OEP initial and subsequent step sizes
real(8) tau0oep,tauoep
! exchange potential and magnetic field
real(8), allocatable :: vxmt(:,:),vxir(:),bxmt(:,:,:),bxir(:,:)
! OEP residual functions
real(8), allocatable :: dvxmt(:,:),dvxir(:),dbxmt(:,:,:),dbxir(:,:)
! magnitude of the OEP residual
real(8) resoep
! hybrid is .true. if a hybrid functional is to be used
logical hybrid,hybrid0
! hybrid functional mixing coefficient
real(8) hybridc

!-------------------------------------------------------------!
!     response function and perturbation theory variables     !
!-------------------------------------------------------------!
! |G| cut-off for response functions
real(8) gmaxrf
! energy cut-off for response functions
real(8) emaxrf
! number of G-vectors for response functions
integer ngrf
! matrix bandwidth of response functions in the G-vector basis
integer mbwgrf
! number of response function frequencies
integer nwrf
! complex response function frequencies
complex(8), allocatable :: wrf(:)
! maximum number of spherical Bessel functions on the coarse radial mesh over
! all species
integer njcmax

!-------------------------------------------------!
!     Bethe-Salpeter equation (BSE) variables     !
!-------------------------------------------------!
! number of valence and conduction states for transitions
integer nvbse,ncbse
! default number of valence and conduction states
integer nvbse0,ncbse0
! maximum number of extra valence and conduction states
integer, parameter :: maxxbse=20
! number of extra valence and conduction states
integer nvxbse,ncxbse
! extra valence and conduction states
integer istxbse(maxxbse),jstxbse(maxxbse)
! total number of transitions
integer nvcbse
! size of blocks in BSE Hamiltonian matrix
integer nbbse
! size of BSE matrix (= 2*nbbse)
integer nmbse
! index from BSE valence states to second-variational states
integer, allocatable :: istbse(:,:)
! index from BSE conduction states to second-variational states
integer, allocatable :: jstbse(:,:)
! index from BSE valence-conduction pair and k-point to location in BSE matrix
integer, allocatable :: ijkbse(:,:,:)
! BSE Hamiltonian
complex(8), allocatable :: hmlbse(:,:)
! BSE Hamiltonian eigenvalues
real(8), allocatable :: evalbse(:)
! if bsefull is .true. then the full BSE Hamiltonian is calculated, otherwise
! only the Hermitian block
logical bsefull
! if hxbse/hdbse is .true. then the exchange/direct term is included in the BSE
! Hamiltonian
logical hxbse,hdbse

!--------------------------!
!     timing variables     !
!--------------------------!
! initialisation
real(8) timeinit
! Hamiltonian and overlap matrix set up
real(8) timemat
! first-variational calculation
real(8) timefv
! second-variational calculation
real(8) timesv
! charge density calculation
real(8) timerho
! potential calculation
real(8) timepot
! force calculation
real(8) timefor

!-----------------------------!
!     numerical constants     !
!-----------------------------!
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: twopi=6.2831853071795864769d0
real(8), parameter :: fourpi=12.566370614359172954d0
! spherical harmonic Y₀₀ = 1/√4π and its inverse
real(8), parameter :: y00=0.28209479177387814347d0
real(8), parameter :: y00i=3.54490770181103205460d0
! complex constants
complex(4), parameter :: czero=(0.e0,0.e0),cone=(1.e0,0.e0)
complex(4), parameter :: ci=(0.e0,1.e0),cmi=(0.e0,-1.e0)
complex(8), parameter :: zzero=(0.d0,0.d0),zone=(1.d0,0.d0)
complex(8), parameter :: zi=(0.d0,1.d0),zmi=(0.d0,-1.d0)
! Pauli spin matrices:
! σ_x = ⎛0  1⎞   σ_y = ⎛0 -i⎞   σ_z = ⎛1  0⎞
!       ⎝1  0⎠         ⎝i  0⎠         ⎝0 -1⎠
! Planck constant in SI units (exact, CODATA 2018)
real(8), parameter :: h_si=6.62607015d-34
! reduced Planck constant ℏ in SI units
real(8), parameter :: hbar_si=h_si/twopi
! speed of light in SI units (exact, CODATA 2018)
real(8), parameter :: sol_si=299792458d0
! speed of light in atomic units (1/α) (CODATA 2018)
real(8), parameter :: sol=137.035999084d0
! scaled speed of light
real(8) solsc
! Hartree in SI units (CODATA 2018)
real(8), parameter :: ha_si=4.3597447222071d-18
! Hartree in eV (CODATA 2018)
real(8), parameter :: ha_ev=27.211386245988d0
! Hartree in inverse meters
real(8), parameter :: ha_im=ha_si/(h_si*sol_si)
! Boltzmann constant in SI units (exact, CODATA 2018)
real(8), parameter :: kb_si=1.380649d-23
! Boltzmann constant in Hartree/kelvin
real(8), parameter :: kboltz=kb_si/ha_si
! electron charge in SI units (exact, CODATA 2018)
real(8), parameter :: e_si=1.602176634d-19
! Bohr radius in SI units (CODATA 2018)
real(8), parameter :: br_si=0.529177210903d-10
! Bohr radius in Angstroms
real(8), parameter :: br_ang=br_si*1.d10
! atomic unit of magnetic flux density in SI
real(8), parameter :: b_si=hbar_si/(e_si*br_si**2)
! atomic unit of electric field in SI
real(8), parameter :: ef_si=ha_si/(e_si*br_si)
! atomic unit of time in SI
real(8), parameter :: t_si=hbar_si/ha_si
! electron g-factor (CODATA 2018)
real(8), parameter :: gfacte=2.00231930436256d0
! electron mass in SI (CODATA 2018)
real(8), parameter :: em_si=9.1093837015d-31
! atomic mass unit in SI (CODATA 2018)
real(8), parameter :: amu_si=1.66053906660d-27
! atomic mass unit in electron masses
real(8), parameter :: amu=amu_si/em_si

!---------------------------------!
!     miscellaneous variables     !
!---------------------------------!
! code version
integer, parameter :: version(3)=[10,6,2]
! maximum number of tasks
integer, parameter :: maxtasks=40
! number of tasks
integer ntasks
! task index
integer itask
! task array
integer tasks(maxtasks)
! current task
integer task
! filename extension for files generated by gndstate
character(256) :: filext='.OUT'
! scratch space path
character(256) scrpath
! number of note lines
integer notelns
! notes to include in INFO.OUT
character(256), allocatable :: notes(:)

end module

