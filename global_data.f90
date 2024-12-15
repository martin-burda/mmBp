MODULE global_data
 use MPI_f08
#if defined (_OPENACC) 
 use OpenACC
#endif
 implicit none

 !--MPI info
 integer :: Nranks, rank ! number of ranks, rank number
 integer :: ierror
 logical :: root 

 !--GPU info
 integer :: Ngpu, gpu ! number of GPUS, GPU number
#if defined (_OPENACC) 
 integer(acc_device_kind) :: devicetype
#endif

 !--Parameters
 integer, parameter :: dp = selected_real_kind(p=12,r=60) 
 integer, parameter :: sp = selected_real_kind(p=6,r=37)  
 integer, parameter :: wp = sp 
 character(len=2), parameter :: dirw = './' 
  
 !--Application parameters
 integer, parameter   :: N = 1191     ! sample size
 integer, parameter   :: Smax = 5     ! max tree depth 
 integer, parameter   :: Dim = 2      ! dimensions of tree
 integer, parameter   :: ndn = 2**Dim ! number of daughter nodes
 integer, parameter   :: eps = tiny(1._wp)

 !--level nodes
 integer, parameter :: lv(0:10) = (/ 1, 2, 2**2, 2**3, 2**4, 2**5, &
                                   & 2**6, 2**7, 2**8, 2**9, 2**10 /)

 !--end of level nodes                                   
 !--needs to be declared with the 'parameter' attribute, hence no loop                                     
 integer, parameter :: elv(0:10) = (/ 1, &
  & 1 + 2, &
  & 1 + 2 + 2**2, &
  & 1 + 2 + 2**2 + 2**3, &
  & 1 + 2 + 2**2 + 2**3 + 2**4, &
  & 1 + 2 + 2**2 + 2**3 + 2**4 + 2**5, &
  & 1 + 2 + 2**2 + 2**3 + 2**4 + 2**5 + 2**6, &
  & 1 + 2 + 2**2 + 2**3 + 2**4 + 2**5 + 2**6 + 2**7, &
  & 1 + 2 + 2**2 + 2**3 + 2**4 + 2**5 + 2**6 + 2**7 + 2**8, &
  & 1 + 2 + 2**2 + 2**3 + 2**4 + 2**5 + 2**6 + 2**7 + 2**8 + 2**9, &
  & 1 + 2 + 2**2 + 2**3 + 2**4 + 2**5 + 2**6 + 2**7 + 2**8 + 2**9 + 2**10 /)

 !--total nodes, up to s
 !--needs to be declared with the 'parameter' attribute, hence no loop  
 integer, parameter :: tn(0:10) = (/ 1, &
  & 1 + 2**(Dim*1), &
  & 1 + 2**(Dim*1) + 2**(Dim*2), &
  & 1 + 2**(Dim*1) + 2**(Dim*2) + 2**(Dim*3), &
  & 1 + 2**(Dim*1) + 2**(Dim*2) + 2**(Dim*3) + 2**(Dim*4),  &
  & 1 + 2**(Dim*1) + 2**(Dim*2) + 2**(Dim*3) + 2**(Dim*4) + &
  & 2**(Dim*5), &
  & 1 + 2**(Dim*1) + 2**(Dim*2) + 2**(Dim*3) + 2**(Dim*4) + &
  & 2**(Dim*5) + 2**(Dim*6), &
  & 1 + 2**(Dim*1) + 2**(Dim*2) + 2**(Dim*3) + 2**(Dim*4) + &
  & 2**(Dim*5) + 2**(Dim*6) + 2**(Dim*7), &
  & 1 + 2**(Dim*1) + 2**(Dim*2) + 2**(Dim*3) + 2**(Dim*4) + &
  & 2**(Dim*5) + 2**(Dim*6) + 2**(Dim*7) + 2**(Dim*8), &
  & 1 + 2**(Dim*1) + 2**(Dim*2) + 2**(Dim*3) + 2**(Dim*4) + &
  & 2**(Dim*5) + 2**(Dim*6) + 2**(Dim*7) + 2**(Dim*8) + 2**(Dim*9), &
  & 1 + 2**(Dim*1) + 2**(Dim*2) + 2**(Dim*3) + 2**(Dim*4) + &
  & 2**(Dim*5) + 2**(Dim*6) + 2**(Dim*7) + 2**(Dim*8) + 2**(Dim*9) + 2**(Dim*10) /)

 !--Application arrays  
 integer, allocatable  :: sst(:), sh(:), ptr(:), ptr_path(:,:), sb(:,:), hst(:), shi(:,:,:)  
 real(wp), allocatable :: prior_ag(:), udat(:,:), ydat(:,:), ttps(:), shr(:), shir(:,:,:)
 real(wp), allocatable :: lnGam1(:), lnGam2(:), lnGam3(:), ln_udat(:,:), ln_1_udat(:,:)
 real(wp), allocatable :: Bdn(:,:,:), Bp(:,:), Ssh(:), Hsh(:)
 real(wp), allocatable :: Nc(:), VNc(:), Wc(:), params(:), prior_bv(:), draw(:,:)
 real(wp), allocatable :: psh(:), psha(:), fn(:,:), prs(:,:), un(:), prh(:,:)
 
 !--MC part
 integer, parameter    :: MCiter = 21000, burnin = 1000
 real(wp), allocatable :: fnj(:,:), fni(:), fna(:)

 !--VaR
 real(wp) :: mn(Dim), sdev(Dim) 
 integer, parameter :: NSim = 10000
 real(wp), parameter    :: NSim_wp = real(NSim)
 real(wp), parameter    :: prior_b = 1. 
 real(wp) :: prior_a, Csim(NSim,dim), Msim(NSim,dim)
 real(wp) :: lprsim(NSim), VaR, wght
 real(wp), allocatable :: VaRg(:), wghtg(:), prm(:)

 !--Index of the VaR quantile
 real(wp), parameter :: alp = 0.05_wp
 real(wp), parameter :: qid_wp = NSim_wp*alp 
 integer, parameter  :: qid = nint(qid_wp) ! nearest integer 

 CONTAINS
 !--------------------------------------------------------------------
 !                    SUBROUTINE allocate_data
 !--------------------------------------------------------------------  
 SUBROUTINE allocate_data

  !--Allocate arrays on host 
  allocate(udat(N,Dim), ydat(N,Dim), ttps(elv(Smax)), sst(N))
  allocate(shir(0:Smax,N,Dim), shi(0:Smax,N,Dim), shr(elv(Smax)), sh(elv(Smax)))
  allocate(lnGam1(elv(Smax)), lnGam2(elv(Smax)), lnGam3(elv(Smax)))
  allocate(ln_udat(N,Dim), ln_1_udat(N,Dim), Bdn(elv(Smax),N,Dim), Bp(tn(Smax),N))
  allocate(sb(tn(Smax),Dim), Ssh(tn(Smax)), Hsh(tn(Smax)))
  allocate(Nc(tn(Smax)), VNc(tn(Smax)), Wc(tn(Smax)), hst(N), params(ndn))
  allocate(prior_bv(ndn), ptr(tn(Smax)), ptr_path(tn(Smax),Smax), draw(tn(Smax-1),ndn))
  allocate(psh(tn(Smax)), psha(tn(Smax)), fn(tn(Smax),N), prs(N,0:Smax), un(N))
  allocate(prh(N,2**(Dim*Smax)), fnj(tn(Smax),N), fni(N), fna(N))
  allocate(prior_ag(Nranks), VaRg(Nranks), wghtg(Nranks), prm(Nranks))

  !--Allocate arrays on device: each MPI rank creates its own instance of the variables
  !$ACC ENTER DATA CREATE(udat, ydat, ttps, sst, shr, sh, shir, shi)
  !$ACC ENTER DATA CREATE(lnGam1, lnGam2, lnGam3, ln_udat, ln_1_udat, Bdn, Bp, sb)
  !$ACC ENTER DATA CREATE(Ssh, Hsh, Nc, VNc, Wc, ptr, ptr_path, draw, psh, psha)
  !$ACC ENTER DATA CREATE(fn, prs, un, prh, hst, fnj, fni, fna)

 END SUBROUTINE allocate_data

 !--------------------------------------------------------------------
 !                    SUBROUTINE deallocate_data
 !--------------------------------------------------------------------  
 SUBROUTINE deallocate_data

  !--Deallocate arrays on device
  !$ACC EXIT DATA DELETE(udat, ydat, ttps, sst, shr, sh, shir, shi)
  !$ACC EXIT DATA DELETE(lnGam1, lnGam2, lnGam3, ln_udat, ln_1_udat, Bdn, Bp, sb)
  !$ACC EXIT DATA DELETE(Ssh, Hsh, Nc, VNc, Wc, ptr, ptr_path, draw, psh, psha)
  !$ACC EXIT DATA DELETE(fn, prs, un, prh, hst, fnj, fni, fna)

  !--Deallocate arrays on device   
  deallocate(udat, ydat, ttps, sst, shr, sh, shir, shi, lnGam1, lnGam2, lnGam3, ln_udat, ln_1_udat, Bdn, Bp, sb)
  deallocate(Ssh, Hsh, Nc, VNc, Wc, params, prior_bv, ptr, ptr_path, draw, psh, psha)
  deallocate(fn, prs, un, prh, hst, prior_ag, VaRg, wghtg)

 END SUBROUTINE deallocate_data

END MODULE global_data
