MODULE functions
 use global_data; use timers; use random; use normals; use sorts
 implicit none
 CONTAINS

 !--------------------------------------------------------------------
 !                    SUBROUTINE eval_mcdf
 !-------------------------------------------------------------------- 
 SUBROUTINE eval_mcdf
 !--Evaluates marginal cdfs
  integer :: d, i
 
   do concurrent(d=1:Dim)
      mn(d) = sum(ydat(:,d))/real(N)
      sdev(d) = sqrt(sum((ydat(:,d)-mn(d))**2)/real(N))      
   end do

   do concurrent(i=1:N,d=1:Dim)
      udat(i,d) = normal_01_cdf((ydat(i,d)-mn(d))/sdev(d))
   end do
 
 END SUBROUTINE eval_mcdf

 !--------------------------------------------------------------------
 !                    SUBROUTINE invert_marginals
 !--------------------------------------------------------------------   
 SUBROUTINE invert_marginals
  integer :: j, d

   do concurrent(j=1:Nsim, d=1:dim)
      Msim(j,d) = normal_cdf_inv(Csim(j,d), mn(d), sdev(d))
   end do

  END SUBROUTINE invert_marginals

 !--------------------------------------------------------------------
 !                    SUBROUTINE eval_lprsim
 !--------------------------------------------------------------------   
 SUBROUTINE eval_lprsim
  !--log portfolio returns simulated
  integer :: j
    
   do j=1,Nsim 
      lprsim(j) = 0.5_wp*Msim(j,1) + 0.5_wp*Msim(j,2) 
   end do
 
   END SUBROUTINE eval_lprsim   

 !--------------------------------------------------------------------
 !                    SUBROUTINE map_sb
 !--------------------------------------------------------------------   
 SUBROUTINE map_sb
  !--sb stores all combinations of dimensions. It is an index vector
  !--for summing Bdn into Bp in eval_Bernstein.
  integer :: j, k, s, t

   sb(1,:) = 1   
   do j=1,Dim
      do s=1,Smax
         t = 1
         do k=1,tn(s)-tn(s-1) 
            sb(tn(s-1)+k,j) = elv(s-1) + t
            if (mod(k,2**((j-1)*s))==0) t = t + 1
            if (t>2**s) t = 1 ! reset counter
         end do
      end do
   end do
   !$ACC UPDATE DEVICE(sb)

 END SUBROUTINE map_sb

 !--------------------------------------------------------------------
 !                    RECURSIVE FUNCTION point_ptr
 !--------------------------------------------------------------------
 PURE RECURSIVE FUNCTION point_ptr(k, levels) result(r)
  integer, intent(in) :: k, levels
  integer :: r, j

   jloop: do j=1,levels
      r = ptr(k)        
      if (levels-1==0) EXIT jloop
      r = point_ptr(r, levels-1)
   end do jloop

 END FUNCTION point_ptr

 !--------------------------------------------------------------------
 !                    SUBROUTINE map_ptr
 !--------------------------------------------------------------------   
 SUBROUTINE map_ptr
  !--ptr stores parent row path for each node of level Smax
  integer :: k, s, j

   ptr(1) = 0
   ptr(2:1+2**Dim) = 1
   do s=1,Smax-1
      do k=1,tn(s)-tn(s-1)
         ptr(tn(s) + (k-1)*2**Dim+1 : tn(s) + k*2**Dim) = tn(s-1) +  k 
      end do
   end do
   !$ACC UPDATE DEVICE(ptr)

   ptr_path = 0
   do s=1,Smax
      do k=tn(s-1)+1,tn(s)
         do j=1,s      
            ptr_path(k,j) = point_ptr(k,j)
         end do
      end do
   end do
   !$ACC UPDATE DEVICE(ptr_path)    

 END SUBROUTINE map_ptr

 !--------------------------------------------------------------------
 !                    SUBROUTINE update_sh_ttps
 !--------------------------------------------------------------------   
 SUBROUTINE update_sh_ttps
  integer :: s, d, k

   sh(1) = 1
   do s=1,Smax
      do k=1,2**s
         sh(elv(s-1)+k) = k
      end do
   end do
   shr = real(sh)

   ttps(1) = 2.
   do s=1,Smax
      ttps(elv(s-1)+1:elv(s)) = real(2**s) 
   end do
   !$ACC UPDATE DEVICE(sh, shr, ttps)
 
 END SUBROUTINE update_sh_ttps

 !--------------------------------------------------------------------
 !                    SUBROUTINE eval_gammas
 !--------------------------------------------------------------------   
 SUBROUTINE eval_gammas
  integer :: i, s, d, k

   !$ACC KERNELS DEFAULT(PRESENT)
   lnGam1 = log_gamma(ttps + 1.)
   lnGam2 = log_gamma(shr)
   lnGam3 = log_gamma(ttps - shr + 1.)
   !$ACC END KERNELS

 END SUBROUTINE eval_gammas

 !--------------------------------------------------------------------
 !                    SUBROUTINE update_shir
 !--------------------------------------------------------------------   
 SUBROUTINE update_shir
  !--shir stores the h index at Smax, along each dimension
  !--shir is cumulative along each dimension, 1,...,elv(Smax)
  integer :: i, s, d, k

   !$ACC KERNELS DEFAULT(PRESENT)
   shir(0,:,:) = 1.
   do s=1,Smax
      do d=1,Dim
         do i=1,N         
            shir(s,i,d) = udat(i,d)*2**s     
         end do
      end do
   end do
   shi = ceiling(shir)

   do s=1,Smax
      shi(s,:,:) = shi(s,:,:) + elv(s-1)
   end do
   shir = real(shi)
   !$ACC END KERNELS

 END SUBROUTINE update_shir

 !--------------------------------------------------------------------
 !                    SUBROUTINE eval_Bernstein
 !--------------------------------------------------------------------   
 SUBROUTINE eval_Bernstein
  integer :: i, s, d, k, j, t

   !$ACC KERNELS DEFAULT(PRESENT)
   ln_udat = log(udat)
   where((ln_udat+1.0)==ln_udat) ln_udat = eps ! bind -inf 
   ln_1_udat = log(1. - udat)
   where((ln_1_udat+1.0)==ln_1_udat) ln_1_udat = eps ! bind -inf   

   Bdn(1,:,:) = 0. !(=log(1.))
   do i=1,N
      do d=1,Dim
         do j=2,elv(Smax)
               Bdn(j,i,d) = lnGam1(j) - lnGam2(j) - lnGam3(j) + &
                          & (shr(j) - 1.)*ln_udat(i,d) + &
                          & (ttps(j) - shr(j))*ln_1_udat(i,d)
         end do
      end do
   end do

   Bp = 0.   
   do k=1,tn(Smax)
      do i=1,N         
         do j=1,Dim         
            Bp(k,i) = Bp(k,i) + Bdn(sb(k,j),i,j) !Bdn is in logs so we sum to get the log of the product   
         end do
      end do                     
   end do
   !$ACC END KERNELS

 END SUBROUTINE eval_Bernstein

 !--------------------------------------------------------------------
 !                    SUBROUTINE draw_SH
 !--------------------------------------------------------------------   
 SUBROUTINE draw_SH
  integer :: i, s, d, j, k

   !--Draw from the posterior
   !$ACC UPDATE HOST(Nc, VNc)
   Ssh = 0.
   do k=1,tn(Smax-1) 
      Ssh(k) = random_beta(1. + Nc(k), prior_a + VNc(k), .true.) 
   end do
   Ssh(tn(Smax-1)+1:tn(Smax)) = 1.
   !$ACC UPDATE DEVICE(Ssh)

   !--Draw from the posterior, stored at the parent nodes 
   !$ACC UPDATE HOST(Wc)
   draw = 0.
   do k=1,tn(Smax-1)
      params = prior_bv + Wc(k)
      draw(k,:) = random_dirichlet(ndn, params)
   end do
   !$ACC UPDATE DEVICE(draw)   

   !--Copy out the draws to the daughters
   !$ACC KERNELS DEFAULT(PRESENT)
   Hsh(1) = 1.
   !$ACC LOOP INDEPENDENT
   do k=1,tn(Smax-1)
      do j=1,ndn
         Hsh(1+(k-1)*ndn+j) = draw(k, j)
      end do
   end do
   !$ACC END KERNELS

 END SUBROUTINE draw_SH

 !--------------------------------------------------------------------
 !                    SUBROUTINE update_psh
 !--------------------------------------------------------------------   
 SUBROUTINE update_psh
  integer :: k, s, r
   
   !$ACC UPDATE HOST(Ssh, Hsh)
   psh = Hsh*Ssh
   
   do s=1,Smax
      do k=tn(s-1)+1,tn(s)
         do r=1,s
            psh(k) = psh(k)*Hsh(ptr_path(k,r))*(1. - Ssh(ptr_path(k,r)))
         end do
      end do
   end do
   !$ACC UPDATE DEVICE(psh)

 END SUBROUTINE update_psh

 !--------------------------------------------------------------------
 !                    SUBROUTINE update_sst
 !--------------------------------------------------------------------   
 SUBROUTINE update_sst
  integer :: s, i

   !$ACC KERNELS DEFAULT(PRESENT)
   do i=1,N
      fn(:,i) = Bp(:,i) + log(psh) ! Bp is in logs
   end do
   fn = exp(fn) ! Pr(s) unnormalized

   prs(:,0) = fn(1,:)
   do i=1,N
      do s=1,Smax
         prs(i,s) = sum(fn(tn(s-1)+1:tn(s),i))
      end do
   end do
   !$ACC END KERNELS
   !$ACC UPDATE HOST(prs)
   do i=1,N
      prs(i,:) = prs(i,:)/sum(prs(i,:)) ! Pr(s) normalized
   end do

   do i=1,N
      do s=Smax,1,-1 ! go backwards for accumulation
         prs(i,s) = sum(prs(i,0:s)) ! cumulative Pr(s)
      end do
   end do
   !$ACC UPDATE DEVICE(prs)

   call random_number(un)
   !$ACC UPDATE DEVICE(un)

   !$ACC KERNELS DEFAULT(PRESENT)
   do i=1,N
      prs(i,:) = prs(i,:) - un(i)
   end do
   where (prs<0.) prs = 1.
   !$ACC END KERNELS      

   !$ACC UPDATE HOST(prs)
   do i=1,N
      sst(i:i) = minloc(prs(i,:)) - 1 ! include level 0, minloc picks position in 0:Smax 
   end do
   !$ACC UPDATE DEVICE(sst) 

 END SUBROUTINE update_sst

 !--------------------------------------------------------------------
 !                    SUBROUTINE update_hst
 !--------------------------------------------------------------------   
 SUBROUTINE update_hst
  !--hst stores the multivariate node index, 1,...,tn(Smax), 
  !--into which i is allocated
  integer :: s, i, k

   !$ACC UPDATE HOST(fn)
   prh = 0.
   do i=1,N
      if (sst(i)>0) then
         prh(i,1:tn(sst(i))-tn(sst(i)-1)) = fn(tn(sst(i)-1)+1:tn(sst(i)),i)/ &
                                          & sum(fn(tn(sst(i)-1)+1:tn(sst(i)),i)) ! normalized
      end if
   end do

   do i=1,N
      if (sst(i)>0) then
         do k=2,tn(sst(i))-tn(sst(i)-1)
            prh(i,k) = prh(i,k-1) + prh(i,k) ! cumulative Pr(s)
         end do
      end if
   end do

   call random_number(un)

   do i=1,N
      if (sst(i)>0) prh(i,1:tn(sst(i))-tn(sst(i)-1)) = prh(i,1:tn(sst(i))-tn(sst(i)-1)) - un(i)
   end do
   where (prh<0.) prh = 1.
   where(sst<0) sst = 0

   do i=1,N
      if (sst(i)==0) then
         hst(i) = 1
      else
         hst(i:i) = minloc(prh(i,1:tn(sst(i))-tn(sst(i)-1))) 
      end if
   end do

   do i=1,N
      if (sst(i)>0) hst(i) = hst(i) + tn(sst(i)-1)
   end do
   !$ACC UPDATE DEVICE(hst) 

 END SUBROUTINE update_hst

 !--------------------------------------------------------------------
 !                    SUBROUTINE update_allocation_counts
 !--------------------------------------------------------------------   
 SUBROUTINE update_allocation_counts
  integer :: i, k, s

   !--Count of subjects that stopped at a given node  
   Nc = 0. 
   do i=1,N
      Nc(hst(i)) = Nc(hst(i)) + 1
   end do
   !$ACC UPDATE DEVICE(Nc)

   !--Count of subjects that reached and continued beyond a given node   
   VNc = 0. 
   do i=1,N
      do s=1,sst(i)
         VNc(ptr_path(hst(i),s)) = VNc(ptr_path(hst(i),s)) + 1
      end do
   end do
   !$ACC UPDATE DEVICE(VNc)

   !--Count of subjects that reached a given node and either stopped 
   !--or continued further
   Wc = Nc + VNc 
   !$ACC UPDATE DEVICE(Wc)

 END SUBROUTINE update_allocation_counts

 !--------------------------------------------------------------------
 !                    SUBROUTINE copula_simulate
 !--------------------------------------------------------------------   
 SUBROUTINE copula_simulate
  real(wp) :: uni, a, b
  real(wp), allocatable :: sPrMn(:) ! cumul multinom probs
  integer, allocatable  :: NCd(:)   ! # of copula draws for each node
  integer  :: Tx, i, d1, d2, nid, q, s, j
  real(wp) :: d1r, d2r
 
   Tx = tn(Smax)
   allocate(sPrMn(Tx), NCd(Tx))

   !--Draw from the multinomial w/ probs psha
   do i=1,Tx
      sPrMn(i) = sum(psha(1:i)) ! cumulative sum
   end do   
 
   NCd = 0
   do j=1,NSim
      call random_number(uni)
      iloop: do i=1,Tx
         if (uni<sPrMn(i)) then
            NCd(i) = NCd(i) + 1 ! # of copula draws for each node
            EXIT iloop
         end if
      end do iloop
   end do
 
   !--Simulate from copula
   if (Dim/=2) then
      write(*,*) 'only Dim=2 supported for this routine'
      STOP
   end if

   q = 0
   do s=1,Smax
      do d1=1,lv(s)
         do d2=1,lv(s)
            nid = tn(s-1) + (d1-1)*lv(s) + d2 ! node index in NCd
            if (NCd(nid)>0) then ! draw from node
               do j=1,NCd(nid)

                  q = q + 1
                  d1r = real(d1)
                  d2r = real(d2)
                  Csim(q,1) = random_beta(d1r, ttps(elv(s))-d1r+1., .true.)  
                  Csim(q,2) = random_beta(d2r, ttps(elv(s))-d2r+1., .true.)  

               end do
            end if
         end do
      end do
   end do

   deallocate(sPrMn, NCd)
 
 END SUBROUTINE copula_simulate

 !--------------------------------------------------------------------
 !                    SUBROUTINE eval_post_weight
 !--------------------------------------------------------------------   
 SUBROUTINE eval_post_weight
  !--Array Bp is evaluated in subroutine eval_Bernstein called previously
  integer :: k, i
  
   !$ACC KERNELS DEFAULT(PRESENT)
   do k=1,tn(Smax)
      do i=1,N
         fnj(k,i) = Bp(k,i) + log(psh(k)) ! Bp is in logs, sum for log of product
      end do
   end do
   fnj = exp(fnj) 

   do i=1,N
      fni(i) = sum(fnj(:,i))
   end do

   !$ACC END KERNELS  
  
 END SUBROUTINE eval_post_weight  

END MODULE functions

