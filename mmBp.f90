PROGRAM mmBp
 use global_data; use timers; use functions
 implicit none
 real(wp) :: eval, agmin, agmax, incr
 integer  :: i, j, d, s, mc, k
 logical  :: run_gpu, run_mpi, run_mpi_gpu
 character(len=3) :: st, ver
  
   !--Initialize MPI
   call MPI_init(ierror)
   call MPI_comm_size(MPI_comm_world, Nranks, ierror)
   call MPI_comm_rank(MPI_comm_world, rank, ierror)
   if (rank==0) then
      root = .true.
   else
      root = .false.
   end if
   if (root) write(*,'(A17,I3,A11)') "This code runs on", Nranks, "  MPI ranks"

   !--Initialize GPUs, link each to an MPI rank
   Ngpu = 0 
#if defined (_OPENACC) 
   devicetype = ACC_get_device_type()
   call ACC_init(devicetype) 
   Ngpu = ACC_get_num_devices(devicetype)
   gpu = rank ! each MPI rank links to a GPU 
   call ACC_set_device_num(gpu, devicetype)
   gpu = ACC_get_device_num(devicetype)
   if (root) then
      write(*,*)
      write(*,'(A13,I3)') 'Device type: ', devicetype
      write(*,'(A25,I3)') 'Number of GPUs on 1 node:', Ngpu      
      write(*,*)
   end if

   do j=0,Nranks-1
      call MPI_barrier(MPI_comm_world, ierror)
      if (rank==j) write(*,'(A4,I3,A17,I3)') "rank", rank, " is linked to GPU", gpu
   end do 
#endif

   !--If compiled without GPU offload or no GPU found
   call MPI_barrier(MPI_comm_world, ierror)
   if (Ngpu==0) then
      if (root) write(*,*) 'No GPUs detected'
   end if
   call MPI_barrier(MPI_comm_world, ierror)     

   !--Initialize --------------------------------------------------------------
   call allocate_data

   agmin = 1.
   agmax = 10.
   incr  = (agmax-agmin)/real(Nranks-1)
   prior_ag(1) = agmin
   do j=2,Nranks
      prior_ag(j) = agmin + real(j-1)*incr
   end do
   prior_a = prior_ag(rank+1)

   if (root) then
      write(*,*)       
      do j=0,Nranks-1
         write(*,'(A4,I3,A13,F6.2)') "rank", j, " has prior a ", prior_ag(j+1)
      end do
      write(*,*)
   end if

   prior_bv = prior_b
   Ssh  = 0.
   fna  = 0. ! average fni over MC run
   eval = 0.
   psha = 0. ! average psh over MC run

   !--Main --------------------------------------------------------------------
   call map_sb
   call map_ptr
   call update_sh_ttps
   call eval_gammas
   if (root) write(*,*) "MCMC iterations:"

   call timing_start

   do mc=1,MCiter
      if ((root).and.(mod(mc,1000)==0)) write(*,*) mc

      !--Data
      if (mc==1) then
         if (root) then 
            !--Read in data
            open(8,file=dirw//'log_returns.txt',action='read')
            do i=1,N
               read(8,*) ydat(i,:)
            end do               
            close(8)            
         end if
         call MPI_bcast(buffer = ydat, count = N*Dim, datatype = MPI_real, root = 0, &
                       & comm = MPI_comm_world, ierror = ierror)
         call eval_mcdf         
         !$ACC UPDATE DEVICE(udat) 

         call update_shir
         call eval_Bernstein
      end if

      if (mc==1) then
         !--Draw from the prior by setting allocation counts to zero
         !$ACC KERNELS DEFAULT(PRESENT)
         Nc  = 0.
         VNc = 0.
         Wc  = 0.
         !$ACC END KERNELS
         call draw_SH     
         call update_psh
         call update_sst
         call update_hst
         call update_allocation_counts

      else

         call draw_SH 
         call update_psh

         if ((Dim<=2).and.(mc>burnin)) then
            eval = eval + 1.
            call eval_post_weight
            !$ACC UPDATE HOST(fni)
            fna = fna*(eval-1.) ! convert to overall sum
            fna = fna + fni     ! add current mc output
            fna = fna/eval      ! convert back to average
         end if

         call update_sst
         call update_hst         
         call update_allocation_counts

      end if
      if (mc>burnin) psha = psha + psh

   end do

   call timing_finish

   psha = psha/real(MCiter-burnin)

   call copula_simulate
   call invert_marginals
   call eval_lprsim 
   call quicksort_nr(lprsim)

   VaR = lprsim(qid)
   wght = sum(fna)/real(N)

   !--root collects VaR values and weights from MPI ranks
   call MPI_gather(recvbuf = VaRg, recvcount = 1, recvtype = MPI_real, &
                  & sendbuf = VaR, sendcount = 1, sendtype = MPI_real, root = 0, &
                  & comm = MPI_comm_world, ierror = ierror)

   call MPI_gather(recvbuf = wghtg, recvcount = 1, recvtype = MPI_real, &
                  & sendbuf = wght, sendcount = 1, sendtype = MPI_real, root = 0, &
                  & comm = MPI_comm_world, ierror = ierror)

   if (root) then

      !--multiply weights by model prior
      do j=1,Nranks
         prm(j) = 0.8**j
      end do
      prm = prm/sum(prm)
      wghtg = wghtg*prm

      !--normalize weights (the marginals will have the same density evaluated 
      !--at the data, so we only need the copula density kernel for weights)
      wghtg = wghtg/sum(wghtg)   

      write(*,*) 'VaR'
      write(*,'(8ES12.3)') VaRg
      write(*,*) 'BMA weights'
      write(*,'(8ES12.3)') wghtg
      write(*,*) 'BMA VaR = ', sum(VaRg*wghtg)

      open(18,file=dirw//'BMA_weights_a.out',action='write',status='replace')      
      do j=1,Nranks
         write(18,'(2F15.8)') prior_ag(j), wghtg(j)
      end do
      close(18)

   end if

   !--Cleanup
10 if (root) then
   end if   
   call deallocate_data
#if defined (_OPENACC)      
   call ACC_shutdown(devicetype)
#endif
   call MPI_barrier(MPI_comm_world, ierror)
   call MPI_finalize(ierror)
 
END PROGRAM mmBp
