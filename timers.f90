MODULE timers
 use global_data
 implicit none
 real(sp) :: start_time, finish_time, time_cpu ! for function cpu_time
 real(sp) :: time_wallclock                    ! for function system_clock
 integer  :: start_time2, finish_time2, crate  ! for function system_clock
 real(dp) :: MPI_time_wallclock                ! for function MPI_Wtime
 real(dp) :: start_time3, finish_time3         ! for function MPI_Wtime
 CONTAINS

 !--------------------------------------------------------------------
 !                    SUBROUTINE timing_start
 !--------------------------------------------------------------------  
 SUBROUTINE timing_start

   call cpu_time(start_time)
   call system_clock(count_rate=crate)
   call system_clock(count=start_time2)
   start_time3 = MPI_Wtime() 

 END SUBROUTINE timing_start

 !--------------------------------------------------------------------
 !                    SUBROUTINE timing_finish
 !--------------------------------------------------------------------  
 SUBROUTINE timing_finish

   call cpu_time(finish_time)
   call system_clock(count=finish_time2)
   finish_time3 = MPI_Wtime()
   time_cpu = finish_time - start_time ! in seconds
   time_wallclock = real(finish_time2 - start_time2)/real(crate) ! in seconds
   MPI_time_wallclock = finish_time3 - start_time3

   if (root) then
      !write(*,*) "rank", rank
      write(*,'(A20,F12.6)') "time_cpu =          ", time_cpu
      write(*,'(A20,F12.6)') "time_wallclock =    ", time_wallclock
      write(*,'(A20,F12.6)') "MPI_time_wallclock =", MPI_time_wallclock
      write(*,*)
   end if

 END SUBROUTINE timing_finish
 
END MODULE timers
