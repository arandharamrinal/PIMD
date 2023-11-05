program predict
use sgdmlvars 
use allvars
implicit none
integer::i,j,k,iostatus
real(8)::En =0.d0,Eai = 0.d0,ti,tf
character(len=50)::ffextfile
real(8),dimension(:),allocatable::X,force,force_ai
call loadmodel()
allocate(X(0:Ncarts-1),force(0:Ncarts-1),force_ai(0:Ncarts-1))
!Read file containing ffext files
open(unit=13,file='forcesfit.dat',status='unknown')
open(unit=11,file='filelist.txt',status='old')
do 
	read(11,*,iostat=iostatus)ffextfile
	if (iostatus.eq.0) then
		print*,ffextfile 
		open(unit=12,file="ffext_files/"//trim(adjustl(ffextfile)),status='old')
		do i = 1,4
			read(12,*)
		enddo
		read(12,*)X(0:Ncarts-1)
		do i = 1,3
    	    read(12,*)
    	enddo
		read(12,*)Eai
		read(12,*)
		read(12,*)force_ai(0:Ncarts-1)
		close(12)
		call cpu_time(ti)
		do i = 1,10000
			call predict_EF(X,En,force)
		enddo
		call cpu_time(tf)
		print*,'Took a total of ',(tf-ti)*0.0001,' seconds per force calculation'
		write(13,'(F16.8)',advance='no')En 
		do i = 0,Ncarts-1
			write(13,'(F16.8)',advance='no')force(i)
 		enddo
		write(13,*)
		print*,'done'
	else
		exit
	endif
enddo 


endprogram 
