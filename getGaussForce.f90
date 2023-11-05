subroutine getForce(ibead,Xcarts,GaussForce,GaussE,writeforce)
use allvars, only:Ncarts,iounitGaussForce
use Gaussparam
use omp_lib
implicit none 
integer,intent(in)::ibead,writeforce
real(8),intent(in)::Xcarts(0:Ncarts-1)
real(8),intent(out)::GaussE, GaussForce(0:Ncarts-1)
real(8):: GaussDip(0:2)
integer::IOunit,readerr,IOstatus,i,j,k,nchar 
character(len=200)::tmpChar,Fullline 
character(len=20)::ibeadChar,dumChar,gaussOutDir
character(len=100)::xyzfile,gaussInpfile,gausschkfile,gaussfchkfile,beadTitle
211 format(49X,F24.16)

!Initialize parameters 
GaussE = 0.d0
GaussDip=0.d0
GaussForce=0.d0

!Directory to which gaussian calculations are run
GaussOutDir = 'g16'
write(ibeadChar,*)ibead

!xyzfile = 'bead'//trim((adjustl(ibeadChar)))//'.xyz'
!gaussian input file name 
gaussInpfile = trim(adjustl(GaussOutDir))//'/bead'//trim((adjustl(ibeadChar)))//'.com'

!Title of gaussian input file 
beadTitle =  trim((adjustl(gaussTitle)))//trim((adjustl(ibeadChar)))

!Gaussian checkpoint file and fchk file names
gausschkfile = 'bead'//trim((adjustl(ibeadChar)))//'.chk' 
gaussfchkfile = 'bead'//trim((adjustl(ibeadChar)))//'.fchk' 

!Write the gaussian input file to gaussian out directory
call writeGaussinp(iounitGaussb+ibead,gaussInpfile,gausschkfile,beadTitle,Xcarts,gaussOutDir)!,atomname)

!Run gaussian and direct the shell output to a log file
call system('g16 < '//trim(adjustl(gaussInpfile))//' > '//trim(adjustl(GaussOutDir))//'/bead'//trim((adjustl(ibeadChar)))//'.log')

!If gaussian calculations are finished extract the forces from the logfile and
!write them to forces${ibead}.out file
tmpCHar = 'grep "Forces (Hartrees/Bohr)" -A 7 '//trim(adjustl(GaussOutDir))//'/bead'//trim((adjustl(ibeadChar)))//'.log  | tail -n 5  | awk '//"'{printf "//'"  %s   %s   %s"'//",$3,$4,$5}'"//'> '//trim(adjustl(GaussOutDir))//'/forces'//trim((adjustl(ibeadChar)))//'.out'
call system(trim(adjustl(tmpCHar)))

!Read forces from forces${ibead}.out file 
IOunit = iounitGaussa+ibead
!log filereading reading 
open(unit=IOunit,file=trim(adjustl(trim(adjustl(GaussOutDir))//'/forces'//trim((adjustl(ibeadChar)))))//'.out',status='old',iostat=IOstatus)
if (IOstatus==0) then
	read(IOunit,*)(GaussForce(j),j=0,Ncarts-1)
else
	write(*,'(A)')'error reading forces.out. sleeping for 1 sec.'
	call system('sleep 1')
	close(IOunit)
	!Try reading again.
	open(unit=IOunit,file=trim(adjustl(trim(adjustl(GaussOutDir))//'/forces'//trim((adjustl(ibeadChar)))))//'.out',status='old',iostat=IOstatus)
	read(IOunit,*)(GaussForce(j),j=0,Ncarts-1)
	if (IOstatus.ne.0) then
		write(*,'(A)')'Error in reading'//trim(adjustl(trim(adjustl(GaussOutDir))//'/forces'//trim((adjustl(ibeadChar)))))//'.out'
	endif
endif
close(IOunit)
!Convert forces to gradient 
GaussForce = -1.d0*GaussForce
!print*,ibead,GaussForce;read(*,*)
!
!!Formchk reading
!call system('formchk -3 '//trim(adjustl(gausschkfile))//'   '//trim(adjustl(gaussfchkfile)) //' >dum'//trim((adjustl(ibeadChar))))
!open(unit=IOunit,file=trim(adjustl(gaussfchkfile)),status='old')
!readerr=0
!do
!	read(IOunit,'(A)',IOSTAT=IOstatus)FullLine
!	if ( IOstatus/=0 ) then
!		close(IOunit)
!		call system('sleep 1s')
!		call system('formchk -3 '//trim(adjustl(gausschkfile))//'   '//trim(adjustl(gaussfchkfile))//' >dum'//trim((adjustl(ibeadChar)))) 
!		call system('sleep 1s')
!		open(unit=IOunit,file=trim(adjustl(gaussfchkfile)),status='old')
!		read(IOunit,'(A)',IOSTAT=IOstatus)FullLine
!		if (IOstatus/=0) then
!			print*,'Error in fchk file reading'
!			stop
!		endif	
!	endif
!	if (FullLine(1:20)=='Total Energy        ') then
!		readerr = readerr + 1
!		read(FullLine,211)GaussE
!	elseif (FullLine(1:20)=='Cartesian Gradient  ') then 
!		readerr = readerr + 1
!		do i = 0,Ncarts/5-1
!			read(IOunit,*)(GaussForce(5*i+j),j=0,4)
!		enddo
!	elseif (FullLine(1:20)=='Dipole Moment       ') then
!		readerr = readerr + 1
!		read(IOunit,*)(GaussDip(j),j=0,2)
!	endif 
!	if (readerr==3) then
!		close(IOunit)
!		exit
!	endif
!enddo
!print*,'gaussF=',GaussForce
if (writeforce==1) then
	write(iounitGaussForce)Xcarts,GaussE,GaussForce,GaussDip
endif
!print*,'gaussE=',GaussE,GaussDip;read(*,*)
end subroutine 


subroutine writeGaussinp(iounit,gaussInpfile,gausschkfile,beadTitle,Xcarts,gaussOutDir)!,atomname)
use Gaussparam
use constants, only:BohrToAng
use allvars,only:Ncarts,Natoms,atomname
implicit none 
real(8),intent(in):: Xcarts(0:Ncarts-1)
integer,intent(in):: iounit
!character(len=5),intent(in)::atomname(0:Natoms-1)
character(len=20),intent(in)::gaussOutDir
character(len=100),intent(in)::gaussInpfile,gausschkfile,beadTitle 
integer::j,k
real(8)::Xtmp(0:Ncarts-1)
Xtmp =  Xcarts * BohrToAng
open(unit=iounit,file=trim(adjustl(gaussInpfile)),status='unknown')
write(iounit,'(A)')'%nosave'
write(iounit,'(A)')trim(adjustl(gaussNproc))
write(iounit,'(A)')trim((gaussMem))
write(iounit,'(A)')'%chk='//trim(adjustl(GaussOutDir))//'/'//trim(adjustl(gausschkfile))
write(iounit,'(A)')trim(adjustl(gaussMethod))
write(iounit,*)
write(iounit,'(A)')trim(adjustl(gaussTitle))
write(iounit,*)
write(iounit,'(A)')trim(adjustl(gaussChargeSpin))
do j = 0,Natoms-1
	write(iounit,'(A,3F16.8)')adjustl(atomname(j)),(Xtmp(3*j+k),k=0,2)
enddo
write(iounit,*)
write(iounit,*)
close(iounit)
end subroutine
