!****************************************************************************************************#  
subroutine writelogfile()
use constants 
use allvars
implicit none
integer::j
open(unit=40, file = "output.log", status = "unknown")
!Write Input parameters
write(40,*),"initRandSeed        ="  ,initRandSeed
write(40,*),"restart          ="  ,restart
write(40,*),"nys              ="  ,nys
write(40,*),"method           ="  ,method                    
write(40,*),"TrpmdIntScheme =" ,  TrpmdIntScheme  
write(40,*),"Cayley           =" , Cayley            
write(40,*),"lambdaTRPMD       =" , lambdaTRPMD 
write(40,*),"thermostat       =" , thermostat                  
write(40,*),"thermostattype   =" , thermostattype                  
write(40,*),"us               =" , us
write(40,*),"Plumed           =" , plumed
write(40,*),"Natoms           =" , Natoms                       
write(40,*),"Nbeads           =" , Nbeads                          
write(40,*),"Ndim             =" , Ndim                           
write(40,*),"Nbonds           =" , Nbonds                     
write(40,*),"Nangles          =" , Nang                   
write(40,*),"Nref            =" , Nref                 
write(40,*),"Nchains          =" , Nchains                  
write(40,*),"write binary     =" , writebinary                  
write(40,*),"JinitFrame      =" , JinitFrame                 
write(40,*),"skipFrames      =" , skipFrames                 
write(40,*),"NtotFrames     =" , NtotFrames             
write(40,*),"NstepSave      =" , NstepSave               
write(40,*),"NrestartOut   =" , NrestartOut          
write(40,*),"NstepRemoveCoM =" , NstepRemoveCoM          
write(40,*),"dt               =" , dt            
write(40,*),"tau              =" , tau            
write(40,*),"simlen         =" , simlen              
write(40,*),"Temp             =" , Temp             
write(40,*),"sys restart file =" , trim(adjustl(inprestartTrajFile))
write(40,*),"mnhc_restart_file=" , trim(adjustl(inprestartmnhcTrajFile))
!Write Input geometry and mass 
write(40,*),'!=================================================================!'
write(40,*),'Atom Name, X,Y,Z, Mass'
!write coordinates and momenta to files
do j=0,Natoms-1
   write(40,'(A6,4E16.8)')currStep,Xinit(3*j:3*j+2)*BohrToAng,mass(j)*AuToAmu
enddo
end subroutine 
!****************************************************************************************************#  
subroutine writeRestartfiles()
    use allvars, only:method,MNHC,CLASSICAL,Nbeads,Natoms,Nchains,currStep,thermostattype,outrestartTrajFile,outrestartmnhcTrajFile,q,p,VelNose,QNose, &
					iounitRestOutTraj,iounitRestOutMNHC
    implicit none
    integer::i,j,k
    111 format (I12,6E16.8)
    open( unit=iounitRestOutTraj, file= outrestartTrajFile,status = "unknown" )
	if (method.ne.CLASSICAL) then
    	call NCtrans()
    	call NCtransP()
	endif
    !write coordinates and momenta to files
    do k =0,Nbeads-1
       do j=0,Natoms-1
             write(iounitRestOutTraj,111)currStep,q(3*j:3*j+2,k), p(3*j:3*j+2,k)
       enddo
    enddo

    !Write bath coordinates and velocities to files
  	if (thermostattype==MNHC) then
     	open( unit=iounitRestOutMNHC, file= outrestartmnhcTrajFile ,status ="unknown" )
	 	do i = 0,Nchains-1
     	   	do k = 0,Nbeads-1
     	       	do j = 0,Natoms-1
     	        	write(iounitRestOutMNHC,111)currStep,QNose(3*j:3*j+2,k,i),VelNose(3*j:3*j+2,k,i)
     	       	enddo
     	   	enddo
     	enddo
	endif
    close(iounitRestOutTraj)
    close(iounitRestOutMNHC)
    !write coordinates in xyz format to (for visualisation) files
end subroutine
!****************************************************************************************************#  
subroutine writeToFiles()
    use constants, only:BohrToAng
    use allvars,only:writebinary,method,CLASSICAL,Nbeads,Natoms,Nchains,currStep,q,p,centroid,centVel,Etot,Epot,Ering,Ek,instTemp,ENoseNoncent,ENoseCent,& 
					iounitTraj,iounitThermo,iounitXYZ,iounitCentTraj,usPotCent,iounitCentXYZ,atomname
    implicit none 
    integer::i,j,k
    111 format (I12,6E16.8)
    222 format (A6,3E24.16)
    333 format (I12,9E16.8)   
    call removeCoMCart()
    call getTotalEnergy()  
	if (method.ne.CLASSICAL) then
    	call getCentroid()
    	call getCentP()
    	call getCentVel()
	endif
    !write system trajectory to file
    if (writebinary==1) then
        write(iounitTraj)currStep,q,p
    else
        do k =0,Nbeads-1
           do j=0,Natoms-1
              write(iounitTraj,111),currStep,q(3*j:3*j+2,k),p(3*j:3*j+2,k)
           enddo
        enddo
    endif
    !write total energy and instantenous temperature to file 
    write(iounitThermo,333)currStep,Etot,Epot,Ering,Ek,instTemp,ENoseNoncent,ENoseCent,usPotCent

    !write centroid trjectory file 
    !do j = 0,Natoms-1
    !   write(iounitCentTraj,111)currStep,centroid(3*j:3*j+2),centroid_p(3*j:3*j+2)
    !enddo 

    !write coordinates in xyz format to (for visualisation) files
    !write(iounitXYZ,*) Natoms * Nbeads
    !write(iounitXYZ,*)
    !do j=0,Natoms-1
    !   do k =0,Nbeads-1
    !      write(iounitXYZ,222),atomname(j),q(3*j:3*j+2,k)*length_au
    !   enddo
    !enddo
     
    !!write coordinates in xyz format to (for visualisation) files
    !write(iounitCentXYZ,*) Natoms !* Nbeads
    !write(iounitCentXYZ,*)
    !do j=0,Natoms-1
    !      write(iounitCentXYZ,222)atomname(j),centroid(3*j:3*j+2)*BohrToAng
    !enddo
end subroutine 
!****************************************************************************************************#  
