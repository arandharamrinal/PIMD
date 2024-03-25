subroutine  setupClassicalNVE()
use allvars!, only:
use brngvars 
!Nbeads,t2d_nm_factor,t2d_cart_factor,d2t_nm_factor,d2t_cart_factor,sqrt_Nbeads_inv,sqrt_Nbeads,Ndim,Natoms,restart,q,p 
implicit none
integer :: i,j,k,l,ii,dum,i_dim,frameNum
integer :: atmCtr,atmNum
real(8) :: rn,meanv(0:2),sigmav(0:2,0:2),randomv(0:2)
character(len=20)::framenamechar
initStep = 0
if (restart==1) then
    open( unit=iounitRestInTraj,file=inprestartTrajFile,status ="old" )
    do j = 0 , Natoms-1
        read(iounitRestInTraj,*)initStep,(q(3*j+i,k),i=0,Ndim-1),(p(3*j+ii,k),ii=0,Ndim-1)         
    enddo
else
    do frameNum= 1,NtotFrames
        open(unit=iounitTrajReadfile,file=trim(adjustl(inptrajFile)),status="unknown")
        NlinesPerFrame=Natoms*Nbeads
        if (frameNum>=JinitFrame) then
            write(framenamechar,*)frameNum
            if (writebinary==1) then
            open(unit=iounitTraj,file=trim(adjustl(trajFile))//trim(adjustl(framenamechar)),status="unknown",form='unformatted')
            else
              open(unit=iounitTraj,file=trim(adjustl(trajFile))//trim(adjustl(framenamechar)),status="unknown")
            endif
            open(unit=iounitThermo,file=trim(adjustl(thermoFile))//trim(adjustl(framenamechar)),status="unknown")
            open(unit=iounitXYZ,file=trim(adjustl(XYZFile))//trim(adjustl(framenamechar)),status="unknown")
            atmNum = 0
            atmCtr = 0
            do j=0,NlinesPerFrame-1
                read(iounitTrajReadfile,*),dum,(q(3*atmNum+i,0),i=0,Ndim-1),(p(3*atmNum,0),ii=0,Ndim-1)
                if (atmNum == Natoms-1) then
                     atmNum = 0
                else
                     atmNum = atmNum + 1
                endif
            enddo
            call getCoM()
            call removeCoMCart()
            call removeCoMVel()
            call classicalNVE() 
        endif
    enddo
endif
end subroutine 

