!******************************************************************************************************!
subroutine getSYparams()
    use allvars,only: nys, YSweight, dt, dtRef, NresCent, NresNonCent, wdtYS2Cent, wdtYS4Cent, & 
         wdtYS8Cent, wdtYS16Cent, wdtYS2NonCent, wdtYS4NonCent, wdtYS8NonCent, wdtYS16NonCent
    implicit none 
    integer :: iys
    real(8) :: dtYS
    !Nose-Hoover Thermostat MTS parameters  
    if ( nys .eq. 1 ) then
        YSweight(0)  = 1.d0
    else if ( nys .eq. 3 ) then
       YSweight(0)   = 1.d0/(2.d0 - 2.d0**(1.d0/3.d0))
       YSweight(1)   = 1.d0 - 2.d0*YSweight(0)
       YSweight(2)   = YSweight(0)
    else if ( nys .eq. 5 ) then
       YSweight(0)   = 1.d0/(4.d0 - 4.d0**(1.d0/3.d0))
       YSweight(1)   = YSweight(0)
       YSweight(2)   = 1.d0 - 4.d0*YSweight(0)
       YSweight(3)   = YSweight(0)
       YSweight(4)   = YSweight(0)
    endif
    do iys = 0,nys-1
       dtYS = (dt * YSweight(iys)) / dble(NresCent)
       wdtYS2Cent(iys)  = 0.5d0     * dtYS
       wdtYS4Cent(iys)  = 0.25d0    * dtYS
       wdtYS8Cent(iys)  = 0.125d0   * dtYS
       wdtYS16Cent(iys) = 0.0625d0  * dtYS
    enddo
    do iys = 0,nys-1
       dtYS = dtRef * YSweight(iys)/dble(NresNonCent)
       wdtYS2NonCent(iys)  = 0.5d0    * dtYS
       wdtYS4NonCent(iys)  = 0.25d0   * dtYS
       wdtYS8NonCent(iys)  = 0.125d0  * dtYS
       wdtYS16NonCent(iys) = 0.0625d0 * dtYS
    enddo
endsubroutine getSYparams
!******************************************************************************************************!
subroutine getMNHCmass()
  use allvars,only: method, Natoms, Ncarts, Ndim, Nbeads, Nchains, qMass, qMassInv, betaNHC, omegaK, & 
      qMassCent, qMassCentInv, tau, omegaP2
  use constants,only : pi

  implicit none
  integer:: i,j,k 
  real(8) :: omega_bath,omega_bath2
  call getOmegaK()
  !Assign masses for the Nose-Hoover chains. Note that the masses are different for Pre-Conditioned 
  !(PC) and Non-pre-conditioned (NPC) PIMD algorithm.
  !If NPC algorithm  is used
  if (method == 1) then 
     do j = 0,Natoms-1
        do i =0,Ndim-1
           do k = 1,Nbeads-1
              qMass(3*j+i,k)= 1.d0 / ( betaNHC * omegaK(k) * omegaK(k) )
              qMassInv(3*j+i,k)= ( betaNHC * omegaK(k) * omegaK(k) )
           enddo
        enddo
     enddo
     do j=0,Ncarts-1
        do k=0,Nchains-1
             qMassCent(j,k) = ( 4.d0 * tau * tau) / betaNHC 
             qMassCentInv(j,k) = 1.d0/qMassCent(j,k)
        enddo
     enddo
  !If PC is used 
  else if (method == 2) then 
     do j = 0,Natoms-1
        do i =0,Ndim-1
           do k = 1,Nbeads-1
              qMass(3*j+i,k)= 1.d0 / betaNHC / omegaP2
              qMassInv(3*j+i,k)= betaNHC * omegaP2
           enddo
        enddo
     enddo
     omega_bath = 2.d0 * pi / tau
     omega_bath2 = omega_bath * omega_bath
     do j=0,Ncarts-1
        do k=0,Nchains-1
             qMassCent(j,k) = 1.d0 / betaNHC / omega_bath2
             qMassCentInv(j,k) = betaNHC * omega_bath2
        enddo
     enddo
  endif
endsubroutine     
!******************************************************************************************************!
subroutine getPIGLETparam()
  use allvars, only : dt,method, Ncarts, Nbeads, cartmass, betaP, sqrtMbyBetan, omegaK, tau, & 
        gammaPIGLET, c1PIGLET, c2PIGLET, c1PIGLETFull, c2PIGLETFull
  implicit none 
  integer :: j,k
  c1PIGLET = 0.d0;c2PIGLET=0.d0
  c1PIGLETFull = 0.d0;c2PIGLETFull=0.d0
  !Get the wK values
  call getOmegaK()
  !Check whether pimd(manolopoulos) or RPMD is used.
  if ((method == 2).or.(method==0)) then 
    print*,"method must be 1 or 4 for piglet"
    call exit()
  endif
  !Check ceriotti paper for reference  
  do j = 0,Ncarts-1
    sqrtMbyBetan(j) = dsqrt(cartmass(j)/betaP) 
      do k = 0,Nbeads-1 
        if (k==0)then 
            gammaPIGLET(k) =  1.d0 / tau  
            c1PIGLET(k)  =  dexp(-dt/2.d0 * gammaPIGLET(k))
            c2PIGLET(k)  =  dsqrt(1-c1PIGLET(k)*c1PIGLET(k))
            c1PIGLETFull(k)  =  dexp(-dt * gammaPIGLET(k))
            c2PIGLETFull(k)  =  dsqrt(1.d0-c1PIGLETFull(k)*c1PIGLETFull(k))
        else 
            gammaPIGLET(k) =  2.d0 * omegaK(k) 
            c1PIGLET(k)  =  dexp(-dt/2.d0 * gammaPIGLET(k))
            c2PIGLET(k)  =  dsqrt(1.d0-c1PIGLET(k)*c1PIGLET(k))
            c1PIGLETFull(k)  =  dexp(-dt * gammaPIGLET(k))
            c2PIGLETFull(k)  =  dsqrt(1.d0-c1PIGLETFull(k)*c1PIGLETFull(k))
        endif        
      enddo     
  enddo     
endsubroutine  
!******************************************************************************************************!
subroutine getLangevinParam()
  use allvars, only : dt,method, Ncarts, cartmass, betaInv, sqrtMbyBetan,tau, & 
        gammaLT, CLT1, CLT2
  implicit none
  integer :: j,k
  !Check ceriotti paper for reference  
  do j = 0,Ncarts-1
    		sqrtMbyBetan(j) = dsqrt(cartmass(j)*betaInv) 
            gammaLT   =  1.d0 / tau  
            CLT1      =  dexp(-dt * 0.5d0 * gammaLT)
            CLT2      =  dsqrt(1.d0-CLT1*CLT1)
  enddo     
endsubroutine  
!******************************************************************************************************!
subroutine trpmdVars()
use allvars, only : dt,Nbeads,Ncarts,lambdaTRPMD,betaPInv,cartmass,Cmat,CmatInv,omegaK,c1PIGLET,c2PIGLET,c1PIGLETFull,c2PIGLETFull,sqrtMbyBetan,gammaCmatbar
implicit none 
integer::j,k,l
!Initialize gamma mat
gammaCmatbar = 0.d0
c1PIGLET = 0.d0;c2PIGLET=0.d0
c1PIGLETFull = 0.d0;c2PIGLETFull=0.d0
do j = 0,Ncarts-1
     sqrtMbyBetan(j) = dsqrt(cartmass(j)*betaPInv) 
enddo
do k = 0,Nbeads-1
     gammaCmatbar(k) = lambdaTRPMD * 2.d0 * omegaK(k)
     c1PIGLET(k)  =  dexp( -0.5d0*dt*gammaCmatbar(k))
     c2PIGLET(k)  =  dsqrt(1.d0-c1PIGLET(k)*c1PIGLET(k))
     c1PIGLETFull(k)  =  dexp( -dt*gammaCmatbar(k))
     c2PIGLETFull(k)  =  dsqrt(1.d0-c1PIGLETFull(k)*c1PIGLETFull(k))
enddo
endsubroutine 
!******************************************************************************************************!
subroutine initializeMNHC()
      use allvars, only : Ncarts, Nbeads, Nchains,  beta, qMassCent, & 
       velNose, qNose, qMass
      implicit none 
      real(8) :: rn, vsigma
      integer::i,j,k,l      
      do l = 0,Nchains-1
         do k = 0,Nbeads-1      
            do j = 0,Ncarts-1
               if (k==0) then
                  qNose(j,k,l) = 0.d0
                  vsigma = dsqrt(1.d0/beta/qMassCent(j,l))
                  call getNormalRand(rn)
                  velNose(j,k,l) = vsigma*rn
               else
                  qNose(j,k,l) = 0.d0
                  vsigma = dsqrt(1.d0/beta/qMass(j,k))
                  call getNormalRand(rn)
                  velNose(j,k,l) = vsigma*rn
               endif
            enddo
         enddo
      enddo
endsubroutine

