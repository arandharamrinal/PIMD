!****************************************************************************************************#
!Thermostats
!****************************************************************************************************#
subroutine  NHthermostatCent()
use allvars, only : Ncarts, nresCent, nys, Nchains, up, EkeNose, massKprimeInv, gNose, qMassCent, & 
         qMassCentInv, qNose, velNose, Gke, wdtYS2Cent, wdtYS4Cent, wdtYS8Cent, wdtYS16Cent
implicit none
integer ::chain,y,iys,nt
integer:: j,ibead
real(8)::scl,AA
do j = 0,Ncarts-1
   EkeNose = up(j,0) * up(j,0) * massKprimeInv(j,0) 
   scl = 1.d0
   gNose(j,0,0) = ( EkeNose - Gke ) * qMassCentInv(j,0)
   do chain = 1 , Nchains-1
        gNose(j,0,chain)   = (qMassCent(j,chain-1) * velNose(j,0,chain-1) * velNose(j,0,chain-1) - Gke) * qMassCentInv(j,chain)
   enddo
   do nt = 0,nresCent-1
      do iys = 0,nys-1
         !Update thermostat Velocities
         velNose(j,0,Nchains-1) = velNose(j,0,Nchains-1) + wdtYS4Cent(iys) * gNose(j,0,Nchains-1)

         !Loop over Nchains chain -update force and velocity
         do chain = 0,Nchains-2
             AA   = dexp( -wdtYS8Cent(iys) * velNose(j,0,Nchains-chain-1) )
             velNose(j,0,Nchains-chain-2) = velNose(j,0,Nchains-chain-2) * AA * AA + wdtYS4Cent(iys) * gNose(j,0,Nchains-chain-2) * AA
         enddo
         !Compute scale factor
         AA = dexp( - wdtYS2Cent(iys) * velNose(j,0,0) )
         scl = scl * AA
         gNose(j,0,0) = ( scl * scl * EkeNose - Gke ) * qMassCentInv(j,0)
         !update Thermostat position
         do chain =0,Nchains-1
             qNose(j,0,chain) = qNose(j,0,chain)  + velNose(j,0,chain) * wdtYS2Cent(iys)
         enddo
         !update thermostat velocity and force 
         do chain = 0 ,Nchains-2
              AA = dexp(-wdtYS8Cent(iys) * velNose(j,0,chain+1))
              velNose(j,0,chain)   = velNose(j,0,chain) * AA * AA + wdtYS4Cent(iys) * gNose(j,0,chain) * AA
              gNose(j,0,chain+1)   = (qMassCent(j,chain) * velNose(j,0,chain) * velNose(j,0,chain)-Gke)*(qMassCentInv(j,chain))
         enddo
         velNose(j,0,Nchains-1) = velNose(j,0,Nchains-1)  + wdtYS4Cent(iys) * gNose(j,0,Nchains-1)
      enddo
   enddo
   !update particle momentum
   up(j,0) = up(j,0) * scl
enddo
end subroutine

!****************************************************************************************************#  
subroutine  NHthermostat(ibead)
use allvars, only : Ncarts, NresNonCent, nys, Nchains, up, EkeNose, massKprimeInv, gNose, qMass, & 
         qMassInv, qNose, velNose, Gke, wdtYS2NonCent, wdtYS4NonCent, wdtYS8NonCent, wdtYS16NonCent
    implicit none
    integer,intent(in)::ibead
    integer :: chain,y,iys,nt,n_c
    integer:: j
    real(8)::scl,AA
    do j = 0,Ncarts-1
        EkeNose = up(j,ibead) * up(j,ibead) * massKprimeInv(j,ibead) 
        scl=1.d0
        gNose(j,ibead,0) =  ( EkeNose-Gke ) * qMassInv(j,ibead)
        do chain = 1 , Nchains-1
           gNose(j,ibead,chain)   = (qMass(j,ibead) * velNose(j,ibead,chain-1) * velNose(j,ibead,chain-1) - Gke)*qMassInv(j,ibead)
        enddo
        do  nt = 0,NresNonCent-1
            do   iys = 0,nys-1
                !UPdate thermostat Velocities
                velNose(j,ibead,Nchains-1) = velNose(j,ibead,Nchains-1) +  wdtYS4NonCent(iys) * gNose(j,ibead,Nchains-1)
                !Loop over Nchains chain -update force and velocity
                do chain = 0,Nchains-2
                     AA = dexp(-wdtYS8NonCent(iys) * velNose(j,ibead,Nchains-chain-1))
                     velNose(j,ibead,Nchains-chain-2) = velNose(j,ibead,Nchains-chain-2) * AA * AA + wdtYS4NonCent(iys) * gNose(j,ibead,Nchains-chain-2) * AA
                enddo
                !Compute s factor
                AA = dexp(-wdtYS2NonCent(iys) * velNose(j,ibead,0))
                scl = scl * AA
                gNose(j,ibead,0) =  ( scl * scl * EkeNose - Gke ) * qMassInv(j,ibead)
                !update Thermostat position
                do chain =0,Nchains-1
                     qNose(j,ibead,chain) = qNose(j,ibead,chain)  +  velNose(j,ibead,chain) * wdtYS2NonCent(iys) 
                enddo
                !update thermostat velocity and force 
                do chain =0,Nchains-2
                     AA = dexp( - wdtYS8NonCent(iys) * velNose(j,ibead,chain+1) )
                     velNose(j,ibead,chain) = velNose(j,ibead,chain) * AA * AA  +  wdtYS4NonCent(iys) * gNose(j,ibead,chain) * AA
                     gNose(j,ibead,chain+1)   = ( qMass(j,ibead) * velNose(j,ibead,chain) * velNose(j,ibead,chain) - Gke ) * qMassInv(j,ibead)
                enddo
                velNose(j,ibead,Nchains-1) = velNose(j,ibead,Nchains-1)  +   wdtYS4NonCent(iys) * gNose(j,ibead,Nchains-1)
            enddo
        enddo
        !update particle momentum
        up(j,ibead) = up(j,ibead) * scl
    enddo
end subroutine
