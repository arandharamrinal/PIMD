!****************************************************************************************************#  
subroutine getMNHCenergy()
     use allvars, only: ENoseNonCent, Nchains, Nbeads, Ncarts, velNose, qMass, betaNHCinv, &
         ENoseCent, QNose, ENoseLocal, massKprimeInv, up, qMassCent
     implicit none 
     integer::i,j,k,chain
     ENoseNonCent=0.d0
     do chain = 0,Nchains-1
        do k = 1,Nbeads-1
           do j = 0,Ncarts-1
              ENoseNonCent = ENoseNonCent + (0.5d0*velNose(j,k,chain)*velNose(j,k,chain)*qMass(j,k) + qNose(j,k,chain)* betaNHCinv)   
           enddo
        enddo
     enddo
     ENoseCent = 0.d0
     do chain = 0,Nchains-1
        do j = 0,Ncarts-1
                 ENoseCent = ENoseCent + (0.5d0 * velNose(j,0,chain) * velNose(j,0,chain) * qMassCent(j,chain) + qNose(j,0,chain) * betaNHCinv)   
        enddo
     enddo
     !do k = 0 , Nbeads-1
     !   do j = 0,Ncarts-1
     !         ENoseLocal (j,k)= 0.d0 
     !         do chain = 0,Nchains-1
     !            ENoseLocal(j,k) = ENoseLocal(j,k) + ( 0.5d0 * velNose(j,k,chain) * velNose(j,k,chain) * qMass(j,k) + qNose(j,k,chain)*betaNHCinv)   
     !         enddo
     !         ENoseLocal(j,k) = ENoseLocal(j,k) + 0.5d0 * up(j,k) * up(j,k) *  massKprimeInv(j,k)
     !   enddo
     !enddo
end subroutine 
