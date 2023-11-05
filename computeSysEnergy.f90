!*****************************************************************************************************!
!ENERGIES
!*****************************************************************************************************!
subroutine getRPkinE()
        use allvars, only : Nbeads,Ncarts,up,Ek,massKprimeInv
        implicit none
        integer ::  j, k
        Ek = 0.0d0
        do k = 0, Nbeads-1
            do j = 0, Ncarts-1
                    Ek = Ek + up(j,k) * up(j,k) *  massKprimeInv(j,k)
            end do
        end do
        Ek = 0.5d0 * Ek
    end subroutine getRPkinE
!*****************************************************************************************************!
subroutine getClassicalKE()
        use allvars, only : Ncarts,p,Ek,massKprimeInv
        implicit none
        integer ::  j, k
        Ek = 0.0d0
        do j = 0, Ncarts-1
                Ek = Ek + p(j,0) * p(j,0) *  massKprimeInv(j,0)
        end do
        Ek = 0.5d0 * Ek
end subroutine getClassicalKE  
!************************************************************************************************!
subroutine getVirialE()
use allvars, only : Evir, Nbeads, Ncarts, q, dVdq, NbeadsInv, extV
implicit none 
integer :: i,j,k
     Evir=0.d0
     do k =0,Nbeads-1
        do j=0,Ncarts-1
                 Evir =Evir +  0.5d0 * q(j,k) * dVdq(j,k) 
        enddo
        Evir =Evir +  extV(k) 
     enddo
     Evir = NbeadsInv*Evir
endsubroutine
!************************************************************************************************!
subroutine getRPenergyCart()
use allvars, only : Ering, Nbeads, Ndim, Natoms, massK, OmegaK, u 
implicit none
integer :: i,j, k
Ering = 0.d0
do k = 0, Nbeads-1
    do j = 0, Natoms-1
        do i = 0,Ndim-1  
            Ering = Ering +  massK(3*j+i,k)  *  OmegaK(k)  *  OmegaK(k)  *  u(3*j+i,k) * u(3*j+i,k)
        end do
    end do
end do
Ering = 0.5d0 * Ering
end subroutine
!************************************************************************************************!
subroutine getRPenergyNM()
use allvars, only : Ering, Nbeads, Ndim, Natoms, massK, omegaP2, u 
implicit none
integer :: i,j, k
Ering = 0.d0
do k = 0, Nbeads-1
    do j = 0, Natoms-1
        do i = 0,Ndim-1             
            Ering = Ering +  massK(3*j+i,k)  *  omegaP2  *  u(3*j+i,k) * u(3*j+i,k)
        end do
    end do
end do
Ering = 0.5d0 * Ering
end subroutine
!****************************************************************************************************#  
subroutine getTotalEnergy()
    use Gaussparam!use eg_pot_force
    use allvars, only : plumed,Ering, Epot, q, NbeadsInv,  method, ENoseNonCent, ENoseNonCent, thermostattype,& 
        instTemp, extV, Etot, Ek, Nbeads, Ndim, Natoms, massK, OmegaK, ENoseCent, totdofInv, &
		PIMDPC,CLASSICAL,MNHC,us,usPotCent
    use constants
    use plumedvars, only:plmdPot
    implicit none
    integer::k
    real(8),dimension(0:Nbeads-1)::Vt
    Etot = 0.d0;Epot=0.d0;Ek=0.d0
    ENoseNonCent = 0.d0
    ENoseCent = 0.d0
    if (method==PIMDPC) then
        Epot = sum(extV) * NbeadsInv
        call getRPenergyNM()
    else
        Epot = sum(extV)
        call getRPenergyCart()
    endif
	if (method==CLASSICAL) then
    	call getClassicalKE()
	else
    	call getRPkinE()
	endif
    if (((method.eq.PIMDPC).or.((method.eq.CLASSICAL))).and.(thermostattype==MNHC)) then  
       call getMNHCenergy()
    endif
    Etot = Epot + Ering  + Ek + ENoseNonCent + ENoseCent 
	if (plumed==1) then
		Etot = Etot + plmdPot*NbeadsInv
	endif
	if (us==1) then
		Etot = Etot + usPotCent
	endif
    instTemp =  2.d0 * Ek * totdofInv * AuToKelvin
end subroutine 
!****************************************************************************************************#  
