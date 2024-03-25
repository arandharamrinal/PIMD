!******************************************************************************************************!
subroutine   initializeAllvars()     
    use allvars; 
    use constants
    !use simvars   
    implicit none
    integer :: i, j, k
    character(len=25) :: tmp_var 
	iounitParams = 10
	iounitMassparam = 11
	iounitInitstr = 12
	iounitTraj = 20
	iounitThermo = 21 
	iounitXYZ = 23
	iounitCentTraj = 24
	iounitCentXYZ = 25
	iounitRestInTraj = 26
	iounitRestInMNHC = 27
	iounitTrajReadfile = 28
	iounitRestOutTraj = 29
	iounitRestOutMNHC = 30 
    open(unit=iounitParams, file = "input.param", status = "old")
    read(iounitParams,*)sysName
    read(iounitParams,*)massFile
    read(iounitParams,*)coordFile
    read(iounitParams,*)initRandSeed
    read(iounitParams,*)restart
    read(iounitParams,*)us                    
    read(iounitParams,*)plumed                    
    read(iounitParams,*)method                    
    read(iounitParams,*)dt            
    read(iounitParams,*)simlen              
    read(iounitParams,*)Temp             
    read(iounitParams,*)Tbeads
    read(iounitParams,*)Natoms                       
    read(iounitParams,*)Nbeads                          
    read(iounitParams,*)Ndim                           
    read(iounitParams,*)Nref                  
    read(iounitParams,*)thermostat                    
    read(iounitParams,*)thermostattype                    
    read(iounitParams,*)tau            
    read(iounitParams,*)nys
    read(iounitParams,*)Nchains                  
    read(iounitParams,*)TrpmdIntScheme
    read(iounitParams,*)Cayley
    read(iounitParams,*)lambdaTRPMD 
    read(iounitParams,*)writebinary                 
    read(iounitParams,*)JinitFrame                 
    read(iounitParams,*)skipFrames             
    read(iounitParams,*)NtotFrames             
    read(iounitParams,*)NstepSave               
    read(iounitParams,*)NrestartOut          
    read(iounitParams,*)inprestartTrajFile
    read(iounitParams,*)inprestartmnhcTrajFile
    N_int = 3*Natoms-6
    Ndih  = Natoms-3
    Ncarts = 3*Natoms
    call AllocateArrays()
    call ReadInpParams()
    totSteps = int(simlen/dt)
    tempKelvin = Temp
	CLASSICAL = 0
	PIMDNPC   = 1
	PIMDPC    = 2 
	RPMD      = 3
	TRPMD     = 4
	MNHC      = 0
	PILE      = 1
	WNLangevin = 2
	if  (Cayley==1) then
		if (trpmdIntscheme==0) then
			BCOCB = 1; BAOAB =0; OBCBO=0; OBABO = 0
		elseif (trpmdIntscheme==1) then
			BCOCB = 0; BAOAB =0; OBCBO=1; OBABO = 0
		endif
	elseif (Cayley==0) then
		if (trpmdIntscheme==0) then
			BCOCB = 0; BAOAB =1; OBCBO=0; OBABO = 0
		elseif (trpmdIntscheme==1) then
			BCOCB = 0; BAOAB =0; OBCBO=0; OBABO = 1
		endif
	endif
    !************************************************************************!
	call checkvars()
	call createAllfiles()
    !************************************************************************!
    !SI units to atomic unit conversion
    call SIToAuTrans()
    !Assign parameters  
    piN             =  pi/dble(Nbeads) 
    beta             =  1.d0/Temp
    omegaP          =  dble(Nbeads)/(beta*hbarAu)!mrinal
    omegaP2         =  omegaP*omegaP
    betaP           =  beta/dble(Nbeads)
    betaP2          =  betaP*betaP
    betaPInv       =  1.d0 / betaP
    betaP2Inv      =  1.d0/betaP2
    betaInv         = 1.d0/beta
    if ((method==PIMDNPC).or.(method==TRPMD)) then
            tau              =  25.d0*betaP*hbarAu 
    endif
    tau2             =  tau*tau
    NatomsInv       =  1.d0/dble(Natoms)
    NbeadsInv       =  1.d0/dble(Nbeads)
    sqrtNbeads      =  sqrt(dble(Nbeads))
    sqrtNbeadsInv  =  1.d0/ (sqrtNbeads)
    dtRef           =  dt /dble(Nref)
    wn               =  sqrtNbeads/(beta*hbarAu)!mrinal
    wn2              =  wn*wn
    betaNHC         = betaP
    Gke            = Nbeads * KbAu  * Temp
    TbeadsInv  = 1.d0 / Tbeads
    !If Usual PIMD method is used.

    if (method == 2) then 
        omegaP = wn
        omegaP2 = wn2
        betaNHC = beta
        Gke    =  KbAu * Temp
    endif 
    betaNHCinv         = 1.d0/betaNHC
    if ((method==RPMD).or.(method==TRPMD)) then
            totDoF          =  dble(Nbeads*Natoms*Ndim)-6
    elseif (method==CLASSICAL) then
            totDoF          =  dble(Natoms*Ndim)
    else
            totDoF          =  dble(Nbeads*Natoms*Ndim)
    endif
    totDoFInv      =  1.d0/totDoF
    totalMass       =  sum(mass)
    totalMassInv   =  1.d0 / totalMass
    Ethermal   =  KbAu * Temp
    !******************************************************************************************************!
    if (thermostattype == 0) then 
            call getMNHCmass()
            !Nose-Hoover
            nresCent = 2
            nresNonCent = 1
            call getSYparams()
    else if (thermostattype == 1) then
            call getPIGLETparam()
    else if (thermostattype == 2 ) then
        call getLangevinParam()
    endif
	open(unit=iounitGaussForce,file='AIenforcedip.dat',form='unformatted',status='unknown') 
endsubroutine     
!******************************************************************************************************!
subroutine SIToAuTrans()
use allvars,only:Temp,Tbeads,dt,tau,simlen,cartmass,mass,cartmassInv,massInv,Xinit,Natoms
use constants
implicit none
integer :: i,j
    !kelvin to a.u.(hartee/kB)
    Temp = Temp * KelvinToAu
    Tbeads = Tbeads * KelvinToAu
    !fs to a.u 
    dt   = dt * FemtoSecToAu
    tau  = tau * FemtoSecToAu    
    simlen   =  simlen * FemtoSecToAu
    cartmass = 0.d0
    do i = 0,Natoms-1
             mass(i) = mass(i) * AmuToAu
             cartmass(3*i:3*i+2) = mass(i)
             cartmassInv(3*i:3*i+2) = 1.d0/mass(i)
             massInv(i) = 1.d0/mass(i) 
    enddo
    Xinit = Xinit * AngToBohr
end subroutine 
!******************************************************************************************************!
subroutine ReadInpParams()
    use allvars,only:Natoms,atomname,Xinit,mass,identityMat,Ndim,massFile,coordFile,iounitMassparam,iounitInitstr
    implicit none
    integer:: i,j 
    open(unit=iounitMassparam, file = ""//trim(adjustl(massFile))//"", status= "old")
    open(unit=iounitInitstr,file = ""//trim(adjustl(coordFile))//"",status="old")
    !atom symbols and read xyz coordinates 
    read(iounitInitstr,*)
    read(iounitInitstr,*)
    do j = 0,Natoms-1
       read(iounitInitstr,*) atomname(j),(Xinit(3*j+i),i=0,Ndim-1)      
    enddo        
    !read atom masses
    do i = 0,Natoms-1
       read(iounitMassparam,*)mass(i)
    enddo

    !!read bond and angle pairs
    identityMat = 0.d0
    do i = 0,Ndim-1
       do j = 0,Ndim-1
          if (i==j) then
             identityMat(i,j) = 1.d0
          endif
       enddo
    enddo
endsubroutine 
!******************************************************************************************************!
subroutine getOmegaK()
    use allvars,only:Ncarts,Nbeads,omegaK,omegaP
	use constants, only:pi
    implicit none
    integer:: i,j,k
    omegaK = 0.d0 
    do k = 0,Nbeads-1
          omegaK(k) = 2.d0 * omegaP * dsin ( k * ( pi / dble(Nbeads ) ) )
    enddo
end subroutine 
!******************************************************************************************************!
