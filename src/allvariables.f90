!******************************************************************************************************!
module allvars
implicit none        
integer             :: initRandSeed,writebinary,Natoms,Ncarts,N_int,Nbeads, Ndim,plumed,us,NstepSave,method,currStep,totSteps
integer             :: BCOCB, BAOAB , OBCBO, OBABO
integer				:: NumAtomtype,CLASSICAL,PIMDPC,PIMDNPC,RPMD,TRPMD,NstepRemoveCoM,Nref,thermostat,thermostattype,initStep
integer             :: Cayley,trpmdIntScheme,skipFrames,MNHC,PILE,WNLangevin
integer				:: iounitParams,iounitInitstr,iounitMassparam,iounitTraj,iounitThermo,iounitXYZ,iounitCentTraj,iounitCentXYZ
integer             :: iounitRestInTraj,iounitRestInMNHC,iounitTrajReadfile,iounitRestOutTraj,iounitRestOutMNHC,iounitGaussForce


real(8)             :: omegaP, omegaP2,  totalEnergy, tau2,  piN,betaNHC,betaNHCinv,Tbeads,TbeadsInv
real(8)             :: Evir, dt, dtRef, Temp, tempKelvin, beta, wn, wn2, tau, NatomsInv, NbeadsInv, sqrtNbeads, sqrtNbeadsInv
real(8),allocatable :: cartmass(:),cartmassInv(:) 
real(8)             :: spring_const, t, betaInv,betaP, betaPInv, betaP2, betaP2Inv, Ek, Ering, simlen,usPotCent   
real(8)             :: instTemp,Ethermal,totalMass,totalMassInv,sysNMmass,sysNMmassInv,totDof,totDofInv
real(8),allocatable :: A(:,:), lamda(:), u(:,:), p(:,:), q(:,:), CoM(:,:),CoMCart(:), CoMP(:,:), CoMVel(:,:)
real(8),allocatable :: extNMpot(:,:), extV(:), dVdu(:,:), dVdq(:,:), usGradCent(:), mass(:), massInv(:), spring_mat(:,:), Xinit(:) 
real(8),allocatable :: massK(:,:),pimdTorpmdMass(:,:), massKprime(:,:),massKprimeInv(:,:)
real(8),allocatable :: centVel(:), centP(:), centroid(:), Rg(:), up(:,:)  
real(8),allocatable :: Cmat(:,:), CmatInv(:,:),CmatA(:),CmatB(:,:)
real(8),allocatable :: poly(:,:,:), polyHalf(:,:,:), CayleyPoly(:,:,:), sqrtCayleyPoly(:,:,:),identityMat(:,:)
real(8),allocatable :: nmKEbead(:,:),nmPotBead(:,:),nmExtPotBead(:,:),nmEtotBead(:,:),OmegaK(:)
integer             :: restart,NtotFrames,JinitFrame,NrestartOut,NlinesPerFrame
real(8)             :: Epot,phi
character(len=25)   :: sysName,testfile,tempChar,nbeadsChar,methodChar
character(len=200)  :: massFile,fnamePrefix,coordFile,outrestartTrajFile,outrestartmnhcTrajFile,inprestartTrajFile,inprestartmnhcTrajFile 
character(len=200)  :: inptrajFile,trajFile,centroidTrajFile,thermoFile,centroidVelFile,XYZFile,centroidXYZFile        


character(len=25),allocatable   :: atomname(:)
!Constants
!Nose-Hoover
real(8),allocatable :: gNose(:,:,:), velNose(:,:,:), qNose(:,:,:), qMass(:,:), qMassCent(:,:),qMassInv(:,:), qMassCentInv(:,:), ysweight(:)
real(8),allocatable :: ENoseLocal(:,:)
real(8)             :: instE,ENose ,Etot,ENoseNonCent, ENoseCent
real(8)             :: EkeNose,Gke 
integer             :: NresCent,NresNonCent,Nchains, nys
real(8),allocatable :: wdtys2Cent(:),wdtys4Cent(:),wdtys8Cent(:),wdtys16Cent(:),wdtys2NonCent(:),wdtys4NonCent(:),wdtys8NonCent(:),wdtys16NonCent(:)
!PIGLET
real(8),allocatable :: sqrtMbyBetan(:),c1PIGLET(:),c2PIGLET(:),c1PIGLETFull(:),c2PIGLETFull(:),gammaPIGLET(:)
real(8) ::zetaPIGLET 
!LANGEVIN
real(8)				:: gammaLT, CLT1, CLT2
!TRPMD
real(8),allocatable :: psiMat(:,:),gammaMat(:),gammaCmat(:),gammaCmatbar(:), sqrtgammaMat(:)
real(8)             :: lambdaTRPMD 
!ammonia
integer             :: Nbonds, Nang,Ndih
real(8),allocatable :: atomType(:)

!Rotational correction
real(8),allocatable :: torqueCent(:), torquePA(:), delqDot(:), delFcent(:),  eigMoI(:), eigMoIinv(:), AngACent(:)
real(8),allocatable :: forceDoublePrime(:), forcePrime(:),MoI(:,:),MoIinv(:,:), Jcent(:), JinPA(:), qDotPrime(:), qDotDoublePrime(:)
real(8),allocatable :: AngVelcent(:),Fcent(:)
real(8)             :: detMoI
!Pre to Non-preconditioned transformation factor 
real(8) ::PreToNonpreNM,PreToNonpreCart,NonpreToPreNM,NonpreToPreCart
end module 
