subroutine AllocateArrays()
use allvars
implicit none 
        allocate(p(0:Ncarts-1,0:Nbeads-1),q(0:Ncarts-1,0:Nbeads-1))
        allocate(u(0:Ncarts-1,0:Nbeads-1),up(0:Ncarts-1,0:Nbeads-1))
        allocate(extV(0:Nbeads-1),extNMpot(0:Natoms-1,0:Nbeads-1))
        allocate(dVdq(0:Ncarts-1,0:Nbeads-1),dVdu(0:Ncarts-1,0:Nbeads-1),usGradCent(0:Ncarts-1))
        allocate(mass(0:Natoms-1),massInv(0:Natoms-1))
        allocate(massK(0:Ncarts-1,0:Nbeads-1),massKprime(0:Ncarts-1,0:Nbeads-1),massKprimeInv(0:Ncarts-1,0:Nbeads-1))
        allocate(CoMCart(0:Ndim-1),CoM(0:Ndim-1,0:Nbeads-1),CoMP(0:Ndim-1,0:Nbeads-1),CoMVel(0:Ndim-1,0:Nbeads-1))
        allocate(centroid(0:Ncarts-1),centVel(0:Ncarts-1),centP(0:Ncarts-1))
        allocate(A(0:Nbeads-1,0:Nbeads-1))
        allocate(Cmat(0:Nbeads-1,0:Nbeads-1),CmatInv(0:Nbeads-1,0:Nbeads-1),lamda(0:Nbeads-1))
        allocate(velNose(0:Ncarts-1,0:Nbeads-1,0:Nchains-1),qNose(0:Ncarts-1,0:Nbeads-1,0:Nchains-1))
        allocate(gNose(0:Ncarts-1,0:Nbeads-1,0:Nchains-1),qMassInv(0:Ncarts-1,0:Nbeads-1))
        allocate(qMass(0:Ncarts-1,0:Nbeads-1),qMassCent(0:Ncarts-1,0:Nchains-1),qMassCentInv(0:Ncarts-1,0:Nchains-1))
        allocate(Rg(0:Natoms-1))
        allocate(atomType(0:NumAtomType-1),Xinit(0:Ncarts-1))
        allocate(ysweight(0:nys-1))
        allocate(poly(0:Natoms-1,0:3,0:Nbeads-1),polyHalf(0:Natoms-1,0:3,0:Nbeads-1),Cayleypoly(0:Natoms-1,0:3,0:Nbeads-1),sqrtCayleypoly(0:Natoms-1,0:3,0:Nbeads-1))
        allocate(atomname(0:Natoms-1))
        allocate(MoI(0:Ndim-1,0:Ndim-1),MoIInv(0:Ndim-1,0:Ndim-1))
        allocate(Jcent(0:Ndim-1),torqueCent(0:Ndim-1),delQdot(0:Ndim-1))
        allocate(AngVelcent(0:Ndim-1), AngACent(0:Ndim-1))
        allocate(qDotDoublePrime(0:Ncarts-1),QdotPrime(0:Ndim-1))
        allocate(delFCent(0:Ndim-1), Fcent(0:Ncarts-1), forcePrime(0:Ndim-1), forceDoublePrime(0:Ncarts-1))
        allocate(identityMat(0:Ndim-1,0:Ndim-1))
        allocate(nmKEbead(0:Ncarts-1,0:Nbeads-1))
        allocate(nmExtPotBead(0:Ncarts-1,0:Nbeads-1))
        allocate(nmPotBead(0:Ncarts-1,0:Nbeads-1))
        allocate(nmEtotBead(0:Ncarts-1,0:Nbeads-1))
        allocate(OmegaK(0:Nbeads-1))
        allocate(ENoseLocal(0:Ncarts-1,0:Nbeads-1))
        allocate(JinPA(0:Ndim-1),eigMoI(0:Ndim-1),eigMoIInv(0:Ndim-1),torquePA(0:Ndim-1))
        allocate(sqrtMbyBetan(0:Ncarts-1),c1PIGLET(0:Nbeads-1),c2PIGLET(0:Nbeads-1),c1PIGLETFull(0:Nbeads-1),c2PIGLETFull(0:Nbeads-1),gammaPIGLET(0:Nbeads-1)) 
        allocate(cartmass(0:Ncarts-1),cartmassInv(0:Ncarts-1))
        allocate(pimdTorpmdMass(0:Ncarts-1,0:Nbeads-1))
        allocate(wdtys2Cent(0:nys-1),wdtys4Cent(0:nys-1),wdtys8Cent(0:nys-1),wdtys16Cent(0:nys-1))
        allocate(wdtys2NonCent(0:nys-1),wdtys4NonCent(0:nys-1),wdtys8NonCent(0:nys-1),wdtys16NonCent(0:nys-1))
        !TRPMD
        allocate( psiMat(0:Ncarts-1,0:Nbeads-1),gammaMat(0:Nbeads-1),gammaCmat(0:Nbeads-1))
        allocate( gammaCmatbar(0:Nbeads-1), sqrtgammaMat(0:Nbeads-1))
        !Initialisation
        Evir            =  0.d0
        instTemp        =  0.d0
        gNose           =  0.d0
        Rg                =  0.d0
        p                =  0.d0
        q                =  0.d0
        u                =  0.d0
        up               =  0.d0
        dVdu             =  0.d0
        Epot              =  0.d0
        extV            =  0.d0
        extNMpot           =  0.d0
        dVdq             =  0.d0
        CoM               =  0.d0
        CoMVel             =  0.d0
        CoMP             =  0.d0
        centroid         =  0.d0
        A                =  0.d0
        Cmat            =  0.d0
        CmatInv        =  0.d0
        velNose         =  0.01d0 
        qMass           =  0.d0
        qMassInv       =  0.d0
        qMassCent      =  0.d0
        qMassCentInv  =  0.d0
        lamda            =  0.d0
        c1PIGLET        =  0.d0
        c2PIGLET        =  0.d0
        c1PIGLETFull   =  0.d0
        c2PIGLETFull   =  0.d0
        gammaPIGLET     =  0.d0
        t                =  0.d0
        instE           =  0.d0
        ENose           =  0.d0
        EkeNose       =  0.d0
        Etot            =  0.d0
        ENoseCent      =  0.d0
        ENoseNonCent      =  0.d0
        JinPA       =  0.d0
        eigMoI          =  0.d0
        torquePA        =  0.d0
        pimdTorpmdMass =  0.d0
		gammaLT			 =  0.d0 
		CLT1			 =  0.d0
		CLT2			 =  0.d0
		usPotcent 		 =  0.d0 
		usGradCent 		 =  0.d0
endsubroutine 
!******************************************************************************************************!
