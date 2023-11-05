subroutine initializeUSparam()
use usvars
use allvars, only :methodChar,sysName,tempChar,NbeadsChar,method,CLASSICAL,PIMDPC,PIMDNPC
implicit none 
integer::dumi,i,j,itype,NuniqIdx,found
character(len=100)::UCoutFile
integer,allocatable::allidx(:)
iounitUSparam = 50;iounitUSucdef=51;iounitUCout=52
open(unit=iounitUSparam,file='us.param',status='old')
read(iounitUSparam,*)numuc
allocate(uctype(0:numuc-1),usK(0:numuc-1),ucEq(0:numuc-1))
read(iounitUSparam,*)(uctype(i),i=0,numuc-1)
read(iounitUSparam,*)(usK(i),i=0,numuc-1)
read(iounitUSparam,*)(ucEq(i),i=0,numuc-1)
read(iounitUSparam,*)NstepSaveUC
if ((method==PIMDPC).or.(method==PIMDNPC)) then
	UCoutfile = "colvars_"//trim(adjustl(methodChar))//"_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"K_"//trim((adjustl(NbeadsChar)))//"beads.out"
elseif (method==CLASSICAL) then
	UCoutfile = "colvars_"//trim(adjustl(methodChar))//"_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"K.out"
endif
open(unit=iounitUCout,file=trim(adjustl(UCoutfile)),status='unknown',action='write')

open(unit=iounitUSucdef,file='usCoord.dat',status='old')
five_pt_1st_deriv = [1.d0,-8.d0,0.d0,8.d0,-1.d0]
five_pt_1st_deriv = five_pt_1st_deriv/12.d0
allocate(ucIdx(itype:numuc-1,0:7))
ucIdx = 0
do itype = 0,numuc-1
	if (uctype(itype) >2) then
		read(iounitUSucdef,*)(ucIdx(itype,i),i=0,7)
	elseif ((uctype(itype) <=2)) then
		read(iounitUSucdef,*)(ucIdx(itype,i),i=0,3)
	endif
		!allocate(allidx(0:7))
		!NuniqIdx = 0
		!allidx = ucIdx(itype,:)
		!do i = 0,7
		!	if (allidx(i)/=0) then
		!		found = 0
		!		do j = i+1,7
		!			if (allidx(j)/=0) then
		!				if (allidx(j)==allidx(i)) then
		!					found = 1
		!				endif
		!			endif
		!		enddo
		!		if (found == 0) NuniqIdx = NuniqIdx + 1
		!	endif
		!enddo
		!allocate(UCuniqIdxarr(0:NuniqIdx-1))
		!NuniqIdx = 0
		!do i = 0,7
		!	if (allidx(i)/=0) then
		!		found = 0
		!		do j = i+1,7
		!			if (allidx(j)/=0) then
		!				if (allidx(j)==allidx(i)) then
		!					found = 1
		!				endif
		!			endif
		!		enddo
		!		if (found == 0) then
		!			UCuniqIdxarr(Nuniqidx) = allidx(i)
		!			NuniqIdx  = Nuniqidx + 1
		!		endif 
		!	endif 
		!enddo
		!elseif ((uctype <=2)) then
		!read(iounitUSucdef,*)(ucIdx(itype,i),i=0,3)
		!allocate(allidx(0:3))
		!allidx(0:3) = ucIdx(itype,:)
		!do i = 0,3
		!	if (allidx(i)/=0) then
		!		found = 0
		!		do j = i+1,7
		!			if (allidx(j)/=0) then
		!				if (allidx(j)==allidx(i)) then
		!					found = 1
		!				endif
		!			endif
		!		enddo
		!		if (found == 0) NuniqIdx = NuniqIdx + 1
		!	endif
		!enddo
		!allocate(UCuniqIdxarr(0:NuniqIdx-1))
		!NuniqIdx = 0
		!do i = 0,3
		!	if (allidx(i)/=0) then
		!		found = 0
		!		do j = i+1,7
		!			if (allidx(j)/=0) then
		!				if (allidx(j)==allidx(i)) then
		!					found = 1
		!				endif
		!			endif
		!		enddo
		!		if (found == 0) then
		!			UCuniqIdxarr(Nuniqidx) = allidx(i)
		!			NuniqIdx  = Nuniqidx + 1
		!		endif 
		!	endif
		!enddo
		!endif
enddo
ucIdx = ucIdx -1
!UCuniqIdxarr = UCuniqIdxarr -1 
end subroutine 
!***********************************************************************************************!
subroutine getBond(r1,r2,lenr)
	use allvars, only:Ndim
    implicit none 
    real(8),intent(in)::r1(0:Ndim-1),r2(0:Ndim-1)
    real(8),intent(out)::lenr
    lenr = norm2(r2-r1)
end subroutine 
!***********************************************************************************************!
subroutine computeAng(r1,r2,r3,Ang)
	use allvars, only:Ndim
	use constants, only : RadianToDegree
    implicit none 
    real(8),intent(in)::r1(0:Ndim-1),r2(0:Ndim-1),r3(0:Ndim-1)
    integer::i,ctr
    real(8),dimension(0:Ndim-1)::v1,v2,u1,u2
    real(8)::nr1,nr2
    real(8),intent(out)::Ang
    !compute all the three
    v1 = r1-r2
    nr1  = norm2(v1)
    v2   = r3-r2
    nr2  = norm2(v2)
    u1 = v1/nr1
    u2 = v2/nr2
    Ang = dacos(dot_product(u1,u2)) * radianTodegree
	Ang =  mod(Ang,360.d0)
end subroutine 
!***********************************************************************************************!
subroutine ComputeDih(r1,r2,r3,r4,dih)
	use allvars, only:Ndim
	use constants, only : RadianToDegree
    implicit none 
    real(8),intent(in)::r1(0:Ndim-1),r2(0:Ndim-1),r3(0:Ndim-1),r4(0:Ndim-1) 
    integer::i,d_ctr
    real(8),dimension(0:Ndim-1)::a1,a2,v1,v2,v3,u1,u2
    real(8)::v1_dot,v2_dot
    real(8)::nr1,nr2
    real(8),intent(out)::dih
    dih = 0.d0
    !compute all the three
    v1 = r2-r1
    nr1  = norm2(v1)
    v2   = r3-r2
    nr2  = norm2(v2)
    v3   = r4-r3
    u1 = v1/nr1
    u2 = v2/nr2
    a1(0) = v1(1)*v2(2) - v2(1)*v1(2)
    a1(1) = v2(0)*v1(2) - v2(2)*v1(0)
    a1(2) = v1(0)*v2(1) - v1(1)*v2(0)
    a2(0) = v2(1)*v3(2) - v3(1)*v2(2)
    a2(1) = v3(0)*v2(2) - v3(2)*v2(0)
    a2(2) = v2(0)*v3(1) - v2(1)*v3(0)
    v1_dot = norm2(v2) * dot_product(v1,a2)
    v2_dot = dot_product(a1,a2)
    dih = datan2(v1_dot,v2_dot) * RadianToDegree 
    dih = mod(dih,360.d0)
    if (dih>180.d0) then
        dih = dih - 360.d0
    endif
end subroutine 
!***********************************************************************************************!
subroutine computeGradDih(atom_idx,X,GradDih)
use allvars, only:Ncarts,Ndim
use usvars, only : ucIdx,five_pt_1st_deriv,dx,dx_inv 
implicit none 
integer,intent(in)::atom_idx
real(8),intent(in)::X(0:Ncarts-1)
real(8),intent(out)::GradDih(0:Ndim-1)
integer::i,k,itype,istep,id1,id2,id3,id4
real(8)::dS_ar(0:4),d1,d2,d3,d4,r1(0:Ndim-1),r2(0:Ndim-1),r3(0:Ndim-1),r4(0:Ndim-1),Xtmp(0:Ncarts-1)
Xtmp = 0.d0;GradDih=0.d0;r1=0.d0;r2=0.d0;r3=0.d0;r4=0.d0
d1=0.d0;d2=0.d0;d3=0.d0;d4=0.d0
id1 = ucIdx(itype,0)
id2 = ucIdx(itype,1)
id3 = ucIdx(itype,2)
id4 = ucIdx(itype,3)
do i = 0,2
    dS_ar = 0.d0
    Xtmp = X
    do istep = -2,2
        Xtmp(3*atom_idx+i) = X(3*atom_idx+i) + istep * dx 
        r1 = Xtmp(3*id1:3*id1+2)
        r2 = Xtmp(3*id2:3*id2+2)
        r3 = Xtmp(3*id3:3*id3+2)
        r4 = Xtmp(3*id4:3*id4+2)
        call computeDih(r1,r2,r3,r4,dS_ar(istep+2))
    enddo
    if (abs(abs(dS_ar(2))-180.d0)<5.d0) then
        d1 = dS_ar(1)-dS_ar(0)
        d2 = dS_ar(2)-dS_ar(1)
        d3 = dS_ar(3)-dS_ar(2)
        d4 = dS_ar(4)-dS_ar(3)
        if (abs(d1)>10.d0) then
           if (d1>0.d0) then
                do k = 1,4
                   dS_ar(k) = dS_ar(k) - 360.d0
                enddo
           else
                do k = 1,4
                   dS_ar(k) = dS_ar(k) + 360.d0
                enddo
           endif
        else if (abs(d2)>10.d0) then
           if (d2>0.d0) then
                do k = 2,4
                   dS_ar(k) = dS_ar(k) - 360.d0
                enddo
           else
                do k = 2,4
                   dS_ar(k) = dS_ar(k) + 360.d0
                enddo
           endif

        else if (abs(d3)>10.d0) then
           if (d3>0.d0) then
                do k = 3,4
                   dS_ar(k) = dS_ar(k) - 360.d0
                enddo
           else
                do k = 3,4
                   dS_ar(k) = dS_ar(k) + 360.d0
                enddo
           endif

        else if (abs(d4)>10.d0) then
           k = 4
           if (d4>0.d0) then
                   dS_ar(k) = dS_ar(k) - 360.d0
           else
                   dS_ar(k) = dS_ar(k) + 360.d0
           endif
        endif
    endif
    GradDih(i) = dx_inv * dot_product(five_pt_1st_deriv(:),dS_ar(:))
enddo
end subroutine 
!***********************************************************************************************!
subroutine computeAnalyticGradAng(a0,a1,a2,Xcarts,theta,gradAng)
use allvars, only:Ncarts,Ndim,currStep
implicit none 
integer,intent(in):: a0,a1,a2
real(8),intent(in):: Xcarts(0:Ncarts-1)
real(8),intent(out):: theta,gradAng(0:Ncarts-1) 
real(8)::bond_vec1(0:Ndim-1),bond_vec2(0:Ndim-1),bond_uvec1(0:Ndim-1),bond_uvec2(0:Ndim-1),bond_len1,bond_len2,cos_theta,bond_dot_prod,sinth_inv
real(8)::rtmp1(0:Ndim-1),rtmp2(0:Ndim-1)
	bond_vec1 = 0.d0;bond_vec2=0.d0;bond_uvec1=0.d0;bond_uvec2=0.d0
	bond_len1 = 0.d0;bond_len2 = 0.d0;cos_theta=0.d0;
	bond_dot_prod=0.d0;sinth_inv=0.d0;theta= 0.d0
	gradAng = 0.d0
    bond_vec1    = Xcarts(3*a1:3*a1+2) - Xcarts(3*a0:3*a0+2)
    bond_vec2    = Xcarts(3*a1:3*a1+2) - Xcarts(3*a2:3*a2+2)
    bond_len1       = dsqrt(dot_product(bond_vec1,bond_vec1))
    bond_len2       = dsqrt(dot_product(bond_vec2,bond_vec2))
    bond_uvec1      = bond_vec1/bond_len1
    bond_uvec2      = bond_vec2/bond_len2
    bond_dot_prod   = dot_product(bond_vec1,bond_vec2)
    cos_theta       = bond_dot_prod/(bond_len1*bond_len2)
    theta           = dacos(cos_theta)
    sinth_inv = 1.d0/dsin(theta)
    rtmp1 = - ( bond_dot_prod * bond_vec1 / (bond_len1**3 * bond_len2 ) - bond_vec2 / (bond_len1*bond_len2) ) * sinth_inv
    rtmp2 = - ( bond_dot_prod * bond_vec2 / (bond_len1 * bond_len2**3 ) - bond_vec1 / (bond_len1*bond_len2) ) * sinth_inv
	!if (currStep>101800) then
	!	print*,theta,bond_len1,bond_len2
	!endif
	!if (currStep>102100) then
	!	stop
	!endif
    gradAng(3*a0:3*a0+2) = gradAng(3*a0:3*a0+2) + rtmp1
    gradAng(3*a2:3*a2+2) = gradAng(3*a2:3*a2+2) + rtmp2
    gradAng(3*a1:3*a1+2) = gradAng(3*a1:3*a1+2) - rtmp1 - rtmp2
endsubroutine 
!***********************************************************************************************!
subroutine computeAnalyticGradBond(b0,b1,Xcarts,bondlen,GradBond)
use allvars, only:Ncarts,Ndim
implicit none
integer,intent(in)::b0,b1 
real(8),intent(in)::Xcarts(0:Ncarts-1)
real(8),intent(out)::GradBond(0:Ncarts-1),bondlen
real(8)::rvec(0:Ndim-1),urvec(0:Ndim-1)
GradBond=0.d0;rvec=0.d0;urvec=0.d0
rvec  = Xcarts(3*b1:3*b1+2) - Xcarts(3*b0:3*b0+2)
bondlen = norm2(rvec)
urvec = rvec/bondlen
GradBond(3*b0:3*b0+2) = -urvec(:)
GradBond(3*b1:3*b1+2) = urvec(:)
end subroutine 
!***********************************************************************************************!
subroutine computeUSgrad(X,dVdxUS,usPot)
	use allvars, only:Ncarts,Ndim,CurrStep,dt
	use constants, only : AuToFemtoSec,BohrToAng,RadianToDegree
	use usvars
    implicit none 
    real(8),intent(in)::X(0:Ncarts-1)
    real(8),intent(out)::usPot
    real(8),intent(out)::dVdxUS(0:Ncarts-1)
    real(8)::dVdxTmp1(0:Ncarts-1),dVdxTmp2(0:Ncarts-1)
    real(8)::Duc,uc(0:numuc-1),Graduc(0:Ndim-1),a1,a2
    integer::i,j,itype,id1,id2,id3,id4,id5,id6,id7,id8
    real(8)::Xtmp(0:Ncarts-1),r1(0:Ndim-1),r2(0:Ndim-1),r3(0:Ndim-1),r4(0:Ndim-1),r5(0:Ndim-1),r6(0:Ndim-1),r7(0:Ndim-1),r8(0:Ndim-1)
    uc=0.d0;
	r1=0.d0;r2=0.d0;r3=0.d0;r4=0.d0
	r5=0.d0;r6=0.d0;r7=0.d0;r8=0.d0
	Graduc = 0.d0
    Xtmp   = 0.d0
    !Convert input coordinates unit from Bohr to angstrom
    Xtmp = X * BohrToAng
	usPot = 0.d0
	dVdxUS=0.d0; 	
	do itype = 0,numuc-1
		if (uctype(itype)==0) then
    		!index of the atoms forming the bond 
    		id1 = ucIdx(itype,0);id2 = ucIdx(itype,1)
			dVdxTmp1 = 0.d0
			call computeAnalyticGradBond(id1,id2,Xtmp,uc(itype),dVdxTmp1)
    		!Compute deviation from the equillibrium
    		Duc = uc(itype)-ucEq(itype)
    		!Compute umbrella potenial 
    		usPot = usPot +  0.5d0 * Duc * Duc * usK(itype)
    		dVdxUS = dVdxUS + usK(itype) * Duc * dVdxTmp1 * BohrToAng
		elseif (uctype(itype)==1) then
    		!index of the atoms forming the Angle 
    		id1 = ucIdx(itype,0);id2 = ucIdx(itype,1);id3 = ucIdx(itype,2);
    		!Compute umbrella potenial and gradient 
			dVdxTmp1 = 0.d0
			call computeAnalyticGradAng(id1,id2,id3,Xtmp,uc(itype),dVdxTmp1)
    		Duc = uc(itype)*RadianToDegree - ucEq(itype)  
    		usPot = usPot + 0.5d0 * Duc * Duc * usK(itype)  
			dVdxUS = dVdxUS + usK(itype) * Duc  * RadianToDegree * BohrToAng * dVdxTmp1
		!elseif (uctype(itype)==2) then 
    	!	!index of the atoms forming the dihedral
    	!	id1 = ucIdx(itype,0);id2 = ucIdx(itype,1);id3 = ucIdx(itype,2);id4 = ucIdx(itype,3)
    	!	!cartesian coordinates(in angstrom) of the dihedral atoms 
		!	r1 = Xtmp(3*id1:3*id1+2);r2 = Xtmp(3*id2:3*id2+2);r3 = Xtmp(3*id3:3*id3+2);r4 = Xtmp(3*id4:3*id4+2)
    	!	!Compute the dihedral 
    	!	call computeDih(r1,r2,r3,r4,uc(itype))
    	!	!Compute deviation from the equillibrium
		!	if (ucEq(itype)>180.d0) then
		!	   ucEq(itype) = ucEq -360.d0
		!	endif
    	!	Duc = uc(itype)-ucEq(itype)
    	!	!Compute umbrella potenial 
    	!	usPot = 0.5d0 * Duc * Duc * usK(itype)
		!	do j = 0,size(UCuniqIdxarr)-1
		!		i = UCuniqIdxarr(j)
    	!		!Compute gradient of uc and the cartesian force(hartee/angstrom)  w.r.t cartesian coordinate of atom 1 
    	!		call computeGradDih(i,Xtmp,Graduc)
    	!		dVdxUS(3*i:3*i+2) = dVdxUS(3*i:3*i+2) + usK(itype)*Duc*Graduc
		!	enddo
    	!	!Convert cartesian force from hartee/angstrom to hartee/bohr
    	!	dVdxUS = dVdxUS + dVdxUS*BohrToAng
		elseif (uctype(itype)==4) then 
    		!index of the atoms forming the umbrella coordinate
    		id1 = ucIdx(itype,0);id2 = ucIdx(itype,1);id3 = ucIdx(itype,2);
			id4 = ucIdx(itype,4);id5 = ucIdx(itype,5);id6 = ucIdx(itype,6);
    		!Compute umbrella potenial and gradient 
			call computeAnalyticGradAng(id1,id2,id3,Xtmp,a1,dVdxTmp1)
			call computeAnalyticGradAng(id4,id5,id6,Xtmp,a2,dVdxTmp2)
			uc(itype)=(a2-a1) * RadianToDegree 
    		Duc = uc(itype)-ucEq(itype)  
    		usPot =  usPot + 0.5d0 * Duc * Duc * usK(itype)  
			dVdxUS = dVdxUS + usK(itype) * Duc *(dVdxTmp2-dVdxTmp1) * RadianToDegree* BohrToAng;
		else 
			print*,'Umbrella coordinate type has not been implemented yet.' 
		endif
	enddo
	if (mod(CurrStep,NstepSaveUC)==0) then
		write(iounitUCout,'(F16.8)',advance='no')currStep*dt*AuToFemtoSec
		do itype = 0,numuc-1
			write(iounitUCout,'(F16.8)',advance='no')uc(itype)
		enddo
		do itype = 0,numuc-1
			write(iounitUCout,'(F16.8)',advance='no')usK(itype)
		enddo
		write(iounitUCout,'(F16.8)')usPot
	endif
    end subroutine 
!***********************************************************************************************!
