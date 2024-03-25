!******************************************************************************************************!
subroutine   createAllfiles()     
    use allvars, only: method,Nbeads,methodChar,nbeadsChar,tempChar,sysName,Temp,thermostat,fnamePrefix,trajFile,centroidTrajFile,thermoFile,& 
					centroidVelFile,XYZFile,centroidXYZFile,outrestartTrajFile,outrestartmnhcTrajFile,thermostattype,inptrajFile,CLASSICAL,&
                    PIMDNPC,PIMDPC,RPMD,TRPMD,MNHC,WNLangevin 
	implicit none 
    if (method==CLASSICAL) then
        methodChar = "classical"
    elseif (method==PIMDNPC) then
        methodChar = "pimd"
    elseif (method==PIMDPC) then
        methodChar = "pimd"
    elseif (method==RPMD) then
        methodChar = "rpmd"
    else if (method==TRPMD) then
        methodChar = "trpmd"
    else
        print*,'Wrong method!'
        stop
    endif
    write(tempChar,*)int(Temp)
    if (method==CLASSICAL) then
        if (thermostat==0) then
                fnamePrefix                 = trim(adjustl(methodChar))//"_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_NVE.out" 
                trajFile                    = trim(adjustl(methodChar))//"_Traj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_NVE.out" 
                centroidTrajFile            = trim(adjustl(methodChar))//"_CentroidTraj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_NVE.out" 
                thermoFile                  = trim(adjustl(methodChar))//"_Thermo_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_NVE.out" 
                centroidVelFile             = trim(adjustl(methodChar))//"_CentroidVel_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_NVE.out" 
                XYZFile                     = trim(adjustl(methodChar))//"_XYZ_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_NVE.xyz" 
                centroidXYZFile             = trim(adjustl(methodChar))//"_centroidXYZ_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_NVE.xyz" 
                outrestartTrajFile          = trim(adjustl(methodChar))//"_Traj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_restart_NVE.out"
        else
                fnamePrefix                 = trim(adjustl(methodChar))//"_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k.out" 
                trajFile                    = trim(adjustl(methodChar))//"_Traj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k.out" 
                centroidTrajFile            = trim(adjustl(methodChar))//"_CentroidTraj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k.out" 
                thermoFile                  = trim(adjustl(methodChar))//"_Thermo_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k.out" 
                centroidVelFile             = trim(adjustl(methodChar))//"_CentroidVel_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k.out" 
                XYZFile                     = trim(adjustl(methodChar))//"_XYZ_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k.xyz" 
                centroidXYZFile             = trim(adjustl(methodChar))//"_centroidXYZ_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k.xyz" 
                outrestartTrajFile          = trim(adjustl(methodChar))//"_Traj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_restart.out"
        endif
    else
        write(nbeadsChar,*)Nbeads
        fnamePrefix                 = trim(adjustl(methodChar))//"_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_"//trim(adjustl(nbeadsChar))//"beads.out" 
        trajFile                    = trim(adjustl(methodChar))//"_Traj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_"//trim(adjustl(nbeadsChar))//"beads.out" 
        centroidTrajFile            = trim(adjustl(methodChar))//"_CentroidTraj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_"//trim(adjustl(nbeadsChar))//"beads.out" 
        thermoFile                  = trim(adjustl(methodChar))//"_Thermo_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_"//trim(adjustl(nbeadsChar))//"beads.out" 
        centroidVelFile             = trim(adjustl(methodChar))//"_CentroidVel_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_"//trim(adjustl(nbeadsChar))//"beads.out" 
        XYZFile                     = trim(adjustl(methodChar))//"_XYZ_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_"//trim(adjustl(nbeadsChar))//"beads.xyz" 
        centroidXYZFile             = trim(adjustl(methodChar))//"_centroidXYZ_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_"//trim(adjustl(nbeadsChar))//"beads.xyz" 
        outrestartTrajFile          = trim(adjustl(methodChar))//"_Traj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_"//trim(adjustl(nbeadsChar))//"beads_restart.out"
    endif
    if (method==CLASSICAL) then
        if (thermostattype==MNHC) then
            outrestartmnhcTrajFile = trim(adjustl(methodChar))//"_Traj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_mnhc_restart.out"
        endif
    else
        if (thermostattype==MNHC) then
            outrestartmnhcTrajFile = trim(adjustl(methodChar))//"_Traj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_"//trim(adjustl(nbeadsChar))//"beads_mnhc_restart.out"
        endif
    endif
    if ((method==RPMD).or.(method==TRPMD)) then
               inptrajFile  = "pimd_Traj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k_"//trim(adjustl(nbeadsChar))//"beads.out" 
    else if (method==CLASSICAL) then
         if (thermostat==0) then
               inptrajFile  = "classical_Traj_"//trim(adjustl(sysName))//"_"//trim(adjustl(tempChar))//"k.out" 
         endif
    endif
endsubroutine 

