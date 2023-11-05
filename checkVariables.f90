subroutine checkVars()
use allvars, only: method,thermostat,Nbeads,cayley,thermostattype,CLASSICAL,PIMDNPC,PIMDPC,RPMD,TRPMD,MNHC,PILE,WNLangevin
implicit none 
if (method==CLASSICAL) then
	if ((thermostattype.ne.WNLangevin).and.(thermostattype.ne.MNHC)) then
		write(*,*)'ERROR!! Wrong thermostat is chosen.'
		write(*,*)'Program supports only Massive Nose-Hoover chains and White-noise Langevin thermostat for classical simulation.'
		stop
	else 
        print*,'Using Langevin thermostattype for temperature control'
	endif
	if (Nbeads.ne.1) then
       print*,'Wrong number beads is used. Nbeads must be 1 for classical Dynamics.'
       stop
    endif

else if  (METHOD==PIMDNPC) then
	if (thermostattype.ne.PILE) then
		write(*,*)'ERROR!! Wrong thermostat is chosen.'
		write(*,*)'Program supports only Path integral Langevin thermostat for Non-preconditioned path integral simulations.'
		stop
	else 
        print*,'Using PILE thermostattype for temperature control'
	endif

   	if (Cayley == 1) then
       	print*,'Using Cayley modification for Free ring polymer evolution'
   	endif 
else if  (METHOD==PIMDPC) then
	if (thermostattype.ne.MNHC) then
		write(*,*)'ERROR!! Wrong thermostat is chosen.'
		write(*,*)'Program supports only Massive Nose-Hoover chains for preconditioned path integral simulations.'
		stop
	endif
else if  (METHOD==TRPMD) then
	if (thermostattype.ne.PILE) then
		write(*,*)'ERROR!! Wrong thermostat is chosen.'
		write(*,*)'Program supports only Path integral Langevin thermostat for Non-preconditioned path integral simulations.'
		stop
	else 
        print*,'Using PILE thermostattype for temperature control'
	endif

   	if (Cayley == 1) then
       	print*,'Using Cayley modification for Free ring polymer evolution'
   	endif 
endif
endsubroutine 
