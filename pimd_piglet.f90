!************************************************************************************************!
!Integration Steps  
!************************************************************************************************!
!************************************************************************************************!
subroutine verlet_step()
    use eg_pot_force
    use allvars, only: q,p,dVdq,dt,cartmass_inv,ext_V,Nbeads_inv,Nbeads,Ncarts
    use omp_lib
    !Updates the velocity and momentums of PIMD simulation using verlet
    !algorithm. 
    !Everything is done in the cartesian space and does not remove the external
    !torque  
    !since thermostat is attached. 
    implicit none
    integer :: k, j
    !print*,q;read(*,*)
    !Get the external potential and the cartesian gradient(dVdq) 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
               call getForce(q(:,k),dVdq(:,k),ext_V(k),0)
    enddo
    !$OMP END PARALLEL DO
    !print*,dVdq;read(*,*)
    !call get_OH_pot()
    !call get_H2O_pot()
    !call grad_nh3()

    !update momentum(Half time step) due to external force.
    p = p - 0.5d0 * dt * dVdq

    ! Update position (full time step)
    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD

    if (Nbeads .eq. 1) then
        ! For a single bead, there are no free ring polymer terms to add,
        ! so we simply update the positions using the momentum, as in
        ! classical trajectories
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmass_inv(j)
        end do
    else
        ! For multiple beads, we update the positions and momenta for the
        ! harmonic free ring term in the Hamiltonian by transforming to
        ! and from normal mode space
        call free_ring_polymer_step()
    end if

    !Get the external potential and the gradient(cartesian) dVdq 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
               call getForce(q(:,k),dVdq(:,k),ext_V(k),0)
    enddo
    !$OMP END PARALLEL DO
    !call get_OH_pot()
    !call get_H2O_pot()
    !call grad_nh3()

    !Update momentum (half time step)
    p = p - 0.5d0 * dt * dVdq

end subroutine verlet_step
!****************************************************************************************************#  
subroutine pimd_piglet()
    use allvars, only : t,ctr_prod,sim_time,Nbeads,n_step_save,dt,psi_mat
    use brngvars
    use mkl_vsl_type
    use mkl_vsl
    implicit none
    integer :: i,j,k,l,chain   
    integer :: atom_ctr,bead_ctr,atom_num
    t=0.d0
    ctr_prod = 0
    do  while (t<sim_time)
         call CN_trans_p()
         errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psi_mat,Nbeads,VSL_MATRIX_STORAGE_DIAGONAL,mean_arr,covmat)
         do k = 0,Nbeads-1 
            call PIGLET(k)
         enddo 
         call NC_trans_p()
         call verlet_step() 
         call CN_trans_p()
         errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psi_mat,Nbeads,VSL_MATRIX_STORAGE_DIAGONAL,mean_arr,covmat)
         do k = 0,Nbeads-1 
            call PIGLET(k)
         enddo
         call NC_trans_p()
         if (mod(ctr_prod,n_step_save)==0) then
            call write_to_files()
         endif 
         ctr_prod = ctr_prod + 1    
         t = t + dt
         enddo    
    Close(6666)
    Close(7777)
    Close(8888)
    Close(9999)
end subroutine
!****************************************************************************************************#  
