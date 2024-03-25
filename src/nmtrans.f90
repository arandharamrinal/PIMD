!****************************************************************************************************#
!cartesian to Normal mode transformation
subroutine CNTrans()
    use allvars, only: Ncarts, u, q, Cmat, NonpreToPreNM
    implicit none
    integer:: j
    u = 0.d0
    do j = 0,Ncarts-1
          u(j,:) = u(j,:) + matmul(Cmat,q(j,:))
    enddo
    u = NonpreToPreNM * u
end subroutine
!****************************************************************************************************#
!NormalMode to Cartesian Transformation
subroutine NCTrans()
    use allvars, only : q, Ncarts, CmatInv, u, NonpreToPreCart
    implicit none
    integer:: j
    q = 0.d0
    do j = 0,Ncarts-1
          q(j,:) = q(j,:)  + matmul(CmatInv,u(j,:))
    enddo
    q = NonpreToPreCart * q 
end subroutine
!****************************************************************************************************#
subroutine CNTransP()
    use allvars, only: up, Ncarts, Cmat, p, NonpreToPreNM
    implicit none     
    integer:: j
    up=0.d0
    do j = 0,Ncarts-1
          up(j,:) = up(j,:) + matmul(Cmat,p(j,:))
    enddo
    up = NonpreToPreNM * up 
    end subroutine
!NormalMode to Cartesian Transformation
!****************************************************************************************************#
subroutine  NCTransP()
   use allvars, only : p, Ncarts, up, CmatInv, NonpreToPreCart
   implicit none     
   integer:: j
   p = 0.d0
   do j = 0,Ncarts-1
        p(j,:) = p(j,:)  + matmul(CmatInv,up(j,:))
   enddo
   p = NonpreToPreCart * p 
end subroutine
!****************************************************************************************************#
subroutine  CNPotTrans()
    use allvars, only : Ncarts, dVdu, Cmat, dVdq, NonpreToPreCart
    implicit none     
    integer :: j
    dVdu = 0.d0
    do j = 0,Ncarts-1
          dVdu(j,:) = dVdu(j,:) + matmul(Cmat,dVdq(j,:))
    enddo
    dVdu = NonpreToPreCart * dVdu
end subroutine
!****************************************************************************************************#
subroutine  NCPotTrans()
    use allvars, only : dVdq, Ncarts, CmatInv, dVdu, NonpreToPreNM
    implicit none     
    integer :: j
    dVdq = 0.d0
    do j = 0,Ncarts-1
          dVdq(j,:) = dVdq(j,:)  + matmul(CmatInv,dVdu(j,:))
    enddo
    dVdq = NonpreToPreNM * dVdq
end subroutine
!****************************************************************************************************#
