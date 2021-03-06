!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_Temp.F90                     !
!    CONTAINS: subroutine SolveImpEqnUpdate_Temp          !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for           !
!     temperature, and updates it to time t+dt            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SolveImpEqnUpdate_Temp_ibm
      use param
      use local_arrays, only : temp,rhs
      use decomp_2d, only: xstart,xend
      use ibm_param, only: forclo
      implicit none
      real, dimension(nx) :: amkl,apkl,ackl, fkl
      integer :: jc,kc,info,ipkv(nx),ic
      real :: betadx,ackl_b
      real :: amkT(nx-1),ackT(nx),apkT(nx-1),appk(nx-2)

!     Calculate the coefficients of the tridiagonal matrix
!     The coefficients are normalized to prevent floating
!     point errors.

      betadx=0.5d0*al*dt/pec

      amkl(1)=0.d0
      apkl(1)=0.d0
      ackl(:)=1.d0
      amkl(nx)=0.d0
      apkl(nx)=0.d0


!     Call to LAPACK library to factor tridiagonal matrix.
!     No solving is done in this call.


      do ic=xstart(3),xend(3)
       do jc=xstart(2),xend(2)

!     Normalize RHS of equation

        fkl(1)= 0.d0
        do kc=2,nxm
         ackl_b=1.0d0/(1.-ac3ssk(kc)*forclo(kc,jc,ic)*betadx)
         amkl(kc)=-am3ssk(kc)*forclo(kc,jc,ic)*betadx*ackl_b
         ackl(kc)=1.d0
         apkl(kc)=-ap3ssk(kc)*forclo(kc,jc,ic)*betadx*ackl_b
        fkl(kc)=rhs(kc,jc,ic)*ackl_b 
        end do
        fkl(nx)= 0.d0
        amkT=amkl(2:nx)
        apkT=apkl(1:nxm)
        ackT=ackl(1:nx) 
!     Solve equation using LAPACK library
        call dgttrf(nx,amkT,ackT,apkT,appk,ipkv,info)
        call dgttrs('N',nx,1,amkT,ackT,apkT,appk,ipkv,fkl,nx,info)
          
!      Update temperature field

        do kc=2,nxm
          temp(kc,jc,ic) = temp(kc,jc,ic) + fkl(kc)
        end do

       enddo
      end do

      return
      end
