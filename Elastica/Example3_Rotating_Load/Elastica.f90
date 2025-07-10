
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION K

      !Equilibrium equations of Elastica. For details see SPC Dhanakoti (2025), "Stability analysis through folds: An end-loaded elastica with a lever arm"
      K=1.0d0 ! Stiffness
      F(1) = -SIN(U(5))  ! d/ds (x) 
      F(2) = COS(U(5))   ! d/ds (y)
      F(3) = 0.0d0       ! d/ds (nx)
      F(4) = 0.0d0       ! d/ds (ny)
      F(5) = U(6)/K      ! d/ds (theta) 
      F(6) = U(3)*COS(U(5)) + U(4)*SIN(U(5))  !d/ds (m)

 	
      RETURN
      END

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION pi
      pi = 4.d0*ATAN(1.d0)
      PAR(1)=0.0d0
      PAR(2)=0.0d0
      PAR(3)=0.0d0
      PAR(4)=0.0d0
      PAR(5)=0.0d0
  
      U= 0.0d0 
      U(2)=1.0d0*T 
 	 
!
      RETURN
      END
!                                                                      
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
      DOUBLE PRECISION pi,Arm,Arm_angle,P,P_ang,Bas_ang,fx,fy,del_r,del_t

	  !Boundary Conditions 
      pi = 4.d0*ATAN(1.d0)
      Arm=PAR(1)       ! Length of the arm (epsilon)
      Arm_angle=PAR(2) ! psi
      P=PAR(3)         ! Magnitude of the force P
      P_ang=PAR(4)     ! Alpha
      Bas_ang= PAR(5)  !Theta_o
      
      fx= P*SIN(P_ang)
      fy=-P*COS(P_ang)
      
      del_r=Arm*SIN(Arm_angle)
      del_t=Arm*COS(Arm_angle)
      
      FB(1)=U0(1)   ! x(0)
      FB(2)=U0(2)   ! y(0)
	 
      FB(3)=U1(3) - fx !nx(1)
      FB(4)=U1(4) - fy !ny(1)
	  
      FB(5)=U0(5) - Bas_ang !theta(0)
      FB(6)=U1(6) - fx*(-del_r*SIN(U1(5)) - del_t*COS(U1(5)) ) - fy*(del_r*COS(U1(5)) - del_t*SIN(U1(5))) !m(l)

         
      RETURN
      END
!                                                                      

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
