
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION epsilon,entry11,entry12,entry13,entry14,k21,k22,k23,k11,k12,k13
      DOUBLE PRECISION entry21,entry22,entry23,entry24,A1,A2,C1,C2
      DOUBLE PRECISION m1(3),d21qn(4),d22qn(4),d23qn(4),pi,L1,L2
      DOUBLE PRECISION m2(3),d11qn(4),d12qn(4),d13qn(4),dt23(3)
      DOUBLE PRECISION d11(3),d12(3),d13(3),strain1(3),uhat11,uhat12,uhat13
      DOUBLE PRECISION d21(3),d22(3),d23(3),strain2(3)
      DOUBLE PRECISION th

	  INTEGER I,I1
!       Right-hand side from Eqs. (4.2.9-12) in Dichmann, Li, and Maddocks,
!         "Hamiltonian Formulations and Symmetries in Rod Mechanics", in the
!         form to be fed to AUTO.


!------------------>  PARAMETRS: 1:Angle between tubes  2: Length of double overlap 3: Length of single overlap. 4,5,6:nx,ny,nz   7: Length of arm    <-------------------------

      uhat11=PAR(1)!1.0d0
      uhat12=0.0d0
      uhat13=PAR(2)
      
  
     A1=1.0
     A2=1000.0    ! To make it planar elasticas
     C1=1000.0/1.3 ! To make it planar elasticas
     
     
     L1=PAR(3)
      
     k21=A1
     k22=A2
     k23=C1


      m2(1) = (u(7)*u(8) + u(6)*u(9) - u(5)*u(10) - u(4)*u(11))/2.0d0
      m2(2) = (-u(6)*u(8) + u(7)*u(9) + u(4)*u(10) - u(5)*u(11))/2.0d0
      m2(3) = (u(5)*u(8) - u(4)*u(9) + u(7)*u(10) - u(6)*u(11))/2.0d0

!       The 4-vectors D_1[q] n, D_2[q] n, and D_3[q] n contain
!         the same 4 entries, but permuted and possibly negated.
!         Here are those entries

      entry21 =  2.0D0*(u(4)*u(12) + u(5)*u(13) + U(6)*u(14))
      entry22 =  2.0D0*(u(5)*u(12) - u(4)*u(13) + U(7)*u(14))
      entry23 =  2.0D0*(u(6)*u(12) - u(7)*u(13) - U(4)*u(14))
      entry24 =  2.0D0*(u(7)*u(12) + u(6)*u(13) - U(5)*u(14))

!       Components of D_1[q] n

      d21qn(1) =  entry21
      d21qn(2) =  -entry22
      d21qn(3) =  -entry23
      d21qn(4) =  entry24

!       Components of D_2[q] n

      d22qn(1) =  entry22
      d22qn(2) =  entry21
      d22qn(3) =  -entry24
      d22qn(4) =  -entry23

!       Components of D_3[q] n

      d23qn(1) =  entry23
      d23qn(2) =  entry24
      d23qn(3) =  entry21
      d23qn(4) =  entry22

!       Components of d1

      d21(1)= u(4)*u(4)-u(5)*u(5)-u(6)*u(6)+u(7)*u(7)
      d21(2)= 2.0d0*u(4)*u(5)+2.0d0*u(6)*u(7)
      d21(3)= 2.0d0*u(4)*u(6)-2.0d0*u(5)*u(7)

!       Components of d2

      d22(1)= 2.0d0*u(4)*u(5)-2.0d0*u(6)*u(7)
      d22(2)= -u(4)*u(4)+u(5)*u(5)-u(6)*u(6)+u(7)*u(7)
      d22(3)= 2.0d0*u(5)*u(6)+2.0d0*u(4)*u(7)

!       Components of d3

      d23(1)= 2.0d0*u(4)*u(6)+2.0d0*u(5)*u(7)
      d23(2)= 2.0d0*u(5)*u(6)-2.0d0*u(4)*u(7)
      d23(3)= -u(4)*u(4)-u(5)*u(5)+u(6)*u(6)+u(7)*u(7)



      strain2(1) = m2(1)/k21 + uhat11
      strain2(2) = m2(2)/k22 + uhat12
      strain2(3) = m2(3)/k23 + uhat13  !

!       Additional terms      

      F(1) = L1*( d23(1))
      F(2) = L1*(d23(2))
      F(3) = L1*(d23(3))
      F(4) = L1*(strain2(1)*u(7)/2.0d0 - strain2(2)*u(6)/2.0d0 + strain2(3)*u(5)/2.0d0)
      F(5) = L1*(strain2(1)*u(6)/2.0d0 + strain2(2)*u(7)/2.0d0 - strain2(3)*u(4)/2.0d0)
      F(6) = L1*(-strain2(1)*u(5)/2.0d0 + strain2(2)*u(4)/2.0d0 + strain2(3)*u(7)/2.0d0)
      F(7) = L1*(-strain2(1)*u(4)/2.0d0 - strain2(2)*u(5)/2.0d0 - strain2(3)*u(6)/2.0d0)
	F(8) =L1*(strain2(1)*u(11)/2.0d0 - strain2(2)*u(10)/2.0d0 + strain2(3)*u(9)/2.0d0 -d23qn(1))
	F(9) =L1*(strain2(1)*u(10)/2.0d0 + strain2(2)*u(11)/2.0d0 - strain2(3)*u(8)/2.0d0 -d23qn(2))
	F(10) =L1*(-strain2(1)*u(9)/2.0d0 + strain2(2)*u(8)/2.0d0 + strain2(3)*u(11)/2.0d0-d23qn(3))
	F(11) =L1*(-strain2(1)*u(8)/2.0d0 - strain2(2)*u(9)/2.0d0 - strain2(3)*u(10)/2.0d0-d23qn(4))
      F(12) = L1*0.0D0
      F(13) = L1*0.0D0
      F(14) = L1*0.0D0

  
 	
      RETURN
      END

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION k1,k2,k3,ubar1,ubar2,pi
      pi = 4.d0*ATAN(1.d0)
      PAR(1)=0.0d0 		!uhat11
      PAR(2)=0.0d0 		!uhat13
      PAR(3)=0.00001d0 
      PAR(4)=0.0d0 
      PAR(5)=0.0d0 
      PAR(6)=0.0d0 
      PAR(7)=0.00d0
      PAR(8)=0.00d0
      PAR(9)=0.00d0
	  
      U= 0.0d0
      
	 
      U(3)=PAR(1)*T! SIN(ubar1*T)/ubar1
 	  
      U(7)=1.0d0 !COS(ubar1*T/2.0d0)

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
      DOUBLE PRECISION pi,nx,ny,nz,e1,d111,d112,d113,d121,d122,d123,d131,d132,d133,e2,e3,Ndir1,Ndir2,Ndir3,Del1,Del2,Del3,bas
      INTEGER I1,I

      pi = 4.d0*ATAN(1.d0)
      
      nz=PAR(4)*COS(-PAR(9))
      ny=PAR(4)*SIN(-PAR(9))
      nx=PAR(4)*0.0d0

     
   	  e1=0.00d0
   	  e2=PAR(7)*SIN( PAR(8)) 	
   	  e3=PAR(7)*COS( PAR(8))   !3 interesting

   
     ! q1,q2,q3,q4, u20,u21,u22,u23   
      d111= u1(4)*u1(4) - u1(5)*u1(5) - u1(6)*u1(6) + u1(7)*u1(7)
      d112= 2.0d0*u1(4)*u1(5) + 2.0d0*u1(6)*u1(7)
      d113= 2.0d0*u1(4)*u1(6) - 2.0d0*u1(5)*u1(7)
   
   	  d121= 2.0d0*u1(4)*u1(5) - 2.0d0*u1(6)*u1(7)
      d122= -u1(4)*u1(4) + u1(5)*u1(5)- u1(6)*u1(6) + u1(7)*u1(7)
      d123= 2.0d0*u1(5)*u1(6) + 2.0d0*u1(4)*u1(7)

      d131= 2.0d0*u1(4)*u1(6) + 2.0d0*u1(5)*u1(7)
      d132= 2.0d0*u1(5)*u1(6) - 2.0d0*u1(4)*u1(7)
      d133= -u1(4)*u1(4) - u1(5)*u1(5) + u1(6)*u1(6) + u1(7)*u1(7)
   
     ! r X N
      Ndir1=nx*d111 + ny*d112 + nz*d113  
      Ndir2=nx*d121 + ny*d122 + nz*d123  
      Ndir3=nx*d131 + ny*d132 + nz*d133  
      
      
      
	  FB(1)=U0(1)
	  FB(2)=U0(2)
	  FB(3)=U0(3)
	  	  
	  FB(4)=U0(4) 
	  FB(5)=U0(5)
	  FB(6)=U0(6) - SIN(PAR(5)/2)
	  FB(7)=U0(7) - COS(PAR(5)/2)
	       
  	  !Moment on 1 tube ! q1,q2,q3,q4, u20,u21,u22,u23..... No (27 + #) or (#)  


      FB(8)= ( U1(7)*U1(8) + U1(6)*U1(9) - U1(5)*U1(10) - U1(4)*U1(11))/2.0d0 +(e2*Ndir3 - e3*Ndir2)
      FB(9)= (-U1(6)*U1(8) + U1(7)*U1(9) + U1(4)*U1(10) - U1(5)*U1(11))/2.0d0 +(e3*Ndir1 - e1*Ndir3)!
      FB(10)=( U1(5)*U1(8) - U1(4)*U1(9) + U1(7)*U1(10) -U1(6)*U1(11))/2.0d0 +(e1*Ndir2 - e2*Ndir1)
           
      FB(11)=  U1(8)*U1(4) + U1(9)*U1(5)+ U1(10)*U1(6)+ U1(11)*U1(7) + 2.0d0*(U1(1)*U1(12) + U1(2)*U1(13) + U1(3)*U1(14)) 
            
         
       ! mu.q + 2 r.n =0
   	  !FB(11) = U0(11)
     
      FB(12)=U1(12) + nx
      FB(13)=U1(13) + ny
      FB(14)=U1(14) + nz
        
   
      
      RETURN
      END
!                                                                      

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
