      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1     RPL,DDSDDT,DRPLDE,DRPLDT,
     2     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3     NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C     
      INCLUDE 'ABA_PARAM.INC'
C     
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4     JSTEP(4)

      DIMENSION STRAND(6),C(6,6),CD(6,6),DCDE(6,6),STRESS0(6)
      PARAMETER (ALPHA=1000.0,LAMBDA=0.0,DMAX=0.9999,DRND=3)
C****************************
C     STRAND.....STRAIN AT THE END OF THE INCREMENT
C     C.....6X6 STIFFNESS MATRIX
C     CD.....6X6 DAMAGED STIFFNESS MATRIX
C     DCDE.....D CD/D E 
C     STRESS0.....STRESS AT THE BEGINNING OF THE INCREMENT
C     STATEV(1).....DAMAGE VARIABLE D
C****************************
C     GET THE MATERIAL PROPERTIES 输入材料属性
      E1 = PROPS(1)              !E1,YOUNG'S MODULUS IN DIRECTION 1 
      E2 = PROPS(2)              !E2=E3,YOUNG'S MODULUS IN DIRECTION 2 & 3
      XNU12 = PROPS(3)           !POISON'S RATIO POI_12,XNU13=XNU12
      XNU21 = XNU12*E2/E1        !POISON'S RATIO POI_21,XNU31=XNU21  
      XNU23 = PROPS(4)           !POISON'S RATIO POI_23,XNU32=XNU23
      G12 = PROPS(5)             !G12=G13,SHEAR MODULUS IN 12 & 13 PLANE
      G23 = PROPS(6)             !G23,SHEAR MODULUS IN 23 PLANE
      STH = PROPS(7)             !FAILURE STRESS IN 1 DIRECTION IN TENSION
C****************************
C     STIFFNESS MATRIX C(6,6)    材料初始的刚度矩阵
      RNU = 1/(1-2*XNU12*XNU21-XNU23**2-2*XNU12*XNU21*XNU23)
      C = 0
      C(1,1) = E1*(1-XNU23**2)*RNU
      C(2,2) = E2*(1-XNU12*XNU21)*RNU
      C(3,3) = C(2,2)
      C(4,4) = G12
      C(5,5) = G12     
      C(6,6) = G23
      C(1,2) = E1*(XNU21+XNU21*XNU23)*RNU
      C(2,1) = C(1,2)
      C(1,3) = C(1,2)
      C(3,1) = C(1,2)
      C(2,3) = E2*(XNU23+XNU12*XNU21)*RNU
      C(3,2) = C(2,3)
C     CALCULATE THE STRAIN AT THE END OF THE INCREMENT 计算增量步结束时的应变
      DO I = 1, 6
          STRAND(I) = STRAN(I) + DSTRAN(I)
      ENDDO
C****************************
C     CALCULATE THE FAILURE COEFFICIENT      
c     损伤判据的计算
      STRANF = STH/E1 
      IF (STRAND(1)>0) THEN 
          F = STRAND(1)/STRANF
      ELSE
          F=0
      ENDIF
C     CALCULATE D,DAMAGE VARIABLE 损伤状态变量的计算
      D = STATEV(1)
      DDDE=0
c     损伤演化的过程
      IF (F>1) THEN
          DV=1-EXP(ALPHA*(1-F))
          IF (DV>D) THEN
              D = D*LAMBDA/(LAMBDA+1)+DV/(LAMBDA+1)
              DDDE=ALPHA*(1-DV)/STRANF/(1+LAMBDA)
              D=ANINT(D*10**DRND)/10**DRND
              DDDE=ANINT(DDDE*10**DRND)/10**DRND
          ENDIF
          IF (D>DMAX) THEN
              D=DMAX
          ENDIF
      ENDIF
      STATEV(1) = D   !UPDATE D
C     DAMAGED STIFFNESS MATRIX CD(6,6) 更新损伤后的刚度矩阵
      CD=C      
      CD(1,1)=(1-D)*C(1,1)
      CD(1,2)=(1-D)*C(1,2)
      CD(1,3)=CD(1,2)
      CD(2,1)=CD(1,2)
      CD(3,1)=CD(1,2)
      CD(4,4)=(1-D)*C(4,4)
      CD(5,5)=(1-D)*C(5,5)            
C****************************      
C     CALCULATE STRESS  更新损伤后材料的应力    
      DO I = 1, 6
          STRESS0(I) = STRESS(I)
          STRESS(I)=0.0
          DO J = 1, 6
              STRESS(I) = STRESS(I)+CD(I,J)*STRAND(J)
          ENDDO
      ENDDO
C     CALCULATE SSE,ELASTIC STRAIN ENERGY 计算应变能
      DO I = 1, NTENS
          SSE = SSE+0.5*(STRESS0(I)+STRESS(I))*DSTRAN(I)
      ENDDO
C****************************      
C     UPDATE THE JACOBIAN 更新雅可比矩阵
      DCDE=0      
      DCDE(1,1)=-DDDE*CD(1,1)
      DCDE(1,2)=-DDDE*CD(1,2)
      DCDE(1,3)=DCDE(1,2)
      DCDE(2,1)=DCDE(1,2)
      DCDE(3,1)=DCDE(1,2)
      DCDE(4,4)=-DDDE*CD(4,4)
      DCDE(5,5)=-DDDE*CD(5,5)
      DDSDDE=CD
      DO I = 1,6
          DO J=1,6
              ATEMP=DCDE(I,J)*STRAND(J)
          ENDDO
          DDSDDE(I,1)=DDSDDE(I,1)+ATEMP
      ENDDO
      RETURN
      END