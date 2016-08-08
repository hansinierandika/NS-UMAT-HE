***********************************************************************
C                    USER SUBROUTINE FOR NORSAND-2D VERSION
***********************************************************************

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
	  
      DOUBLE PRECISION  P,Q,PI,J2,J3,J3_J2,ETA,sin3teta,cos3teta,teta,
     &         XM_tc,emax,emin,XN,H,A,power,chi_tc,xnu,AJ2eta,AJ2eta_new
     &         CMN,FM,FDM,M_teta,M_teta_new,YIELD,G,BULK,DLAMDA,
     &          DF_DPI,DPI_DEFS,DF_DPI_DEFS,e_c,psi,e_c_i,psi_i,
     &          PI_max,DIL_max,DF_DP,DF_DQ,P_V_DSTRAN,P_S_DSTRAN
     &          DF_DSIGMA_q,DF_DSIG_v,DF_DSIG,DF_DM,DM_DTETA,DTETA_DSIG,
     &          V_DSTRAN,S_DSTRAN,DM_DTETA_BOTTOM,DM_DTETA_UP,
     &          DTETA_DSIG_PART1,DTETA_DSIG_PART2
C
C
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C
  
C
      DIMENSION  STRE(6),DDS(6,6),STRA(6),DSTRA(6),DPSTRAN(6),
     +           DESTRAN(6), DSTRESS(6),S(6),DF_DS(6),TDF_DS(6),
     +           DELASTIC(6,6),DPLASTIC(6,6),P_STRAN(6),P_E(6),
     +           DV1(6),DV2(6),DV3(6),DV4(6,6),DV5(6,6),TERM1(6,6),
     +           DF_DSIG(6),TDF_DSIG(6),DP_DSIG(6),DQ_DSIG(6),
     +           DTETA_DSIG(6),DJ2_DSIG(6),DJ3_DSIG(6),DF_DPC_DEF(6),
     +           DTETA_DSIG_PART3(6),DSUM(6,6),DF_DSIG_S(6),
     +           subdstran(6),subdstran_s(6),subdstran_s_t(6),
     +           XNDF_DSIG_S(6),DSTRESS_TANG(6),vector_ETA(6)

      WRITE(6,*) 'Entered the subroutine*********************'

C    ASSIGN VALUES OBTAINED FROM ABAQUS INTO INTERNAL VARIABLES

C    STRESS, STRAN AND DDSDDE
      DO I=1,6
       STRE(I)=STRESS(I)
       STRA(I)=STRAN(I)
       DSTRA(I)=DSTRAN(I)
      DO J=1,6
      DDS(I,J)=DDSDDE(I,J)
      END DO
      END DO
	  
      IF (NTENS .EQ. 4) THEN
      DO I=5,6
       STRE(I)=0
       STRA(I)=0
       DSTRA(I)=0
      DO J=5,6
      DDS(I,J)=0
      END DO
      END DO
	  
      END IF	   
	   
       WRITE(6,*) 'STRESS1', STRE(1),STRE(2),STRE(3),
     +                          STRE(4),STRE(5),STRE(6)
       WRITE(6,*) 'DSTRAN', DSTRA(1),DSTRA(2),DSTRA(3),
     +                          DSTRA(4),DSTRA(5),DSTRA(6)
       WRITE (6,*) 'STRAN1', STRA(1),STRA(2),STRA(3),
     +                          STRA(4),STRA(5),STRA(6)
	   
 
C   SPECIFY MATERIAL PROPERTIES
          
      XM_tc=PROPS(1)
      emax=PROPS(2)
      emin=PROPS(3)	
      XN=PROPS(4)
      HARD=PROPS(5)
      A=PROPS(6)
      power=PROPS(7)
      xnu=PROPS(8)
      chi_tc=PROPS(9)
      switch=PROPS(10)
      tol=PROPS(11)
          
       WRITE(6,*)'XM_tc',XM_tc,'emax',emax,'emin',emin ,'XN',XN,
     &     'HARD',HARD,'A',A,'POWER',power,'xnu',xnu,'chi_tc',chi_tc
       
C     DEFINE STATE VARIABLES
 
      e         =  STATEV(1)
      	  
C	  DEFINE ZERO MATRIX

      CALL ZERO1(S)
      CALL ZERO1(DSTRESS)
	  
      CALL ZERO1(DV1)
      CALL ZERO1(DV2)
      CALL ZERO1(DV3)
	  
      CALL ZERO2(DV4)
      CALL ZERO2(DV5)
	  
      CALL ZERO2(DELASTIC)
      CALL ZERO2(DPLASTIC)
      CALL ZERO2(TERM1)
      CALL ZERO2(DSUM)
	  
      CALL PQ(STRE,NTENS,P,Q,ETA,S)
	   
      WRITE(6,*) 'S1', S

      WRITE(6,*) 'STRESS invariants 1',  P,Q,ETA
      
      WRITE(6,*) 'void ratio 1', e
	  
      CALL  KTHETA(STRE,NTENS,STATEV,NSTATV,PROPS,NPROPS,XM_tc,
     * S,teta,J3,J2,CMN,cos3teta,M_teta)
	  
      
      WRITE(6,*)'M_teta 1', M_teta
	  
C    CALCULATE IMAGE PRESSURE

      CALL PIMAGE(P,ETA,XN,M_teta,PI)
	  
      WRITE(6,*)'Image pressure 1',PI
	  
C      SUBSTEPPING BEGIN
      WRITE(6,*)'Begin sub stepping******************'
	  	  
       T=0.0
      DT=1.0/5.0
      ICOUNT=0
	  
110   DO WHILE(T .LT. 1.0)
      ICOUNT=ICOUNT+1
	  
      WRITE(6,*)'ICOUNT********',ICOUNT
	  
      DO I=1,6
      subdstran(I)=DSTRA(I)*DT
      END DO
	  
C    SET UP ELASTICITY MATRIX

      CALL ELASTIC(P,A,power,xnu,NTENS,G,BULK,DELASTIC)
       
       WRITE(6,*)'Elastic stiffness matrix',DELASTIC
 	   
C     DETERMINE YIELD FUNCTION
      YIELD = Q-P*(M_teta/XN)*(1+(XN-1)*(P/PI)**(XN/(1-XN)))
	  
      
C      CHECK THE YIELD CRITERION	  
      WRITE(6,*) 'Yield',YIELD
	  
      IF (YIELD .LT. -1.0E-10) THEN
      WRITE(6,*)'ELASTIC RANGE'
       DLAMDA=0
      DO I=1,6
       DO J=1,6
        DDS(I,J) = DELASTIC(I,J)
        DSUM(I,J) = DSUM(I,J)+DDS(I,J)
       END DO
      END DO
      END IF
	  
	  
      IF (YIELD .GE. -1.0E-10) THEN
	  
      WRITE(6,*)'PLASTIC RANGE'

C     CALCULATE ELASTO PLASTIC STIFFNESS MATRIX

C     CALCULATE DIRECTION OF FLOW (DF_DSIG)

      CALL DFDSIG(XN,P,Q,ETA,PI,S, 
     &             M_teta,teta,J2,J3,cos3teta,CMN,DF_DSIG)
	  
      WRITE (6,*)'Plastic flow direction',DF_DSIG
	  

	  
C    CALCULATE CONTRIBUTION OF HARDENING

      CALL HARD_TERM(XN,XM_tc,HARD,emax,emin,chi_tc,e,
     &        STATEV,NSTATEV,P,Q,ETA,PI,S,M_teta,DF_DSIG,PI_max,TERM3)
	   
C      CALCULATION OF PLASTIC PART OF STIFFNESS MATRIX	   
	  
      CALL KMLT1(DELASTIC,DF_DSIG,DV1)
      CALL KMLT4(DF_DSIG,DELASTIC,DV2)
      CALL KMLT3(DV1,DV2,TERM1)
      
      CALL KMLT4(DF_DSIG,DELASTIC,DV2)
      CALL KMLT2(DV2,DF_DSIG,TERM2)
      
	  
      WRITE(6,*)'TERM1',TERM1
      WRITE(6,*)'TERM2',TERM2
      WRITE(6,*)'TERM3',TERM3
      
	  
      DO I=1,6
      DO J=1,6
      DPLASTIC(I,J)=TERM1(I,J)/(TERM2-TERM3)
      END DO
      END DO
	  
      WRITE(6,*)'DPLASTIC',DPLASTIC
	  
C     CALCULATE TOTAL STIFFNESS MATRIX (from elastic and plastic normal) 
	  
      DO I=1,6
      DO J=1,6
      DDS(I,J)=DELASTIC(I,J)-DPLASTIC(I,J)
      DSUM(I,J)=DSUM(I,J)+DDS(I,J)
      END DO
      END DO
	  
      WRITE(6,*)'DDS',DDS
       
C    CALULATE STRESS INCREMENT(from elastic and plastic normal)
      CALL KMLT1(DDS,subdstran,DSTRESS,6)
	  
      WRITE(6,*)'DSTRESS a',DSTRESS

C     CALCULATE TANGENTIAL STRAIN INCREMENT	  
	  
      CALL TANGENTIAL_STRAIN(DF_DSIG,subdstran,STATEV,NSTATEV, 
     &                                  subdstran_s_t)
C     CALCULATE TANGENTIAL STRAIN INCREMENT	
	  
       CALL TANGENTIAL_STRESS(P,PI,PI_max,S,G,subdstran_s_t,
     &                                     STATEV,NSTATEV,DSTRESS_TANG)
	 
C    REDUCE DEVIOTORIC TANGENTIAL STRESS INCREMENT	  
C      DO I=1,6
C      DSTRESS(I)=DSTRESS(I)-DSTRESS_TANG(I)
C      END DO
C      WRITE(6,*)'DSTRESS b',DSTRESS
	  
C    UPDATE STRESS 

      DO I=1,6
      STRE(I)=STRE(I)+DSTRESS(I)
      END DO
      WRITE(6,*)'STRESS2',STRE
	  
C    UPDATE STRAN
      DO I=1,6
	     STRA(I) = STRA(I)+ subdstran(I)  
      END DO
      WRITE(6,*)'STRAN2',STRA     
	  
C    CALCULATE STRAIN PARAMETERS
      CALL STRAN_PARA(subdstran,STRA,V_DSTRAN,S_DSTRAN,
     &                          V_STRAN,S_STRAN)
	 

C    CALCULATE VOID RATIO
      e=e+V_DSTRAN*(1.0+e)
        
      WRITE(6,*)'VOID',e
       
	   
C   UPDATE P,Q,PI
      CALL PQ(STRE,NTENS,P,Q,ETA,S)

      WRITE(6,*) 'S2', S

      WRITE(6,*) 'STRESS invariants 2',  P,Q,ETA
      
      WRITE(6,*) 'void ratio 2', e

C    CALCULATE M_teta

      CALL  KTHETA(STRE,NTENS,STATEV,NSTATV,PROPS,NPROPS,XM_tc,
     * S,teta,J3,J2,CMN,cos3teta,M_teta)
	  
      WRITE(6,*)'M_teta 2', M_teta	  
	  
C    CALCULATE IMAGE PRESSURE

      CALL PIMAGE(P,ETA,XN,M_teta,PI)
	  
      WRITE(6,*)'Image pressure 2',PI
	  
      END IF
	  
       T=T+DT
      END DO
	  
C     END SUBSTEPPING	

      WRITE(6,*)'End sub stepping******************'
	  
      DO I=1,6
      DO J=1,6
       DDS(I,J)=DSUM(I,J)/ICOUNT
      END DO
      END DO  

C   UPDATE P,Q,PI       

      CALL PQ(STRE,NTENS,P,Q,ETA,S)
	   
      WRITE(6,*) 'S3', S

      WRITE(6,*) 'STRESS invariants 3',  P,Q,ETA
      
      WRITE(6,*) 'void ratio 3', e
C    CALCULATE M_teta

      CALL  KTHETA(STRE,NTENS,STATEV,NSTATV,PROPS,NPROPS,XM_tc,
     &               S,teta,J3,J2,CMN,cos3teta,M_teta)
	 
      WRITE(6,*)'M_teta 3', M_teta	  
	  
C    CALCULATE IMAGE PRESSURE

      CALL PIMAGE(P,ETA,XN,M_teta,PI)
	  
      WRITE(6,*)'Image pressure 3',PI
      
	  
C    UPDATE STATE VARIABLES
      STATEV(1)=e
      STATEV(2)=P
      STATEV(3)=Q
      STATEV(4)=PI
      
      
      STATEV(7)=V_STRAN
      STATEV(8)=S_STRAN
	  
      STATEV(9)=PI_max
      STATEV(10)=M_teta
      STATEV(11)=psi
	  
      STATEV(18)=ETA
      
	  
  
	  
C    ASSIGN INTERNAL VARIABLES INTO ABAQUS VARIABLES

C    STRESS, STRAN AND DDSDDE
      DO I=1,6
       STRESS(I)=STRE(I)
       STRAN(I)=STRA(I)
       
      DO J=1,6
      DDSDDE(I,J)=DDS(I,J)
      END DO
      END DO
	  
      IF (NTENS .EQ. 4) THEN
      DO I=5,6
       STRESS(I)=0
       STRAN(I)=0
       DSTRAN(I)=0
      DO J=5,6
      DDSDDE(I,J)=0
      END DO
      END DO
      END IF	  

      WRITE(6,*)'***********************************END'
      RETURN
      END

****************************************************
C              INTERNAL SUBROUTINES
****************************************************

************************************************************
C               CALCULATE STRESS INVARIANTS
**************************************************
      SUBROUTINE PQ(STRE1,NTENS,P1,Q1,ETA1,S1)
		  
      INCLUDE 'ABA_PARAM.INC'
      DOUBLE PRECISION P1,Q1,ETA1,S1
	
      DIMENSION STRE1(NTENS)
      DIMENSION S1(6)
	  
      CALL ZERO1(S1)
	 
C	  CALCULATE HYDROSTATIC STRESS
      P1=(STRE1(1)+STRE1(2)+STRE1(3))/(-3.0)
	  
      IF (P1 .LT. 0.0)P1=1.0
	  
C	  CALCULATE EFFECTIVE STRESS
       S1(1)=-STRE1(1)-P1
       S1(2)=-STRE1(2)-P1
       S1(3)=-STRE1(3)-P1
       S1(4)=-STRE1(4)
       S1(5)=-STRE1(5)
       S1(6)=-STRE1(6)
	   
C     CALCULATE DEVIOTORIC STRESS
       Q1 = SQRT(STRE1(1)*(STRE1(1)-STRE1(2))+
     &    STRE1(2)*(STRE1(2)-STRE1(3))+ STRE1(3)*(STRE1(3)-STRE1(1))+ 
     &   3*STRE1(4)*STRE1(4)+3*STRE1(5)*STRE1(5)+3*STRE1(6)*STRE1(6))
	 
       IF (Q1.LT.0.001) Q1 = 0.001 
	  
       ETA1=Q1/P1 
      RETURN
      END
	  
****************************************************************
C      CALCULATE IMAGE PRESSURE
**************************************************	  
      SUBROUTINE PIMAGE(P1,ETA1,XN1,M_teta1,PI1)
      INCLUDE 'ABA_PARAM.INC'
	  
      DOUBLE PRECISION P1,ETA1,M_teta1,PI1,XN1
	  
C    CALCULATE IMAGE PRESSURE
      

      PI1 = P1*((1.0/(1.0-XN1))-
     &  (XN1/(1.0-XN1))*ETA1/M_teta1)**((XN1-1.0)/XN1)
	  
      IF((1.0/(1.0-XN1))-(XN1/(1.0-XN1))*ETA1/M_teta1 .LT. 0.0) PI1=P1
	  
      RETURN
      END
	  
******************************************************************
C       CLACULATE LODE ANGLE
*****************************************************	  
      SUBROUTINE KTHETA(STRE1,NTENS,STATEV,NSTATV,PROPS,NPROPS,XM_tc1,
     * S1,teta1,J3_1,J2_1,CMN1,cos3teta1,M_teta1)
	 
      INCLUDE 'ABA_PARAM.INC'
	  
      DOUBLE PRECISION teta1,J3_1,J2_1,CMN1,cos3teta1,M_teta1,J3_J2,
     *            sin3teta,FM,FDM,AJ2eta,AJ2eta_new,XM_tc1
      DIMENSION STATEV(NSTATV),STRE1(NTENS),PROPS(NPROPS)
      DIMENSION S1(6)
 
C    CALCULATE M_teta
      

C    CALCULATE LODE ANGLE TETA

      J3_1=S1(1)*S1(2)*S1(3)-S1(1)*S1(6)*S1(6)-S1(2)*S1(5)*S1(5)-
     &                   S1(3)*S1(4)*S1(4)+2*S1(4)*S1(5)*S1(6)
	 
      J2_1=(1/SQRT(6.0))*(SQRT((STRE1(1)-STRE1(2))**2+
     &             (STRE1(2)-STRE1(3))**2+(STRE1(3)-STRE1(1))**2+
     &               STRE1(4)**2+STRE1(5)**2+STRE1(6)**2))
	 
      IF (J2_1 .GT. -1.E-8 .AND. J2_1.LT.1.E-8) THEN
      J2_1=0
      J3_J2=0
      ELSE
      J3_J2=J3_1/J2_1**3
      END IF
	  
      sin3teta = 3.0*SQRT(3.0)*J3_J2/2.0
      IF (sin3teta .GT. 0.99) sin3teta =1.0 
      IF (sin3teta .LT. -0.99) sin3teta= -1.0
	  
      WRITE(6,*)'J2',J2_1
      WRITE(6,*)'J3',J3_1
      WRITE(6,*)'J3_J2',J3_J2
	  
      teta1 = -ASIN(sin3teta)/3.0
      IF (teta1 .GT. 0.523598) teta1 =0.523598 
      IF (teta1 .LT. -0.523598) teta1= -0.523598
	  
      cos3teta1=COS(3.0*teta1)
      sinteta=SIN(teta1)
	  
      IF (sinteta .GT. 0.499) sinteta=0.5
      IF (sinteta .LT. -0.499) sinteta=-0.5
	  
      IF (cos3teta1 .GT. -1.E-8 .AND. cos3teta1 .LT. 1.E-8) cos3teta=0
	  
      WRITE(6,*)'sin3teta',sin3teta
      WRITE(6,*)'cos3teta',cos3teta1
      WRITE(6,*)'teta',teta1
      WRITE(6,*)'sinteta',sinteta
	  
C    CALCULATE CMN
	  
      CMN1 = (27.0-3.0*XM_tc1**2.0)/(3.0-XM_tc1**2.0+
     * (2.0*XM_tc1**3.0/9.0))
      WRITE(6,*)'CMN',CMN1
	  
C    CALCULATE MATSOKO NAKAI FUNCTION OF M_teta
  
C    USE NEWTON METHOD TO FIND ROOTS OF M_teta
      AJ2eta=0.1
      tol=1.0e-5
      l=0.0
	  
 250   IF( l .LT. 50.0) THEN
 
         FM=(CMN1-3.0)*AJ2eta+2.0/SQRT(27.0)*CMN1*SIN(3.0*teta1)*
     &     AJ2eta**1.5-(CMN1-9.0)
	     FDM=(CMN1-3.0)+2.0/SQRT(27.0)*CMN1*3.0/2.0*SIN(3.0*teta1)*
     &    AJ2eta**0.5
	  
      AJ2eta_new = AJ2eta-(FM/FDM)
      END IF
      IF (ABS(AJ2eta_new-AJ2eta) .GT. tol .AND. l .LT. 50.0) THEN
      WRITE(6,*)'FM,FDM,AJ2eta_new',FM,FDM,AJ2eta_new
      AJ2eta = AJ2eta_new
      l = l+1 
      GO TO 250
      ELSE
      AJ2eta = AJ2eta_new
      END IF
      WRITE(6,*)'l',l
	  
      WRITE(6,*)'J2eta', AJ2eta
	  
      M_teta1=sqrt(3.0*AJ2eta)
	  
      RETURN
      END
	  
****************************************************************
C           CALCULATE ELASTIC STIFFNESS MATRIX
*******************************************************
      SUBROUTINE ELASTIC(P1,A1,power1,xnu1,NTENS,G,BULK,DELASTIC1)
C
      INCLUDE 'ABA_PARAM.INC'
	  
      DOUBLE PRECISION P1,A1,power1,xnu1,G,BULK
C
      DIMENSION DELASTIC1(6,6)
	  
      CALL ZERO2(DELASTIC1)

C    SET UP ELASTICITY MATRIX

      G = A1*ABS(P1)**power1
C   
      BULK = 2.0*(1.0+xnu1)*G/(3.0*(1-2.0*xnu1))
      
      DO K1 = 1, 3
         DO K2 = 1, 3
           DELASTIC1(K2,K1) = BULK-2.0*G/3.0
         END DO
        DELASTIC1(K1,K1) = BULK+4.0*G/3.0
      END DO
C
      IF (NTENS .EQ. 6.0) THEN
      DO K3 = 4,6
      DELASTIC1(K3,K3) = G
      END DO
      ELSE
      DELASTIC1(4,4) = G
      END IF
    
      RETURN
      END
****************************************************************
C         CALCULATE NORMAL PLASTIC FLOW DIRECTION
******************************************************
      SUBROUTINE DFDSIG(XN1,P1,Q1,ETA1,PI1,S1, 
     &             M_teta1,teta1,J2_1,J3_1,cos3teta1,CMN1,DF_DSIG)
C
      INCLUDE 'ABA_PARAM.INC'
	  
      DOUBLE PRECISION XN1,P1,Q1,ETA1,PI1,S1,
     &                 M_teta1,teta1,J2_1,J3_1,cos3teta1,CMN1,
     &                 DF_DP,DF_DQ,DF_DM,DM_DTETA_UP,DM_DTETA_BOTTOM,
     &                 DM_DTETA,DTETA_DSIG_PART1,DTETA_DSIG_PART2

      DIMENSION     S1(6),DP_DSIG(6),DQ_DSIG(6),DJ2_DSIG(6),DJ3_DSIG(6),
     &               DTETA_DSIG(6),DTETA_DSIG_PART3(6),DF_DSIG(6)
	 
C     CALCULATE DIRECTION OF FLOW (DF_DSIG)


      CALL ZERO1(DP_DSIG)
      CALL ZERO1(DQ_DSIG)
      CALL ZERO1(DJ2_DSIG)
      CALL ZERO1(DJ3_DSIG)
      CALL ZERO1(DTETA_DSIG)
      CALL ZERO1(DTETA_DSIG_PART3)
      CALL ZERO1(DF_DSIG)

C      Calculate DF_DP

      DF_DP = -(M_teta1/XN1)*(1.0+(XN1-1.0)*(1.0+(XN1/(1.0-XN1)))*
     &           (P1/PI1)**(XN1/(1.0-XN1)))
	 
      WRITE (6,*)'DF_DP',DF_DP

C      Calculate DF_DQ
      DF_DQ = 1.0
	  
C      Calculate DF_DM
      DF_DM = -P1/XN1*(1.0+(XN1-1.0)*(P1/PI1)**(XN1/(1.0-XN1)))
      WRITE (6,*)'DF_DM',DF_DM

C      Calculate DP_DSIG	  
      DO K= 1,3
       DP_DSIG(K) = 1.0/3.0
      END DO
      WRITE (6,*)'DP_DSIG',DP_DSIG

C      Calculate DQ_DSIG	  
      DO K= 1,3
       DQ_DSIG(K) = (3.0/(2.0*Q1))*S1(K)
      END DO
	  
      DO K= 4,6
       DQ_DSIG(K) = (3.0/(2.0*Q1))*2.0*S1(K)
      END DO
	  
      WRITE (6,*)'DQ_DSIG',DQ_DSIG
	  
C      Calculate DM_DTETA    
      DM_DTETA_UP = (-2.0/9.0)*cos3teta1*CMN1*M_teta1**3
      DM_DTETA_BOTTOM = ((2.0/3.0)*(CMN1-3.0)*M_teta1+
     &           (2.0/9.0)*(CMN1*sin(3*teta1)*M_teta1**2))
     
      DM_DTETA = DM_DTETA_UP/DM_DTETA_BOTTOM
               
	 
      WRITE (6,*)'DM_DTETA',DM_DTETA_UP,DM_DTETA_BOTTOM,DM_DTETA
	 
C    Calculate DTETA_DSIG

C    Calculate DJ2_DSIG
      IF(ABS(J2_1) .LT. 1.E-17) THEN
	  
      DO K= 1,3
       DJ2_DSIG(K) = (1.0/(2.0))*S1(K)
      END DO
	  
      DO K= 4,6
       DJ2_DSIG(K) = (1.0/(2.0))*2.0*S1(K)
      END DO
	  
      ELSE
	  
      DO K= 1,3
       DJ2_DSIG(K) = (1.0/(2.0*J2_1))*S1(K)
      END DO
	  
      DO K= 4,6
       DJ2_DSIG(K) = (1.0/(2.0*J2_1))*2.0*S1(K)
      END DO
	  
      END IF
  
      WRITE (6,*)'DJ2_DSIG',DJ2_DSIG
	  
C    Calculate DJ3_DSIG
      
       DJ3_DSIG(1) = (2.0/3.0)*S1(2)*S1(3)-(1.0/3.0)*S1(1)*S1(3)-
     &   (1.0/3.0)*S1(1)*S1(2)-(2.0/3.0)*S1(6)**2+
     &    (1.0/3.0)*S1(4)**2+(1.0/3.0)*S1(5)**2
       DJ3_DSIG(2) = (2.0/3.0)*S1(1)*S1(3)-(1.0/3.0)*S1(2)*S1(3)-
     &   (1.0/3.0)*S1(1)*S1(2)-(2.0/3.0)*S1(5)**2+
     &    (1.0/3.0)*S1(4)**2+(1.0/3.0)*S1(6)**2
      DJ3_DSIG(3) = (2.0/3.0)*S1(1)*S1(2)-(1.0/3.0)*S1(1)*S1(3)-
     &   (1.0/3.0)*S1(3)*S1(2)-(2.0/3.0)*S1(4)**2+
     &    (1.0/3.0)*S1(5)**2+(1.0/3.0)*S1(6)**2
      DJ3_DSIG(4) = -2.0*S1(3)*S1(4)+2.0*S1(5)*S1(6)
      DJ3_DSIG(5) = -2.0*S1(2)*S1(5)+2.0*S1(6)*S1(4)
      DJ3_DSIG(6) = -2.0*S1(1)*S1(6)+2.0*S1(4)*S1(5)
	  
      WRITE (6,*)'DJ3_DSIG',DJ3_DSIG
	  
      IF (ABS(J2_1).LT. 1.E-17 .OR. ABS(cos3teta1) .LT. 1.E-17) THEN
	  
      DO K= 1,6
      DTETA_DSIG(K) = 0.0
      END DO
	  
      ELSE
      DTETA_DSIG_PART1=SQRT(3.0)/(2.0*cos3teta1*J2_1**3)
       DTETA_DSIG_PART2 = 3.0*J3_1/J2_1
      DO K= 1,6
      DTETA_DSIG_PART3(K)=3.0*J3_1/J2_1*DJ2_DSIG(K)
      DTETA_DSIG(K) = SQRT(3.0)/(2.0*cos3teta1*J2_1**3)*
     &                   ((3.0*J3_1*DJ2_DSIG(K)/J2_1)-DJ3_DSIG(K))
      END DO
	  
      END IF 
      
      WRITE (6,*)'DTETA_DSIG_PART1',DTETA_DSIG_PART1
      WRITE (6,*)'DTETA_DSIG_PART2',DTETA_DSIG_PART2
      WRITE (6,*)'DTETA_DSIG_PART3',DTETA_DSIG_PART3
      WRITE (6,*)'DJ3_DSIG',DJ3_DSIG
      WRITE (6,*)'DTETA_DSIG',DTETA_DSIG
	  
C    CALCULATE DF_DSIG (Normal Vector to the yield surface)

      DO K= 1,6
      DF_DSIG(K) = DF_DP*DP_DSIG(K)+DF_DQ*DQ_DSIG(K)+
     &                   DF_DM*DM_DTETA*DTETA_DSIG(K)             
      END DO
	  

      RETURN
      END
***************************************************************
C             CALCULATE HARDENING TERM
******************************************************
      SUBROUTINE HARD_TERM(XN1,XM_tc1,HARD1,emax1,emin1,chi_tc1,e,
     & STATEV,NSTATEV,P1,Q1,ETA1,PI1,S1,M_teta1,DF_DSIG1,PI_max1,TERM3)
    
C
      INCLUDE 'ABA_PARAM.INC'
	  
      DOUBLE PRECISION XN1,XM_tc1,HARD1,emax1,emin1,chi_tc1,e,
     &                 P1,Q1,ETA1,PI1,S1,M_teta1,
     &                 DF_DPI,e_c,psi,e_c_i,psi_i,DIL_max,PI_max,
     &                 DPI_DEFS,DF_DPI_DEFS,DF_DSIG_v,DF_DSIG_q,TERM3
C
      DIMENSION        STATEV(NSTATEV),DF_DSIG1(6)

C    CALCULATE CONTRIBUTION OF HARDENING

C    Calculate DF_DPI
     
      DF_DPI = -M_teta1*(P1/PI1)**(1+XN1/(1-XN1))
	  
      WRITE (6,*)'DF_DPI',DF_DPI

C    Calculate DPI_DEFS

      e_c = emax1-(emax1-emin1)/(10-log(P1))
      psi = e- e_c
	  
	  
      WRITE (6,*)'e_critical',e_c
      WRITE (6,*)'psi',psi
      e_c_i = emax1-(emax1-emin1)/(10-log(PI1))
      psi_i = e - e_c_i
	  
      WRITE (6,*)'e_critical_image',e_c_i
      WRITE (6,*)'psi_i',psi_i
	  
      DIL_max=chi_tc1*psi_i*M_teta1/XM_tc1
      PI_max1= P1*(1+(DIL_max*XN1/XM_tc1))**((XN1-1)/XN1)
	  
      WRITE (6,*)'Maximum dilation',DIL_max
      WRITE (6,*)'Maximum image pressure',PI_max1
	  
      DPI_DEFS = HARD1*(M_teta1/XM_tc1) *(PI_max1-PI1)
	  
      WRITE (6,*)'DPI_DEFS',DPI_DEFS
	  
C    Calculate DF_DPI_DEFS

      DF_DPI_DEFS = DF_DPI*DPI_DEFS
    
      WRITE (6,*)'DF_DPI_DEFS',DF_DPI_DEFS
	  
C    Calculate DF_DSIG_q
      DF_DSIG_v = (DF_DSIG1(1)+DF_DSIG1(2)+DF_DSIG1(3))/3.0
      DF_DSIG_q = SQRT(2.0/3.0)*SQRT((DF_DSIG1(1)-DF_DSIG_v)**2+
     &            (DF_DSIG1(2)-DF_DSIG_v)**2+(DF_DSIG1(3)-DF_DSIG_v)**2+
     &              2*DF_DSIG1(4)**2+2*DF_DSIG1(5)**2+2*DF_DSIG1(6)**2)
	 
      WRITE (6,*)'DF_DSIGMA_q',DF_DSIG_q
	  
      TERM3 = DF_DPI_DEFS*DF_DSIG_q
	  
      STATEV(44)=e_c
      STATEV(45)=psi
      STATEV(46)=e_c_i
      STATEV(47)=psi_i
      STATEV(48)=DIL_max
      
    
      RETURN
      END
	  
***************************************************************  
C           CALCULATE STRAIN INVARIANTS
******************************************************	  
      SUBROUTINE STRAN_PARA(subdstran1,STRA1,V_DSTRAN,S_DSTRAN,
     &                          V_STRAN,S_STRAN)
C
      INCLUDE 'ABA_PARAM.INC'

      DOUBLE PRECISION V_DSTRAN,S_DSTRAN,V_STRAN,S_STRAN
      DIMENSION subdstran1(6),STRA1(6)

C    CALCULATE TOTAL VOLUMETRIC STRAIN INCREMENT
       V_DSTRAN=0.0
      DO I=1,3
	     V_DSTRAN=V_DSTRAN+subdstran1(I)
      END DO
	  
C    CALCULATE TOTAL DEVIOTORIC STRAIN INCREMENT
            S_DSTRAN=SQRT(2.0/3.0)*SQRT((subdstran1(1)-V_DSTRAN/3.0)**2+
     &  (subdstran1(2)-V_DSTRAN/3.0)**2+(subdstran1(3)-V_DSTRAN/3.0)**2+
     &	  subdstran1(4)**2/2.0+subdstran1(5)**2/2.0+subdstran1(6)**2/2.0)
	 
C    CALCULATE TOTAL VOLUMETRIC STRAIN 
       
      V_STRAN=STRA1(1)+STRA1(2)+STRA1(3)
     
	  
C    CALCULATE TOTAL DEVIOTORIC STRAIN 
      S_STRAN=SQRT(2.0/3.0)*SQRT((STRA1(1)-V_STRAN/3.0)**2+
     &          (STRA1(2)-V_STRAN/3.0)**2+(STRA1(3)-V_STRAN/3.0)**2+
     &	         STRA1(4)**2/2.0+STRA1(5)**2/2.0+STRA1(6)**2/2.0)

    
      RETURN
      END
	  
***************************************************************
C     CALCULATE TANGENTIAL STRAIN INCTREMENT
******************************************************	  
       SUBROUTINE TANGENTIAL_STRAIN(DF_DSIG1,subdstran1,
     &                            STATEV, NSTATEV, subdstran_s_t)
C
      INCLUDE 'ABA_PARAM.INC'
	  
      DOUBLE PRECISION DF_DSIG_v,XMDF_DSIG_S,XND
	  
      DIMENSION        STATEV(NSTATEV),
     &                 DF_DSIG1(6),DF_DSIG_S(6),XNDF_DSIG_S(6),
     &                 subdstran1(6),subdstran_s(6) ,subdstran_s_t(6)
	 
      CALL ZERO1(DF_DSIG_S)
      CALL ZERO1(XNDF_DSIG_S)
      CALL ZERO1(subdstran_s)
      CALL ZERO1(subdstran_s_t)

C    Calculate DF_DSIG_S(Deviotoric part of normal vector to the yield surface)

      DF_DSIG_v = (DF_DSIG1(1)+DF_DSIG1(2)+DF_DSIG1(3))/3.0

       DO I=1,3
        DF_DSIG_S(I) = DF_DSIG1(I)-DF_DSIG_v
       END DO
       DO I=4,6
        DF_DSIG_S(I) = DF_DSIG1(I)
       END DO
      WRITE(6,*)'deviotoric vector',DF_DSIG_S
	  
C    Calculate modulus of DF_DSIG_S
	   
      XMDF_DSIG_S = SQRT(DF_DSIG_S(1)**2+DF_DSIG_S(2)**2+DF_DSIG_S(3)**2
     &                +DF_DSIG_S(4)**2+DF_DSIG_S(5)**2+DF_DSIG_S(6)**2)
	 
      IF (XMDF_DSIG_S .EQ. 0.0) XMDF_DSIG_S = 0.00001
	  
      WRITE(6,*)'Mod of deviotoric vector',XMDF_DSIG_S
	  
C       Calculate the unit vector along deviotoric part of normal vector

       DO I=1,6
        XNDF_DSIG_S(I) = DF_DSIG_S(I)/XMDF_DSIG_S
       END DO
	   
      WRITE(6,*)'deviotoric unit vector',XNDF_DSIG_S
	   
C     Transform strain increment into deviatoric part
       DO I=1,3
        subdstran_s(I) = subdstran1(I)-
     &         (subdstran1(1)+subdstran1(2)+subdstran1(3))/3.0
       END DO
       DO I=4,6
        subdstran_s(I) = subdstran1(I)
       END DO
      WRITE(6,*)'subdstran_s',subdstran_s
	  
C    Transform  deviatoric part of strain increment into tangential part
      CALL KMLT2(XNDF_DSIG_S,subdstran1,XND)
	  
       DO I=1,6
        subdstran_s_t(I) = subdstran_s(I)-XND*XNDF_DSIG_S(I)
       END DO
	   
       WRITE(6,*)'subdstran_s_t',subdstran_s_t
	   
      STATEV(12)=DF_DSIG1(4)
      STATEV(13)=DF_DSIG_S(4)
      STATEV(14)=subdstran1(4)
      STATEV(15)=subdstran_s(4)
      STATEV(16)=subdstran_s_t(4)
      STATEV(17)=XNDF_DSIG_S(4)
	  
      STATEV(21)=DF_DSIG1(1)
      STATEV(22)=DF_DSIG_S(1)
      STATEV(23)=subdstran1(1)
      STATEV(24)=subdstran_s(1)
      STATEV(25)=subdstran_s_t(1)
      STATEV(26)=XNDF_DSIG_S(1)
	  
      STATEV(27)=DF_DSIG1(2)
      STATEV(28)=DF_DSIG_S(2)
      STATEV(29)=subdstran1(2)
      STATEV(30)=subdstran_s(2)
      STATEV(31)=subdstran_s_t(2)
      STATEV(32)=XNDF_DSIG_S(2)
	  
      STATEV(33)=DF_DSIG1(3)
      STATEV(34)=DF_DSIG_S(3)
      STATEV(35)=subdstran1(3)
      STATEV(36)=subdstran_s(3)
      STATEV(37)=subdstran_s_t(3)
      STATEV(38)=XNDF_DSIG_S(3)
    
      RETURN
      END	
	  
***************************************************************
C     CALCULATE TANGENTIAL STRESS INCREMENT
******************************************************
      SUBROUTINE TANGENTIAL_STRESS(P1,PI1,PI_max1,S1,G1,subdstran_s_t1,
     &                                       STATEV,NSTATEV,DSTRESS_TANG)
C
      INCLUDE 'ABA_PARAM.INC'

      DOUBLE PRECISION P1,PI1,PI_max1,S1,G1,
     &                R,xmod_vector_ETA,ratio,T
	 
      DIMENSION  STATEV(NSTATEV),
     &           S1(6),subdstran_s_t1(6),vector_ETA(6),DSTRESS_TANG(6)
	 
      CALL ZERO1(vector_ETA)
      CALL ZERO1(DSTRESS_TANG)
	  
C      CALCULATION OF TANGENTIAL STRESS RATE 

      R = PI1/PI_max1
      WRITE(6,*)'R',R
	  
      DO I=1,6
	      vector_ETA(I) = S1(I)/P1
      END DO
	  
       xmod_vector_ETA = SQRT(vector_ETA(1)**2+vector_ETA(2)**2+
     &                       vector_ETA(3)**2+vector_ETA(4)**2+
     &                       vector_ETA(5)**2+vector_ETA(6)**2)
	 
      WRITE(6,*)'mod_vector_ETA',xmod_vector_ETA
	  
      ratio = xmod_vector_ETA/1.17
      T = 2.0*G1*0.02*ratio**1.0/P1*R**2
	  
      WRITE(6,*)'T',T
       
       DO I=1,6
        DSTRESS_TANG(I) = (4*G1**2/(10000+2*G1))*subdstran_s_t1(I)
       END DO
      WRITE(6,*)'DSTRESS_tang',DSTRESS_TANG
      
      STATEV(5)=R
      STATEV(6)=T
	  
      STATEV(19)=xmod_vector_ETA
      STATEV(20)=ratio
    
      STATEV(40)=DSTRESS_TANG(1)
      STATEV(41)=DSTRESS_TANG(2)
      STATEV(42)=DSTRESS_TANG(3)
      STATEV(43)=DSTRESS_TANG(4)
	
      RETURN
      END
	  
	  
***************************************************************
C     INITIALISE VOID RATIO
******************************************************
      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)

      STATEV(1)=0.55
    
      RETURN
      END
	  
	  
	  
***************************************************
C              MATH SUBROUTINES
***************************************************
**         MULTIPLY 6X6 MATRIX WITH 6X1 VECTOR    *
***************************************************
*USER SUBROUTINE
      SUBROUTINE KMLT1(DM1,DM2,DM)
C      
      INCLUDE 'ABA_PARAM.INC'

      DIMENSION DM1(6,6),DM2(6),DM(6)
C
      DO I=1,6
      X= 0.0
      DO K=1,6 
      Y=DM1(I,K)*DM2(K)
      X=X+Y
      END DO
      DM(I)=X
      END DO
      RETURN
      END
	  
***************************************************
**         MULTIPLY 1X6 MATRIX WITH 6X1 VECTOR (DOT PRODUCT)   *
***************************************************
*USER SUBROUTINE
      SUBROUTINE KMLT2(CM1,CM2,CM)
C      
      INCLUDE 'ABA_PARAM.INC'

      DIMENSION CM1(6),CM2(6)
C
      X= 0.0
      DO K=1,6 
      Y=CM1(K)*CM2(K)
      X=X+Y
      END DO
      CM =X
      RETURN
      END
	  
***************************************************
**         MULTIPLY 6X1 VECTOR WITH 1X6 VECTOR (DYDE PRODUCT)   *
***************************************************
*USER SUBROUTINE
      SUBROUTINE KMLT3(BM1,BM2,BM)
C      
      INCLUDE 'ABA_PARAM.INC'

      DIMENSION BM1(6),BM2(6),BM(6,6)
C
     
      DO I=1,6
      DO J=1,6
      BM(I,J)=BM1(I)*BM2(J)
      END DO
      END DO
      
      RETURN
      END
	  
***************************************************
**         MULTIPLY 1X6 MATRIX WITH 6X6 MATRIX    *
***************************************************
*USER SUBROUTINE
      SUBROUTINE KMLT4(AM1,AM2,AM)
C      
      INCLUDE 'ABA_PARAM.INC'

      DIMENSION AM1(6),AM2(6,6),AM(6)
C
      DO I=1,6
      X= 0.0
      DO K=1,6 
      Y=AM1(K)*AM2(K,I)
      X=X+Y
      END DO
      AM(I)=X
      END DO
      RETURN
      END
	  
	  
C	  DEFINE ZERO MATRIX

      SUBROUTINE ZERO1(Z1)     
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION Z1(6)

      DO I=1,6
      Z1(I)=0.0
      END DO
    
      RETURN
      END
	  
      SUBROUTINE ZERO2(Z2)     
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION Z2(6,6)

      DO I=1,6
      DO J=1,6
      Z2(I,J)=0.0
      END DO
      END DO
    
      RETURN
      END

