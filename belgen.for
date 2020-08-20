	PROGRAM BELGEN
	IMPLICIT REAL*8(A-H,O-Z)

!	Program to calculate isoscalar transition rates for various multipoles
!	according to internal report KVI77i(1981).
!	Written by M.N. Harakeh

	PARAMETER	(N_CASES =	1)
	COMMON R(3),AR(3),SR(3),L_MULT,AMASS,Z,A13,EPSILON
	COMMON /DEF/ BETA(3),LDFRM(3),VBETA,GAMMA,BETACOR,XK,BETA0_DEF,BETA3_DEF
	COMMON /INT/ DR,KR
	REAL*8 ME(3),MECOR(3)
	CHARACTER TYPE(2)*7,DEFDEN*1,BAND*2,KMD*2,APPR*1,TRAN_DEN*4
	CHARACTER DEFP*3
	DIMENSION R_I(3),R1(3),R2(3),R4(3),RL(3),BE(3),PS(3),XV(3)
! 	DIMENSION EX(N_CASES), BETA_FERMI(N_CASES),DEF_LEN(N_CASES)
    DIMENSION EX(100), BETA_FERMI(100), DEF_LEN(100)
	DATA TYPE(1),TYPE(2) /'VOLUME ','SURFACE'/
	DATA LDFRM /2,4,6/
	PARAMETER (FOURPI=12.56637062)

	OPEN(UNIT=9,FILE='BEL',STATUS='NEW')

!	Read input parameters.

  100	CALL CIOAX (' PR[oceed] or EX[it] :',KMD,'CD','PR')
	IF(KMD.EQ.'EX') STOP
	CALL CIOF  (' MASS OF THE NUCLEUS : ',AMASS)
	CALL CIOF  (' CHARGE OF THE NUCLEUS : ',Z)
	CALL CIOWRT(' WOODS-SAXON POTENTIAL;  1=VOLUME   :')
	CALL CIOI  ('                         2=SURFACE  :',INX)
	WRITE(6,1)
  1	FORMAT('  RADIUS AND DIFFUSENESS OF THE WOODS-SAXON OPTICAL MODEL POTENTIAL')
	CALL CIOF  (' RR/(A**1/3) :',RR)
	CALL CIOF  (' AR :',ARR)
	CALL CIOF  (' MAXIMUM INTEGRATION RADIUS :',RMAX)
	CALL CIOF  (' INTEGRATION STEP :',DR)
	CALL CIOI  (' L = THE MULTIPOLARITY OF THE TRANSITION :',L_MULT)
	KR=NINT(RMAX/DR)
	A13=AMASS**0.33333333
	R(1)=1.2*A13
	AR(1)=0.D0
	R(2)=RR*A13
	AR(2)=ARR
	R(3)= 1.115*A13-0.53/A13
	T=2.5D0
	IF(AMASS.EQ.16.)T=2.D0
	AR(3)=T/4.4
	L_TEMP=L_MULT
	IF(L_MULT.EQ.0.OR.L_MULT.EQ.1) L_TEMP=L_MULT+2
	XL=(3+L_TEMP)*(3+L_TEMP)
	BSP=9.*FLOAT(2*L_MULT+1)*R(1)**(2*L_TEMP)/(FOURPI*XL)
	IF(L_MULT.EQ.0) BSP= FOURPI*BSP/4.
	IF(L_MULT.EQ.1) BSP=BSP/4.

!	If DWBA or Coupled channel analysis was performed in the symmetric
!	or asymmetric rotational model special attention should be given
!	to the calculation of transition rates.

	CALL CIOWRT(' Were the deformation parameters obtained in the')  
	CALL CIOWRT(' DWBA or CC analysis using deformed transition')
	CALL CIOAX (' densities (Y[es] or N[o]) :',DEFDEN,'C',' ')

	IF(DEFDEN.EQ.'Y')THEN
	  XK=0.D0
	  CALL CIOWRT(' To which band does the state for which the')
	  CALL CIOWRT(' transition rate is to be calculated belongs')
	  CALL CIOAX (' GR[ound],GA[mma],BE[ta] or OC[tupole] band :',BAND,'C',' ')

	  IF (BAND.EQ.'GR'.AND.(L_MULT.EQ.2.OR.L_MULT.EQ.4)) THEN
	    CALL CIOWRT(' Should the approximate formula to calculate the')
	    CALL CIOAX (' transition rate be used (Y[es] or N[o]) :',APPR,'C',' ')
	    IF(INX.EQ.2) APPR='N'
	  ENDIF

	  CALL CIOWRT(' EXCITATION ENERGY OF THE STATE IN QUESTION')
	  CALL CIOF  (' EX :',EX_DEF)
	  CALL CIOWRT(' GROUND STATE DEFORMATIONS')
	  CALL CIOF  (' QUADRUPOLE DEFORMATION;BETA2 : ',BETA2_DEF)
	  CALL CIOF  (' HEXADECAPOLE DEFORMATION;BETA4 :',BETA4_DEF)
	  CALL CIOF  (' 6thPOLE DEFORMATION;BETA6 :',BETA6_DEF)
	  DEF_LEN2=BETA2_DEF*R(2)
	  DEF_LEN4=BETA4_DEF*R(2)
	  DEF_LEN6=BETA6_DEF*R(2)

!	For the states of the Gamma band.

	  IF(BAND.EQ.'GA')THEN
	    XK=2.D0
	    CALL CIOF  (' AMPLITUDE OF GAMMA VIBRATION;GAMMA :',GAMMA)
	    IF(L_MULT.EQ.4) THEN
	      CALL CIOWRT(' Hexadecapole vibrational amplitude for the')
	      CALL CIOF  (' coupling from g.s. to 4+(gamma);BETA4 :',BETA4_GAM)
	      DEF_LEN4G=BETA4_GAM*R(2)
	    ENDIF

!	For the Beta band states. For a Beta vibration a correction term
!	to take care of volume conservation should be included. Moreover,
!	for a monople transition a breathing mode coupling parameter
!	in addition to the Beta vibration amplitude should be included
!	in case this was also done in the DWBA or CC calculation.

	  ELSE IF(BAND.EQ.'BE')THEN
	    XK=0.D0
	    CALL CIOWRT(' BETA VIBRATIONAL AMPLITUDE FOR THE BETA BAND')
	    CALL CIOF  (' BETA-BETA2 :',BETA0_DEF)
	    CALL CIOWRT(' BETA VIBRATION CORRECTION TERM COUPLING')
	    CALL CIOF  (' PARAMETER; BETA0_COR :',BETA0COR)
	    IF(ABS(BETA0_DEF).LT.1.0E-6) BETA0_DEF=1.D0
	    DEF_LEN0=BETA0_DEF*R(2)
	    BETACOR=BETA0COR/BETA0_DEF
	    IF(L_MULT.EQ.0)THEN
	      CALL CIOWRT(' If in addition to the beta vibrational form')
	      CALL CIOWRT(' factor a breathing mode form factor [Tassie')
	      CALL CIOWRT(' model,i.e. Satchler form factor I] was used,')
	      CALL CIOF  (' then enter the amplitude BETA0_MON :',BETA0_MON)
	      DEF_LENM=BETA0_MON*R(2)
	    ENDIF

!	For the Octupole band states. For an Octupole vibration a correction
!	term to take care of center of mass motion should be included. For a
!	dipole transition an isoscalar dipole coupling parameter in addition
!	to the Octupole vibration amplitude should be included if indeed used
!	in the DWBA or CC calculation.

	  ELSE IF(BAND.EQ.'OC') THEN
	    CALL CIOWRT(' OCTUPOLE VIBRATIONAL AMPLITUDE FOR THE BAND')
	    CALL CIOF  (' BETA3 :',BETA3_DEF)
	    CALL CIOWRT(' OCTUPOLE VIBRATION CORRECTION TERM COUPLING')
	    CALL CIOF  (' PARAMETER;BETA1_COR :',BETA1COR)
	    IF(ABS(BETA3_DEF).LT.1.0E-6) BETA3_DEF=1.D0
	    DEF_LEN3=BETA3_DEF*R(2)
	    BETACOR=BETA1COR/BETA3_DEF
	    CALL CIOI  (' THE PROJECTION ON THE SYMMETRY AXIS K :',KPROJ)
	    XK=FLOAT(KPROJ)
	    IF(L_MULT.EQ.1) THEN
	      CALL CIOWRT(' If in addition to the octupole vibration form')
	      CALL CIOWRT(' factor a Harakeh-Dieperink isoscalar dipole form')
	      CALL CIOF  (' factor was used then enter amplitude BETA1 :',BETA1_GR)
	      DEF_LEN1=BETA1_GR*R(2)
	    ELSE IF(L_MULT.EQ.5) THEN
	      CALL CIOWRT(' If in addition to the octupole vibration form')
	      CALL CIOWRT(' factor a l=5 surface vibration form factor was')
	      CALL CIOF  (' used then enter amplitude BETA5 :',BETA5_GR)
	      DEF_LEN5=BETA5_GR*R(2)
	    ENDIF
	  ENDIF
	  NX=INX
	  IF(BAND.NE.'GR') NX=INX+1
	ELSE

!	For deformation parameters obtained using non-statically deformed
!	form factors.

	  CALL CIOWRT(' Were the deformation parameters obtained using')
	  IF(L_MULT.EQ.0) THEN
	    CALL CIOWRT(' Satchler''s VER1(version I) or VER2(version II) transition')
	  ELSEIF(L_MULT.EQ.1) THEN
	    CALL CIOWRT(' HARA[keh-Dieperink] or ORLA[ndini et al.] transition')
	ELSE
	    CALL CIOWRT(' SURF[ace] vibration or TASS[ie] type transition')
	ENDIF
	  CALL CIOAX (' density :',TRAN_DEN,'C',' ')
	  CALL CIOWRT(' Is deformation parameter BET[a] or deformation')
	  CALL CIOAX (' length BR read as input? :',DEFP,'C',' ')
	  IF(DEFP.EQ.'BR') THEN
	    CALL CIOWRT(' EXCITATION ENERGY , DEFORMATION LENGTH')
	    DO K=1,N_CASES
! 	    WRITE(*,*)K,EX(K)
! 	    CALL CIOFX (' EX :',EX(K),'>',*3)
        CALL CIOF  (' Ex :',EX(K))
	    CALL CIOF  (' Deformation length B*R :',DEF_LEN(K))
	    BETA_FERMI(K)=DEF_LEN(K)/R(2)
	    ENDDO
	  ELSE
	    CALL CIOWRT(' EXCITATION ENERGY , DEFORMATION PARAMETER')
	    DO K=1,N_CASES
! 	    CALL CIOFX (' EX :',EX(K),'>',*3)
        CALL CIOF  (' Ex :',EX(K))
	    CALL CIOF  (' BETA :',BETA_FERMI(K))
	    DEF_LEN(K)=BETA_FERMI(K)*R(2)
	    ENDDO
	  ENDIF
	  GOTO 4
  3	  K=K-1
  4	  CONTINUE
	ENDIF

	CALL SUMRUL

!	Set sum rule value obtained from Woods-Saxon geometry to that obtained
!	from Bernstein mass density.

! 	SR(2)=SR(3)

!	Write out input parameters.

	WRITE(9,101)AMASS,Z
  101	FORMAT('1','                 NUCLEUS',' MASS :',F8.1,',          CHARGE :',F8.1)
	WRITE(9,102)TYPE(INX),RR,ARR
  102	FORMAT(/1X,A7,' WOODS-SAXON POTENTIAL PARAMETERS',' RR/(A**1/3) :',F8.3,' Fm,   AR :',F5.3,' Fm')

!	Write out some calculated parameters.

	WRITE(9,103) R(1),AR(1),R(2),AR(2),R(3),AR(3),RMAX,DR
  103	FORMAT(//' SHAPE OF THE MASS DENSITY USED IN CALCULATING &
	1 THE TRANSITION RATES',&
	2      /27X,'RADIUS',10X,'DIFFUSENESS'&
	3      /28X,'(Fm)',14X,'(Fm)'&
	4      //' UNIFORM DENSITY',9X,F8.3,10X,F0.0&
	5      /' WOODS-SAXON POTENTIAL',3X,F8.3,10X,F0.0&
	6      /' BERNSTEIN MASS DENSITY',2X,F8.3,10X,F0.0&
	7     //' MAX. INTEG. RADIUS :'4X,F8.3,'   INTEG. STEP :',F5.3)
	WRITE(9,104)BSP,2*L_TEMP
  104	FORMAT(//' SINGLE PARTICLE UNIT (s.p.u.)   :',E15.6,'  e**2*Fm** &
	& ',I2)
	WRITE(9,105)2*L_TEMP,SR(1),SR(3)
  105	FORMAT(/' ENERGY WEIGHTED SUM RULE (EWSR) :   MeV*e**2*Fm**',I2, &
	&     //' UNIFORM DENSITY                 :',E15.6 &
	&     /' BERNSTEIN MASS DENSITY          :',E15.6)
	IF(DEFDEN.EQ.'Y') THEN
	  WRITE(9,106)
  106	  FORMAT(//' THE FOLLOWING RESULTS HAVE BEEN CALCULATED UNDER THE' &
	1       /' ASSUMPTION THAT ALL DEFORMATION PARAMETERS HAVE BEEN' &
	2       /' OBTAINED FROM DEFORMED TYPE FORM FACTORS')
	  IF(APPR.EQ.'Y') THEN
	    WRITE(9,1080) L_MULT
 1080	    FORMAT(' For this case of L=',I1,' however, the transition rates &
	1 '/' are calculated from an approximate formula (B&M, vol II)')
	  ENDIF
	  IF(BAND.EQ.'GR') THEN
	    WRITE(9,107)
	  ELSE IF(BAND.EQ.'GA') THEN
	    WRITE(9,108)
	  ELSE IF(BAND.EQ.'BE') THEN
	    WRITE(9,109)
	  ELSE IF(BAND.EQ.'OC') THEN
	    WRITE(9,110) KPROJ
	  ENDIF
  107	  FORMAT(//25X,'GROUND STATE BAND')
  108	  FORMAT(//25X,'GAMMA  BAND')
  109	  FORMAT(//25X,'BETA  BAND')
  110	  FORMAT(//25X,'K = ',I1,'  OCTUPOLE BAND')
	  WRITE(9,1090)BETA2_DEF,BETA4_DEF,BETA6_DEF
 1090	  FORMAT(' QUADRUPOLE DEFORMATION :',F8.4,' HEXADECAPOLE &
	1 DEFORMATION :',F8.4,' 6th POLE DEFORMATION :',F8.4)
	  IF(BAND.EQ.'GA') THEN
	    WRITE(9,111) GAMMA
	    IF(L_MULT.EQ.4) WRITE(9,112) BETA4_GAM
	  ELSE IF(BAND.EQ.'BE') THEN
	    WRITE(9,113) BETA0_DEF,BETA0COR
	    IF(L_MULT.EQ.0) WRITE(9,114) BETA0_MON
	  ELSE IF(BAND.EQ.'OC') THEN
	    WRITE(9,115) BETA3_DEF,BETA1COR
	    IF(L_MULT.EQ.1) WRITE(9,116) BETA1_GR
	    IF(L_MULT.EQ.5) WRITE(9,117) BETA5_GR
	  ENDIF
  111	  FORMAT(' GAMMA VIBRATION AMPLITUDE :',F8.4)
  112	  FORMAT(' BETA4 OF 4+ STATE OF GAMMA BAND :',F8.4)
  113	  FORMAT(' BETA-BETA2 :',F8.4,'  CORRECTION TERM COUPLING &
	1 PARAMETER :',F8.4)
  114	  FORMAT(' BREATHING MODE COUPLING PARAMETER :',F8.4)
  115	  FORMAT(' BETA3  :',F8.4,'  C.M. CORRECTION TERM COUPLING &
	1 PARAMETER :',F8.4)
  116	  FORMAT(' ISOSCALAR GIANT DIPOLE COUPLING PARAMETER :',F8.4)
  117     FORMAT(' L=5 SURFACE VIBRATION COUPLING PARAMETER :',F8.4)
	ELSE
	  IF(L_MULT.EQ.0) THEN
	    IF(TRAN_DEN.EQ.'VER1') THEN
	      WRITE(9,118)
  118	      FORMAT('/ Monopole transition rates are calculated &
	& assuming SATCHLER''s VERSION I transition density')
	    ELSEIF(TRAN_DEN.EQ.'VER2') THEN
	      WRITE(9,119)
  119	      FORMAT('/ Monopole transition rates are calculated &
	& assuming SATCHLER''S VERSION II transition density')
	    ENDIF
	  ELSEIF(L_MULT.EQ.1) THEN
	    IF(TRAN_DEN.EQ.'HARA') THEN
	      WRITE(9,120)
  120	      FORMAT('/ Dipole transition rates are calculated &
	& assuming HARAKEH-DIEPERINK transition density')
	    ELSEIF(TRAN_DEN.EQ.'ORLA') THEN
	      WRITE(9,121)
  121	      FORMAT('/ Dipole transition rates are calculated &
	& assuming Orlandini et al. transition density')
	    ENDIF
	  ELSE
	    IF(TRAN_DEN.EQ.'SURF') THEN
	      WRITE(9,122)
  122	      FORMAT('/ Multipole transition rates are calculated &
	& assuming SURFACE vibration transition density')
	    ELSEIF(TRAN_DEN.EQ.'TASS') THEN
	      WRITE(9,123)
  123	      FORMAT('/ Multipole transition rates are calculated &
	& assuming TASSIE type transition density')
	    ENDIF
	  ENDIF
	ENDIF
	WRITE(9,124)L_MULT,L_MULT,L_MULT,L_MULT,L_MULT
  124	FORMAT(//3X,'EX',5X,'BETA'4X,'M(E',I1,')[UNIFM]',4X,'B(E',I1,')',4X,'M(E',I1,')[FERMI]',4X,'B(E',I1,')[FERMI]' &
	&      4X,'M(E',I1,')[BERNS]',5X,'%EWSR',5X,'%EWSR',5X,'%EWSR')
	WRITE(9,125)L_TEMP,L_TEMP,2*L_TEMP,L_TEMP
  125	FORMAT(2X,'(MeV)',11X,' (e-Fm**',I2,')',5x,'(spu)',5X,'(e-Fm**',I2,')',6X,'(e2-Fm**',I2,')',5X,'(e-Fm**',I2,')' &
	&      4X,'(UNIFORM)',2X,'(FERMI)',1X,'(BERNSTEIN)')

!	Go here for monopole transitions.

	IF(L_MULT.EQ.0)THEN

!	The monople matrix element is defined w.r.t. to the operator (1/2)*r**2

	  TEMP0=Z*SQRT(FOURPI)/FOURPI
	  DO 6 I=1,3
	  IF(INX.EQ.1) THEN
	    CALL RFERMI(R(I),AR(I),2,R2(I))
	  ELSEIF(INX.EQ.2) THEN
	    CALL RFERMI(R(I),AR(I),1,R2_1)
	    CALL RFERMI(R(I),AR(I),-1,R_1)
	    R2(I)=2*R2_1/R_1
	  ENDIF
  6	  CONTINUE

!	For the bandhead of the Beta band (deformed transition density).

	  IF(DEFDEN.EQ.'Y')THEN
	    TEMP=100.*EX_DEF
	    DO 7 I=1,3
	    BETA(1)=DEF_LEN2/R(I)
	    BETA(2)=DEF_LEN4/R(I)
	    BETA(3)=DEF_LEN6/R(I)
	    VBETA=DEF_LEN0/R(I)
	    CALL INTEG(BAND,0,NX,R(I),AR(I),ME(I),MECOR(I))
	    ME(I)=Z*ME(I)
	    MECOR(I)=Z*MECOR(I)

!	For breathing mode addition.

	    IF (BETA0_MON.NE.0.) THEN
	      BETA0=DEF_LENM/R(I)
	      ME(I)=ME(I)+TEMP0*BETA0*R2(I)
	    ENDIF

	    BE(I)=ME(I)*ME(I)
  7	    PS(I)=TEMP*BE(I)/SR(I)
	    BSP_UNIF=BE(1)/BSP
	    WRITE(9,126)EX_DEF,ME(1),BSP_UNIF,ME(2),BE(2),ME(3),PS(1),PS(2),PS(3)
	    WRITE(9,128) (MECOR(I),I=1,3)
	    WRITE(6,128) (MECOR(I),I=1,3)
	  ELSE

!	For states analysed using breathing mode form factors.

	    IF(TRAN_DEN.EQ.'VER2') THEN
	      DO I=1,3
		IF(INX.EQ.1) THEN
		  CALL RFERMI(R(I),AR(I),1,R1(I))
		ELSEIF(INX.EQ.2) THEN
		  CALL RFERMI(R(I),AR(I),0,R1_1)
		  CALL RFERMI(R(I),AR(I),-1,R_1)
		  R1(I)=(3/2)*R1_1/R_1
		ENDIF
		XV(I)=(FOURPI/4.)*AR(I)/R(I)
		XV(I)=XV(I)**2
		XV(I)=(3+XV(I))/(1+XV(I))
	      ENDDO
	    ENDIF
	    DO 8 I=1,K
	    EX_TEMP=100.*EX(I)
	    DO 9 J=1,3
	    BETA_TEMP=DEF_LEN(I)/R(J)
	    IF(TRAN_DEN.EQ.'VER1') THEN
	      ME(J)=TEMP0*BETA_TEMP*R2(J)
	    ELSEIF(TRAN_DEN.EQ.'VER2') THEN
	      ME(J)=TEMP0*BETA_TEMP*(4*R(J)*R1(J)-XV(J)*R2(J))/2
	    ENDIF
	    BE(J)=ME(J)*ME(J)
  9	    PS(J)=EX_TEMP*BE(J)/SR(J)
	    BSP_UNIF=BE(1)/BSP
	    WRITE(9,127)EX(I),BETA_FERMI(I),ME(1),BSP_UNIF,ME(2),BE(2),ME(3),PS(1),PS(2),PS(3)
  8	    CONTINUE
	  ENDIF

!	Go here for dipole transitions.

	ELSEIF(L_MULT.EQ.1)THEN

!	The dipole matrix element is defined w.r.t. to the operator
!	(1/2)r**3*Y1.

	  TEMP1=Z/(2.*FOURPI)
	  DO 10 I=1,3
	  IF(INX.EQ.1) THEN
	    CALL RFERMI(R(I),AR(I),2,R2(I))
	    CALL RFERMI(R(I),AR(I),4,R4(I))
	  ELSEIF(INX.EQ.2) THEN
	    CALL RFERMI(R(I),AR(I),1,R2_1)
	    CALL RFERMI(R(I),AR(I),3,R4_1)
	    CALL RFERMI(R(I),AR(I),-1,R_1)
	    R2(I)=2*R2_1/R_1
	    R4(I)=3*R4_1/R_1
	  ENDIF
  10	  CONTINUE

!	For the 1- state of the Octupole band (deformed transition density).

	  IF(DEFDEN.EQ.'Y') THEN
	    TEMP=100.*EX_DEF
	    DO 11 I=1,3
	    BETA(1)=DEF_LEN2/R(I)
	    BETA(2)=DEF_LEN4/R(I)
	    BETA(3)=DEF_LEN6/R(I)
	    VBETA=DEF_LEN3/R(I)
	    CALL INTEG(BAND,1,NX,R(I),AR(I),ME(I),MECOR(I))
	    ME(I)=Z*ME(I)
	    MECOR(I)=Z*MECOR(I)

!	For isoscalar GR mode addition.
        WRITE(*,*)'BETA1_GR:',BETA1_GR
	    IF(BETA1_GR.NE.0.) THEN
	    WRITE(*,*)'IS GR Mode'
	      BETA1=DEF_LEN1/R(I)
	      ME(I)=ME(I)+TEMP1*BETA1*(11.*R4(I)-(25./3.)*R2(I)*R2(I)-10.*EPSILON*R2(I))/R(I)
	    ENDIF
	    BE(I)=ME(I)*ME(I)
  11	    PS(I)=TEMP*BE(I)/SR(I)
	    BSP_UNIF=BE(1)/BSP
	    WRITE(9,126)EX_DEF,ME(1),BSP_UNIF,ME(2),BE(2),ME(3),PS(1),PS(2),PS(3)
	    WRITE(9,129) (MECOR(I),I=1,3)
	    WRITE(6,129) (MECOR(I),I=1,3)
	  ELSE

!	For states analyzed using isoscalar dipole GR form factor.

	    DO 12 I=1,K
	    EX_TEMP=100.*EX(I)
	    DO 13 J=1,3
	    BETA_TEMP=DEF_LEN(I)/R(J)
	    IF(TRAN_DEN.EQ.'HARA') THEN
	      ME(J)=11.*R4(J)-(25./3.)*R2(J)*R2(J)-10.*EPSILON*R2(J)
	      WRITE(*,*)'J,R4(J),R2(J),EPSILON',J,R4(J),R2(J),EPSILON
! 	      WRITE(*,*)'J:',J,'ME(J):',ME(J)
	    ELSEIF(TRAN_DEN.EQ.'ORLA') THEN
	      ME(J)=26.*R4(J)-(100./3.)*R2(J)*R2(J)
	    ENDIF
	    ME(J)=TEMP1*BETA_TEMP*ME(J)/R(J)
	    BE(J)=ME(J)*ME(J)
  13	    PS(J)=EX_TEMP*BE(J)/SR(J)
	    BSP_UNIF=BE(1)/BSP
	    WRITE(9,127)EX(I),BETA_FERMI(I),ME(1),BSP_UNIF,ME(2),BE(2),ME(3),PS(1),PS(2),PS(3)
  12	    CONTINUE
	  ENDIF

!	Go here for L_MULT.GE.2.
!	The matrix element for L>1 is defined w.r.t. the operator   r**L*YL.

	ELSE
	  TEMPL=FLOAT(L_MULT+2)*Z/FOURPI
	  L_1=L_MULT-1
	  X=1./FLOAT(L_MULT)
	  DO 14 I=1,3
	  IF(INX.EQ.1) THEN
	    CALL RFERMI(R(I),AR(I),L_1,RL(I))
	  ELSEIF(INX.EQ.2) THEN
	    CALL RFERMI(R(I),AR(I),L_1-1,RL_2)
	    CALL RFERMI(R(I),AR(I),-1 ,R_1)
	    RL(I)=(L_1+2)*RL_2/(2*R_1)
	  ENDIF
  14	  CONTINUE
	  RUL=(RL(2)*R(2)*FLOAT(L_MULT+2)/3.)**X
	  RSUL=RUL/A13
	  WRITE(6,15) R(1),R(2),R(3)
  15	  FORMAT(1H0,' R(UNIFORM)=',F8.3,', R(FERMI)=',F8.3,', R(BERNSTEIN)=',F8.3)
	  WRITE(6,16)RUL,RSUL
  16	  FORMAT(1HO,' RTR=',F8.4,', RTR/(A**1/3)=',F8.4)

!	For states analysed using deformed transition densities.

	  IF(DEFDEN.EQ.'Y')THEN
	    TEMP=100.*EX_DEF

!	For the 2+ states  analysed using deformed form factor.

	    IF(L_MULT.EQ.2)THEN

!	The following option is available only for the 2+ of the
!	ground state band.

	      IF(BAND.EQ.'GR'.AND.APPR.EQ.'Y')THEN
	        TEMP1=0.08605059D0*Z*DEF_LEN2*DEF_LEN2
	        DO 17 I=1,3
	        ME(I)=TEMPL*RL(I)*DEF_LEN2
	        IF(I.EQ.1) THEN
	          ME(I)=ME(I)+TEMP1
	        ELSE
	          CALL RFERMI(R(I),AR(I),-1,R_I(I))
	          ME(I)=ME(I)+TEMP1*(1.+0.75*(1.-2./3.*R_I(I)*R(I)))
	        ENDIF
  17	        CONTINUE
	      ELSE

!	Here it is again for a general case.

	        DO 18 I=1,3
	        BETA(1)=DEF_LEN2/R(I)
	        BETA(2)=DEF_LEN4/R(I)
	        BETA(3)=DEF_LEN6/R(I)
	        VBETA=DEF_LEN0/R(I)
	        CALL INTEG(BAND,2,NX,R(I),AR(I),ME(I),MECOR(I))
	        ME(I)=Z*ME(I)
  18	        CONTINUE
	      ENDIF

!	For the 3- state of an Octupole band.

	    ELSE IF(L_MULT.EQ.3) THEN
	      DO 19 I=1,3
	      BETA(1)=DEF_LEN2/R(I)
	      BETA(2)=DEF_LEN4/R(I)
	      BETA(3)=DEF_LEN6/R(I)
	      VBETA=DEF_LEN3/R(I)
	      CALL INTEG(BAND,3,NX,R(I),AR(I),ME(I),MECOR(I))
	      ME(I)=Z*ME(I)
  19	      CONTINUE

!	For the 4+ states of the g.s., Gamma and Beta bands.

	    ELSE IF(L_MULT.EQ.4) THEN

!	This option is available only for the 4+ of the Ground state band.

	      IF (BAND.EQ.'GR'.AND.APPR.EQ.'Y') THEN
	        TEMP2=0.28862209D0*Z*DEF_LEN2*DEF_LEN2
	        DO 20 I=1,3
	        ME(I)=TEMPL*RL(I)*DEF_LEN4
	        CALL RFERMI(R(I),AR(I),1,R1(I))
	        CALL RFERMI(R(I),AR(I),2,R2(I))
	        IF(I.EQ.1) THEN
	          ME(I)=ME(I)+TEMP2*R2(I)
	        ELSE
	          ME(I)=ME(I)+TEMP2*R2(I)*(1.-2./3.*(1.-0.8*R1(I)*R(I)/R2(I)))
	        ENDIF
  20	        CONTINUE
	      ELSE

!	Here it is again for the general case.

	        DO 21 I=1,3
	        BETA(1)=DEF_LEN2/R(I)
	        BETA(2)=DEF_LEN4/R(I)
	        BETA(3)=DEF_LEN6/R(I)
	        VBETA=DEF_LEN0/R(I)
	        CALL INTEG(BAND,4,NX,R(I),AR(I),ME(I),MECOR(I))
	        ME(I)=Z*ME(I)
	        IF (BAND.EQ.'GA') ME(I)=ME(I)+TEMPL*RL(I)*DEF_LEN4G
  21	        CONTINUE
	      ENDIF

!	For the 5- state of an Octupole band.

	    ELSE IF(L_MULT.EQ.5) THEN
	      DO 22 I=1,3
	      BETA(1)=DEF_LEN2/R(I)
	      BETA(2)=DEF_LEN4/R(I)
	      BETA(3)=DEF_LEN6/R(I)
	      VBETA=DEF_LEN3/R(I)
	      CALL INTEG(BAND,5,NX,R(I),AR(I),ME(I),MECOR(I))
	      ME(I)=Z*ME(I)
	      IF(BETA5_GR.NE.0.) ME(I)=ME(I)+TEMPL*RL(I)*DEF_LEN5
  22	      CONTINUE

!	For the 6+ state of the g.s., Gamma and Beta bands.

	    ELSE IF(L_MULT.EQ.6) THEN
	      DO 23 I=1,3
	      BETA(1)=DEF_LEN2/R(I)
	      BETA(2)=DEF_LEN4/R(I)
	      BETA(3)=DEF_LEN6/R(I)
	      VBETA=DEF_LEN0/R(I)
	      CALL INTEG(BAND,6,NX,R(I),AR(I),ME(I),MECOR(I))
	      ME(I)=Z*ME(I)
  23	      CONTINUE
	    ENDIF
	    DO 24 I=1,3
	    BE(I)=ME(I)*ME(I)
  24	    PS(I)=TEMP*BE(I)/SR(I)
	    BSP_UNIF=BE(1)/BSP
	    WRITE(9,126)EX_DEF,ME(1),BSP_UNIF,ME(2),BE(2),ME(3),PS(1),PS(2),PS(3)
	    WRITE(6,126)EX_DEF,ME(1),BSP_UNIF,ME(2),BE(2),ME(3),PS(1),PS(2),PS(3)
	    WRITE(6,126)EX_DEF,BE(1),BSP_UNIF,BE(2),BE(2),BE(3),PS(1),PS(2),PS(3)
  126	    FORMAT(/F7.3,8X,E15.6,3X,F6.2,3(1X,E15.6),3X,3(F8.3,2X))
  128	    FORMAT(' VOL. CORR. TERM ',E13.6,12X,E13.6,19X,E13.6)
	  ELSE

!	For states analysed using Surface or Tassie vibration densities.

	    IF(TRAN_DEN.EQ.'TASS') THEN
	      TEMPL=FLOAT(2*L_MULT+1)*Z/FOURPI
	      L_1=2*L_1
	      DO 25 I=1,3
	      IF(INX.EQ.1) THEN
		CALL RFERMI(R(I),AR(I),L_1,RL(I))
	      ELSEIF(INX.EQ.2) THEN
		CALL RFERMI(R(I),AR(I),L_1-1,RL_2)
		CALL RFERMI(R(I),AR(I),-1 ,R_1)
		RL(I)=(L_1+2)*RL_2/(2*R_1)
	      ENDIF
  25	      CONTINUE
	      L_1=L_1/2
	      DO 26 I=1,3
  26	      RL(I)=RL(I)/R(I)**L_1
	    ENDIF
	    DO 27 I=1,K
	    EX_TEMP=100.*EX(I)
	    DO 28 J=1,3
	    ME(J)=TEMPL*RL(J)*DEF_LEN(I)
	    BE(J)=ME(J)*ME(J)
  28	    PS(J)=EX_TEMP*BE(J)/SR(J)
	    BSP_UNIF=BE(1)/BSP
	    WRITE(9,127)EX(I),BETA_FERMI(I),ME(1),BSP_UNIF,ME(2),BE(2),ME(3),PS(1),PS(2),PS(3)
  127	    FORMAT(/F7.3,F7.3,1X,E15.6,3X,F6.2,3(1X,E15.6),3X,3(F8.3,2X))
  129	    FORMAT(' COM CORR. TERM  ',E13.6,12X,E13.6,19X,E13.6)
  27	    CONTINUE
	  ENDIF
	ENDIF
	DO 29 I=1,3
  29	BETA(I) = 0.D0
	VBETA=0.D0
	BETACOR=0.D0
	XK=0.D0
	IF(BETA0_COR.EQ.1.0) BETA0_COR=0.D0
	IF(BETA3_COR.EQ.1.0) BETA3_COR=0.D0
	GOTO 100
	END
