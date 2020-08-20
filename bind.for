	PROGRAM MAIN
	IMPLICIT REAL*8(A-H,O-Z)

C MNH	This program was written to calculate width of single particle
C MNH	resonances in a real Woods-Saxon potential including Coulomb and
C MNH	angular momentum barriers. It also calculates the radial matrix
C MNH	element for r**L for a transition between two single particle states.
C MNH	In this case the first particle state can be a resonant wavefunction.
C MNH	Some subroutines from the extended version of program DWUCK4 have
C MNH	been slightly modified to be used here.

C MNH	Written by: M.N. Harakeh, NPL, University of Washington, December'83.
*****************************************************************************

	PARAMETER (LR = 3200)
	COMMON ALPHA(15),BETA(2),ICON(20),ANGLE(6),HBARC,AMU,AMASS,CHGSQ
	1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
	2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
	3,BYTES, K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,IFF,LTRT(8),JTRT(8)
	4,ISTRT(8),JSAV(8),IBF(8),GDR
	COMMON /LOCAL/  U(2*LR,2),V(2*LR,2),SPACE(8*LR),UB(2*LR)
	COMMON /LDEPEN/ ULD(2*LR,2),ILDSAV(8),ILD,IDSOSV(8),IDSO
	COMMON /ISOSPN/ TWBF(2*LR),TT,TR,CANAL,CCORE
	COMMON /UBOUND/ DY,YMAX,EISIG0(2)
	COMMON /BIND2/ E,FZ,VRR,RR,AR,VSOR,GAM,VS,RS,AS,FSPERM,FISW,VD
	COMMON /CHAR/ CON,SEARCH,NOWF,CHOICE
	CHARACTER CON*2,SEARCH*1,NOWF*1,CHOICE*1,WRITE*1

	FSCONT = 137.03602
	HBARC = 197.329
	AMU = 931.502
	CHGSQ = HBARC/FSCONT
	AC(3)=0.65
	STEP =0.1
	DR(3)=0.1
	RMAX =20.
	RZ   =1.25
	RY   =1.25
	AR   =0.65
	VSOR =25.
	PI = 4*ATAN(1.)
1	CALL CIOA  (' PR[oceed] or EX[it] : ',CON)
	IF(CON.EQ.'EX') STOP
! 	Output files for resonance phase and
! 	1 radial matrix element as a function of     excitation energy
! 	2 are written to BINDSP.DAO and BINDME.DAO, respectively'
	CALL CIOWRT(' Are you interested in calculating radial moment')
	CALL CIOA  (' for a transition Y[es] or N[o]: ',NOWF)
	IF(NOWF.EQ.'Y') THEN
	  CALL CIOI  (' Multipole order of radial moment: ',L_MULT)
	  NPART=2
	ELSEIF(NOWF.EQ.'N') THEN
	  NPART=1
	ENDIF
	CALL CIOF  (' Integration step in fm : ',STEP)
	CALL CIOF  (' Max integration radius in fm : ',RMAX)
	IF(STEP.NE.0) DR(3)=STEP
	K=NINT(RMAX/DR(3))

	IF(NPART.EQ.2) THEN
! 	  CALL CIOWRT(' The final wavefunction has to be a bound wavefunction')
      WRITE(*,*)' The final wavefunction has to be a bound wavefunction'
! 	  CALL CIOWRT(' -----------------------------------------------------')
      WRITE(*,*)' -----------------------------------------------------'
	ENDIF
	DO 2 I=NPART,1,-1
	IF(I.EQ.1) THEN
	  CALL CIOWRT('Read parameters for initial wavefunction')
	ELSEIF(I.EQ.2) THEN
	  CALL CIOWRT('Read parameters for final wavefunction')
	ENDIF
	FZ =-1.0
	VRR =-1.0
	CALL CIOWRT(' Vary resonance/binding  E[nergy]')
	CALL CIOA (' or potential D[epth] : ',SEARCH)
3	IF(SEARCH.EQ.'E') THEN
	  FISW=1.0
	  CALL CIOF  (' Depth of potential in MeV (+ve): ',VD)
	  VD=ABS(VD)
	  CALL CIOWRT(' Appr. resonance (+ve)/binding (-ve) energy')
	  CALL CIOF  (' in C.O.M. in MeV : ',E)
	ELSEIF(SEARCH.EQ.'D') THEN
	  FISW=0.
	  CALL CIOWRT(' Resonance (+ve)/binding (-ve) energy in C.O.M.')
	  CALL CIOF  (' in MeV : ',E)
	  E_TEMP=E
	  FISW_TEMP=FISW
	ENDIF

	IF(I.EQ.2.AND.E.GT.0.) THEN
	  CALL CIOWRT(' FINAL STATE MUST BE BOUND; Try again!')
	  GO TO 3
	ENDIF

	CALL CIOF (' Particle mass in amu              : ',FM(3))
	CALL CIOF (' Particle charge                   : ',Z(3))
	CALL CIOF (' Target mass in amu                : ',FMA(3))
	AFACT = FMA(3)**.333333333
	CALL CIOF (' Target charge                     : ',ZA(3))
	CALL CIOF (' Charge radius in A**1/3           : ',RZ)
	RC(3) = RZ*AFACT
CMN	CALL CIOF (' Charge diffuseness in fm          : ',AC(3))
	CALL CIOF (' Mass radius in A**1/3 fm          : ',RY)
	RR = RY*AFACT
	CALL CIOF (' Mass diffuseness in fm            : ',AR)
	CALL CIOWRT('Thomas spin-orbit coupling term if entered negative')
	CALL CIOWRT('more input for a fixed S-O potential will be asked ')
	CALL CIOF (' for; VSOR                          : ',VSOR)
	IF (VSOR.LT.0.) THEN
	  CALL CIOF (' S-O potential in DWUCK form       : ',VS)
	  CALL CIOF (' S-O pot radius in a**1/3 fm       : ',RS)
	  CALL CIOF (' S-O pot diffuseness in fm         : ',AS)
	ENDIF
	CALL CIOWRT(
	1'Note - nr of nodes excludes nodes at zero and infinity')
	CALL CIOF (' Nr of nodes of wavefunction       : ',FN)
	CALL CIOF (' Orbital angular momentum          : ',FL)
	CALL CIOF (' 2*Total angular momentum          : ',FJ2)
	CALL CIOF (' 2*Spin of particle                : ',FSPERM)
	CALL FORMF(U,V(1,I),3,0)
2	CONTINUE
	IF(NPART.EQ.1) THEN
	  CALL CIOWRT('Do you want to write out the wavefunction and its')
	  CALL CIOA  (' derivative Y[es] or N[o]: ',WRITE)
	  IF(WRITE.EQ.'Y') THEN
	    R=0.
	    DO M=2,K-1
	      R=R+DR(3)
	      I=LR+M
	      DER=(V(I+1,1)-V(I,1))/DR(3)
	      WRITE(8,12) R, V(I,1), DER
12	      FORMAT(1X,F6.2,2(3X,E12.4))
	    ENDDO
	    WRITE='N'
	  ENDIF  
	  IF(E.GT.0.)THEN
	    CALL CIOWRT(' Do you want the phase in fine energy steps')
	    CALL CIOA  (' over the resonance Y[es] or N[o]: ',CHOICE)
	    IF(CHOICE.EQ.'N') GO TO 1
	    OPEN (UNIT=8,STATUS='NEW',FILE='BINDSP.DAO')
	    WRITE(8,4)
4	    FORMAT(' Resonance energy in MeV',4X,'  DELTA in Deg.')
	    FISW=2.
	    IF(GAM.GT.0.) THEN
	      ESTEP=GAM/7
	      E=E-3*GAM
	      IF(E.LT.0.) E=0.
	      NTEMP=43
	    ELSEIF(GAM.LT.0.) THEN
	      ESTEP=E/20
	      E=E/2
	      NTEMP=21
	    ENDIF
	    EMIN=E
	    EMAX=E+(NTEMP-1)*ESTEP
	    ELOW=EMIN
	    EHIGH=EMAX
	    E_STEP=ESTEP
	    CALL CIOF  (' Lowest energy  ; ELOW  : ',ElOW)
	    CALL CIOF  (' Highest energy ; EHIGH : ',EHIGH)
	    CALL CIOF  (' Energy step    ; E_STEP: ',E_STEP)
	    IF(ELOW.NE.EMIN.OR.EHIGH.NE.EMAX.OR.E_STEP.NE.ESTEP) THEN
	      EMIN=ELOW
	      EMAX=EHIGH
	      ESTEP=E_STEP
	      NTEMP=NINT((EMAX-EMIN)/ESTEP)+1
	      E=EMIN
	    ENDIF
	    DO 5 KK=1,NTEMP
	    CALL FORMF(U,V(1,1),3,0)
	    E=E+ESTEP
5	    CONTINUE
	    FISW=FISW_TEMP
	    E=E_TEMP
	    CLOSE (UNIT=8)
	  ENDIF
	  GO TO 1
	ELSEIF(NPART.EQ.2) THEN
	  IF(E.GT.0.)THEN
	    CALL CIOWRT(' Do you want the m.e. in fine energy steps')
	    CALL CIOA  (' over the resonance Y[es] or N[o]: ',CHOICE)
	    IF(CHOICE.EQ.'N') GO TO 1
	    OPEN (UNIT=8,STATUS='NEW',FILE='BINDSP.DAO')
	    WRITE(8,4)
	    OPEN (UNIT=9,STATUS='NEW',FILE='BINDME.DAO')
	    WRITE(9,6) L_MULT
6	    FORMAT(' Resonance energy in MeV',4X,'  M.E. in fm**',I1)
	    FISW=2.
	    IF(GAM.GT.0.) THEN
	      ESTEP=GAM/7
	      E=E-3*GAM
	      IF(E.LT.0.) E=0.
	      NTEMP=43
	    ELSEIF(GAM.LT.0.) THEN
	      ESTEP=E/20
	      E=E/2
	      NTEMP=21
	    ENDIF
	    EMIN=E
	    EMAX=E+(NTEMP-1)*ESTEP
	    ELOW=EMIN
	    EHIGH=EMAX
	    E_STEP=ESTEP
	    CALL CIOF  (' Lowest energy  ; ELOW  : ',ElOW)
	    CALL CIOF  (' Highest energy ; EHIGH : ',EHIGH)
	    CALL CIOF  (' Energy step    ; E_STEP: ',E_STEP)
	    IF(ELOW.NE.EMIN.OR.EHIGH.NE.EMAX.OR.E_STEP.NE.ESTEP) THEN
	      EMIN=ELOW
	      EMAX=EHIGH
	      ESTEP=E_STEP
	      NTEMP=NINT((EMAX-EMIN)/ESTEP)+1
	      E=EMIN
	    ENDIF
	    DO 7 KK=1,NTEMP
	    CALL FORMF(U,V(1,1),3,0)
	    SUM=0.
	    R=0.
	    DO 8 M=1,K
	    R=R+DR(3)
	    I=LR+M
8	    SUM=SUM+V(I,1)*V(I,2)*R**(L_MULT+2)
	    SUM=SUM*DR(3)
	    WRITE(9,9) E,SUM
9	    FORMAT(2(9X,1PE15.7))
	    WRITE(6,10) L_MULT,SUM,L_MULT
10	    FORMAT(/' Radial Matrix Element <2|R**',I1,'|1> = ',F12.5,
	1   '  fm**',I1/)
	    E=E+ESTEP
7	    CONTINUE
	    FISW=FISW_TEMP
	    E=E_TEMP
	    CLOSE (UNIT=8)
	    CLOSE (UNIT=9)
	    GO TO  1
	  ENDIF
	  SUM=0.
	  R=0.
	  DO 11 M=1,K
	  R=R+DR(3)
	  I=LR+M
11	  SUM=SUM+V(I,1)*V(I,2)*R**(L_MULT+2)
	  SUM=SUM*DR(3)
	  WRITE(6,10) L_MULT,SUM,L_MULT
	  GO TO  1
	ENDIF
	END

      SUBROUTINE FORMF(U,V,N,LAM)
      IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER LR = 3200
      COMMON ALPHA(15),BETA(2),ICON(20),ANGLE(6),HBARC,AMU,AMASS,CHGSQ
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,BYTES, K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,IFF,LTRT(8),JTRT(8)
     4,ISTRT(8),JSAV(8),IBF(8),GDR
      COMMON /LDEPEN/ ULD(2*LR,2),ILDSAV(8),ILD,IDSOSV(8),IDSO
      COMMON /ISOSPN/ TWBF(2*LR),TT,TR,CANAL,CCORE
      COMMON /BIND2/ E,FZ,VRR,RR,AR,VSOR,GAM,VS,RS,AS,FSPERM,FISW,VD
	COMMON /CHAR/ CON,SEARCH,NOWF,CHOICE
	CHARACTER CON*2,SEARCH*1,NOWF*1,CHOICE*1
      DIMENSION U(2*LR),V(2*LR),G(4)
      EQUIVALENCE (G(1),FN)
      DATA ETA6,PI/60.D0,3.141593/

C     READ IN CARD SET 5,6,OR 7   ENERGY CARD
C
CMN   READ (5,9000)E,FM(N),Z(N),FMA(N),ZA(N),RY,AC(N),PNLOC(N),FS(N),QCD
      IF(N.NE.2) GO TO 6
      Q=E
      IF(QCD.GT.0.0) GO TO 6
      IF(QCD.LT.0.0) Q=Q+QCD
      E=(FM(2)+FMA(2))*(ECM(1)+Q)/FMA(2)
    6 CONTINUE
      IF(N.NE.3) GO TO 8
      E=E+QCD
    8 CONTINUE
      IS(N)=FS(N)
      NS(N)=IS(N)+1
CMN   IF(AMASS.EQ.0.0) AMASS=FMA(1)
      IF(AMASS.EQ.0.0) AMASS=FMA(3)
      AFACT=FMA(N)**.333333333
      BFACT=FM (N)**.333333333
CMN   RC(N)=ABS(RY)*AFACT
CMN   IF(RY.LT.0.0) RC(N)=RC(N)+ABS(RY)*BFACT
CMN   DR(N)=DRF*AMASS/FMA(N)
      RY=RC(N)/AFACT
      DO 12 M=1,2*LR
      U(M)=0.D0
      V(M)=0.D0
      IF(N.GE.3) GO TO 12
      ULD(M,N)=0.D0
   12 CONTINUE
      IF(E.EQ.0.0) GO TO 66
      IF(ICON(10).EQ.0) GO TO 30
C
C     ICON(10).NE.0 GIVES RELATIVISTIC KINEMATICS
C
      FM1=FM(N)*AMU
      FM2=FMA(N)*AMU
      FMT=FM1+FM2
      IF(N.NE.2) GO TO 26
      IF(QCD.GT.0.0) GO TO 26
      E=E+(ECM(1)+Q)**2/(2.0*FM2)
   26 CONTINUE
      IF(N.EQ.3) E=E*(0.5*E+FMT)/FM2
      T1=SQRT(2.0*E*FM2+FMT*FMT)
      W1=(FMT*FM1+FM2*E)/T1
      W2=FM2*(E+FMT)/T1
      FMU(N)=W1*W2/(W1+W2)
      ECM(N)=T1-FMT
      FACT=2.0*FMU(N)/(HBARC*HBARC)
      FMU(N)=FMU(N)/AMU
      FK2(N)=(E*E+2.0*E*FM1)*(FM2/T1)**2/(HBARC*HBARC)
      GO TO 36
   30 CONTINUE
      FMU(N)=FM(N)*FMA(N)/(FM(N)+FMA(N))
      ECM(N)=FMU(N)*E/FM(N)
      IF(N.EQ.3) ECM(N)=E
      FACT=2.0*FMU(N)*AMU/(HBARC*HBARC)
      FK2(N)=FACT*ECM(N)
   36 CONTINUE
      FK(N)=SQRT(ABS(FK2(N)))
      ETAK=CHGSQ*Z(N)*ZA(N)*FACT
      ETA(N)=ETAK*0.5/FK(N)
C
C     ADD COULOMB AND KINETIC ENERGIES TO U
C
      RCX=RC(N)
      IF(RCX.EQ.0.0) RCX=DR(N)
      F1=0.5*ETAK/RCX
      RC2=RCX*RCX
      R=0.D0
      DO 50 M=1,LR
      R=R+DR(N)
      MK=M+M-1
      IF(R.GE.RCX) GO TO 40
      F2=F1*(3.0-R*R/RC2)
      GO TO 43
   40 CONTINUE
      F2=ETAK/R
   43 CONTINUE
      IF(N.NE.3) GO TO 46
      U(MK+1)=FK2(N)-F2
      GO TO 50
   46 CONTINUE
      U(MK  )=FK2(N)-F2
   50 CONTINUE
      GO TO 67
   66 CONTINUE
      FK(N)=0.D0
      ETA(N)=0.D0
      ECM(N)=0.D0
   67 CONTINUE
CMN   WRITE(6,9010)N
	IF(E.GE.0.AND.CHOICE.EQ.'Y') GO TO 69
      WRITE(6,9503)E,RY,AC(N),FS(N)
      WRITE(6,9504)FM(N),FMA(N)
      WRITE(6,9505)Z(N),ZA(N),PNLOC(N)
      WRITE(6,9500)
      RHO=FK(N)*RC(N)
      WRITE(6,9506)ECM(N),RC(N),RHO
      WRITE(6,9507)FK(N),ETA(N),DR(N)
      WRITE(6,9008)
   69 CONTINUE
      FS(N)=FS(N)/2.0
      IF(N.NE.3) GO TO 80
      IF(ICON(19).LT.4) GO TO 80
      DO 75 M=1,2*LR
      U(M)=0.D0
      TWBF(M)=0.D0
   75 CONTINUE
C     READ(9) K,CANAL,CCORE,(TWBF(2*M-1),M=1,K),(U(2*M-1),M=1,K)
C     WRITE(6,9600) K,CANAL,CCORE
      GO TO 3000
   80 CONTINUE
      IBF(1)=0
      CALL POTEN(U,V,AFACT,BFACT,FACT,RM,N,LAM)
      IF(N.NE.3) GO TO 3000
      IF(E.EQ.0.0) GO TO 3000
      IF(IBF(1).NE.0) GO TO 3000
      IBF(4)=0
C
C     SINGL PARTICLE ORBITAL
C
C
C     READ IN QUANTUM NUMBERS FOR SINGLE PARTICLE ORBITAL
C
C     READ(5,9000)G,VTRIAL,FISW,DAMP
      ISW=FISW
	IF(E.GE.0.AND.CHOICE.EQ.'Y') GO TO 81
	  WRITE(6,9500)
81	CONTINUE
      FJ=FJ2/2.0
      	IF(ISW.EQ.0) THEN
	  VTRIAL=ETA6
	ELSEIF(ISW.EQ.1) THEN
	  VTRIAL=VD
	ENDIF
	FSS=FSPERM
      	IF(E.GE.0.AND.CHOICE.EQ.'Y') GO TO 82
	WRITE(6,9508)G,FISW,VTRIAL
	IF(DAMP.NE.0.0) WRITE(6,9511) DAMP
82	CONTINUE
      FSS=FSS/2.0
      FACT=(FJ*FJ+FJ-FL*FL-FL-FSS*FSS-FSS)*0.5
      DO 2028 M=1,LR
      MK=M+M-1
      U(MK+1)=U(MK+1)+V(MK+1)*FACT
      V(M)=V(MK  )*FACT
      U(MK  )=U(MK  )+V(M)
 2028 CONTINUE
	IF(E.GE.0.AND.CHOICE.EQ.'Y') GO TO 83
	  WRITE(6,9500)
83	CONTINUE
      CALL BIND(U,  V(LR+1),DR(3),RM,FN,FL,K,FK(3),ETA(3),VTRIAL,ECM(3)
     1,FK2(3)    ,ISW,IBF(3),GAM)
      IBF(2)=RM/DR(3)
      FACT=PNLOC(3)**2/4.0
C
C     NON-LOCAL CORRECTION FOR SINGLE PARTICLE FUNCTION
C
      SUM=0.D0
      R=0.D0
      DO 2075 M=1,K
      MK=M+M-1
      R=R+DR(3)
      U(MK)=U(MK)*VTRIAL+U(MK+1)
      TEMP=FACT*(U(MK)-FK2(3))
      IF(PNLOC(3)) 2072,2075,2070
 2070 CONTINUE
      V(M+LR)=V(M+LR)/SQRT(1.0+TEMP)
      GO TO 2074
 2072 CONTINUE
      V(M+LR)=V(M+LR)*EXP(-TEMP/2.0)
 2074 CONTINUE
      SUM=SUM+(V(M+LR)*R)**2
 2075 CONTINUE
      IF(FACT.EQ.0.0) SUM=1.D0
      IF(FACT.EQ.0.0) GO TO 2076
      SUM=1.0/SQRT(SUM*DR(3))
 2076 CONTINUE
      IF(DAMP.EQ.0.0) GO TO 2090
C
C     APPLY DAMPING FACTOR, EXP(-DAMP*R*R) TO FORM FACTOR
C
      R=0.D0
      H=ABS(DAMP)
      DO 2085 M=1,K
      R=R+DR(3)
      F2=H*R*R
      F2=-MIN(F2,30.0D0)
      F1=EXP(F2)
      V(M+LR)=V(M+LR)*F1
 2085 CONTINUE
 2090 CONTINUE
      DO 2100 M=1,LR
      MK=M+M-1
      V(M)=U(MK)
      IF(M.LE.K) GO TO 2095
      V(M)=0.D0
      V(M+LR)=0.D0
 2095 CONTINUE
      U(MK  )=V(M+LR)*SUM
      U(MK+1)=0.D0
 2100 CONTINUE
 3000 CONTINUE
      RETURN
 9000 FORMAT(10F8.4)
 9008 FORMAT(21H0POTENTIAL PARAMETERS)
 9010 FORMAT( 9H0PARTICLE,I2,120(1H*))
 9500 FORMAT(1H )
 9503 FORMAT(18H INPUT DATA       ,9H   ECM  =,F9.4,9H   RC0  =,F9.4,
     19H   AC   =,F9.4,9H   2*STR=,F9.4)
 9504 FORMAT(18X,9H   MASSP=,F9.4,9H   MASST=,F9.4)
 9505 FORMAT(18X,9H   ZP   =,F9.4,9H   ZT   =,F9.4,9H   PNLOC=,F9.4)
 9506 FORMAT(18H DERIVED DATA     ,9H   ECM  =,F9.4,9H   RC   =,F9.4,
     19H   RHO  =,F9.4)
 9507 FORMAT(18X,9H   K    =,F9.4,9H   ETA  =,F9.4,9H   DR   =,F9.4)
 9508 FORMAT(18X,9H   NODES=,F9.4,9H   L    =,F9.4,9H   2*J  =,F9.4,
     19H   2*S  =,F9.4,9H   FISW =,F9.4,9H   VTRL =,F9.4)
 9511 FORMAT(18X,9H   DAMP = ,F9.4)
 9600 FORMAT(1H0,10X,4HK = ,I3,5X,12HCANALOGUE = ,F7.5,5X,8HCCORE = ,F7.
     15)
      END

      SUBROUTINE POTEN(U,V,AFACT,BFACT,FACT,RM,N,LAM)
      IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER LR = 3200
      COMMON ALPHA(15),BETA(2),ICON(20),ANGLE(6),HBARC,AMU,AMASS,CHGSQ
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,BYTES, K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,IFF,LTRT(8),JTRT(8)
     4,ISTRT(8),JSAV(8),IBF(8),GDR
      COMMON /LDEPEN/ ULD(2*LR,2),ILDSAV(8),ILD,IDSOSV(8),IDSO
      COMMON /BIND2/ E,FZ,VRR,RR,AR,VSOR,GAM,VS,RS,AS,FSPERM,FISW,VD
	COMMON /CHAR/ CON,SEARCH,NOWF,CHOICE
	CHARACTER CON*2,SEARCH*1,NOWF*1,CHOICE*1
      DIMENSION U(2*LR),V(2*LR),UT(5),G(4),ABETA(3),FLDF(3),LDFRM(3)
      DIMENSION CN(16),YLAM(16),X(8),W(8),CP(16),UB(2*LR),VB(2*LR)
      REAL*8 B(3,19)
      EQUIVALENCE (G(1),FN)
      DATA B(1,1),B(2,1),B(3,1)/6H NX=1 ,6HVOLUME,6H W-S  /
      DATA B(1,2),B(2,2),B(3,2)/6H NX=2 ,6HSURFAC,6HE W-S /
      DATA B(1,3),B(2,3),B(3,3)/6H NX=3 ,6H2ND DE,6HRIV   /
      DATA B(1,4),B(2,4),B(3,4)/6H NX=4 ,6HL.S VO,6HLUME  /
      DATA B(1,5),B(2,5),B(3,5)/6H NX=5 ,6HL.S SU,6HRFACE /
      DATA B(1,6),B(2,6),B(3,6)/6H NX= O,6HUT OF ,6HRANGE /
      DATA F4PI / 12.56637061 /

C     POINTS FOR GAUSS LEGENDRE INTEGRATION
C
      DATA X/0.0950125098,0.2816035507,0.4580167776,0.6178762444,
     1       0.7554044083,0.8656312023,0.9445750230,0.9894009439/
C
C     WEIGHTS FOR GAUSS LEGENDRE INTEGRATION
C
      DATA W/0.1894506104,0.1826034150,0.1691565193,0.1495959888,
     1       0.1246289712,0.0951585116,0.0622535239,0.0271524594/
      DATA ETA4,ETA5/4.D0,10.D0/
      DATA SQRPI/1.772453851/
C     THE NEXT 1 STATEMENT FAKES OUT A BUG IN THE IBM 360 COMPILER,
C       FORTRAN IV, LEVEL H, OPTIMIZER=2.  N IS ALWAYS LESS THAN 4.
      IF (N.GT.10) GO TO 83
      RM=0.D0
      ILD=0
      ITM=0
	LFLAG=0
	VR=VRR
	VI=0.D0
C
C     READ IN CARD SET 5,6,OR 7   POTENTIAL CARDS
C
   70 CONTINUE
      ITM=ITM+1
CMN   READ (5,9000)FZ,VR,RY,AR,VSOR,VI,RZ,AI,VSOI,PWR
      NX=ABS(FZ)
      IF(NX.LE.0) NX=6
      IF(NX.GT.6) NX=6
CMN   RR=ABS(RY)*AFACT
      RY=RR/AFACT
      RI=ABS(RZ)*AFACT
CMN   IF(RY.LT.0.0) RR=RR+ABS(RY)*BFACT
CMN   IF(RY.LT.0.0) RI=RI+ABS(RZ)*BFACT
	IF(E.GE.0.AND.CHOICE.EQ.'Y') GO TO 75
	WRITE(6,9509)(B(J,NX),J=1,3),VR,RY,AR,RR,VSOR
	WRITE(6,9510)                VI,RZ,AI,RI,VSOI,PWR
   75 CONTINUE
CMN   RY=ABS(RY)
      RY=ABS(RR)
      RZ=ABS(RZ)
C      NX=1 VOLUME WOODS-SAXON
C      NX=2 SURFACE WOODS-SAXON
C      NX=3 SECOND DERIVATIVE WOODS-SAXON
C      NX=4 L.S POTENTIAL FOR WOODS-SAXON VOLUME
C      NX=5 L.S POTENTIAL FOR WOODS-SAXON SURFACE
      KFLAG=0
      IF(N-3)79,76,79
   76 CONTINUE
      IF(E)81,80,81
   79 CONTINUE
      VR=VR*FACT
      VI=VI*FACT
      KT=FK(N)*MAX(RR,RI)+ETA5
      LPLUS=MAX0(LPLUS,KT)
   80 KT=(2.3*ETA4*MAX(AR,AI)+MAX(RR,RI))/DR(N)
      GO TO 82
   81 CONTINUE
      RM=MAX(RR,RM)
      VR=VR*FACT
      VI=VI*FACT
      IF(E.GT.0.0) GO TO 80
      ALPHB=2.3*ETA4
      TEMP=1.0+ALPHB/10.0
      KT=((ALPHB-ETA(N)*TEMP)/FK(N)+RR)/DR(N)
   82 K=MIN0(MAX0(K,KT),LR)
   83 CONTINUE
      IF(AR.NE.0.0) GO TO 85
      F1=0.D0
      F2=0.D0
      GO TO 86
   85 F2=EXP(-DR(N)/AR)
      F1=EXP( RR/AR)
   86 CONTINUE
      IF(AI.NE.0.0) GO TO 95
      F3=0.D0
      F4=0.D0
      GO TO 96
   95 F4=EXP(-DR(N)/AI)
      F3=EXP( RI/AI)
   96 CONTINUE
      IF(NX.GT.5) GO TO 2009
      GO TO (100,200,300,400,500),NX
C
C     VOLUME WOODS SAXON
C
  100 CONTINUE
      DO 160 M=1,K
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      U(MK  )=U(MK  )-VR*F1/(1.0+F1)
      U(MK+1)=U(MK+1)-VI*F3/(1.0+F3)
  160 CONTINUE
      GO TO 2000
C
C     1ST DER WOODS SAXON
C
  200 CONTINUE
      DO 260 M=1,K
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      U(MK  )=U(MK  )+VR*F1/(1.0+F1)**2
      U(MK+1)=U(MK+1)+VI*F3/(1.0+F3)**2
  260 CONTINUE
      GO TO 2000
C
C     2ND DER WOODS SAXON
C
  300 CONTINUE
      DO 360 M=1,K
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      U(MK  )=U(MK  )-VR*F1*(1.0-F1)/(1.0+F1)**3
      U(MK+1)=U(MK+1)-VI*F3*(1.0-F3)/(1.0+F3)**3
  360 CONTINUE
      GO TO 2000
C
C     L.S VOLUME WOODS SAXON
C
  400 CONTINUE
      IF(AR.NE.0.0) VR=VR/AR
      IF(AI.NE.0.0) VI=VI/AI
      IF((VR.NE.0.0).OR.(VI.NE.0.0)) IBF(4)=1
      R=0.D0
      DO 460 M=1,K
      R=R+DR(N)
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      V(MK  )=V(MK  )-VR*F1/(R*(1.0+F1)**2)
      V(MK+1)=V(MK+1)-VI*F3/(R*(1.0+F3)**2)
  460 CONTINUE
      GO TO 2000
C
C     L.S 1ST DER WOODS SAXON
C
  500 CONTINUE
      IF(AR.NE.0.0) VR=VR/AR
      IF(AI.NE.0.0) VI=VI/AI
      IF((VR.NE.0.0).OR.(VI.NE.0.0)) IBF(4)=1
      SPO=0.D0
      IF(VSOR.GT.0.0) SPO=1.D0
      R=0.D0
      DO 560 M=1,K
      R=R+DR(N)
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      TEMP=1.0+F1
      V(MK  )=V(MK  )+VR*F1*((1.0-F1)/TEMP+SPO*AR/R)/(R*TEMP*TEMP)
      TEMP=1.0+F3
      V(MK+1)=V(MK+1)+VI*F3*((1.0-F3)/TEMP+SPO*AI/R)/(R*TEMP*TEMP)
  560 CONTINUE
 2000 CONTINUE
      IF(KFLAG.NE.0) GO TO 2009
	IF (VSOR.LT.0) THEN
	  IF(LFLAG.NE.0) GO TO 2009
	  FZ=-4
	  VI=VS
	  VR=0.D0
	  RZ=RS
	  AI=AS
	  LFLAG=1
	  GO TO 70
	ENDIF
      IF(ABS(VSOR)+ABS(VSOI).EQ.0.0) GO TO 2009
      NX=NX+3
      IF(NX.GT.5) GO TO 2009
      KFLAG=1
      VR=VR*VSOR/45.2
      VI=VI*VSOI/45.2
      GO TO 83
 2009 CONTINUE
      IF(FZ.GT.0.0) GO TO 70
      RETURN
 9000 FORMAT(10F8.4)
 9001 FORMAT(18H0PARAMETERS       ,9H   QCODE=,F9.4,9H   RANGE=,F9.4,
     17H   V0 =,F11.4)
 9100 FORMAT(5E16.7)
 9508 FORMAT(18X,9H   NODES=,F9.4,9H   L    =,F9.4,9H   2*J  =,F9.4, 
     19H   2*S  =,F9.4,9H   FISW =,F9.4)
 9509 FORMAT(3A6,9H   V RL =,F9.4,9H   R0RL =,F9.4,9H   A RL =,F9.4, 
     19H   R RL =,F9.4,9H   VSOR =,F9.4)
 9510 FORMAT(18X,9H   V IM =,F9.4,9H   R0IM =,F9.4,9H   A IM =,F9.4, 
     19H   R IM =,F9.4,9H   VSOI =,F9.4,9H   POWR =,F9.4)
 9512 FORMAT(18X,9H   BETA1=,F9.4,9H   LDFR1=,F9.4,9H   BETA2=,F9.4, 
     19H   LDFR2=,F9.4,9H   BETA3=,F9.4,9H   LDFR3=,F9.4)
 9513 FORMAT(18X,9H   F1   =,F9.4,9H   F2   =,F9.4,9H   F3   =,F9.4)
      END

      SUBROUTINE BIND(U,F,DR,RM,FNODE,FL,K,FK,ETA,V,E,FK2,ISW,IERR,GAM)
      IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER LR = 3200
	COMMON /CHAR/ CON,SEARCH,NOWF,CHOICE
	CHARACTER CON*2,SEARCH*1,NOWF*1,CHOICE*1
      DIMENSION U(2*LR),F(LR)
      DATA ICMAX/16/
      DATA EPS,FLAMX,DETX/1.0D-06,1.0D-04,1.0D-03/
      DATA PI,UMAX,UMAXSQ/3.141593,1.0D+16,1.0D+32/

      ICNT=0
      FLP=FL*(FL+1.0)
      DR2=DR*DR/12.0
      MT=K+K-1
      LL=FL+1.D0
C
C     CALCULATE OUTER BOUNDARY CONDITION
C
   10 CONTINUE
      ICNT=ICNT+1
      RNORM=0.D0
      R=DR*FLOAT(K)
      IF(FK2.LT.0.0) GO TO 23
      M=MAX0(LL+1,3)
      T3=FK*(R-DR)
      T4=FK*R
      R1=T3*1.01
      R2=T4*1.01
      ETAP=ETA/1.01
      FLL=FL+1.0
      CALL COU(R1,R2,ETAP,M,DR,F( 1),F(21),F(41),F(61),F(81))
      SURF=1.01*FK*(SQRT(FLL*FLL+ETAP*ETAP)*F(LL+61)/F(LL+60))/FLL
      R1=0.99*T3
      R2=0.99*T4
      ETAP=ETA/0.99
      CALL COU(R1,R2,ETAP,M,DR,F( 1),F(21),F(41),F(61),F(81))
      SURF=0.99*FK*(SQRT(FLL*FLL+ETAP*ETAP)*F(LL+61)/F(LL+60))/FLL-SURF
      R1=T3
      R2=T4
      CALL COU(R1,R2,ETA ,M,DR,F( 1),F(21),F(41),F(61),F(81))
      T6=F(LL+60)
      T5=F(LL+40)
      T4=F(LL+20)
      T3=F(LL   )
      SURF=SURF*25.0*T6*T6/FK2
      IF(ISW.NE.2) GO TO 40
      KM=K
      GO TO 120
   23 CONTINUE
      T3=FLP/(R-DR)**2-U(MT-1)-V*U(MT-2)
      T4=FLP/ R    **2-U(MT+1)-V*U(MT  )
      T3=SQRT(T3)
      T4=SQRT(T4)
      T5=DR*(T3+T4)/2.0
      T5=SQRT(T4/T3)*EXP(T5)
      T5=T5*EPS
      T6=EPS
   40 FNCT=0.D0
      KT=K-2
C
C     INTEGRATE FROM INFINITY TO CLASSICAL TURNING POINT
C
      IFLAG=0
   41 CONTINUE
      R=DR*FLOAT(K-1)
      G6=U(MT+1)-FLP/(R+DR)**2
      G5=U(MT-1)-FLP/ R    **2
      IF(FK2.LT.0.0) GO TO 46
      IF(G5.LT.0.0) IFLAG=IFLAG+1
   46 CONTINUE
      Q6=1.0+DR2*G6
      Q5=1.0+DR2*G5
      F6=T6
      F5=T5
      W2=0.D0
      FNORM2=0.D0
      F(K  )=F6
      F(K-1)=F5
      DO 100 M=1,KT
      MM=K-M-1
      MK=MM+MM-1
      R=R-DR
      G4=U(MK+1)+V*U(MK  )-FLP/R**2
      Q4=1.0+DR2*G4
      F4=((12.0-10.0*Q5)*F5-Q6*F6)/Q4
      Q6=Q5
      Q5=Q4
      F6=F5
      F5=F4
      G6=G5
      G5=G4
      F(MM)=F4
      IF(G6*G5.GT.0.0) GO TO 90
      IFLAG=IFLAG+1
      IF(FK2.LT.0.0) GO TO 110
      IF(IFLAG.GE.2) GO TO 110
   90 CONTINUE
      TEMP=F6*F6
      IF(TEMP.LT.1.E-18) GO TO 100
      FNORM2=FNORM2+TEMP
      W2=W2+U(MK+2)*TEMP
      IF(TEMP.LE.UMAXSQ) GO TO 100
      RNORM=RNORM+1.D0
      FNORM2=FNORM2/UMAXSQ
      W2=W2/UMAXSQ
      F5=F5/UMAX
      F6=F6/UMAX
      MP=M+2
      DO 95 I=1,MP
      MPP=K+1-I
      IF(ABS(F(MPP)).GE.1.0/UMAX) GO TO 92
      F(MPP)=0.D0
      GO TO 95
   92 CONTINUE
      F(MPP)=F(MPP)/UMAX
   95 CONTINUE
  100 CONTINUE
      IF(IFLAG.GE.2) GO TO 110
      IFLAG=2
      KT=NINT(RM/DR)+1
      KT=K-KT-1
      GO TO 41
  110 CONTINUE
C
C     INTEGRATE FROM ORIGIN   TO CLASSICAL TURNING POINT
C
      KM=MM+1
  120 CONTINUE
      KS=NINT(FL/3.4)+2
      W1=0.D0
      FNORM1=0.D0
      F2=0.D0
      Q2=0.D0
      R =0.D0
      DO 200 M=1,KM
      MK=M+M-1
      R=R+DR
      Q3=1.0+DR2*(U(MK+1)+V*U(MK  )-FLP/R**2)
      IF(M.GT.KS) GO TO 150
      F3=R**LL
      GO TO 151
  150 CONTINUE
      F3=((12.0-10.0*Q2)*F2-Q1*F1)/Q3
  151 CONTINUE
      Q1=Q2
      Q2=Q3
      F1=F2
      F2=F3
      F(M)=F3
      TEMP=F2*F2
      FNORM1=FNORM1+TEMP
      W1=W1+U(MK)*TEMP
      IF(F1*F2.LT.0.0) FNCT=FNCT+1.D0
      IF(TEMP.LE.UMAXSQ) GO TO 200
      FNORM1=FNORM1/UMAXSQ
      W1=W1/UMAXSQ
      F1=F1/UMAX
      F2=F2/UMAX
      DO 195 I=1,M
      IF(ABS(F(I)).GE.1.0/UMAX) GO TO 192
      F(I)=0.D0
      GO TO 195
  192 CONTINUE
      F(I)=F(I)/UMAX
  195 CONTINUE
  200 CONTINUE
      IF(ISW.EQ.2) GO TO 650
      DET=(F1*F6-F5*F2)/(F2*F6*DR)
      FN=FNODE-FNCT
      R=FLOAT(KM)*DR
      FNORM1=FNORM1/(F2*F2)
      FNORM2=FNORM2/(F6*F6)
      FNORM=FNORM1+FNORM2
      IF(ICNT.EQ.ICMAX) GO TO 600
      IF(ISW.EQ.1) GO TO 451
C
C     CHOOSE NEXT GUESS ON WELL DEPTH
C
      IF(FN.EQ.0.0) GO TO 280
      FLAM=1.0+ABS(E)*3.0*FN/(V*FK*RM)
      GO TO 382
  280 FLAM=1.0-DET/(V*DR*(W1/(F2*F2)+W2/(F6*F6)))
  382 CONTINUE
      IF((ABS(DET).LE.DETX).AND.(ABS(FLAM-1.0).LE.FLAMX)) GO TO 500
      IF(FLAM.GT.1.2) FLAM=1.2D0
      IF(FLAM.LT.0.8) FLAM=0.8D0
      V=V*FLAM
      GO TO 10
C
C     CHOOSE NEXT GUESS ON BINDING ENERGY
C
  451 CONTINUE
      IF(FN.EQ.0.0) GO TO 480
      FLAM=1.0+ABS(E)*3.0*FN/(E*FK*RM)
      GO TO 482
  480 CONTINUE
      FLAM=1.0-DET/(DR*FK2*FNORM1)
  482 CONTINUE
      IF((ABS(DET).LE.DETX).AND.(ABS(FLAM-1.0).LE.FLAMX)) GO TO 500
      IF(FLAM.GT.1.2) FLAM=1.2D0
      IF(FLAM.LT.0.8) FLAM=0.8D0
      TEMP=SQRT(FLAM)
      FK=FK*TEMP
      ETA=ETA/TEMP
      TEMP=FK2*FLAM-FK2
      FK2=FK2+TEMP
      E=E*FLAM
      DO 485 M=1,K
      MK=M+M-1
      U(MK+1)=U(MK+1)+TEMP
  485 CONTINUE
      GO TO 10
C
C     NORMALIZE FUNCTION
C     FOR E .GT. 0., NORMALIZE TO VINCENT AND FORTUNE EQ 7 OR 30
C
  500 CONTINUE
	IF(E.GE.0.AND.CHOICE.EQ.'Y') GO TO 501
	WRITE(6,9500) V,DET,FNCT,R,E,ICNT
501	CONTINUE
      VOL=FNORM*DR*F6*F6
      FNORM=SQRT(FNORM*DR)
      IF(FK2.LT.0.0) GO TO 505
      RNORM=UMAX**NINT(RNORM)
      TEMP=VOL+SURF
      IF(RNORM.GT.2.0) TEMP=VOL
      GAM=2.0*E/(FK*TEMP)
      VOL=VOL*RNORM*RNORM
      WRITE(6,9502) GAM,VOL,SURF,K*DR
      FNORM=SQRT(ABS(TEMP))/F6
  505 CONTINUE
      TEMP=1.0/(F2*FNORM)
      R=0.D0
      DO 510 M=1,KM
      R=R+DR
      F(M)=F(M)*TEMP/R
  510 CONTINUE
      KM=KM+1
      TEMP=1.0/(F6*FNORM)
  515 CONTINUE
      DO 520 M=KM,K
      R=R+DR
      F(M)=F(M)*TEMP/R
  520 CONTINUE
      RM=FLOAT(KM)*DR
      RETURN
C
C     ERROR MESSAGE
C
  600 WRITE(6,9001) ICMAX
      IERR=1
      GO TO 500
C
C     NORM TO SIN(K*R+DELTA)/(K*R)
C
  650 CONTINUE
      DET=T3*T6-T4*T5
      A1=(F1*T6-F2*T5)/DET
      B1=(F2*T3-F1*T4)/DET
      DET=1.0/SQRT(A1*A1+B1*B1)
      A1=A1*DET
      B1=B1*DET

      TEMP=DET*SQRT(FK/(PI*E))
! 	DELTA=ACOSD(A1)
      DELTA=180.*ACOS(A1)/PI
	WRITE(8,9504)E,DELTA
	WRITE(6,9503)E,DELTA
      WRITE(6,9501)A1,B1
      KM=1
      R=0.D0
      GO TO 515
 9001 FORMAT(1H0,28HBOUND STATE SEARCH FAILS IN ,I2,11H ITERATIONS)
 9500 FORMAT(21X,6HV    =,F9.4,3X,6HDET  =,F9.4,3X,6HNODES=,F9.4,3X,
     1 6HRM   =,F9.4,3X,6HE    =,F9.4,3X,6HITER.=,I4)
 9501 FORMAT(21X,6HCOSD =,F9.4,9H   SIND =,F9.4)
 9502 FORMAT(31H0UNBOUND FORM FACTOR    GAMMA =,E12.4,15H MEV      VOL =
     1,E12.4,10H    SURF =,E12.4,10H    RMAX =, F7.2)
 9503 FORMAT(' Resonance energy in MeV =',F8.4,'  DELTA in Deg. =',F9.4)
 9504 FORMAT(2(8X,1PE15.7))
      END

      SUBROUTINE COU(R,RP,E,L,H,F,FP,G,GP,S)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(51),FP(51),G(51),GP(51),S(51)

      TE=2.*E
      TF=E*E
      LL=L
      IF(LL-50)20,35,35
   20 ELP=50.D0
      J=50
      GO TO 45
   35 ELP=LL
      J=LL
   45 A=ATAN (E/ELP)
      B=SQRT (TF+ELP*ELP)
      Y=A*(ELP-0.5)+E*(LOG(B)-1.)-SIN (A)/(12.*B)
     1  +SIN (3.*A)/(360.*B**3)-SIN (5.*A)/(1260.*B**5)
     2  +SIN (7.*A)/(1680.*B**7)-SIN (9.*A)/(1188.*B**9)
      K=J-1
      IF(J-LL)65,65,70
   65 S(J)=Y
   70 DO 100 I=1,K
      ELP=ELP-1.
      J=J-1
      Y=Y-ATAN (E/ELP)
      IF(J-LL)95,95,100
   95 S(J)=Y
  100 CONTINUE
      DEL1=R-TE
      RMAX=.41666667*(TF+4.*E+3.)
      IF(RMAX.LT.10.0) RMAX=10.D0
      DEL=R-RMAX
      IF(E-5.)280,130,130
  130 IF(ABS (DEL1)-ABS (DEL))140,140,280
  140 DEL=DEL1
      IF(DEL)147,145,147
  145 I=2
      GO TO 150
  147 I=1
  150 X=TE
      T1=TF
      T2=T1*T1
      T3=E** .666666667
      T4=T3*T3
      T5=T4*T4
      T6=T3*T5
      T7=T4*T6
      T8=T3*T7
      T9=E** .166666667
      Y=1.22340402*T9*(1.+.495957017D-1/T4-.888888889D-2/T1+.245519918
     1D-2/T6-.910895806D-3/T2+.845362D-3/T8)
      Z=-.707881773/T9*(1.-.172826039/T3+.317460317D-3/T1-.358121485
     1D-2/T5+.311782468D-3/T2-.907396643D-3/T7)
      GO TO 665
  280 IF(E)285,290,285
  285 IF(DEL)310,290,290
  290 X=R
      I=2
      GO TO 320
  310 X=RMAX
      I=1
  320 T1=TF
      T2=2.*X
      T3=X-E*LOG(T2)+S(1)
      T4=E/T2
      SS=1.D0
      TS=0.D0
      SL=0.D0
      TL=1.-E/X
      SSS=1.D0
      STS=0.D0
      SSL=0.D0
      STL=TL
      EN=0.D0
      DO 620 K=1,25
      T5=EN+1.
      T6=T5+EN
      T7=EN*T5
      T8=T6*T4/T5
      T9=(T1-T7)/(T2*T5)
      T5=T8*SS-T9*TS
      TS=T8*TS+T9*SS
      SS=T5
      IF(ABS (SS/SSS)-1.0E-10)630,630,540
  540 T5=T8*SL-T9*TL-SS/X
      TL=T8*TL+T9*SL-TS/X
      SL=T5
      SSS=SSS+SS
      STS=STS+TS
      SSL=SSL+SL
      STL=STL+TL
      EN=EN+1.
  620 CONTINUE
  630 T8=SIN (T3)
      T9=COS (T3)
      Y=SSS*T9-STS*T8
      Z=SSL*T9-STL*T8
  665 GO TO (670,810),I
  670 M=1
  671 N=ABS (DEL/H)+1.0
      DX=DEL/FLOAT(N)
      T1=DX/2.
      T2=DX/8.
      T3=TE
      DO 805 I=1,N
      T4=DX*(T3/X-1.)*Y
      X=X+T1
      T5=DX*(T3/X-1.)*(Y+T1*Z+T2*T4)
      X=X+T1
      T6=DX*(T3/X-1.)*(Y+DX*Z+T1*T5)
      Y=Y+DX*(Z+(T4+2.*T5)/6.)
      Z=Z+(T4+4.*T5+T6)/6.
  805 CONTINUE
      GO TO (810,828),M
  810 G(1)=Y
      M=2
      DEL=RP-R
      W=Z
      GO TO 671
  828 GP(1)=Y
      T1=TF
      T8=SQRT (1.+T1)
      G(2)=((1./R+E)*G(1)-W)/T8
      GP(2)=((1./RP+E)*Y-Z)/T8
      T2=1.D0
      T3=2.D0
      DO 910 I=3,LL
      T4=T2+T3
      T5=T2*T3
      T6=T3*SQRT (T2*T2+T1)
      T7=T2*SQRT (T3*T3+T1)
      G (I)=(T4*(E+T5/R )*G (I-1)-T6*G (I-2))/T7
      GP(I)=(T4*(E+T5/RP)*GP(I-1)-T6*GP(I-2))/T7
      T2=T2+1.
      T3=T3+1.
  910 CONTINUE
      I=L+11
      N=2.*R+11.
      IF(I-N)960,960,950
  950 N=I
  960 Y=1.0D-20
      YP=Y
      X=Y
      XP=X
      Z=0.D0
      ZP=Z
      T2=N
 1000 T3=T2+1.
      T4=T2+T3
      T5=T2*T3
      T6=T2*SQRT (T3*T3+T1)
      T7=T3*SQRT (T2*T2+T1)
      Y =(T4*(E+T5/R )*Y -T6*Z )/T7
      YP=(T4*(E+T5/RP)*YP-T6*ZP)/T7
      IF(N-LL)1060,1060,1080
 1060 F(N)=Y
      FP(N)=YP
      GO TO 1120
 1080 IF(1.-ABS (Y))1090,1090,1120
 1090 CONTINUE
      Y =Y *1.0D-20
      YP=YP*1.0D-20
      X =X *1.0D-20
      XP=XP*1.0D-20
 1120 N=N-1
      Z=X
      ZP=XP
      X=Y
      XP=YP
      T2=T2-1.
      IF(N)1150,1150,1000
 1150 Y=F(1)*G(2)-F(2)*G(1)
      YP=FP(1)*GP(2)-FP(2)*GP(1)
      Z=1./(Y*T8)
      ZP=1./(YP*T8)
      DO 1180 I=1,LL
      FP(I)=FP(I)*ZP
 1180 F(I)=F(I)*Z
      RETURN
      END
