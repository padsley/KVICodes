	PROGRAM ALLJ
	IMPLICIT REAL*8(A-H,O-Z)
C++
C TITLE: ALLJ
C ABSTRACT:
C	THIS PROGRAM CALCULATES CLEBSCH GORDON COEFFICIENTS,
C	RACAHS, 3JS, 6JS, AND 9JS SYMBOLS.
C	IT COULD BE USED EASILY AND IS MEANT TO BE USED AS A
C	SUBSTITUTE FOR THE 3J,6J HANDBOOK FOR CALCULATING
C	THESE COEFFICIENTS.
C	CONSULT THE 3J,6J HANDBOOK IN THE KVI LIBRARY 
C	FOR THE FORMULAS AND NOTATIONS USED IN THIS PROGRAM.
C PROCEDURES CALLED:
C	TJ, SIXJ, ANINEJ, W, CG
C	CIO - package.
C AUTHOR: M.N.Harakeh. KVI, Groningen, Netherlands.
C CREATION DATE:
C MODIFIED BY:
C PAK 28-oct-1980.	Modified for VAX.
C--

	PARAMETER	NCMDS = 6
	CHARACTER*5	CMDTBL(NCMDS)
	DATA CMDTBL /'C-G','3J','RACAH','6J','9J','EXIT'/

C Calculate binomial coefficients once and for all

    1	CALL CIOCMD ('ALLJ>', CMDTBL, NCMDS, I)
	GOTO (10,20,30,40,50,900),I
  900	STOP 'End of ALLJ'

C Clebsh Gordon coefficient

   10	CALL CIOWRT ('Input for Clebsh-Gordon coefficient')
	CALL CIOFX ('J1', XJ1, '$', *1)
	CALL CIOFX ('J2', XJ2, '$', *1)
	CALL CIOFX ('M1', XM1, '$', *1)
	CALL CIOFX ('M2', XM2, '$', *1)
	CALL CIOFX ('J3', XJ3, '$', *1)
	CALL CIOFX ('M3', XM3, '$', *1)

	PRINT 9014, CG(XJ1,XM1,XJ2,XM2,XJ3,XM3)
 9014	FORMAT(' Result is ',F14.6)
	GOTO 1

C 3J symbols.

   20	CALL CIOWRT ('Input for 3J symbols')
	CALL CIOFX ('J1', XJ1, '$', *1)
	CALL CIOFX ('J2', XJ2, '$', *1)
	CALL CIOFX ('J3', XJ3, '$', *1)
	CALL CIOFX ('M1', XM1, '$', *1)
	CALL CIOFX ('M2', XM2, '$', *1)
	CALL CIOFX ('M3', XM3, '$', *1)

	PRINT 9014, TJ(XJ1,XJ2,XJ3,XM1,XM2,XM3)
	GOTO 1

C Racah.

   30	CALL CIOWRT ('Input for RACAH')
	CALL CIOFX ('J1', XJ1, '$', *1)
	CALL CIOFX ('J2', XJ2, '$', *1)
	CALL CIOFX ('L2', XL2, '$', *1)
	CALL CIOFX ('L1', XL1, '$', *1)
	CALL CIOFX ('J3', XJ3, '$', *1)
	CALL CIOFX ('L3', XL3, '$', *1)

	PRINT 9014, W(XJ1,XJ2,XL2,XL1,XJ3,XL3)
	GOTO 1

C 6J symbols.

   40	CALL CIOWRT ('Input for 6J symbols')
	CALL CIOFX ('J1', XJ1, '$', *1)
	CALL CIOFX ('J2', XJ2, '$', *1)
	CALL CIOFX ('J3', XJ3, '$', *1)
	CALL CIOFX ('L1', XL1, '$', *1)
	CALL CIOFX ('L2', XL2, '$', *1)
	CALL CIOFX ('L3', XL3, '$', *1)

	PRINT 9014, SIXJ(XJ1,XJ2,XJ3,XL1,XL2,XL3)
	GOTO 1

C 9J symbols.

   50	CALL CIOWRT ('Input for 9J symbols')
	CALL CIOFX ('A1', A1, '$', *1)
	CALL CIOFX ('A2', A2, '$', *1)
	CALL CIOFX ('A3', A3, '$', *1)
	CALL CIOFX ('B1', B1, '$', *1)
	CALL CIOFX ('B2', B2, '$', *1)
	CALL CIOFX ('B3', B3, '$', *1)
	CALL CIOFX ('C1', C1, '$', *1)
	CALL CIOFX ('C2', C2, '$', *1)
	CALL CIOFX ('C3', C3, '$', *1)

	PRINT 9014, ANINEJ(A1,A2,A3,B1,B2,B3,C1,C2,C3)
	GOTO 1
	END
