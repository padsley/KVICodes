	SUBROUTINE CIOX
!++
! TITLE:  CIOX
! FACILITY:  CIO  - conversational I/O package.
! ABSTRACT: This module contains several entries which do free
!	format input with optional prompt, default, lower and
!	upper limits, and special characters which may control
!	the flow of a program.
! ENVIRONMENT:  FORTRAN-77.
! PROCEDURES CALLED:
!	CIODSS, CIOERL, CIOGIT, CIOINV, CIOMSG, CIOYNL,
!	CIOLUNFREE, INDEX, STR__UPCASE
! AUTHOR: Peter A. Kroon, KVI, Groningen, NL.
! CREATION DATE:  21-oct-1980.
! MODIFIED BY:
! PAK 28-nov-1980.	Added CIOS and CIOSX entries.
! PAK 26-jan-1981.	Changed call to LUNFREE to LIB$FREE_LUN.
!			Changed call to CIOCLU to STR__UPCASE.
!			Other F77 related changes.
!			Included HELP option.
! PAK 27-feb-1981.	Added option to specify default radix in specification
!			string for CIOI and CIOIX.
! PAK 13-may-1981.	Implemented qualifiers.
! PAK 25-jan-1984.	Allow TRUE and FALSE as valid input for CIOL(X)
! F.Z. 30-Aug-1990.	Fixed bug in CIOBX, which did not used the BDEFAULT
!			argument but the LDEF argument of CIOLX.
! F.Z. 13-Dec-1993.	Adapt to AXP. Used IARGCOUNT instead of NOARG.
! FUNCTIONAL DESCRIPTION:
! ENTRIES:
!	CIOI	Integer input
!	CIOIX	Idem, extended calling sequence.
!	CIOF	Floating (real) input
!	CIOFX	Idem, extended calling sequence.
!	CIOD	Double precision (real*8) input
!	CIODX	Idem, extended calling sequence.
!	CIOB	Boolean (logical*1) input
!	CIOBX	Idem, extended calling sequence.
!	CIOL	Logical (logical*4) input
!	CIOLX	Idem, extended calling sequence.
!	CIOA	Alphanumeric (character) input
!	CIOAX	Idem, extended calling sequence.
!	CIOS	Alphanumeric (character) input; returns string length.
!	CIOSX	Idem, extended calling sequence.
!
! CALLING SEQUENCE: See separate documentation.
!--

	INCLUDE 'ciocom.inc'
	CHARACTER	FORMAT*25,C*1
	CHARACTER	PROMPT*(*)
	CHARACTER	SPECST*(*)
	CHARACTER	CIOYNL*3
	INTEGER		IARG1, IARG2, IARG3
	REAL*4		FARG1, FARG2, FARG3
	REAL*8		DARG1, DARG2, DARG3
	INTEGER		IVALUE, ITEMP, ILOLIM, IHILIM, IDEFAULT
	REAL*8		FVALUE
	REAL*8		DVALUE, DTEMP, DLOLIM, DHILIM, DDEFAULT
	LOGICAL		LVALUE, LTEMP, LDEFAULT, LDFTMP
	LOGICAL*1	BVALUE, BDEFAULT
	INTEGER		IALEN
	LOGICAL*1	LIALEN
	CHARACTER*(*)	AVALUE,ADEF
	CHARACTER*(*) POSTIV, NEGTIV
	PARAMETER	(POSTIV	= 'YESJATRUE')
	PARAMETER	(NEGTIV	= 'NONEENFALSE')
	INTEGER		NARG
!g77	INTRINSIC	IARGCOUNT
	INTEGER		KEY1, KEY2
	INTEGER		ITEML
	INTEGER		I, J
	INTEGER		KRETUR
	INTEGER		STATUS

	CHARACTER*(*) TXCUR, TXDEF, TXLOW, TXHIGH, TXSIZ
	PARAMETER	(TXCUR	= 'current value (=) ...')
	PARAMETER	(TXDEF	= 'default value (*) ...')
	PARAMETER	(TXLOW	= 'lower limit .........')
	PARAMETER	(TXHIGH	= 'upper limit .........')
	PARAMETER	(TXSIZ	= 'maximum stringsize ..')

	CHARACTER*4	CTEMP,DSFLOW
	DATA		DSFLOW /S_FLOW/


!----------------------------------------------------------------------
! CIOI, CIOIX
!----------------------------------------------------------------------

	ENTRY CIOI (PROMPT,IVALUE)
	CALL CIODSS (' ')
	GOTO 1000

	ENTRY CIOIX (PROMPT,IVALUE,SPECST,IARG1,IARG2,IARG3,*,*,*,*)

	CALL CIODSS (SPECST)
	! KEY1 is the label where we go if an invalid flow control char is typed
 1000	ASSIGN 1050 TO KEY1
! Extract low, high and default if specified.
 1010	GOTO (1011,1012,1013), KSPLOW
	GOTO 1020
 1011	ILOLIM = IARG1
	GOTO 1020
 1012	ILOLIM = IARG2
	GOTO 1020
 1013	ILOLIM = IARG3
 1020	GOTO (1021,1022,1023), KSPHIG
	GOTO 1030
 1021	IHILIM = IARG1
	GOTO 1030
 1022	IHILIM = IARG2
	GOTO 1030
 1023	IHILIM = IARG3
 1030	GOTO (1031,1032,1033), KSPDEF
	GOTO 1040
 1031	IDEFAULT = IARG1
	GOTO 1040
 1032	IDEFAULT = IARG2
	GOTO 1040
 1033	IDEFAULT = IARG3
 1040	CONTINUE

! Extract item from input buffer.
 1050	ITEMP = IVALUE
 1055	CALL CIOGIT (PROMPT,CWORK,ITEML,*1060,*1070,*1080, *100)
!	alternate returns:		 =     *     ?    flow

! Examine first character on 'O' for octal, 'X' for hex or 'D' for
! decimal conversion.
! If none of these, assume radix as specified in specification string.
	CALL STR__UPCASE (C, CWORK(1:1))
	IF (C.EQ.'O') THEN		! user prefixed value by 'O'
	  goto 1056
	ELSE IF (C.EQ.'X') THEN		! user prefixed value by 'X'
	  C = 'Z'
	  goto 1056
	ELSE IF (C.EQ.'D') THEN		! user prefixed value by 'D'
	  C = 'I'
 1056	  CWORK = CWORK(2:ITEML)	! remove the radix indicator
	  ITEML = ITEML - 1
	ELSE				! no prefix
	  C = CSPRDX			! use default radix.
	ENDIF

! Construct the proper format, then decode the item.
	WRITE (FORMAT,91050) C,ITEML
91050	FORMAT('(',A,I2,')')
	READ (CWORK,FORMAT,ERR=300) ITEMP
! Check against limits if specified.
 1060	IF (KSPLOW.GT.0) THEN
	  IF (ITEMP.LT.ILOLIM) GOTO 400
	ENDIF
	IF (KSPHIG.GT.0) THEN
	  IF (ITEMP.GT.IHILIM) GOTO 400
	ENDIF
! Give value to caller and join common return code.
 1065	IVALUE = ITEMP
	GOTO 700
! Use default if specified, then go check agains limits if specified.
 1070	IF (KSPDEF.EQ.0) GOTO 500
	ITEMP = IDEFAULT
	GOTO 1060
! Show current, default and limits as far as specified,
! and radix if not decimal.
! Then join common code to display the defined path's.
 1080	WRITE (FORMAT,91080) CSPRDX		! encode dynamic format
91080	FORMAT ('(1X,:,T8,A,T30,',A,'12)')	! according to default radix.
	WRITE(LUNO,FORMAT) !blank line
	WRITE(LUNO,FORMAT) TXCUR,IVALUE
	IF (KSPDEF.GT.0) WRITE(LUNO,FORMAT) TXDEF,IDEFAULT
	IF (KSPLOW.GT.0) WRITE(LUNO,FORMAT) TXLOW,ILOLIM
	IF (KSPHIG.GT.0) WRITE(LUNO,FORMAT) TXHIGH,IHILIM
	IF (CSPRDX.EQ.'O') THEN
	  WRITE (LUNO,'(T8,''input in radix ......        octal'')')
	ELSEIF (CSPRDX.EQ.'Z') THEN
	  WRITE (LUNO,'(T8,''input in radix ......  hexadecimal'')')
	ENDIF
	GOTO 600

!----------------------------------------------------------------------
! CIOF, CIOFX
!----------------------------------------------------------------------

	ENTRY CIOF (PROMPT,FVALUE)
	CALL CIODSS (' ')
	GOTO 2000

	ENTRY CIOFX (PROMPT,FVALUE,SPECST,FARG1,FARG2,FARG3,*,*,*,*)

	CALL CIOMSG(0,'CIOAX L189')
	CALL CIOMSG(0,SPECST)
	CALL CIODSS(SPECST)
 2000 CONTINUE
	! KEY1 is the label where we proceed on invalid input
	ASSIGN 2050 TO KEY1
	! KEY2 is the label where we proceed on valid input
	ASSIGN 2065 TO KEY2
! Extract low, high and default if specified.
 2010	GOTO (2011,2012,2013), KSPLOW
	GOTO 2020
 2011	DLOLIM = FARG1
	GOTO 2020
 2012	DLOLIM = FARG2
	GOTO 2020
 2013	DLOLIM = FARG3
 2020	GOTO (2021,2022,2023), KSPHIG
	GOTO 2030
 2021	DHILIM = FARG1
	GOTO 2030
 2022	DHILIM = FARG2
	GOTO 2030
 2023	DHILIM = FARG3
 2030	GOTO (2031,2032,2033), KSPDEF
	GOTO 2040
 2031	DDEFAULT = FARG1
	GOTO 2040
 2032	DDEFAULT = FARG2
	GOTO 2040
 2033	DDEFAULT = FARG3
 2040	CONTINUE

! Copy current value to temporary.
! Then join the code for double precision conversion.
 2050	DTEMP = FVALUE
	GOTO 3055

! Give value to caller and join common return code.
 2065	FVALUE = DTEMP
	GOTO 700

!----------------------------------------------------------------------
! CIOD, CIODX
!----------------------------------------------------------------------

	ENTRY CIOD (PROMPT,DVALUE)
	CALL CIODSS (' ')
	GOTO 3000

	ENTRY CIODX (PROMPT,DVALUE,SPECST,DARG1,DARG2,DARG3,*,*,*,*)

	CALL CIODSS (SPECST)
 3000	CONTINUE
	! KEY1 is the label where we proceed on invalid input
	ASSIGN 3050 TO KEY1
	! KEY2 is the label where we proceed on valid input
	ASSIGN 3065 TO KEY2
! Extract low, high and default if specified.
 3010	GOTO (3011,3012,3013), KSPLOW
	GOTO 3020
 3011	DLOLIM = DARG1
	GOTO 3020
 3012	DLOLIM = DARG2
	GOTO 3020
 3013	DLOLIM = DARG3
 3020	GOTO (3021,3022,3023), KSPHIG
	GOTO 3030
 3021	DHILIM = DARG1
	GOTO 3030
 3022	DHILIM = DARG2
	GOTO 3030
 3023	DHILIM = DARG3
 3030	GOTO (3031,3032,3033), KSPDEF
	GOTO 3040
 3031	DDEFAULT = DARG1
	GOTO 3040
 3032	DDEFAULT = DARG2
	GOTO 3040
 3033	DDEFAULT = DARG3
 3040	CONTINUE

! Copy current value to temporary.
 3050	DTEMP = DVALUE
! Extract item from input buffer and decode it.
! The single-floating (CIOF) code joins us here.
 3055	CALL CIOGIT (PROMPT,CWORK,ITEML,*3060,*3070,*3080, *100)
!	alternate returns:		 =     *     ?    flow
	WRITE (FORMAT,93050) ITEML
93050	FORMAT('(F',I2,'.0)')
	READ (CWORK,FORMAT,ERR=300) DTEMP
! Check against limits if specified.
 3060	IF (KSPLOW.GT.0) THEN
	  IF (DTEMP.LT.DLOLIM) GOTO 400
	ENDIF
	IF (KSPHIG.GT.0) THEN
	  IF (DTEMP.GT.DHILIM) GOTO 400
	ENDIF
! Give value to caller and join common return code.
! The single-floating (CIOF) code leaves us here again.
 3062	GOTO KEY2,(2065,3065)
 3065	DVALUE = DTEMP
	GOTO 700
! Use default if specified, then go check agains limits if specified.
 3070	IF (KSPDEF.EQ.0) GOTO 500
	DTEMP = DDEFAULT
	GOTO 3060
! Show current, default and limits as far as specified.
! Then join the common code to display the defined path's.
 3080	WRITE(LUNO,93080) !blank line
	WRITE(LUNO,93080) TXCUR,DVALUE
	IF (KSPDEF.GT.0) WRITE(LUNO,93080) TXDEF,DDEFAULT
	IF (KSPLOW.GT.0) WRITE(LUNO,93080) TXLOW,DLOLIM
	IF (KSPHIG.GT.0) WRITE(LUNO,93080) TXHIGH,DHILIM
93080	FORMAT(1X,:,T8,A,T30,G16.6)
	GOTO 600

!----------------------------------------------------------------------
! CIOB, CIOBX, CIOBXD
!----------------------------------------------------------------------

	ENTRY CIOB (PROMPT,BVALUE)
	CALL CIODSS (' ')
	GOTO 4000

	ENTRY CIOBX  (PROMPT,BVALUE,SPECST,BDEFAULT,*,*,*,*)

	CALL CIODSS (SPECST)
 4000	CONTINUE
	! KEY1 is the label where we proceed on invalid input
	ASSIGN 4050 TO KEY1
	! KEY1 is the label where we proceed on valid input
	ASSIGN 4065 TO KEY2

! Copy current value and default to temporaries.
! Then join code for LOGICAL*4 (CIOL)
 4050	LTEMP  = BVALUE
	LDFTMP = BDEFAULT
	GOTO 5055

! Give value to caller and join common return code.
 4065	BVALUE = LTEMP
	GOTO 700

!----------------------------------------------------------------------
! CIOL and CIOLX
!----------------------------------------------------------------------

	ENTRY CIOL (PROMPT,LVALUE)
	CALL CIODSS (' ')
	GOTO 5000

	ENTRY CIOLX (PROMPT,LVALUE,SPECST,LDEFAULT, *,*,*,*)

	CALL CIODSS (SPECST)
 5000	CONTINUE
	! KEY1 is the label where we go on invalid input
	ASSIGN 5050 TO KEY1
	ASSIGN 5065 TO KEY2

! Copy current value and default to temporaries.
 5050	LTEMP  = lVALUE
	LDFTMP = LDEFAULT

! Extract item from input buffer and decode it.
 5055	CALL CIOGIT (PROMPT,CWORK,ITEML,*5060,*5070,*5080, *100)
!	alternate returns:		 =     *     ?    flow
! Translate CWORK in either .TRUE. or .FALSE.
! First translate CWORK from lower to upper case.
	CALL STR__UPCASE (CWORK(:ITEML), CWORK(:ITEML))	!??? g77
	IF (INDEX(POSTIV,CWORK(:ITEML)).NE.0) THEN
	  LTEMP = .TRUE.
	ELSEIF (INDEX(NEGTIV,CWORK(:ITEML)).NE.0) THEN
	  LTEMP = .FALSE.
	ELSE
	  GOTO 300
	ENDIF
! Give value to caller and join common return code.
 5060	GOTO KEY2,(4065,5065)
 5065	LVALUE = LTEMP
	GOTO 700
! Use default if specified.
 5070	IF (KSPDEF.EQ.0) GOTO 500
	LTEMP = LDFTMP
	GOTO 5060
! Show current and default as far as specified.
! Then join the common code to display the defined path's.
 5080	WRITE(LUNO,95080) !blank line
	WRITE(LUNO,95080) TXCUR,CIOYNL(LVALUE)
	IF (KSPDEF.GT.0) WRITE(LUNO,95080) TXDEF,CIOYNL(LDEFAULT)
95080	FORMAT(1X,:,T8,A,T30,A)
	GOTO 600

!----------------------------------------------------------------------
! CIOS and CIOSX
!----------------------------------------------------------------------

	ENTRY CIOS (PROMPT, AVALUE, IALEN)
	LIALEN = .TRUE.
	GOTO 6001

	ENTRY CIOA (PROMPT, AVALUE)
	LIALEN = .FALSE.
 6001	CALL CIODSS (' ')
	GOTO 6000

	ENTRY CIOSX (PROMPT, AVALUE, IALEN, SPECST, ADEF,*,*,*,*)
	LIALEN = .TRUE.
      GOTO 6002

	ENTRY CIOAX (PROMPT, AVALUE, SPECST, ADEF, *,*,*,*)
	LIALEN = .FALSE.
 6002	CALL CIODSS (SPECST)

 6000	CONTINUE
	! KEY1 is the label where we go on invalid input
	ASSIGN 6050 TO KEY1

! Extract item from input buffer and decode it.
 6050	CALL CIOGIT (PROMPT,CWORK,ITEML,*6090,*6070,*6080, *100)
!	alternate returns:		 =     *     ?    flow
! Check if callers string is long enough to contain CWORK.
	IF (ITEML.GT.LEN(AVALUE)) THEN
	  ITEML = LEN(AVALUE)
	  CALL CIOMSG(0,'warning - string truncated to')
	  CALL CIOMSG(0,CWORK(:LEN(AVALUE)))
	ENDIF
! Give value to caller and join common return code.
 6060	AVALUE = CWORK(:ITEML)
	GOTO 6072
! Use default if specified and join common return code.
 6070	IF (KSPDEF.EQ.0) GOTO 500
	AVALUE = ADEF
	ITEML = LEN (ADEF)
! Convert to upper case if caller requested.
 6072	IF (KSPCNV.NE.0) CALL STR__UPCASE (AVALUE, AVALUE)
! If called via CIOAL or CIOALX return string length to caller.
	IF (LIALEN) IALEN = ITEML
! Join common return code.
	GOTO 700
! Show current and default as far as specified.
! Then join the common code to display the defined path's.
 6080	WRITE(LUNO,96080) !blank line
	WRITE(LUNO,96080) TXCUR,AVALUE
	IF (KSPDEF.GT.0) WRITE(LUNO,96080) TXDEF,ADEF
	WRITE(LUNO,96081) TXSIZ,LEN(AVALUE)
96080	FORMAT(1X,:,T8,A,T30,A)
96081	FORMAT(1X,:,T8,A,T30,I5)
	GOTO 600
! Come here when equal sign entered. For CIOS and CIOSX entries
! we now must find out the length of the string currently in
! AVALUE, since the users expects us to return the length.
 6090	DO 6092 ITEML = LEN(AVALUE), 1, -1
	IF (AVALUE(ITEML:ITEML).NE.' ') GOTO 6072
 6092	CONTINUE
	ITEML = 0
	GOTO 6072

!----------------------------------------------------------------------
! This section is common to all entries in this module.
!----------------------------------------------------------------------

! Come here if CIO control character typed.
! See if alternate return has been specified for this control
! character. If so, take that return. If not complain, and
! invalidate input buffer and alternate files.

  100	I = INDEX (S_FLOW, CWORK(1:1))
	IF (KSPCTL(I).EQ.0) THEN
	  CALL CIOERL('you can''t go that way')
	  GOTO KEY1
	ELSE
	  KRETUR = KSPCTL(I)
	  GOTO 710
	ENDIF

! Come here if decode error
! Complain and invalidate inputbuffer and alternate files.
  300	CALL CIOERL ('bad input')
	GOTO KEY1

! Come here if limits exceeded
! Complain and invalidate inputbuffer and alternate files.
  400	CALL CIOERL ('value outside limits')
	GOTO KEY1

! Come here if default requested while no default specified.
! Complain and invalidate inputbuffer and alternate files.
  500	CALL CIOERL ('no default available')
	GOTO KEY1

! Come here to display the currently defined flow control characters.
! Only display if at least one is defined.
  600	CTEMP = ' '
	J = 0
	DO 610 I=1,4
	IF (KSPCTL(I).NE.0) THEN
	  J = J + 1
	  CTEMP(J:J) = DSFLOW(I:I)
	ENDIF
  610	CONTINUE
! 	IF (J.GT.0)
! 	  WRITE (LUNO,'(T8,''you can go to .......'',T30,4(4X,A))')
! 	  (CTEMP(I:I),I=1,J)
! Display the legal qualifiers, if any.
	CALL CIOQSH
! And a final blank line.
	WRITE (LUNO,'(1X)')
	GOTO KEY1

! Come here just before returning to the caller in order to purge
! the internal promptbuffer, delete the possibly open prompt-
! scratchfile, and deallocate the corresponding logical unit number.
  700	KRETUR = 0
! Code for alternate returns caused by flow-control-char's joins here.
  710	IPBPTR = 0
	IF (LUNPRF .NE. 0) THEN
	  CLOSE (UNIT = LUNPRF)
	  CALL CIOLUNFREE (LUNPRF, STATUS) !???g77 status testen
	  LUNPRF = 0
	ENDIF
	RETURN KRETUR

	END
