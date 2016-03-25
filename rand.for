      function rand (n)

C     PORTABLE RANDOM NUMBER GENERATOR (SINGLE PRECISION VERSION)
C     ***********************************************************
C
C     A PORTABLE, GOOD AND MODERATELY FAST RANDOM NUMBER GENERATOR,
C     PROGRAMMED BY R. BRENT, 26 OCTOBER 1972.
C
C     WHEN CALLED WITH N = 0 RAND RETURNS A RANDOM NUMBER UNIFORMLY
C     DISTRIBUTED IN (0, 1)  (INCLUDING 0 BUT NOT 1).   THE GENERATOR
C     MAY BE INITIALIZED BY CALLING WITH NONZERO N (OTHERWISE A DEFAULT
C     INITIALIZATION WITH N = 8190 IS USED).
c     ******* modified by JRR to replace 0.0 by 1.e-10 **********
C
C     IF THE FLOATING-POINT WORD HAS A FRACTION OF IW BITS, THE
C     ALGORITHM RETURNS X(N)/(2**IW), WHERE
C     X(N) = X(N-1) + X(N-127) MODULO 2**IW.   SINCE 1 + X + X**127 IS
C     PRIMITIVE MODULO 2, THE PERIOD IS AT LEAST 2**127 - 1 > 10**38
C     (SEE KNUTH, THE ART OF COMPUTER PROGRAMMING, VOL. 2, 1969,
C     PP. 26, 32-34 AND 464).   X(N) IS STORED AS R3 = X(N)/2**IW.
C
C     THE NUMBER 127 MAY BE REPLACED BY 22, 60, 63, 153 OR 532
C     THROUGHOUT (SEE ZIEGLER AND BRILLHART, ON PRIMITIVE TRINOMIALS
C     (MOD 2), INFORM. CONTR. 13(1968), 541-554 AND 14(1969), 566-569).
C     SUGGESTED VALUES ARE 60 IF THE TABLE SIZE OR INITIALIZATION TIME
C     IS CRITICAL, AND 532 IF THE DEGREE OF RANDOMNESS IS CRITICAL.
C
C     RESTRICTIONS ARE THAT IW > 12, 17*8190 DOES NOT OVERFLOW IN
C     INTEGER ARITHMETIC, AND A BINARY COMPUTER WITH REASONABLE
C     FLOATING-POINT ARITHMETIC IS USED.
C
c     *************************************************************************

      real*4	r3(127),rand,rm,rmc,s,t,r1
      integer*4	n,i2,ikt,ic,id,iw,i
      common /randsv/ i2,r1,r3
      save	rm,rmc,ikt,id
      data	istrt /0/
      data	iw /-1/
c     **** if this bloc is saved and restored by calling prog ****
c     **** then the sequence will be repeated- even accross machines! *****

      if ((istrt.eq.0 .and. r1.eq.0.0) .or. n.ne.0) then
c	************** INITIALIZATION SEQUENCE ************
 
        if (iw.eq.-1) then
c 	  ************* determine fraction length iw etc. **************
	  s= 0.e0
	  t= 1.e0
   10	  iw = iw + 1
	    t = 0.5e0*t
	    r1 = s
	    s = s + t
	    if (s.gt.r1 .and. s.lt.1.e0) go to 10

	  ikt = (iw-1) / 12
	  ic = iw - 12*ikt
	  id= 1
	  do i= 1,13-ic
	    id = 2 * id
	    end do

c         **** compute rmc = 2**(-ic) exactly. ****
	  rmc= 1.e0
	  do i = 1,ic
	    rmc = 0.5e0*rmc
	    end do
c	  **** the following constant is 2**(-12) and must be exact. ****
	  rm = 0.015625e0*0.015625e0
	  end if

c	**** determine starting value for small generator. ****
	ir = mod (abs(n), 8190) + 1

c	**** initialize r3 ****
	do i2= 127,1,-1
	  r1 = 0.e0
	  do i = 1,ikt
c	    **** generate random integer ir in (0, 8191). ****
	    ir = mod (17*ir, 8191)
c	    **** add twelve random bits to r1. ****
	    r1 = (r1 + float(ir/2)) * rm
	    end do
	  ir = mod (17*ir, 8191)
c	  **** add last ic random bits to r1. ****
	  r1 = (r1 + float(ir/id)) * rmc
	  r3(i2) = r1
	  end do
	end if

c     *************** generate next random number *****************

      istrt= 1
      if (i2.eq.0) i2 = 127
c     **** add last and oldest saved random numbers. ****
      t = r1 + r3(i2)
c     **** reduce modulo one making sure arithmetic is exact. ****
      if (t.ge.1.e0) t = (r1 - 0.5e0) + (r3(i2) - 0.5e0)
      r1 = t
c     **** save for future use. ****
      r3(i2) = r1
c     **** update counter for next call. ****
      i2 = i2 - 1
c     **** return the random number generated. ****
      if (r1.lt.0.E0 .or. r1.ge.1.E0) then
	write (6,*) 'ERROR: RAND NUMBER=',r1
	stop 'RAND'
      else if (r1.eq.0.E0) then
	rand= 1.E-10
      else
        rand = r1
      end if

      return
      end
