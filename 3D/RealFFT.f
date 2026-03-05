      PROGRAM RFFT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Programa para probar la subrutina de transofrmada rapida de Fourier	C
C para una funcion real.						C
C									C
C Septiembre 2008							C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none !toda variable se tiene que declarar
	Integer i,Nk
	parameter(Nk=2**14)
	Double precision pi,dr,dk,R,K,Fr(Nk) 
	pi=4.D0*datan(1.D0)
	dr=50.D0/dble(Nk) !5.D-2
	dk=2.D0*pi/(dble(Nk)*dr)
	do i=1,Nk
	 R=dble(i)*dr
	 Fr(i)=dexp(-R*R)
	 write(20,*)R,Fr(i)
	enddo
	call realft(Fr,Nk,1)
	do i=1,Nk/2
	 K=dble(i)*dk
	 write(30,*)K,Fr(2*i-1)*2.D0*dr
	enddo 
	call realft(Fr,Nk,-1)
	do i=1,Nk
	 R=dble(i)*dr
	 write(40,*)R,Fr(i)*2.D0/dble(Nk)
	enddo 
	stop
      END



      SUBROUTINE sinft(y,n)
  	INTEGER n
  	Double precision y(n)
! C USES realft
! Calculates the sine transform of a set of n real-valued data points stored in array y(1:n
! The number n must be a power of 2. On exit y is replaced by its transform. This program
! without changes, also calculates the inverse sine transform, but in this case the output
! array 
!  Should be multiplied by 2/n.
  	INTEGER j
  	REAL sum,y1,y2
  	DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp   !Double precision in the trigonometric
						! recurrences.
  	theta=3.141592653589793d0/dble(n)  !Initialize the recurrence.
  	wr=1.0d0
  	wi=0.0d0
  	wpr=-2.0d0*sin(0.5d0*theta)**2
  	wpi=sin(theta)
  	y(1)=0.0
  	do j=1,n/2
       	 wtemp=wr
       	 wr=wr*wpr-wi*wpi+wr                  ! Calculate the sine for the auxiliary array.
       	 wi=wi*wpr+wtemp*wpi+wi             !The cosine is needed to continue the recurrence.
       	 y1=wi*(y(j+1)+y(n-j+1))            !   Construct the auxiliary array.
       	 y2=0.5*(y(j+1)-y(n-j+1))
         y(j+1)=y1+y2                        !  Terms j and N − j are related.
       	 y(n-j+1)=y1-y2
  	enddo 
  	call realft(y,n,+1)                       ! Transform the auxiliary array.
  	sum=0.0
  	y(1)=0.5*y(1)               !Initialize the sum used for odd terms below.
  	y(2)=0.0
  	do j=1,n-1,2
       	 sum=sum+y(j)
       	 y(j)=y(j+1)                    !Even terms in the transform are determined directly.
         y(j+1)=sum                            !Odd terms are determined by this running sum.
 	enddo 
  	return
      END

      SUBROUTINE realft(data,n,isign)
	INTEGER isign,n
	Double precision data(n)
C USES four1
!       Calculates the Fourier transform of a set of n real-valued data points. Replaces this data
!       (which is stored in array data(1:n)) by the positive frequency half of its complex Fourier
!       transform. The real-valued ﬁrst and last components of the complex transform are returned
!       as elements data(1) and data(2), respectively. n must be a power of 2. This routine
!       also calculates the inverse transform of a complex data array if it is the transform of real
!       data. (Result in this case must be multiplied by 2/n.)
	INTEGER i,i1,i2,i3,i4,n2p3
	REAL c1,c2,h1i,h1r,h2i,h2r,wis,wrs
	DOUBLE PRECISION theta,wi,wpi,wpr, wr,wtemp    !Double precision for the trigonometric recurrences.
	theta=3.141592653589793d0/dble(n/2)           !Initialize the recurrence.
	c1=0.5
	if (isign.eq.1) then
	 c2=-0.5
         call four1(data,n/2,+1)                  !The forward transform is here.
  	else
	 c2=0.5                                   !Otherwise set up for an inverse transform.
	 theta=-theta
	endif
	wpr=-2.0d0*sin(0.5d0*theta)**2
	wpi=sin(theta)
	wr=1.0d0+wpr
	wi=wpi
	n2p3=n+3
	do i=2,n/4                                 !Case i=1 done separately below.
	 i1=2*i-1
	 i2=i1+1
	 i3=n2p3-i2
	 i4=i3+1
	 wrs=sngl(wr)
	 wis=sngl(wi)
	 h1r=c1*(data(i1)+data(i3))              ! The two separate transforms are separated out of
	 h1i=c1*(data(i2)-data(i4))               !     data.
	 h2r=-c2*(data(i2)+data(i4))
	 h2i=c2*(data(i1)-data(i3))
	 data(i1)=h1r+wrs*h2r-wis*h2i            ! Here they are recombined to form the true trans-
	 data(i2)=h1i+wrs*h2i+wis*h2r               !   form of the original real data.
	 data(i3)=h1r-wrs*h2r+wis*h2i
	 data(i4)=-h1i+wrs*h2i+wis*h2r
	 wtemp=wr                                 !The recurrence.
	 wr=wr*wpr-wi*wpi+wr
	 wi=wi*wpr+wtemp*wpi+wi
	enddo
	if (isign.eq.1) then
	 h1r=data(1)
	 data(1)=h1r+data(2)
	 data(2)=h1r-data(2)     ! Squeeze the ﬁrst and last data together to get
	else                             !them all within the original array.
	 h1r=data(1)
	 data(1)=c1*(h1r+data(2))
	 data(2)=c1*(h1r-data(2))
	 call four1(data,n/2,-1)  !This is the inverse transform for the case isign=-1.
	endif
	return
      END

      SUBROUTINE four1(data,nn,isign)
	INTEGER isign,nn
	Double precision data(2*nn)
!       Replaces data(1:2*nn) by its discrete Fourier transform, if isign is input as 1; or replaces
!       data(1:2*nn) by nn times its inverse discrete Fourier transform, if isign is input as −1.
!       data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn
!       MUST be an integer power of 2 (this is not checked for!).
	INTEGER i,istep,j,m,mmax,n
	REAL tempi,tempr
	DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp    ! Double precision for the trigonometric recurrences.
	n=2*nn 
	j=1
	do i=1,n,2                 ! This is the bit-reversal section of the routine.
	 if(j.gt.i)then
	  tempr=data(j)        !Exchange the two complex numbers.
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
         endif
         m=n/2
 1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
          goto 1
         endif
         j=j+m
        enddo 
	mmax=2                        ! Here begins the Danielson-Lanczos section of the routine.
 2	if (n.gt.mmax) then           ! Outer loop executed log2 nn times.
	 istep=2*mmax
	 theta=6.28318530717959d0/(isign*mmax)              ! Initialize for the trigonometric recur-
	 wpr=-2.d0*sin(0.5d0*theta)**2
	 wpi=sin(theta)
	 wr=1.d0
	 wi=0.d0
	 do m=1,mmax,2          !Here are the two nested inner loops.
	  do i=m,n,istep
           j=i+mmax                  !This is the Danielson-Lanczos formula:
           tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
           tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
           data(j)=data(i)-tempr
           data(j+1)=data(i+1)-tempi
           data(i)=data(i)+tempr
           data(i+1)=data(i+1)+tempi
          enddo
          wtemp=wr              !Trigonometric recurrence.
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
         enddo 
         mmax=istep
	 goto 2                        !Not yet done.
	endif                         !All done.
	return
      END
