      Program Fourier
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    C
C Programa para probar la subrutina de transformada de Fourier	    C
C								    C
C octubre 2006							    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none
CC
C Numero de puntos de las funciones
CC
	Integer N
	parameter(N=2**10)
CC
C Arreglos para las funciones, original y transformada
CC
	Double precision Sdk(N,2),Gdr(N,2)
CC
C Auxiliares
CC
	integer i
	Double precision Hdk(N,2)
	Double precision R,dK,dR,pi,rho
CC
C funciones externas
CC
	Double precision facdes
	pi=4.D0*DATAN(1.D0)
CC
	write(*,*)'Definicion de la funcion'
CC
	dK=0.1D0
	dR=1.D0/(N*dK)
	rho=6.D0*0.4D0/pi
	DO i=1,N
	 Sdk(i,1)=i*dK
	 Sdk(i,2)=facdes(Sdk(i,1),0.4D0)
	 Hdk(i,1)=Sdk(i,1)
	 R=i*dR
	 Hdk(i,2)=Sdk(i,1)*(Sdk(i,2)-1.D0)/(rho*R)
! 	 write(40,*)Hdk(i,1),Hdk(i,2)
	ENDDO
CC
	write(*,*)'ahora hacemos la transformacion'
CC
	call sfft(Hdk,Gdr,N,2)
CC
	write(*,*)'Salidas'
CC
	open(8,file='Gder.dat')
	open(9,file='Sdek.dat')
	DO i=1,N
	 write(8,*)Gdr(i,1),Gdr(i,2)+1.D0
	 write(9,*)Sdk(i,1),Sdk(i,2)
	ENDDO
	close(8)
	close(9)
      End
	

       Complex function G(x,emp)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Esta funcion calcula la transformada de Lapalce de la funcion de  C
C distibucion radial. Como auxliar a calculo de factor de estructuraC
C de esfera dura                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        Implicit none !aguas toda variable se tiene que declarar

C declaracion de argumentos
        Complex x
        double precision emp

C variables locales
        Complex a,b,c,d,num,den1,den2,den3,den

C calculos auxiliares

         a= cmplx(1.D0+emp/2.D0,0.D0)
         b= cmplx(1.D0+2.D0*emp,0.D0)
         c= cmplx(1.D0-emp,0.D0)
         d= cmplx(18.D0,0.D0)

         num= x*(a*x+b)
         den1= cmplx(12.D0*emp,0.D0)*(a*x+b)
         den2= (c**2*x**3+cmplx(6.D0*emp,0.D0)*c*x**2)*exp(x)
         den3=(d*cmplx(emp,0.D0)**2*x-cmplx(12.D0*emp,0.D0)*b)*exp(x)
         den= den1+den2+den3

C valor que regresa la funcion

         G= num/den

         return

       end
      
       double precision function facdes(z,fv)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Esta funcion calcula el factor de estructura.C
C de esfera dura.Necesario usar la funcion G   C
C auxiliar.                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
        Implicit none !aguas toda variable se tiene que declarar

C decalracion de argumentos

        double precision  z, fv,corte

C decalraciones locales 
        Complex cz, mcz

C declaracion funcion auxiliar      
        complex  G
        
C variables auxiliares

        double precision T0,T2,F1,F2 
        
        parameter(corte=0.5D0) !parametro para fijar k's pequeñas
        
C Evitando problemas a k's chicas

        if(z.lt.corte) then
        
            T0=(2.D0*(-1.D0+fv)**4*(0.5D0+fv))/(1.D0+2.D0*fv)**3
            F1=0.4D0*(-1.D0+fv)**4*fv*(0.5D0+fv)
            F2=(4.D0+(-2.75D0+fv)*fv)/((1.D0+2.D0*fv)**5) 
            T2=F1*F2
            
            facdes=T0+T2*z**2
           
            Return
          
        else

          cz= dcmplx(0.D0,z)
          mcz= dcmplx(0.D0,(-1.D0)*z)

C valor de la funcion
          facdes= 1.D0 + 12.D0*fv*(dble((G(mcz,fv)-G(cz,fv))/cz))
      
          Return
          
        endif
        
       End

      Subroutine sfft(F,TF,N,inv)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    C
C Esta subrutina sirve para sacar la transformada seno de Fourier   C 
C de una funcion. Recuerde que la transformada y la antitransformadaC
C tienen la misma forma.					    C
C Parametros:							    C
C	*F: funcion a transformar,arreglo de 2 columnas, primera    C
C	    contine variable independiente y segunda variable dep.  C
C	*TF: transformada de la funcion, tambien arreglo de 2 colum-C
C	     nas.						    C
C	*N: numero de puntos en el arreglo, debe ser un multiplo    C
C 	    entero de 2.					    C
C	*inv: entero para especificar si es transf.(1) o 	    C
C	      antitransf. (2).					    C
C Nota: Estoy usando una notacion generalizada, t es la variable en C
C el espacio real y f la variable en el espacio reciproco.          C
C								    C
C octubre 2006							    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Implicit none
CC
C Declaracion de parametros
CC
	Integer N,inv
	Double Precision F(N,2),TF(N,2)
CC
C Arreglo auxiliar para transformar
CC
	Real AF(N)
CC
C Variables para construir la malla
CC
	Integer im
	Double precision df,dt
CC
C llenamos el arreglo auxiliar para transformar
CC
	DO im=1,N
	 AF(im)=Real(F(im,2))
	ENDDO
CC
C Ahora usamos la transformada rapida de Fourier
CC
	call sinft(AF,N)
CC
C Llenamos el arreglo de la funcion transformada
CC
	df=1.D0/F(N,1)
	if(inv.eq.1)then
	 DO im=1,N
	  TF(im,1)=im*df
	  TF(im,2)=dble(AF(im))
	 ENDDO
	endif
	if(inv.eq.2)then
	 write(*,*)'inversa'
	 DO im=1,N
	  TF(im,1)=im*df
	  TF(im,2)=2.D0*dble(AF(im))/Dble(N)
	 ENDDO
	endif
	Return
      END
	
      SUBROUTINE sinft(y,n)
	INTEGER n
  	REAL y(n)
! C USES realft
!  Calculates the sine transform of a set of n real-valued data points
!  stored in array y(1:n).
!  The number n must be a power of 2. On exit y is replaced by its
!  transform. This program, without changes, also calculates the
!  inverse sine transform, but in this case the output array should be
!  multiplied by 2/n.

	INTEGER j
	REAL sum,y1,y2
! Double precision in the trigonometric recurrences.
 	DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
! Initialize the recurrence.
	theta=3.141592653589793d0/dble(n) 
  	wr=1.0d0
  	wi=0.0d0
  	wpr=-2.0d0*sin(0.5d0*theta)**2
  	wpi=sin(theta)
  	y(1)=0.0
	do j=1,n/2
         wtemp=wr
! Calculate the sine for the auxiliary array.
         wr=wr*wpr-wi*wpi+wr
!  The cosine is needed to continue the recurrence.
         wi=wi*wpr+wtemp*wpi+wi
! Construct the auxiliary array.
         y1=wi*(y(j+1)+y(n-j+1))
         y2=0.5*(y(j+1)-y(n-j+1))
! Terms j and N − j are related.
         y(j+1)=y1+y2
         y(n-j+1)=y1-y2
  	enddo
! Transform the auxiliary array.
  	call realft(y,n,+1)
  	sum=0.0
! Initialize the sum used for o dd terms below.
  	y(1)=0.5*y(1)
 	y(2)=0.0
  	do  j=1,n-1,2
         sum=sum+y(j)
!  Even terms in the transform are determined directly.
         y(j)=y(j+1)
! Odd terms are determined by this running sum.
         y(j+1)=sum
  	enddo
  	return
      END

      SUBROUTINE realft(data,n,isign)
  	INTEGER isign,n
  	REAL data(n)
C USES four1
! Calculates the Fourier transform of a set of n real-valued data
! points. Replaces this data (which is stored in array data(1:n)) by
! the positive frequency half of its complex Fouriertransform. The
! real-valued ﬁrst and last components of the complex transform are
! returned as elements data(1) and data(2), respectively. n must be a
! power of 2. This routine also calculates the inverse transform of a
! complex data array if it is the transform of real data. (Result in
! this case must be multiplied by 2/n.)

	INTEGER i,i1,i2,i3,i4,n2p3
  	REAL c1,c2,h1i,h1r,h2i,h2r,wis,wrs
! Double precision for the trigonometric recurrences.
  	DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
! Initialize the recurrence.
  	theta=3.141592653589793d0/dble(n/2)
  	c1=0.5
  	if (isign.eq.1) then
       	 c2=-0.5
! The forward transform is here.
         call four1(data,n/2,+1)
  	else
! Otherwise set up for an inverse transform.
         c2=0.5
         theta=-theta
  	endif
  	wpr=-2.0d0*sin(0.5d0*theta)**2
  	wpi=sin(theta)
  	wr=1.0d0+wpr
  	wi=wpi
  	n2p3=n+3
! Case i=1 done separately below.
  	do i=2,n/4
         i1=2*i-1
         i2=i1+1
         i3=n2p3-i2
         i4=i3+1
         wrs=sngl(wr)
         wis=sngl(wi)
! The two separate transforms are separated out of data.
         h1r=c1*(data(i1)+data(i3))
         h1i=c1*(data(i2)-data(i4))
         h2r=-c2*(data(i2)+data(i4))
         h2i=c2*(data(i1)-data(i3))
! Here they are recombined to form the true transform of the original
! real data.
         data(i1)=h1r+wrs*h2r-wis*h2i
         data(i2)=h1i+wrs*h2i+wis*h2r
         data(i3)=h1r-wrs*h2r+wis*h2i
         data(i4)=-h1i+wrs*h2i+wis*h2r
! The recurrence.
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
  	enddo 
	if (isign.eq.1) then
    	 h1r=data(1)
    	 data(1)=h1r+data(2)
! Squeeze the ﬁrst and last data together to get
    	 data(2)=h1r-data(2)
! them all within the original array.
	else
    	 h1r=data(1)
    	 data(1)=c1*(h1r+data(2))
    	 data(2)=c1*(h1r-data(2))
! This is the inverse transform for the case isign=-1.
    	 call four1(data,n/2,-1)
	endif
	return
      END

      SUBROUTINE four1(data,nn,isign)
  	INTEGER isign,nn
  	REAL data(2*nn)

!  Replaces data(1:2*nn) by its discrete Fourier transform, if isign
!  is input as 1; or replaces data(1:2*nn) by nn times its inverse
! discrete Fourier transform, if isign is input as −1.data is a 
! complex array of length nn or, equivalently, a real array of length
! 2*nn. nn MUST be an integer power of 2 (this is not checked for!).

 	INTEGER i,istep,j,m,mmax,n
  	REAL tempi,tempr
! Double precision for the trigonometric recurrences.
  	DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
  	n=2*nn
  	j=1
! This is the bit-reversal section of the routine.
  	do  i=1,n,2
       	 if(j.gt.i)then
! Exchange the two complex numbers.
           tempr=data(j)
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
! Here begins the Danielson-Lanczos section of the routine.
  	mmax=2
! Outer lo op executed log2 nn times.
 2 	if (n.gt.mmax) then
         istep=2*mmax
! Initialize for the trigonometric recurrence
         theta=6.28318530717959d0/(isign*mmax)
         wpr=-2.d0*sin(0.5d0*theta)**2
         wpi=sin(theta)
         wr=1.d0
         wi=0.d0
! Here are the two nested inner lo ops.
         do  m=1,mmax,2
          do  i=m,n,istep
! This is the Danielson-Lanczos formula:
           j=i+mmax
           tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
           tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
           data(j)=data(i)-tempr
           data(j+1)=data(i+1)-tempi
           data(i)=data(i)+tempr
           data(i+1)=data(i+1)+tempi
          enddo 
! Trigonometric recurrence.
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
    	 enddo 
    	 mmax=istep
! Not yet done.
	 goto 2
! All done.
	endif
	return
      END
