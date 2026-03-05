      PROGRAM factornorgodico
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Este programa sirve para calcular los parametros no ergodicos de   C
C las funciones de dispersion intemedia conociendo el valor de gamma C
C el cual proviene del criterio.				     C
C								     C
C Enero 2007							     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none !toda varible se tiene que declarar
	integer Nk,i
	parameter (Nk=2**12)
	Double precision fdek(Nk,2), Sdek(Nk,2),g0(Nk)
	Double precision parametro,gama,kmi,yy,F2,fvw
	double precision yminima,lambda
	character*1 dec
CC
	write(*,*)'¿quieres tambien el self?(s o n)'
	read(*,*)dec
	write(*,*)'¿Valor del parametro de control?'
	read(*,*)parametro
	write(*,*)'¿cuanto vale gamma?'
	read(*,*)gama

CC
	write(*,*)'Calculando S(k)'
CC
C llame aqui la subrutina para el factor de estructura
CC
	fvw=parametro-parametro**2/16.D0
	call Sdeked(Sdek,fvw,Nk,60.D0,10.D0)
CC
	write(*,*)'calculando la k del primer minimo'
	kmi=yminima(Sdek,Nk)
	goto 100
	write(*,*)'definiendo funcion g0'
	do i=1,Nk
	 yy=Sdek(i,1)
	 F2=(4.D0/yy**2)*dcos(yy)+(2.D0/yy-4.D0/yy**3)*dsin(yy)
	 if(yy.eq.0.0D0)then
	  F2=2.D0/3.D0-yy**2/3.D0
	 endif
	 g0(i)=1.D0-(3.D0/2.D0)*F2
CC
C La definicion de gcero depende de la Vineyard que escojamos.
CC
	enddo
	write(*,*)'Calculando f(k)'
	call facnoerg(gama,Sdek,g0,kmi,Nk,fdek)
	open (10,file='facnoergomultfv563.dat')
	do i=1,Nk
	 write(10,*)fdek(i,1),fdek(i,2)
	enddo
	close(10)
 100	continue
	if(dec.eq.'s')then
	 do i=1,Nk
	  Sdek(i,2)=1.D0
	  g0(i)=1.D0
CC
C El caso self no lleva S(k).
CC
	 enddo
	 write(*,*)'Calculando fs(k)'
	 call facnoerg(gama,Sdek,g0,kmi,Nk,fdek)
	 open (20,file='facnoergselfduraM.dat')
	 do i=1,Nk
	  write(20,*)fdek(i,1),fdek(i,2)
	 enddo
	 close(20)
	endif
      END
CC
C Subrutinas de S(k)
CC
      Subroutine  Sdeked(Sdk,fracvol,Np,kmaxima,rmaxima)
!Sdeked(Sdk,Gr,fracvol,Np,kmaxima,rmaxima)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C Subrutina para evaluar el factor de estructura (S(k)), la función  C
C de distribución radial (g(r))  sistema con potencial de interac-   C
C ción tipo esfera dura utilizando la estrucutura de esfera dura     C
C obtenida con la solución exacta propuesta por M.S. Wertheim (PRL   C
C  vol.10,no.8, pp.321,1963) para la cerradura de Perckus-Yevick  y  C
C la expresion para S(k) en funcion de la transformada de Lapalace   C
C de g(r) (Physica 149A (1988),123-163).                             C
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        Implicit none !aguas toda variable se tiene que declarar
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C Parametros: *Sdk:arreglo para el factor de estructura (S(k))       C
C             *Gr: arreglo para la funcion de dist. radial (g(r))    C
C             *Yr: arreglo para la funcion Y(r)                      C
C             *fracvol:fraccion de volumen                           C
C             *Np:numero de puntos a evaluar en S(k) y g(r)          C
C             *kmaxima:valor de k maximo para evaluar S(k)           C
C             *rmaxima:valor de r maximo para evaluar g(r)           C
C             *n: exponente del potencial                            C
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        integer Np,Npg
        double precision fracvol,fv,rmaxima,kmaxima
        double precision Sdk(Np,2)!,Gr(Np,2),Yr(Np,2)
!parametro para extende k para ganar buena g(r)
        parameter(Npg=50)
        integer i
        double precision K,dK
!arreglo auxiliar para buena g(r)
!        double precision Sdkg(Npg*Np,2)
        double precision C1,C2
C funciones auxiliares para evaluar S(k)
        complex G
        double precision facdes
C Solo por comodidad
        fv=fracvol
C primero evaluamos una buena S(k) usado los resultados exactos 
C mediante la función facdes
        K=0.001D0
        dK=kmaxima/dble(Np)
        do i=1,Np
         Sdk(i,1)=K
         Sdk(i,2)=facdes(K,fv)                 
         K=K+dK
        enddo
CC
C definición del espaciado para una buena g(r)
CC
!         dK=1000.D0/Dble(Np)
!         K=0.001D0
CC
C Se usa el arreglo Sdkg para guardar información a k's
C grandes y poder hacer una buena evaluación de g(r)
CC                
!          do i=1,Npg*Np
!           Sdkg(i,1)=K
!           Sdkg(i,2)=facdes(K,fv)
!           K=K+dK
!          enddo
CC
C ya tenemos S(k) en hasta una k grande ahora le sacamos la
C transformada inversa de Fourier para obtener g(r)
CC
!         call  TIF3D(Npg*Np,Sdkg,Np,Gr,fv,rmaxima)
CC        
C ya se calcularon S(k) y g(r)
C Calcualmos Y(r)
CC
!         Do i=1,Np
!          Yr(i,1)=Gr(i,1)
!          if(Gr(i,1)/lam.le.1.D0) then
!                   
!            C1=(1.D0+2.D0*fve)**2-6.D0*fve*(1.D0+0.5D0*fve)**2*Gr(i,1)
!              C2=fve*(1.D0+2.D0*fve)**2*0.5D0*Gr(i,1)**3           
!              Yr(i,2)=(C1+C2)/(1.D0-fve)**4
!              Gr(i,2)=Yr(i,2)*dexp(-(1.D0/(lam*Gr(i,1)))**n)
!            else
!             
!              Yr(i,2)=Gr(i,2)
! 	     Gr(i,2)=Yr(i,2)*dexp(-(1.D0/(lam*Gr(i,1)))**n)           
!            endif
!            
!         Enddo
        Return
       End

       Complex function G(x,emp)
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Esta funcion calcula la transformada de Lapalce de la funcion de  C
C distibuion radial. Como auxliar a calculo de factor de estructura C
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
CC
C Subrutinas y funciones para calcular el factor no ergodico
CC
      Subroutine  facnoerg(gamma,Sk,gcero,kmin,nk,f)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Esta subrutiuna evalua los factores no ergodicos predichos por el  C
C criterio de atrapamiento.					     C
C Parametros: *gamma: valor del parametro gamma del criterio	     C
C	      *Sk:factor de estructura				     C
C	      *gcero:funcion relacionada con Vineyard		     C
C	      *kmin:posicion del primer minimo de S(k)		     C
C	      *nk:numero de puntos en k				     C
C	      *f:arreglo de salida				     C
C								     C
C Agosto de 2006						     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Implicit none !toda variable se tiene que declarar
CC
C Parametros
CC
	Integer nk
	Double precision kmin,gamma
	Double precision Sk(nk,2),gcero(nk),f(nk,2)
CC
C variables internas
CC
	Integer i
	Double precision l,den,num
	Double precision lambda
CC
C evaluando el factor no ergodico
CC
	DO i=1,nk
	 l=lambda(Sk(i,1),kmin)
	 den=Sk(i,2)*gcero(i)*l
	 num=gamma*(Sk(i,1))**2
	 f(i,1)=Sk(i,1)
	 f(i,2)=1.D0/(1.D0+ (num/den))
	ENDDO
	Return
      End

      Double precision function lambda(x,xminima)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								   C
C Funcion interpoladora lambda(k),se usa en la cerradura para      C 
C Cs(k,z) en la Teoria Autoconsistente. Ver detalles en la tesis   C
C de Laura Yeomans Reyna.                                          C
C Parametros: *x: valor de k en que queremos evaluar la funcion    C
C	      *xminima: valor de k en donde esta el primer minimo  C
C			de S(k)					   C
C								   C
C Agosto de 2006						   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       Implicit none

       Double precision x,xminima

       lambda= 1.D0/(1.D0 + (x/xminima)**2)
       Return

      End

      Double precision function yminima(Sk,N)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C							             C
C Esta funcion sirve para calcular el primer minimo de S(k), pero    C
C puede adapatarse para otra funcion solo hay que ponerla en un      C
C arreglo bidimensional                                              C
C Parametros: *Sk: arreglo en que esta S(k)                          C
C                 (primera columna k y segunda S(k))                 C
C	      *N: Tamaño del arreglo                                 C
C								     C
C Agosto del 2006
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       Implicit none      
       Integer N
C       Double precision emp

       Integer h,m,j
       Double precision Sk(N,2)
C       Double precision y,dy
       Double precision smin,smax,kmin,kmax

c       Double precision facdes

C       y=3.D0
C       dy=0.001D0
       
C       Do i=1,N

C          S(i,1)=y
C          S(i,2)= facdes()
C          y= y+ dy

C       enddo

       smax=Sk(1,2)
               
       Do j=2, N
          if (Sk(j,2).GT.smax)then
            
             smax=Sk(j,2)
             m= j
          endif

       enddo
       
       smin= Sk(m,2)

       Do h=m , N-m

          if (Sk(h,2).LT.smin)then

             smin= Sk(h,2)
             kmin= Sk(h,1)

          endif

       enddo

       yminima=kmin
       Return

       End