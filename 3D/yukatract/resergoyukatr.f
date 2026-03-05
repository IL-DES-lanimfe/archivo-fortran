      Program critresergo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C Este programa aplica el criterio de atrapamiento al estilo         C
C "restauracion de la egodicidad" desarrollado por Marco Chavez Rojo,C
C Pedro E. Ramirez Gonzalez y el Dr. Magdaleno Medina, basandose en  C
C la forma en que MCT obtiene el suyo.Esta pensado para dos para-    C
C metros de control,la aplicacion del criterio a sistemas diferentes C
C se da cambiando el factor de  estructura para el caso especifico   C
C que estemos tratando.                                              C
C																   C
C Agosto 2006                                                        C
C 																   C
C He modificado este programa para que calcule el inverso de deltaz  C
C (gamma) en vez de deltaz, esto es para estandarizarme con el arti- C
C culo que publicado en PRE 76, 041504 (2007).					   C
C Adaptado a yukawas atractivas.									   C
C																   C
C Enero 2007														   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        Implicit none !Toda variable se tiene que declarar
CC
C Declaraciones
CC
C Parametros de control
CC
	 Integer Npar2
	 Double precision parametro2
	 Double precision par1max,par2max
	 Double precision par1min,par2min
	 Double precision dpar2
	 Double precision fv
CC
C variables para bisecciones sucesivas
CC
	 Double precision par1post,par1ant,par1med
	 Double precision integralpost,integralant,integralmed
	 Double precision solpost,solant,solmed
	 Double precision solprod,prodant,prodpost
CC
C variables de estatica 
CC
	 Integer Nk
	 parameter(Nk=2**12)
	 Double precision SdeK(Nk,2)
	 Double precision kmin
CC
C variables de integracion
CC
	 Integer iint
	 Double precision y,dy,S,l,gama,num,den,integral
	 Double precision T(Nk)
CC
C variables axiliares
CC
	 Integer i
	 Double precision pi,a,gamam
	 Double precision Smx,mc
	 Double precision fnoer(Nk,2)
CC
C funciones que se usan
CC
	 Double precision lambda,yminima,Smax
CC
C Variables del factor de estructura
CC
	 Double precision z,Rk(Nk),Sk(Nk),xi
CC
C Saludo inicial y lectura de datos
CC
	 write(*,*)' '
	 write(*,*)'Este programa es para generar un diagrama de fases'
	 write(*,*)'de vitrificacion con dos parametros con RE.'
	 write(*,*)'Estamos en el caso de pot. de yukawa atractivo'
	 write(*,*)'parametro 1:fraccion de volumen'
	 write(*,*)'parametro 2:K(inverso de la temperatura)'
 	 write(*,*)'diga el intervalo de busqueda del parametro 1'
	 write(*,*)'maximo?'
	 read(*,*)par1max
	 write(*,*)'minimo?'
	 read(*,*)par1min
 	 write(*,*)'diga el intervalo de busqueda del parametro 2'
	 write(*,*)'maximo?'
	 read(*,*)par2max
	 write(*,*)'minimo?'
	 read(*,*)par2min
! 	write(*,*)'numero de puntos?'
! 	read(*,*)Npar2
CC
C ya entraron todos los datos obligatorios
CC	
	 pi= dacos(-1.D0)
C Si ninguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente
CC
C	 write(*,*)'¿fraccion de volumen?'
C	 read(*,*)fv
CC
	 write(*,*)' '
	 write(*,*)'ya rugistes leon, ahora a trabajar'
	 write(*,*)' '
       open(8,file='vidrioyukatrz10.dat') 
       open(9,file='facnoergoyukatrz10.dat')
       open(10,file='Skmaxyukatrz10.dat')
CC
C se usa el metodo de bisecciones sucesivas. Parametro 2 fijo y 
C bisecciones sucesivas en parametro 1 para localizar la transicion
C luego se cambia el parametro 2 y se repite la operacion
CC
	 z=10.D0
	 dpar2=0.2D0 !(par2max-par2min)/dble(Npar2)
	 DO parametro2=par2min,par2max,dpar2
	  par1post=par1max
        par1ant= par1min
CCCCCCCCCCCCCCCCCCCCCCC Modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  write(*,*)'calculando S(k) cota inferior'
C ponga aqui la subrutina de S(k) evaluandola en par1ant y parametro2
CC
	  call msa(par1ant,parametro2,z,0.01D0,Nk,Rk,Sk,xi)
	  do i=1,Nk
	   Sdek(i,1)=2.D0*Rk(i)
	   Sdek(i,2)=Sk(i)
	  enddo
C inicia el calculo de S(k) y g(r)
CC
	  dy=Sdek(2,1)-Sdek(1,1)
	  write(*,*)'calculando k del primer minimo de S(k)'
	  kmin=yminima(Sdek,Nk)
CCCCCCCCCCCCCCCCCCCC Fin modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCC 

	  write(*,*)'evaluando criterio en cota inferior'
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
	  fv=par1ant
C    	  fv=parametro2
CC
CC inicia la evaluacion de la integral del criterio
CC
	  integralant=0.D0
	  Do gama=0.D0,0.1D0,1.D-5
	   a=gama
         do iint=1, Nk
	    y=Sdek(iint,1)
          S=Sdek(iint,2)
          l= lambda(y,kmin)
	    num=a*((S-1.D0)*l*y**2)**2
	    den=(a*y**2+S*l)*(a*y**2+l)
          T(iint)=num/(36.D0*pi*fv*den)
         enddo
	   call SIMPSON(T,dy,Nk,integral)
	   if(integral.gt.integralant)then
	    integralant=integral
	    gamam=a
	   endif
	  EndDo
	  solant=1.D0 - integralant
CC
C Criterio de convergencia
CC
	  if (dabs(solant).lt.1.D-3)then
         write(8,*)sngl(par1ant),sngl(parametro2),sngl(gamam)
         write(*,*)'se encontro solucion para: ',par1ant,parametro2
         write(*,*)'evaluando factor no ergodico'
         call facnoerg(gamam,Sdek,kmin,Nk,fnoer)
	   write(9,*)'p',par1ant,parametro2
         do i=1,Nk
          write(9,*)fnoer(i,1),fnoer(i,2)
         enddo
         write(*,*)'calculando el maximo de S(k)'
         Smx=Smax(Sdek,Nk)
         write(10,*)sngl(par1ant),sngl(parametro2),sngl(Smx),
     &	        sngl(dsqrt(gamam)*(6.d0*fv/pi)**(1.d0/3.D0))
         goto 200
        endif
CC
C fin criterio de convergencia
CC
	  write(*,*)'la evaluacion de la integral en',par1ant
	  write(*,*)' es',integralant
CC
C Termina la evaluacion de la  integral para valor de parametro 1 mas
C pequeño 
CC
CCCCCCCCCCCCCCCCCCCCCC Modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCCCCC
        write(*,*)'calculando S(k) en la cota superior'
C ponga aqui la subrutina de S(k) evaluandola en par1post y parametro2
	  call msa(par1post,parametro2,z,0.01D0,Nk,Rk,Sk,xi)
	  do i=1,Nk
	   Sdek(i,1)=2.D0*Rk(i)
	   Sdek(i,2)=Sk(i)
	  enddo
	  dy=Sdek(2,1)-Sdek(1,1)
	  write(*,*)'calculando k del primer minimo de S(k)'
	  kmin=yminima(Sdek,Nk)
CCCCCCCCCCCCCCCCCCCC Fin modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCCC 
	  write(*,*)'evaluando criterio en la cota superior'
CC
C Se inicia el calculo del parametro gamma
CC
CC
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
	  fv=par1post
C    	  fv=parametro2
CC
C inicia la evaluacion de la integral del criterio
CC
	  integralpost=0.D0
	  Do gama=0.D0,0.1D0,1.D-5
	   a=gama
         do iint=1, Nk
	    y=Sdek(iint,1)
          S=Sdek(iint,2)
          l= lambda(y,kmin)
	    num=a*((S-1)*l*y**2)**2
	    den=(a*y**2+S*l)*(a*y**2+l)
          T(iint)=num/(36.D0*pi*fv*den)
         enddo
	   call SIMPSON(T,dy,Nk,integral)
	   if(integral.gt.integralpost)then
	    integralpost=integral
	    gamam=a
	   endif
	  EndDo
	  solpost=1.D0 - integralpost
CC
C Criterio de convergencia
CC
	  if (dabs(solpost).lt.1.D-3)then
         write(8,*)sngl(par1post),sngl(parametro2),sngl(gamam)
	   write(*,*)'se encontro solucion para: ',par1post,parametro2
         write(*,*)'evaluando factor no ergodico'
         call facnoerg(gamam,Sdek,kmin,Nk,fnoer)
	   write(9,*)'p',par1post,parametro2
         do i=1,Nk
          write(9,*)fnoer(i,1),fnoer(i,2)
         enddo
         write(*,*)'calculando el maximo de S(k)'
         Smx=Smax(Sdek,Nk)
         write(10,*)sngl(par1post),sngl(parametro2),sngl(Smx),
     &	        sngl(dsqrt(gamam)*(6.d0*fv/pi)**(1.d0/3.D0))
         goto 200
        endif
CC
C Fin criterio de convergencia
CC
	  write(*,*)'la evaluacion de la integral en',par1post
	  write(*,*)' es',integralpost
CC
C Termina la evaluacion de la  integral para valor de parametro 1 mas
C grande 
CC
        solprod=solant*solpost
        IF(solprod.lt.0.D0)then
 100     continue
	   write(*,*)'bisecciones sucesivas'
         par1med=(par1post+par1ant)/2.D0
CCCCCCCCCCCCCCCCCCCCCC Modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCCCCC
	   write(*,*)'calculando S(k) en el valor intermedio'
C llame aqui a la subrutina de S(k) evaluandola en par1med y
C parametro2
	   call msa(par1med,parametro2,z,0.01D0,Nk,Rk,Sk,xi)
	   do i=1,Nk
	    Sdek(i,1)=2.D0*Rk(i)
	    Sdek(i,2)=Sk(i)
	   enddo
	   dy=Sdek(2,1)-Sdek(1,1)
	   write(*,*)'calculando k del primer minimo de S(k)'
	   kmin=yminima(Sdek,Nk)
CCCCCCCCCCCCCCCCCCCC Fin modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCC 
	   write(*,*)'evaluando criterio en valor intermedio'
CC
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
	  fv=par1med
C    	  fv=parametro2
CC
CC inicia la evaluacion de la integral del criterio
CC
	  integralmed=0.D0
	  Do gama=0.D0,0.1D0,1.D-5
	   a=gama
         do iint=1, Nk
	    y=Sdek(iint,1)
          S=Sdek(iint,2)
          l= lambda(y,kmin)
	    num=a*((S-1)*l*y**2)**2
	    den=(a*y**2+S*l)*(a*y**2+l)
          T(iint)=num/(36.D0*pi*fv*den)
         enddo
	   call SIMPSON(T,dy,Nk,integral)
	   if(integral.gt.integralmed)then
	    integralmed=integral
	    gamam=a
	   endif
	  EndDo
	  solmed=1.D0 - integralmed
CC
C Criterio de convergencia
CC
	  if (dabs(solmed).lt.1.D-3)then
         write(8,*)sngl(par1med),sngl(parametro2),sngl(gamam)
         write(*,*)'se encontro solucion para: ',par1med,parametro2
         write(*,*)'evaluando factor no ergodico'
         call facnoerg(gamam,Sdek,kmin,Nk,fnoer)
	   write(9,*)'p',par1med,parametro2
         do i=1,Nk
          write(9,*)fnoer(i,1),fnoer(i,2)
         enddo
         write(*,*)'calculando el maximo de S(k)'
         Smx=Smax(Sdek,Nk)
         write(10,*)sngl(par1med),sngl(parametro2),sngl(Smx),
     &	        sngl(dsqrt(gamam)*(6.d0*fv/pi)**(1.d0/3.D0))
         goto 200
        endif
CC
C Fin criterio de convergencia
CC
	  write(*,*)'la evaluacion de la integral en',par1med
	  write(*,*)' es',integralmed

CC
C Termina la evaluacion de la  integral para valor de parametro 1 
C intermedio 
CC
C Decidiendo el intervalo en donde esta la solucion
CC
         prodpost= solmed*solpost
         prodant=solmed*solant
         if(prodant.lt.0.D0)then
          par1post=par1med
	    solpost=solmed
          goto 100
         endif
         if(prodpost.lt.0.D0)then
          par1ant=par1med
	    solant=solmed
          goto 100
         endif
	  ENDIF
        write(*,*)'la solucion no esta en ese intervalo'
 200    continue
	 ENDDO
       close(8)
	 close(9)
	 close(10)
	 stop
      End

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                  C
C Subrutina del factor de estructura.Para evitar broncas con la    C
C precision de la integracion que se necesita para resolver el     C
C criterio es necesario que el valor maximo de k para el cual      C
C tienes datos de S(k) sea al menos 100, aunque siempre puedes     C
C probar para ver cual es el valor optimo.			             C
C                                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine msa(vf,XDK,b,DQ,N,Q,SQ,chi)
! This program finds the MSA-Yukawa S(Q)
! The solution of Cummings and Smith (Chem. Phys. 42, 241 (79))
! is used. J. Bergenholtz U. Konstanz 01.28.98
! Variables:
! vf   - volume fraction
! XDK  - Yukawa prefactor in units of kT (XDK > 0 for attraction) 
! b    - screening parameter in units of sigma
! DQ   - increment of q vector in units of RADIUS 
! N    - number of q vectors 
! Q    - q vector array in units of RADIUS
! SQ   - structure factor 
      double precision SQ(N), Q(N)
      double precision SMSA,vf,XDK,b,DQ
      double precision aa,bb,dd,beta,XA,XB,dde,chi 
      integer i
! ************************
      CALL QUART(XDK,b,vf,aa,bb,dd,beta,dde,chi)
!       print *,'MSA SOLUTION:'
!       print *,'beta=',beta
!       print *,'d=',dd,' dde=',dde
!       print *,'vf=',vf
      do i=1,N
! DQ given in units of radius      
        Q(i)=2.*DQ*(i-1)
	SQ(i)=SMSA(Q(i),b,vf,aa,bb,dd,beta,XA,XB,dde)
! Q to be converted to units of radius
	Q(i)=Q(i)/2.
      enddo
      return
      end
      
      function SMSA(Q,b,vf,aa,bb,dd,beta,XA,XB,dde)
! Calculates S(Q) where Q=q*sigma for the MSA-Yukawa
! J. Bergenholtz U. Konstanz 01/23/98
      double precision SMSA,Q,b,vf,aa,bb,dd,beta,dde
      double precision XA,XB
      if(Q.eq.0.0d0) then
	SMSA=1./aa/aa
      else
          XA=1.-12.*vf*(aa*(Q*dcos(Q)-dsin(Q))/Q/Q/Q+
     ^       bb*(dcos(Q)-1.)/Q/Q+beta*dd*dsin(Q)/Q+
     ^         beta*b*dde/(b*b+Q*Q)+
     ^     beta*dd*(b*dcos(Q)-Q*dsin(Q))/(b*b+Q*Q))
          XB=-12.*vf*(aa*(Q*dsin(Q)+dcos(Q)-1.-
     ^     .5*Q*Q)/Q/Q/Q+bb*(dsin(Q)-Q)/Q/Q+
     ^       beta*dd*(1.-dcos(Q))/Q+
     ^       beta*Q*dde/(b*b+Q*Q)+
     ^     beta*dd*(b*dsin(Q)+Q*dcos(Q))/(b*b+Q*Q))
	SMSA=1./(XA*XA+XB*XB)
      endif
      return
      end
      
      subroutine QUART(XDK,b,vf,aa,bb,dd,beta,dde,chi)
! This program solves the MSA for
! Yukawa particles and returns the constants
! needed for S(Q)
! J. Bergenholtz U. Konstanz 01/23/98 
      double precision vf,XDK,b
      double precision X,Y,W0,W1,W2,W3,W4
      double precision root1,root2,x1,x2,rtsafe
      double precision xb1(2),xb2(2)
      double precision aa,bb,dd,beta,chi
      double precision T,S,F,E,D,dde
      double precision gamma,rho,tau,efdiff 
      integer i,j,nb
      COMMON/CON/W0,W1,W2,W3,W4
! ****************************** 
      IF(XDK.eq.0.0d0) THEN
! PY LIMIT:      
        aa=(1.+2.*vf)/(1.-vf)/(1.-vf)
	bb=-1.5*vf/(1.-vf)/(1.-vf)
	dd=0.0		! dd = -E/F but dd is always multiplied by beta
	beta=0.0
      ELSE
! same as in Cummings & Smith
      X=6.*vf*(b*dexp(-b)-6.*vf*(2.-2.*b-
     ^        dexp(-b)*(2.-b*b))/b/b/(1.-vf)-
     ^       18.*vf*vf*(2.-b-dexp(-b)*(2.+b))/
     ^          b/b/(1.-vf)/(1.-vf))
! same as in Cummings & Smith
      Y=b-6.*vf*(2.-b*b-2.*dexp(-b)*(1.+b))/b/b/(1.-vf)-
     ^   18.*vf*vf*(2.-b-dexp(-b)*(2.+b))/b/b/(1.-vf)/(1.-vf)
! COMPUTE COEFFICIENTS:
! same as in Cummings & Smith
      W4=36.*vf*vf
      W3=-X
      W2=12.*vf*XDK
      W1=-XDK*Y
      W0=XDK*XDK
!       print *,'W0,W1,W2,W3,W4=',W0,W1,W2,W3,W4
! SOLVE FOR REAL POSITIVE ROOTS: ********
! BRACKET AND USE RTSAFE:
      x1=0.0
      x2=500.
      nb=2
! ***************************************      
      call zbrak(x1,x2,1000,xb1,xb2,nb)
! ROOT1 IN xb1(1)-xb2(1)
! ROOT2 IN xb1(2)-xb2(2)
! FIND ONE ROOT:      
      root1=rtsafe(xb1(1),xb2(1),1.D-10)
!       print *,'root1=',root1
! ATTEMPT AT SECOND ROOT:      
      root2=rtsafe(xb1(2),xb2(2),1.D-10)
!       print *,'root2=',root2
! IDENTIFY PROPER ROOT (LOWER):
      if(root1.lt.root2) then
        beta=root1
      else
        beta=root2
      endif
!       print *,'beta=',beta
! COMPUTE REMAINING CONSTANTS:
! same as in Cummings & Smith 
      T=12.*vf*(1.-b-dexp(-b))/b
! same as in Cummings & Smith
      S=12.*vf*(1.-.5*b*b-dexp(-b)*(1.+b))/b/b
      gamma=12.*vf/b/(1.-vf)/(1.-vf)
      tau=1.+2.*vf-6.*vf/b
      rho=1.5*vf+(1.-4.*vf)/b
! same as in Cummings & Smith
      F=-6.*vf*(1.-dexp(-b))**2+gamma*(tau*(1.-dexp(-b)*(1.+b))+
     ^    3.*vf*b*dexp(-b))*S-gamma*(rho*(1.-dexp(-b)*(1.+b))-
     ^      .5*(1.-4.*vf)*b*dexp(-b))*T
! same as in Cummings & Smith
      E=-6.*vf+gamma*tau*S-gamma*rho*T
! same as in Cummings & Smith 
      D=b-(1.+2.*vf)*S/(1.-vf)/(1.-vf)+3.*vf*T/2./(1.-vf)/(1.-vf)
! same as small d in Cummings & Smith
      dd=((-XDK+beta*D)*dexp(-b)+E*beta*beta)/F/beta/beta
! efdiff=exp(b)(E-F)
      efdiff=-6.*vf*(2.-dexp(-b))+gamma*S*(tau*(1.+b)-3.*vf*b)-
     ^         gamma*T*(rho*(1.+b)+(.5-2.*vf)*b)
      chi=(-XDK+beta*D)/F/beta/beta+efdiff/F
      dde=-efdiff/F+(XDK-beta*D)/F/beta/beta
!       print *,'chi=',chi,' (d-1)/eps=',(dd-1.)*dexp(b),' dde=',dde
! rearrange Cummings & Smith expression for aa
      aa=(1.+2.*vf+12.*vf*beta*((1.+2.*vf-6.*vf/b)*
     ^     (-dde-dd*(1.+b))+3.*dd*vf*b)/b)/(1.-vf)/(1.-vf) 
!      print *,'aa=',aa
      if(aa.lt.0.0d0) print *,'warning, probably inside spinodal'
      bb=-(1.5*vf+12.*vf*beta*((1.5*vf+(1.-4.*vf)/b)*
     ^     (-dde-dd*(1.+b))-.5*(1.-4.*vf)*dd*b)/b)/(1.-vf)/(1.-vf)
!       print *,'bb=',bb
!       print *,'dde=',dde
!       print *,'compr=',1./aa/aa
      ENDIF
      return
      end

      subroutine funcd(x,f,df)
      double precision x,f,df
      double precision W0,W1,W2,W3,W4
      COMMON/CON/W0,W1,W2,W3,W4
      f=W4*X*X*X*X+W3*X*X*X+W2*X*X+W1*X+W0
      df=4.*W4*X*X*X+3.*W3*X*X+2.*W2*X+W1
      return
      end
      
      function fx(x)
      double precision x,fx,W0,W1,W2,W3,W4
      COMMON/CON/W0,W1,W2,W3,W4
      fx=W4*X*X*X*X+W3*X*X*X+W2*X*X+W1*X+W0
      return
      end
      
      function rtsafe(x1,x2,xacc)
      integer maxit
      double precision rtsafe,x1,x2,xacc
      external funcd
      parameter (maxit=3000)
      integer j
      double precision df,dx,dxold,f,fh,fl,temp,xh,xl
      call funcd(x1,fl,df)
      call funcd(x2,fh,df)
      if((fl.gt.0.d0.and.fh.gt.0.d0).or.
     ^       (fl.lt.0.d0.and.fh.lt.0.d0))
     ^    pause 'root must be bracketed in rtsafe'
      if(fl.eq.0.d0) then
        rtsafe=x1
        return
      else if(fh.eq.0.d0) then
        rtsafe=x2
        return
      else if(fl.lt.0.d0) then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
      rtsafe=.5*(x1+x2)
      dxold=dabs(x2-x1)
      dx=dxold
      call funcd(rtsafe,f,df)
      do 11 j=1,maxit
      if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0.d0
     ^        .or.abs(2.*f).gt.dabs(dxold*df) ) then
        dxold=dx
        dx=.5*(xh-xl)
        rtsafe=xl+dx
        if(xl.eq.rtsafe) return
      else
        dxold=dx
        dx=f/df
        temp=rtsafe
        rtsafe=rtsafe-dx
        if(temp.eq.rtsafe) return
      endif
      if(abs(dx).lt.xacc) return
      call funcd(rtsafe,f,df)
      if(f.lt.0.d0) then
        xl=rtsafe
      else
        xh=rtsafe
      endif
   11 enddo
      pause 'rtsafe exceeding maxits'
      return
      end
      
      subroutine zbrak(x1,x2,n,xb1,xb2,nb)
      integer n,nb
      double precision x1,x2,xb1(nb),xb2(nb),fx
      external fx
      integer i,nbb
      double precision dx,fc,fp,x
      nbb=0
      x=x1
      dx=(x2-x1)/n
      fp=fx(x)
      do i=1,n
        x=x+dx
	fc=fx(x)
	if(fc*fp.lt.0.d0) then
	  nbb=nbb+1
	  xb1(nbb)=x-dx
	  xb2(nbb)=x
	  if(nbb.eq.nb) goto 1
	endif
	fp=fc
      enddo
    1 continue
      nb=nbb
      return
      end

        

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                  C
C Hasta aqui todo lo necesario para calcular S(k).                 C
C Las siguentes funciones se utilizan en el criterio               C
C                                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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

       smax=Sk(251,2)
               
       Do j=252, N
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
      
      SUBROUTINE SIMPSON(XIN,DR,M,RES1)
C     ****************************************************
C     SUBRUTINA DE INTEGRACION POR SIMPSON
C     ****************************************************

	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION XIN(M)
	
	RES1=0.D0
	S0=0.D0
	S1=0.D0
	S2=0.D0
   
	DO 100 I=2,M-1,2
	S0=S0+XIN(I-1)
	S1=S1+XIN(I)
	S2=S2+XIN(I+1)
100	CONTINUE
	RES1=DR*(S0+4.D0*S1+S2)/3.D0
        
	IF(MOD(M,2).EQ.0) RES1=RES1+DR*(5.D0*XIN(M)+8.D0*XIN(M-1)-
     #  XIN(M-2))/12.D0
	RETURN
	END
	
      Subroutine  facnoerg(gamma,Sk,kmin,nk,f)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C																	C
C Esta subrutiuna evalua los factores no ergodicos predichos por el	C
C criterio de atrapamiento.											C
C Parametros: *gamma: valor del parametro gamma del criterio			C
C	      *Sk:factor de estructura									C
C	      *kmin:posicion del primer minimo de S(k)					C
C	      *nk:numero de puntos en k									C
C	      *f:arreglo de salida										C
C																	C
C Agosto de 2006														C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	 Implicit none !toda variable se tiene que declarar
CC
C Parametros
CC
	 Integer nk
	 Double precision kmin,gamma
	 Double precision Sk(nk,2),f(nk,2)
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
	  den=Sk(i,2)*l
	  num=gamma*(Sk(i,1))**2
	  f(i,1)=Sk(i,1)
	  f(i,2)=1.D0/(1.D0+ (num/den))
	 ENDDO
	 Return
      End

       Double precision function Smax(Sk,N)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C							             C
C Esta funcion sirve para calcular el maximo de S(k), pero           C
C puede adapatarse para otra funcion solo hay que ponerla en un      C
C arreglo bidimensional                                              C
C Parametros: *Sk: arreglo en que esta S(k)                          C
C                 (primera columna k y segunda S(k))                 C
C	      *N: Tamaño del arreglo                                C
C								     C
C Agosto del 2006                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        Implicit none
       
        Integer N
        Integer j
        Double precision Sk(N,2)
        Smax=Sk(251,2)
        Do j=252, N
         if (Sk(j,2).GT.smax)then
          Smax=Sk(j,2)
         endif
        enddo
        Return
       End
