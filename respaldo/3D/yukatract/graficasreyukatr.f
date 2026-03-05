      Program graficasdz

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C 
C Este programa evalua la integral de la funcional del criterio de  C 
C atrapamiento para el caso de "restauracion de la ergdicidad".Esta C
C diseñado para dos parametros de control (por ejemplo frac. de vol.C 
C y temperatura), pero se puede expandir facilmente para mas para-  C
C metros.Para cambiar de sistema solo es necesario cambiar subrutinaC
C y variables para el factor de estructura. Lo he transformado para C
C que calcule todo en base al parametro gamma que se esta utilizandoC
C en el articulo sobre transicion vitrea que se va a publicar y asi C
C estandarizar nuestra notación. Recordando que gamma=1/dzeta, ya   C
C la dzeta que se usa en este programa ya esta adimensionalizada.   C
C Sin embargo por comodidad se sigue llamando dzeta pero teniendo laC
C reserva que nos referimos a gamma.				    C
C Adaptado a Yukawas atractivos					    C
C								    C
C Enero de 2007                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        Implicit none
C indices para los ciclos en general
        integer iint,ig0,i 
C variables para la integral sobre vectores de onda
        double precision  y,dy,S  
C variables para evaluar la funcional del criterio
        double precision dzeta,ddzeta,dzetamax,dzetamin !gamma inversa
C numero de puntos a evaluar en dzeta
	integer Ndzeta
C variables para las integrales de funcional (temporales)
        double precision num,den
C variables las evaluaciones de la funcional del criterio
        double precision integral
C variables para parametros del criterio (labda,gcero,1er min S(k))
        double precision l,g,kmin
C Numero de puntos en vecs. de onda
	Integer nk
	parameter (nk=2**12)
C arreglos para S(k) y g(r)
        double precision  Sdek(nk,2),T(nk),g0(nk)
C parametros S(k) y g(r)
	Double precision parametro1,parametro2 
C variables para evaluar en los parametros de control
	double precision par2max,par2min,dpar2
C variables auxiliares
        Double precision pi,fv,z,xi
	Double precision Rk(nk),Sk(nk)
CC
CCCCCCCCCCCCCCCC Funciones que utiliza el programa CCCCCCCCCCCCCCCCCCC
CC
        double precision  lambda,kminima

CCCCCCCCCCCCCCCCCCC Archivos de salida CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C archivo de salida funcional criterio 
        open(9,file='funcritre.dat')
	open(10,file='solucionesre.dat')
CC
        pi= dacos(-1.D0)
CC
C saludo de beinvenida y entrada de datos
CC
	write(*,*)' '
	write(*,*)'Bienvenido a graficas, vamos a calcular'
	write(*,*)'las integrales del criterio y de la funcional de'
	write(*,*)'gamma para asi tener una idea de donde esta la'
	write(*,*)'transición vítrea. Necesito unos parametros'
	write(*,*)' '
	write(*,*)'Es el caso de pot. de yukawa atractivo'
	write(*,*)'parametro 1: fraccion de volumen'
	write(*,*)'parametro 2: K (inverso de la temperatura)'
	write(*,*)' '
CC
C El programa te hace las graficas de la funcional del criterio como 
C funcion del parametro 2 manteniendo el parametro 1 fijo, es decir 
C cada corrida del programa te genera una iso-parametro1 por decirlo 
C de alguna manera.
CC
        write(*,*)'¿cuanto vale el parametro 1?'
        read (*,*)parametro1
        write(*,*)'¿entre que intervalo evaluamos el parametro 2?'
        write(*,*)'¿maximo?'
        read (*,*)par2max
	write(*,*)'¿minimo?'
        read (*,*)par2min
	write(*,*)'¿Longitud del incremento en el parametro 2?'
        read (*,*)dpar2
! 	write(*,*)'¿cuantos puntos quieres en la fun. del criterio?'
!         read (*,*)Ndzeta
CC
C Empieza el ciclo para barrer sobre parametro 2
CC
	fv=parametro1
	z=12.5D0
	
	DO parametro2=par2min,par2max,dpar2
          write(*,*)'parametros',parametro1,parametro2
CC
CCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
	  write(*,*)'calculado S(k)'
CC Aqui se inserta la subrutina del factor de estrucutura
CC
	 call msa(parametro1,parametro2,z,0.01D0,nk,Rk,Sk,xi)
	 do i=1,nk
	  Sdek(i,1)=2.D0*Rk(i)
	  Sdek(i,2)=Sk(i)
	  write(50,*)Sdek(i,1),Sdek(i,2)
	 enddo
C se concoce S(k) y g(r)
CC

	  dy=Sdek(2,1)-Sdek(1,1)
          write(*,*)'calculando la k del primer minimo de S(k)'
          kmin=kminima(SdeK,nk)
	  write(*,*)'kmin',kmin
CC
C ahora se conoce la posicion del primer minimo de S(k)
CC
  	  write(*,*)'evaluando la funcion g-cero'
	  do ig0=1,nk
	   g0(ig0)=1.D0
CC
C la definicion de g-cero depende de la vineyard, para variar gcero se
C puede hacer una funcion
CC
C          y=SdeK(ig0)
C          g0(ig0)=gcero(y,parametro1,parametro2)
	  enddo
CC
CCCCCCCCCCCCCCCCCCCCCC Fin modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCC
CC
	  write(*,*)'calculando la funcional del criterio'
CC
C dependiendo del caso
CC
!	    fv=parametro2
	    dzetamin=0.D0
	    dzetamax=0.1D0
	    ddzeta=1.D-5
            Do dzeta=dzetamin,dzetamax,ddzeta
              integral= 0.D0
CC
C evaluacion del integrando de la funcional del criterio
CC
              do iint=1, nk
		 y=Sdek(iint,1)
                 S=SdeK(iint,2)
                 g= g0(iint)
                 l= lambda(y,kmin)
                 num= g*dzeta*((S-1.D0)*l*y**2)**2
                 den= 36.D0*pi*fv*(dzeta*y**2+S*g*l)*(dzeta*y**2+l)
                 T(iint)= num/den
              enddo
CC
C evaluacion de la integral
CC
              call SIMPSON(T,dy,nk,integral)
	      write(9,*)dzeta,1.D0/dzeta,integral
CC
C buscamos soluciones dentro del lado vitreo
CC
	if((dzeta.le.0.01).and.(dabs(1.D0-integral).le.1.D-4))then
	 write(10,*)parametro2,dzeta,1.D0/dzeta
	endif
	if((dzeta.gt.0.01).and.(dabs(1.D0-integral).le.1.D-5))then
	 write(10,*)parametro2,dzeta,1.D0/dzeta
	endif
CC
CC Termina la evaluacion de la  integral de la funcional de gamma
CC
C aumentamos gamma para el siguiente punto
CC
          EndDo
          write(*,*)'funcional del criterio lista',parametro2
CC Termina el ciclo de evaluacion de la funcional
       ENDDO
CC Termina ciclo de parametros
       close(9)
       close(10)
CCC Fin del programa CCCC
       stop
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                  C
C Subrutinas del factor de estructura, proporcionadas 		   C
C 	 							   C
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
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    C
C Hasta aqui factor de estrucutura				    C
C							            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                    C
C Las siguentes funciones se utilizan en el criterio C
C                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double precision function kminima(Sk,N)

       Implicit none
       
       Integer N
C       Double precision emp

       Integer h,m,j
       Double precision Sk(N,2)
C       Double precision y,dy
       Double precision smin,smax,kmin,kmax

       Double precision facdes

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

       Do h=m , N

          if (Sk(h,2).LT.smin)then

             smin= Sk(h,2)
             kmin= Sk(h,1)

          endif

       enddo

       kminima=kmin
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
