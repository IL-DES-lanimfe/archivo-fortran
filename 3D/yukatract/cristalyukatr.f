       Program cristal

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C Este programa aplica el criterio de Verlet para determinar el puntoC
C de cristalizacion en un liquido.Esta diseñado para dos parametros  C
C de control,la aplicacion del criterio a sistemas diferentes se da  C
C cambiando el factor de  estructura para el caso especifico que     C
C estemos tratando.                                                  C
C								     C
C Septiembre 2006                                                    C
C								     C
C Adaptacion a Yukawas atractivas				     C
C								     C
C Enero 2007
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
	Double precision solpost,solant,solmed
	Double precision solprod,prodant,prodpost
	Double precision Smx
CC
C variables de estatica 
CC
	Integer Nk
	parameter(Nk=2**12)
	Double precision SdeK(Nk,2)
CC
C variables axiliares
CC
	Integer i
	Double precision pi,z,xi
	Double precision Rk(Nk),Sk(Nk)
CC
C funciones que se usan
CC
	Double precision Smax
CC
C Saludo inicial y lectura de datos
CC
	write(*,*)' '
	write(*,*)'Este programa es para generar un diagrama de fases'
	write(*,*)'de cristalizacion con dos parametros'
	write(*,*)'Este es el caso del pot. de Yukawa atractivo'
        write(*,*)'parametro1:fraccion de volumen'
	write(*,*)'parametro2:K(inverso de temperatura)'
	write(*,*)' '
 	write(*,*)'establezca el intervalo de busqueda del parametro 1'
	write(*,*)'maximo?'
	read(*,*)par1max
	write(*,*)'minimo?'
	read(*,*)par1min
	write(*,*)'establezca el intervalo de busqueda del parametro 2'
	write(*,*)'maximo?'
	read(*,*)par2max
	write(*,*)'minimo?'
	read(*,*)par2min
! 	write(*,*)'numero de puntos?'
! 	read(*,*)Npar2
CC
C ya entraron todos los datos obligatorios
CC	
C Si ninguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente
CC
C	 write(*,*)'¿fraccion de volumen?'
C	 read(*,*)fv
CC
	write(*,*)' '
	write(*,*)'ya rugistes leon, ahora a trabajar'
	write(*,*)' '
	pi= dacos(-1.D0)
	z=100.D0
	dpar2=0.1D0 !dabs(par2max-par2min)/dble(Npar2)
        open(8,file='cristalyukatrz100.dat') 
      DO parametro2=par2min, par2max,dpar2
CC
C se usa el metodo de bisecciones sucesivas. Parametro 2 fijo y 
C bisecciones sucesivas en parametro 1 para localizar la transicion
C luego se cambia el parametro 2 y se repite la operacion
CC
	 par1post=par1max
       par1ant= par1min
CCCCCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	 write(*,*)'calculando S(k) cota inferior'
CC
C ponga aqui la subrutina de S(k) evaluandola en par1ant y parametro2
CC
	 call msa(par1ant,parametro2,z,0.01D0,Nk,Rk,Sk,xi)
	 do i=1,Nk
	  Sdek(i,1)=Rk(i)
	  Sdek(i,2)=Sk(i)
	 enddo
CC
       write(*,*)'calculando el maximo de S(k)'
       Smx=Smax(Sdek,Nk)
	 solant=2.85D0-Smx
	 write(*,*)'maximo de S(k)',Smx,par1ant
CC
C Criterio de convergencia
CC
 	 if (dabs(solant).lt.1.D-4)then
          write(8,*)par1ant,parametro2,Smx
          write(*,*)'se encontro solucion para: ',par1ant,parametro2
          goto 200
         endif	
CC
C fin criterio de convergencia
CC
CCCCCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         write(*,*)'calculando S(k) en la cota superior'
CC
C ponga aqui la subrutina de S(k) evaluandola en par1post y parametro2
CC
	 call msa(par1post,parametro2,z,0.01D0,Nk,Rk,Sk,xi)
	 do i=1,Nk
	  Sdek(i,1)=Rk(i)
	  Sdek(i,2)=Sk(i)
	 enddo
CC
	write(*,*)'calculando el maximo de S(k)'
        Smx=Smax(Sdek,Nk)
	solpost=2.85D0-Smx
	write(*,*)'maximo de S(k)',Smx,par1post
CC
C Criterio de convergencia
CC
        if (dabs(solpost).lt.1.D-4)then
          write(8,*)par1post,parametro2,Smx
          write(*,*)'se encontro solucion para: ',par1post,parametro2
          goto 200
         endif 
CC
C Termina la evaluacion de la  integral para el valor del parametro1
C mas grande
CC
         solprod=solant*solpost
         IF(solprod.lt.0.D0)then
	  write(*,*)'bisecciones sucesivas'
 100	  continue
          par1med=(par1post+par1ant)/2.D0
CCCCCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  write(*,*)'calculando S(k) en el valor intermedio'
CC
C llame aqui a la subrutina de S(k) evaluandola en par1med y parametro2
CC
	  call msa(par1med,parametro2,z,0.01D0,Nk,Rk,Sk,xi)
	  do i=1,Nk
	   Sdek(i,1)=Rk(i)
	   Sdek(i,2)=Sk(i)
	  enddo
CC
	  write(*,*)'calculando el maximo de S(k)'
          Smx=Smax(Sdek,Nk)
	  solmed=2.85D0 - Smx
	  write(*,*)'maximo de S(k)',Smx,par1med
CC
C Criterio de convergencia
CC
          if (dabs(solmed).lt.1.D-4)then
           write(8,*)par1med,parametro2,Smx
           write(*,*)'se encontro solucion para: ',par1med,parametro2
           goto 200
          endif
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
 200     continue
        ENDDO 
        close(8)
	stop
      End

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                  C
C Subrutina del factor de estructura 			           C
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
CC
	Smax=Sk(251,2)
        Do j=252, N
         if (Sk(j,2).GE.Smax)then
          Smax=Sk(j,2)
         endif
        enddo
        Return
       End

