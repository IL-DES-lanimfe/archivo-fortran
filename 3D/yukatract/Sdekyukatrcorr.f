      PROGRAM yukawaatr
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Este programa sirve para probar al subrutina de S(k) para el poten-C
C cial de yukawa atractivo.					     C
C								     C
C Enero 2007							     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none !toda variable se tiene que declarar
	Integer Nk,i,kpeq
	parameter (Nk=2**12)
	Double precision fv,Ka,z
	Double precision dk,Rk(Nk),Sk(Nk)
	double precision xi,pi,fvw,kdw,DelC
	double precision Spy(Nk,2),Svw(Nk,2)
CC
C Parametros
CC
	pi=dacos(-1.D0)
	fv=0.563D0
	fvw=fv-fv**2/16.D0
	Ka=10.D0
	z=12.5D0
	dk=0.01D0
CC
C Ahora llamamos a la subrutina de S(k)
CC
	call msa(fv,Ka,z,dk,Nk,Rk,Sk,xi)
CC
C correccion de Verlet-Weis y salida
CC
	call Sdeked(Spy,fv,Nk,81.9D0,10.D0)
        call Sdeked(Svw,fvw,Nk,81.9D0,10.D0)
	do i=1,Nk
! 	 kdw=2.D0*Rk(i)*(fvw/fv)**(1.D0/3.D0)
	 DelC=(Svw(i,2)-Spy(i,2))/(Svw(i,2)*Spy(i,2))
! 	 if(1.D0-DelC*Sk(i).lt.3.D0)DelC=0.D0
! 	 Sk(i)=Sk(i)/(1.D0-DelC*Sk(i))
	 write(50,*)2.D0*Rk(i),Sk(i)
	 write(60,*)2.D0*Rk(i),Sk(i)/(1.D0-DelC*Sk(i))
	enddo
      END
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
!        parameter(Npg=50)
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
!       Double precision function DELTAC(y,phi)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C								     C
! C Esta funcion evalua la expresion analítica de la transformada de   C
! C Fourier de la funcion de correlación directa para el caso de esfe- C
! C ras duras.							     C
! C								     C
! C Febrero 2007							     C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 	implicit none !toda variable se tiene que declarar
! 	Double precision y,phi
! 	Double precision T1,T2,numT2,denT2,Tkcuad,Tkcero
! 	Double precision Tcos,Tsen,T1DC,T2DC
! 	Double precision pi
! CC
! 	pi=dacos(-1.D0)
! CC
! C Expansion a k's pequeñas
! CC
! 	if(y.le.0.5D0) then
! CC
! C coeficiente del termino de orden k^2
! CC
! 	 T1=(-1.D0/15.D0+phi/40.D0+phi**3/24.D0)/(-1.D0+phi)**4
! 	 numT2=-16.D0*(-4096.D0/15.D0+512.D0*phi/5.D0-32.D0*phi**2/5.D0
!      #   +512.D0*phi**3/3.D0-32.D0*phi**4+2.D0*phi**5-phi**6/24.D0)
! 	 denT2=(16.D0+(-16.D0+phi)*phi)**4
! 	 T2=numT2/denT2
! 	 Tkcuad=T1+T2
! CC
! C coeficiente termino de orden k^0
! CC
! 	 T1=(2.D0/3.D0-phi/5.D0+phi**2/5.D0-13.D0*phi**3/60.D0)/
!      #      (-1.D0+phi)**4
! 	numT2=-16.D0*(8192.D0/3.D0-4096.D0*phi/5.D0+4352.D0*phi**2/5.D0
!      #   -14848.D0*phi**3/15.D0+848.D0*phi**4/5.D0-52.D0*phi**5/5.D0
!      #   +13.D0*phi**6/60.D0)
! 	 T2=numT2/denT2
! 	 Tkcero=T1+T2
! CC
! C Valor de deltaC
! CC
! 	 DELTAC=(Tkcuad*y**2+Tkcero)/dsqrt(2.D0*pi)
! 	 return
! 	Endif
! CC
! C Para toda k
! CC
! 	T1=-6.D0*phi*(y**2*(2.D0+phi)**2+4.D0*(1.D0+2.D0*phi)**2)
! 	T2=24.D0*phi*(1.D0+2.D0*phi)**2+y**4*(2.D0-3.D0*phi+phi**3)
!      #  -6.D0*y**2*phi*(-2.D0+phi*(4.D0+7.D0*phi))
! 	Tcos=-(T1+T2)*dcos(y)/(-1.D0+phi)**4
! 	Tsen=2.D0*y*(12.D0*phi*(1.D0+2.D0*phi)**2+y**2*(-1.D0+6.D0*phi
!      #  -5.D0*phi**3))*dsin(y)/(-(-1.D0+phi)**4)
! 	T1DC=(Tcos+Tsen)/(y**6*dsqrt(2.D0*pi))
! 
! 	T1=16.D0*(-6.D0*(-16.D0+phi)*phi*(y**2*(-32.D0+(-16.D0+phi)*phi
!      #  )*phi+16.D0*(-8.D0+(-16.D0+phi)*phi)**2))
! 	T2=16.D0*((96.D0*(-16.D0+phi)*phi*(-8.D0+(-16.D0+phi)*phi)**2
!      #  +y**4*(-32.D0+(-16.D0+phi)*phi)*(16.D0+(-16.D0+phi)*phi)**2-
!      #  6.D0*y**2*(-16.D0+phi)*phi*(-512.D0+(-16.D0+phi)*phi*(-64.D0+
!      #  7.D0*(-16.D0+phi)*phi))))
! 	Tcos=-(T1+T2)*dcos(y)/(16.D0+(-16.D0+phi)*phi)**4
!       Tsen=16.D0*(-2.D0*y*(-48.D0*(-16.D0+phi)*phi*(-8.D0+(-16.D0+phi)
!      #*phi)**2+y**2*(16.D0+(-16.D0+phi)*phi)**2+y**2*(16.D0+(-16.D0
!      #+phi)*phi)*(-256.D0+5.D0*(-16.D0+phi)*phi*(-16.D0+(-16.D0+phi)
!      #*phi)))*dsin(y))/(-(16.D0+(-16.D0+phi)*phi)**4)
! 	T2DC=(Tcos+Tsen)/(y**6*dsqrt(2.D0*pi))
! CC
! C valor de deltaC
! CC
! 	DELTAC=T1DC+T2DC
! 	return
!       END	

CC
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
      print *,'MSA SOLUTION:'
      print *,'beta=',beta
      print *,'d=',dd,' dde=',dde
      print *,'vf=',vf
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
      print *,'W0,W1,W2,W3,W4=',W0,W1,W2,W3,W4
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
      print *,'root1=',root1
! ATTEMPT AT SECOND ROOT:      
      root2=rtsafe(xb1(2),xb2(2),1.D-10)
      print *,'root2=',root2
! IDENTIFY PROPER ROOT (LOWER):
      if(root1.lt.root2) then
        beta=root1
      else
        beta=root2
      endif
      print *,'beta=',beta
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
      print *,'chi=',chi,' (d-1)/eps=',(dd-1.)*dexp(b),' dde=',dde
! rearrange Cummings & Smith expression for aa
      aa=(1.+2.*vf+12.*vf*beta*((1.+2.*vf-6.*vf/b)*
     ^     (-dde-dd*(1.+b))+3.*dd*vf*b)/b)/(1.-vf)/(1.-vf) 
      print *,'aa=',aa
      if(aa.lt.0.0d0) print *,'warning, probably inside spinodal'
      bb=-(1.5*vf+12.*vf*beta*((1.5*vf+(1.-4.*vf)/b)*
     ^     (-dde-dd*(1.+b))-.5*(1.-4.*vf)*dd*b)/b)/(1.-vf)/(1.-vf)
      print *,'bb=',bb
      print *,'dde=',dde
      print *,'compr=',1./aa/aa
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
      
