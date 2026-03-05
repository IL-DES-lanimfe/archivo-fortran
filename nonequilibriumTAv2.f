      PROGRAM GENERALRELAJA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Este programa utiliza la solucion para la relajacion temporal del 	C
C factor de estructura estatico usando los conceptos del proceso de 	C
C Ornstein-Uhlenbeck aplicados a la ecuacion de difusion.		C
C									C
C Proceso de enfriamiento con yukawa atractiva				C
C									C
C Febrero 2009								C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none !toda variable se tiene que declarar
	Integer i,Ndk,nsteps,dec,dec2,ikmax,ikshow,ikshow2,Lmaximo
	Integer step,wrt
	parameter(Ndk=2**12)
	Double Precision ntaus,fv,fvi,fvf,Ti,Tf,Ka,t,chhi
	Double Precision Rk(Ndk),Sk(Ndk),Sdeki(Ndk,2),Sdekf(Ndk,2)
	Double precision tau(Ndk),Sdekt(Ndk),D,DLT
	Double precision Tp,fvp,delT,deltem,st,Smx
	Double precision crit,Criterio,delzfin(7000,2)
	write(*,*)'tiempo del experimento?(cuantas veces el tiempo de relax)'
	read(*,*)ntaus
	write(*,*)'numero de pasos?'
	read(*,*)nsteps
	write(*,*)'enfriamiento(1) o compactacion(2)?'
	read(*,*)dec
	write(*,*)'con memoria (1) o sin memoria (2)?'
	read(*,*)dec2
	if (dec.eq.1)then
	 write(*,*)'fraccion de volumen?'
	 read(*,*)fv
	 write(*,*)'temperatura inicial?'
	 read(*,*)Ti
	 write(*,*)'temperatura final?'
	 read(*,*)Tf
	endif
	if (dec.eq.2)then
	 write(*,*)'Temperatura?'
	 read(*,*)Ti
	 write(*,*)'fraccion de volumen inicial?'
	 read(*,*)fvi
	 write(*,*)'fraccion de volumen final?'
	 read(*,*)fvf
	endif

CC
C Relajacion en un solo paso 
CC
C factor de estructura inicial
CC
	if(dec.eq.1) Ka=1.D0/Ti
	if (dec.eq.2)then
	 Ka=1.D0/Ti
	 fv=fvi
	endif
	call msa(fv,Ka,20.D0,0.1D0,Ndk,Rk,Sk,chhi)
	do i=1,Ndk
	 Sdeki(i,1)=2.D0*Rk(i)
	 Sdeki(i,2)=Sk(i)
	enddo
	ikmax=Lmaximo(Sdeki,NdK)
CC
C factor de estructura final
CC
	if(dec.eq.1) Ka=1.D0/Tf
	if (dec.eq.2)then
	 Ka=1.D0/Ti
	 fv=fvf
	endif
	call msa(fv,Ka,20.D0,0.1D0,Ndk,Rk,Sk,chhi)
	do i=1,Ndk
	 Sdekf(i,1)=2.D0*Rk(i)
	 Sdekf(i,2)=Sk(i)
! 	 write(50,*)Sdekf(i,1),Sdekf(i,2)
CC
C tiempo de relajacion total
CC
	 tau(i)=Sdekf(i,2)/(2.D0*Sdekf(i,1)**2)
	enddo
! 	stop
CC
C ciclo temporal en un solo paso
CC
	if (dec.eq.1) then
! 	 open(15,file='Sdekttaus1000.dat')
	 open(20,file='relaxkmaxtaus5000a.dat')
	 ikshow=ikmax
	endif
	if (dec.eq.2)then
! 	 open(15,file='SdektCtaus1000.dat')
	 open(20,file='relaxCkmaxtaus5000a.dat')
	 ikshow=ikmax
	endif
	do t=0.D0,ntaus*tau(ikshow),tau(ikshow)/100.D0
CC
C ciclo en vectores de onda
CC
! 	 if(mod(t,10.D0*tau(ikshow)).eq.0.D0)write(15,*)'tiempo',t
	 ikshow2=0
	 Smx=0.D0
	 do i=1,Ndk
	  Sdekt(i)=Sdeki(i,2)*dexp(-t/tau(i))
     %               +(1.D0-dexp(-t/tau(i)))*Sdekf(i,2)
! 	  if(mod(t,10.D0*tau(ikshow)).eq.0.D0)then
! 	    write(15,*)Sdeki(i,1),Sdekt(i)
! 	  endif
	  if((i.gt.1).and.(Sdekt(i).gt.Smx))then
	   ikshow2=i
	   Smx=Sdekt(i)
	  endif
	 enddo
CC termina ciclo vectores de onda
	 write(*,*)t,ikshow2,ikshow
	 write(20,*)t,Sdekt(ikshow2)
	enddo
	close(20)
	close(15)
CC termina ciclo temporal un solo paso
CC
C Relajacion en pasos
CC
	if(dec.eq.1) delT=(Tf-Ti)/dble(nsteps)
	if(dec.eq.2) delT=(fvf-fvi)/dble(nsteps)
	deltem=ntaus*tau(ikshow)/dble(nsteps)
CC
C Relajacion en pasos
CC
	if(dec2.eq.1)then
	 if(dec.eq.1)then
	  Tp=Ti+delT
	  open(25,file='Sdektau1pasostaus5000TAa.dat')
	  open(30,file='relax1pasoskmaxtaus5000TAa.dat')
	 endif
	 if(dec.eq.2)then
	  fvp=fvi+delT
	  open(25,file='SdektauC1pasostaus5000TAa.dat')
	  open(30,file='relaxC1pasoskmaxtaus5000TAa.dat')
	 endif 
	endif
CC
	if(dec2.eq.2)then
	 if(dec.eq.1)then
	  Tp=Ti+delT
	  open(25,file='Sdektau1pasostaus5000a.dat')
	  open(30,file='relax1pasoskmaxtaus5000a.dat')
	 endif
	 if(dec.eq.2)then
	  fvp=fvi+delT
	  open(25,file='SdektauC1pasostaus5000a.dat')
	  open(30,file='relaxC1pasoskmaxtaus5000a.dat')
	 endif 
	endif
CC
C ciclo de las terrazas
CC
	step=0
CC
C Proceso de enfriameinto
CC
	if(dec.eq.1)then
 200	 continue
	 if(Tp.lt.Tf) goto 100
	  write(*,*)'nueva terraza',Ti,'-',Tp
	  step=step+1
	  if(mod(step,10).eq.0)then
	   wrt=1
	  else
	   wrt=0
	  endif
	  if(step.eq.1)wrt=1
	  if(dec2.eq.1)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C modulo de teoria autoconsistente C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C condicion de arresto dinamico
CC
	   crit=Criterio(Sdeki,Ndk,fv)
	   if(crit.ge.1.D0)then
	    write(*,*)'arresto'
	    goto 1000
	   endif
CC
C solucion completa
CC	
	   call SCGLE(Sdeki,fv,delzfin,wrt)
CC
C Calculo de coeficiente de difusion a tiempos largos
CC
	   D=DLT(delzfin)
	  endif
CC
	  if(dec2.eq.2) D=1.D0
CC
C factor de estructura final
CC
	  Ka=1.D0/Tp
	  call msa(fv,Ka,20.D0,0.1D0,Ndk,Rk,Sk,chhi)
	  do i=1,Ndk
	   Sdekf(i,1)=2.D0*Rk(i)
	   Sdekf(i,2)=Sk(i)
CC
C tiempo de relajacion total
CC
	   tau(i)=Sdekf(i,2)/(2.D0*D*Sdekf(i,1)**2)
	  enddo
CC
C ciclo temporal en un solo paso
CC
	  do t=0.D0,deltem,deltem/100.D0
CC
C ciclo en vectores de onda
CC
	   if(wrt.eq.1)write(25,*)'tiempo',t,st
	   ikshow2=0
	   Smx=0.D0
	   do i=1,Ndk
	    Sdekt(i)=Sdeki(i,2)*dexp(-t/tau(i))
     %               +(1.D0-dexp(-t/tau(i)))*Sdekf(i,2)
	    if(wrt.eq.1)then
	     write(25,*)Sdeki(i,1),Sdekt(i)
	    endif
	    if((i.gt.1).and.(Sdekt(i).gt.Smx))then
	     ikshow2=i
	     Smx=Sdekt(i)
	    endif
	   enddo
CC termina ciclo vectores de onda
	   write(*,*)t,ikshow2,ikshow
	   write(30,*)st+t,(st+t)/(ntaus*D*tau(ikshow)),Sdekt(ikshow2)
	  enddo
CC  termina ciclo temporal en la terraza
	  st=st+deltem
	  do i=1,Ndk
	   Sdeki(i,2)=Sdekt(i)
	  enddo
	  Tp=Tp+delT
	  goto 200
 100	 continue
	endif
CC
C fin proceso de enfriamiento
CC
C Compactacion
CC

	if(dec.eq.2)then
 400	 continue
	 if(fvp.gt.fvf) goto 300
	  write(*,*)'nueva terraza',fvi,'-',fvp
	  step=step+1
	  if(mod(step,10).eq.0)then
	   wrt=1
	  else
	   wrt=0
	  endif
	  if(dec2.eq.1)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C modulo teoria autoconsistente  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C criterio de arresto dinamico
CC
	   crit=Criterio(Sdeki,Ndk,fvi)
	   if(crit.ge.1.D0)then
	    write(*,*)'arresto'
	    goto 1000
	   endif
CC
C solucion completa
CC
	   call SCGLE(Sdeki,fvi,delzfin,wrt)
CC
C calculo del coeficiente de difusion a tiempos largos
CC
	   D=DLT(delzfin)
	  endif
CC
	  if(dec2.eq.2) D=1.D0
CC
C factor de estructura final
CC
	  Ka=1.D0/Ti
	  fv=fvp
	  call msa(fv,Ka,20.D0,0.1D0,Ndk,Rk,Sk,chhi)
	  do i=1,Ndk
	   Sdekf(i,1)=2.D0*Rk(i)
	   Sdekf(i,2)=Sk(i)
CC
C tiempo de relajacion total
CC
	   tau(i)=Sdekf(i,2)/(2.D0*D*Sdekf(i,1)**2)
	  enddo
CC
C ciclo temporal en un solo paso
CC
	  do t=0.D0,deltem,deltem/100.D0
CC
C ciclo en vectores de onda
CC
	   if(wrt.eq.1) write(25,*)'tiempo',t,st
	   do i=1,Ndk
	    Sdekt(i)=Sdeki(i,2)*dexp(-t/tau(i))
     %               +(1.D0-dexp(-t/tau(i)))*Sdekf(i,2)
	    if(wrt.eq.1)then
	     write(25,*)Sdeki(i,1),Sdekt(i)
	    endif
	   enddo
CC termina ciclo vectores de onda
	   write(*,*)t
	   write(30,*)st+t,(st+t)/(ntaus*tau(ikshow)),Sdekt(ikshow)
	  enddo
CC  termina ciclo temporal en la terraza
	  st=st+deltem
	  do i=1,Ndk
	   Sdeki(i,2)=Sdekt(i)
	  enddo
	  fvi=fvp
	  fvp=fvp+delT
	  goto 400
 300	 continue
	endif
	close(30)
 1000	continue
	stop
      END	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C factor de estructura yukawa						C
C									C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    	C
C Hasta aqui factor de estrucutura				    	C
C							            	C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Integral del criterio							C
C									C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double precision Function Criterio(Sdek,nk,fv)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C	 
C Este programa evalua la integral de la funcional del criterio de	C	 
C atrapamiento para el caso de "restauracion de la ergdicidad".Esta	C
C hecho para dos parametros de control (por ejemplo frac. de vol.	C 
C y temperatura), pero se puede expandir facilmente para mas para-	C
C metros.Para cambiar de sistema solo es necesario cambiar subrutina	C
C y variables para el factor de estrucutura.				C	  
C									C
C Agosto de 2006							C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 	Implicit none
C Numero de puntos en vecs. de onda
	Integer nk,IKMI
C arreglos para S(k) y g(r)
	double precision  Sdek(nk,2),fv,T(nk)
C variables para evaluar el criterio
	double precision gama,gamamin,gamamax,dgama,maxcrit
C indices para los ciclos en general
	integer iint 
C variables para la integral sobre vectores de onda
	double precision  y,dy,S  
C variables para las integrales de funcional (temporales)
	double precision num,den
C variables las evaluaciones de la funcional del criterio
	double precision integral
C variables para parametros del criterio (labda,1er min S(k))
	double precision l,kmin
C variables auxiliares
	Double precision pi 
CC
CCCCCCCCCCCCCCCC Funciones que utiliza el programa CCCCCCCCCCCCCCCCCCC
CC
       double precision  xlambda
	Integer Lminima
CC
       pi= dacos(-1.D0)
CC
CCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
	dy=Sdek(2,1)-Sdek(1,1)
!         write(*,*)'calculando la k del primer minimo de S(k)'
       	IKMI=Lminima(Sdek,NK)
       	kmin=Sdek(IKMI,1)
CC
CCCCCCCCCCCCCCCCCCCCCC Fin modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCC
CC
	write(*,*)'calculando la funcional del criterio'
	gamamin=0.D0
	gamamax=0.1D0
	dgama=1.D-5
	maxcrit=0.D0
        Do gama=gamamin,gamamax,dgama
         integral= 0.D0
CC
C evaluacion del integrando de la funcional del criterio
CC
         do iint=1, Nk
	  y=Sdek(iint,1)
          S=SdeK(iint,2)
          l= xlambda(y,kmin)
          num=gama*((S-1.D0)*l*y**2)**2
          den=36.D0*pi*fv*(gama*y**2+S*l)*(gama*y**2+l)
          T(iint)= num/den
         enddo
CC
C evaluacion de la integral
CC
         call SIMPSON(T,dy,Nk,integral)
CC
C determinando el maximo
CC
	 if(integral.gt.maxcrit) maxcrit=integral
        EndDo
        write(*,*)'funcional del criterio lista',maxcrit
	Criterio=maxcrit
	return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Solucion completa SCGLE						C
C									C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine  SCGLE(Sdek,fv,deltazf,NWR)
C     PROGRAMA BASE: TAC3.f
C     ESTE PROGRAMA CALCULA LAS FUNCIONES DE DISPERSION INTERMEDIA 
C     F(k,t) y SU PARTE SELF Fs(k,t), DE LA TEORIA AUTOCONSISTENTE 
C     PROPUESTA EN LA TESIS LYR. ADICIONALMENTE CALCULA LAS 
C     FUNCIONES MEMORIA C(t), C(k,t) y Cs(k,t).
C     REQUIERE COMO INSUMOS: EL FACTOR DE ESTRUCTURA S(k) Y LOS 
C     COEFICIENTES DE LA MEMORIA EXPONENCIAL SIMPLE COLECTIVA Y 
C     SELF.
C     ESTE PROGRAMA CORRESPONDE A CASO DE SISTEMAS COLOIDALES 
C     TRIDIMENSIONALES.
C     workshop 2006, tratando de implementar algoritmo de Fuchs para 
C     resolver a tiempos largotes.
       	IMPLICIT INTEGER*4 (I-N),Double precision(A-H,O-Z)
	Parameter (NK=2**12)
!NUMERO DE PUNTOS EN TIEMPOS CORTOS           
       	PARAMETER (NT2=10)
!tamaño de la malla temporal      
       	PARAMETER (MT=500)   
!NUMERO DE ITERACIONES           
      	 PARAMETER (NITER=5000)
!reescalamientos
	Parameter (NRES=25)
!MALLA K, MALLA T, FACTOR DE ESTRUCTURA            
       	DIMENSION Y(NK),T(0:MT+1),S(NK) 
!DELTAZ EN CADA REESCALAMIENTO      
       	DIMENSION DELTAZ(MT),DELTAZOLD(MT)
!FUNCIONES EN CADA REESCALAMIENTO      
       	DIMENSION FS(NK,MT),XC(NK,MT)
!MEMORIA SELF PARA CADA REESCALAMIENTO      
       	DIMENSION CS(NK,MT)
!MEMORIA COLECTIVA PARA CADA REESCALAMIENTO      
       	DIMENSION CC(NK,MT)
!INTEGRALES DE COLA EN CADA REESCALAMIENTO      
       	DIMENSION DFS(NK,MT),DXC(NK,MT)
       	DIMENSION DCS(NK,MT),DCC(NK,MT)
!SUMAS EN CADA REESCALAMIENTO      
       	DIMENSION A(MT),ASS(MT)
!G VINEYARD, GG CANTIDADES ESTATICAS      
       	DIMENSION G(NK,MT),GGO(NK)
C auxiliares
       	CHARACTER*2 DEC
       	DIMENSION Sdek(NK,2)
	Dimension deltazf(7000,2)
CC Mensaje inicial
       	write(*,*)' '
       	write(*,*)'resolviendo el esquema autoconsistente',NWR
       	write(*,*)' '
! 	IF(NWR.EQ.1) open(24,file='Fdekttau200pasostaus1000TA.dat')
       	PI=4.*ATAN(1.0)
       	IRES=0
C     DATOS DE ENTRADA:
	PHI=fv 
CC
C     ESCALA DE DT
       	DT=1.D-5
       	TOL=0.00001
C      dk=100.D0/dble(NK)
C      Y(1)=0.001d0
CC
CC minimo de S(k) para lambda y maximo para graficar
       	IKMI=Lminima(Sdek,NK)
       	YMIN=Sdek(IKMI,1) !*0.9113D0
       	write(*,*)'el minimo de S(k)',IKMI,YMIN
       	IKM=Lmaximo(Sdek,NK)
       	write(*,*)'el maximo de S(k)',IKM,Sdek(IKM,1),Sdek(IKM,2)
CC
       	DO  IK=1,NK
	 Y(IK)=Sdek(IK,1)
         S(IK)=Sdek(IK,2)
C     FUNCION G(k)
         GGO(IK)=((Y(IK)**4*(S(IK)-1.D0)**2)/S(IK))/(36.D0*PI*PHI)
       	ENDDO 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DIFERENCIAL DE K
        DY=Y(2)-Y(1)
C     calculo de deltaz(t=0)
       	CALL SIMPSON(GGO,DY,NK,DEZ0)
C     MALLA DEL EJE TEMPORAL:
       	T(0)=0.0
       	DO IT=1,MT
         T(IT)=T(IT-1)+DT
 !     READ(13,*)T(IT)!,DELTAZ(IT)
         GRANDO = 0.
CC
C Primeros puntos se aproximan con la aproximacion de tiempos cortos
CC
         DO IK=1,NK
          FS(IK,IT)=EXP(-Y(IK)**2*T(IT))
          XC(IK,IT)=EXP(-Y(IK)**2*T(IT)/S(IK))
          DFS(IK,IT)=FS(IK,IT)
          DXC(IK,IT)=XC(IK,IT)
          GRANDO = GRANDO + GGO(IK)*FS(IK,IT)*XC(IK,IT)
         ENDDO
         DELTAZ(IT)=GRANDO*DY
CC
C Ahora podemos tener las primeras expresiones para las memorias
CC
	 DO IK=1,NK
	  CS(IK,IT)=DELTAZ(IT)*xlambda(Y(IK),YMIN)
	  CC(IK,IT)=CS(IK,IT)
          DCS(IK,IT)=CS(IK,IT)
          DCC(IK,IT)=CC(IK,IT)
	 ENDDO
	 if(IT.le.NT2/2)then
! 	  write(26,*)T(IT),DELTAZ(IT)
! 	  IF(NWR.EQ.1)THEN
! 	   WRITE(24,*)'TIEMPO',T(IT)
! 	   DO IK=1,NK
! 	    WRITE(24,*)Y(IK),XC(IK,IT)
! 	   ENDDO 
! 	  ENDIF
	  deltazf(IT,1)=T(IT)
	  deltazf(IT,2)=DELTAZ(IT)
	 endif
       	ENDDO
C  AQUI EMPIEZA LO SABROSO>
	write(*,*)' '
	write(*,*)'Por fin termina estatica, empieza lo SABROSO'
	H=DT
       	DO IT=NT2/2+1,MT,1  !COMIENZA INTERVALO TEMPORAL
         DELTAZ(IT)=DELTAZ(IT-1) !INPUT INICIAL DE DELTAZ
         DO IA=1,NITER         !CICLO DE ITERACIONES        
          DELTAZOLD(IT)=DELTAZ(IT)
          GRANDO=0.
          DO IK=1,NK         !CICLO SOBRE VECTORES DE ONDA
           CS(IK,IT)=DELTAZ(IT)*xlambda(Y(IK),YMIN)
c              convol=grando2(ik)+CS(IK,IT)exp(-BT(ik)*T(it))*H     !ojo con los reescal
CCCCCCC   vineyard CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	    CC(IK,IT)=CS(IK,IT)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           SUMAS=0.
           SUMAC=0.
           IF(MOD(IT,2).EQ.0.) THEN
            DO IS=2,IT/2    !CICLO DE SUMATORIA
             SUMAS=SUMAS+DFS(IK,IS)*(CS(IK,IT-IS)-CS(IK,IT-IS+1))
     ^          +DCS(IK,IS)*(FS(IK,IT-IS)-FS(IK,IT-IS+1))
             SUMAC=SUMAC+DXC(IK,IS)*(CC(IK,IT-IS)-CC(IK,IT-IS+1))
     ^          +DCC(IK,IS)*(XC(IK,IT-IS)-XC(IK,IT-IS+1))
            ENDDO           !TERMINA CICLO DE SUMATORIA
             ASS(IT)=-CS(IK,IT/2)*FS(IK,IT/2)+DFS(IK,1)*CS(IK,IT-1)
     ^            +(1./H + DCS(IK,1))*FS(IK,IT-1)+SUMAS
             DEN1=(1./H+Y(IK)**2+DCS(IK,1))
             FS(IK,IT)=(CS(IK,IT)*(1.-DFS(IK,1))+ASS(IT))/DEN1
             A(IT)=-CC(IK,IT/2)*XC(IK,IT/2)+DXC(IK,1)*CC(IK,IT-1)
     ^          +(1./H + DCC(IK,1))*XC(IK,IT-1)+SUMAC
             DEN2=(1./H+(Y(IK)**2)/S(IK)+DCC(IK,1))
	     XC(IK,IT)=(CC(IK,IT)*(1.-DXC(IK,1))+A(IT))/DEN2
           ELSE
            DO IS=2,IT/2    !CICLO DE SUMATORIA
             SUMAS=SUMAS+DFS(IK,IS)*(CS(IK,IT-IS)-CS(IK,IT-IS+1))
     ^          +DCS(IK,IS)*(FS(IK,IT-IS)-FS(IK,IT-IS+1))
             SUMAC=SUMAC+DXC(IK,IS)*(CC(IK,IT-IS)-CC(IK,IT-IS+1))
     ^          +DCC(IK,IS)*(XC(IK,IT-IS)-XC(IK,IT-IS+1))
            ENDDO           !TERMINA CICLO DE SUMATORIA
             SUMAS=SUMAS+DFS(IK,IT-IT/2)*(CS(IK,IT-IT/2)-CS(IK,IT/2))/2.
     ^          +DCS(IK,IT-IT/2)*(FS(IK,IT-IT/2)-FS(IK,IT/2))/2.
             SUMAC=SUMAC+DXC(IK,IT-IT/2)*(CC(IK,IT-IT/2)-CC(IK,IT/2))/2.
     ^         +DCC(IK,IT-IT/2)*(XC(IK,IT-IT/2)-XC(IK,IT/2))/2.
             PROM=(CS(IK,IT-IT/2)+CS(IK,IT/2))*(FS(IK,IT-IT/2)+
     ^             FS(IK,IT/2))/4.
             ASS(IT)=-PROM+DFS(IK,1)*CS(IK,IT-1)
     ^            +(1./H + DCS(IK,1))*FS(IK,IT-1)+SUMAS
             DEN1=(1./H+Y(IK)**2+DCS(IK,1))
      	     FS(IK,IT)=(CS(IK,IT)*(1.-DFS(IK,1))+ASS(IT))/DEN1
             PROM=(CC(IK,IT-IT/2)+CC(IK,IT/2))*(XC(IK,IT-IT/2)
     ^             +XC(IK,IT/2))/4.
             A(IT)=-PROM+DXC(IK,1)*CC(IK,IT-1)
     ^          +(1./H + DCC(IK,1))*XC(IK,IT-1)+SUMAC
             DEN2=(1./H+(Y(IK)**2)/S(IK)+DCC(IK,1))
	     XC(IK,IT)=(CC(IK,IT)*(1.-DXC(IK,1))+A(IT))/DEN2
           ENDIF
           DFS(IK,IT)=FS(IK,IT)
           DXC(IK,IT)=XC(IK,IT)
           DCS(IK,IT)=CS(IK,IT)
           DCC(IK,IT)=CC(IK,IT)
	   GRANDO=GRANDO+GGO(IK)*FS(IK,IT)*XC(IK,IT)
          ENDDO             !CICLO SOBRE VECTORES DE ONDA
	  DELTAZ(IT)=GRANDO*DY
          GRANDO=0.
!CRITERIO DE CONVERGENCIA
	  IF((DELTAZ(IT)-DELTAZOLD(IT)).GE.0.D0)THEN
	   CRITERIO=(DELTAZ(IT)-DELTAZOLD(IT))/DELTAZ(IT)
	   IF(CRITERIO.LT.TOL)THEN
! 	    write(26,*)T(IT),DELTAZ(IT)
! 	   IF(NWR.EQ.1)THEN
! 	     WRITE(24,*)'TIEMPO',T(IT)
! 	     DO IK=1,NK
! 	      WRITE(24,*)Y(IK),XC(IK,IT)
! 	     ENDDO 
! 	    ENDIF
	    deltazf(IT,1)=T(IT)
	    deltazf(IT,2)=DELTAZ(IT)
            GOTO 1001
	   ENDIF
          ENDIF
	  IF((DELTAZOLD(IT)-DELTAZ(IT)).GT.0.D0)THEN
	   CRITERIO=(DELTAZOLD(IT)-DELTAZ(IT))/DELTAZOLD(IT)
	   IF(CRITERIO.LT.TOL)THEN
! 	    write(26,*)T(IT),DELTAZ(IT)
! 	    IF(NWR.EQ.1)THEN
! 	     WRITE(24,*)'TIEMPO',T(IT)
! 	     DO IK=1,NK
! 	      WRITE(24,*)Y(IK),XC(IK,IT)
! 	     ENDDO 
! 	    ENDIF
	    deltazf(IT,1)=T(IT)
	    deltazf(IT,2)=DELTAZ(IT)
            GOTO 1001
	   ENDIF
          ENDIF
!FIN CRITERIO DE CONVERGENCIA
          IF(IA.EQ.NITER) then
           WRITE(*,*) "NO CONVIRGIO",T(IT)
	   deltazf(ITAC+irel,2)=0.D0
	   return
          ENDIF
         ENDDO              !CICLO DE ITERACIONES         
 1001    continue		
       	ENDDO                !TERMINA INTERVALO TEMPORAL
CC
	ITAC=MT   !para llevar la cuenta en el tiempo
	write(*,*)'tiempo acumulado: ',ITAC,deltazf(ITAC,1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
COMIENZAN LOS BRINCOS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 201    IRES=IRES+1
        DT=DT+DT
        H=DT
	irel=1
C   REDEFINO FUNCIONES Y DIFERENCIALES ("decimation")
        DO IT=1,MT/2
         DELTAZ(IT)=DELTAZ(2*IT)
         T(IT)=T(2*IT)
         DO IK=1,NK
          FS(IK,IT)=FS(IK,2*IT)
          XC(IK,IT)=XC(IK,2*IT)
          CS(IK,IT)=CS(IK,2*IT)
          CC(IK,IT)=CC(IK,2*IT)
          DFS(IK,IT)=.5*(DFS(IK,2*IT)+DFS(IK,2*IT))
          DXC(IK,IT)=.5*(DXC(IK,2*IT)+DXC(IK,2*IT))
          DCS(IK,IT)=.5*(DCS(IK,2*IT)+DCS(IK,2*IT))
          DCC(IK,IT)=.5*(DCC(IK,2*IT)+DCC(IK,2*IT))
         ENDDO
        ENDDO
        DO IT=MT/2+1,MT,1  !COMIENZA INTERVALO TEMPORAL
         T(IT)=T(IT-1)+DT
CINPUT INICIAL DE DELTAZ
         DELTAZ(IT)=2.*DELTAZ(IT-1)-DELTAZ(IT-2) 
         DO IA=1,NITER         !CICLO DE ITERACIONES        
          DELTAZOLD(IT)=DELTAZ(IT)
          GRANDO=0.
          DO IK=1,NK         !CICLO SOBRE VECTORES DE ONDA
           CS(IK,IT)=DELTAZ(IT)*xlambda(Y(IK),YMIN)
c              convol=grando2(ik)+CS(IK,IT)exp(-BT(ik)*T(it))*H     !ojo con los reescal
CCCCCCC   vineyard CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	   CC(IK,IT)=CS(IK,IT)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           SUMAS=0.
           SUMAC=0.
          IF(MOD(IT,2).EQ.0.) THEN
C OJO! FALTA UTILIZAR LOS VALORES YA CALCULADOS EN LUGAR DE LOS PROMEDIOS
            DO IS=2,IT/2    !CICLO DE SUMATORIA
             SUMAS=SUMAS+DFS(IK,IS)*(CS(IK,IT-IS)-CS(IK,IT-IS+1))
     ^            +DCS(IK,IS)*(FS(IK,IT-IS)-FS(IK,IT-IS+1))
             SUMAC=SUMAC+DXC(IK,IS)*(CC(IK,IT-IS)-CC(IK,IT-IS+1))
     ^            +DCC(IK,IS)*(XC(IK,IT-IS)-XC(IK,IT-IS+1))
            ENDDO           !TERMINA CICLO DE SUMATORIA
            ASS(IT)=-CS(IK,IT/2)*FS(IK,IT/2)+DFS(IK,1)*CS(IK,IT-1)
     ^              +(1./H + DCS(IK,1))*FS(IK,IT-1)+SUMAS
            DEN1=(1./H+Y(IK)**2+DCS(IK,1))
      	    FS(IK,IT)=(CS(IK,IT)*(1.-DFS(IK,1))+ASS(IT))/DEN1
            A(IT)=-CC(IK,IT/2)*XC(IK,IT/2)+DXC(IK,1)*CC(IK,IT-1)
     ^            +(1./H + DCC(IK,1))*XC(IK,IT-1)+SUMAC
            DEN2=(1./H+(Y(IK)**2)/S(IK)+DCC(IK,1))
	    XC(IK,IT)=(CC(IK,IT)*(1.-DXC(IK,1))+A(IT))/DEN2
          ELSE
           DO IS=2,IT/2    !CICLO DE SUMATORIA
            SUMAS=SUMAS+DFS(IK,IS)*(CS(IK,IT-IS)-CS(IK,IT-IS+1))
     ^            +DCS(IK,IS)*(FS(IK,IT-IS)-FS(IK,IT-IS+1))
            SUMAC=SUMAC+DXC(IK,IS)*(CC(IK,IT-IS)-CC(IK,IT-IS+1))
     ^           +DCC(IK,IS)*(XC(IK,IT-IS)-XC(IK,IT-IS+1))
           ENDDO           !TERMINA CICLO DE SUMATORIA
            SUMAS=SUMAS+DFS(IK,IT-IT/2)*(CS(IK,IT-IT/2)-CS(IK,IT/2))/2.
     ^             +DCS(IK,IT-IT/2)*(FS(IK,IT-IT/2)-FS(IK,IT/2))/2.
            SUMAC=SUMAC+DXC(IK,IT-IT/2)*(CC(IK,IT-IT/2)-CC(IK,IT/2))/2.
     ^             +DCC(IK,IT-IT/2)*(XC(IK,IT-IT/2)-XC(IK,IT/2))/2.
            PROM=(CS(IK,IT-IT/2)+CS(IK,IT/2))*(FS(IK,IT-IT/2)
     ^       +FS(IK,IT/2))/4.
            ASS(IT)=-PROM+DFS(IK,1)*CS(IK,IT-1)
     ^             +(1./H + DCS(IK,1))*FS(IK,IT-1)+SUMAS
            DEN1=(1./H+Y(IK)**2+DCS(IK,1))
      	    FS(IK,IT)=(CS(IK,IT)*(1.-DFS(IK,1))+ASS(IT))/DEN1
            PROM=(CC(IK,IT-IT/2)+CC(IK,IT/2))*(XC(IK,IT-IT/2)
     ^        +XC(IK,IT/2))/4.
            A(IT)=-PROM+DXC(IK,1)*CC(IK,IT-1)
     ^            +(1./H + DCC(IK,1))*XC(IK,IT-1)+SUMAC
            DEN2=(1./H+(Y(IK)**2)/S(IK)+DCC(IK,1))
	     XC(IK,IT)=(CC(IK,IT)*(1.-DXC(IK,1))+A(IT))/DEN2
          ENDIF
c              if(IT.EQ.184)WRITE(*,*)'sum',SUMAS,SUMAC,ASS(IT),A(IT)
          DFS(IK,IT)=FS(IK,IT)
          DXC(IK,IT)=XC(IK,IT)
          DCS(IK,IT)=CS(IK,IT)
          DCC(IK,IT)=CC(IK,IT)
	  GRANDO=GRANDO+GGO(IK)*FS(IK,IT)*XC(IK,IT)
         ENDDO             !CICLO SOBRE VECTORES DE ONDA
	 DELTAZ(IT)=GRANDO*DY
         GRANDO=0.
!CRITERIO DE CONVERGENCIA
	 IF((DELTAZ(IT)-DELTAZOLD(IT)).GE.0.D0)THEN
	   CRITERIO=(DELTAZ(IT)-DELTAZOLD(IT))/DELTAZ(IT)
	  IF(CRITERIO.LT.TOL)THEN
! 	   write(26,*)T(IT),DELTAZ(IT)
	   deltazf(ITAC+irel,1)=T(IT)
	   deltazf(ITAC+irel,2)=DELTAZ(IT)
! 	   IF(NWR.EQ.1)THEN
! 	     WRITE(24,*)'TIEMPO',deltazf(ITAC+irel,1)
! 	     DO IK=1,NK
! 	      WRITE(24,*)Y(IK),XC(IK,IT)
! 	     ENDDO 
! 	    ENDIF
           GOTO 2006
	  ENDIF
	 ENDIF
	 IF((DELTAZOLD(IT)-DELTAZ(IT)).GT.0.D0)THEN
	  CRITERIO=(DELTAZOLD(IT)-DELTAZ(IT))/DELTAZOLD(IT)
	  IF(CRITERIO.LT.TOL)THEN
! 	   write(26,*)T(IT),DELTAZ(IT)
           deltazf(ITAC+irel,1)=T(IT)
	   deltazf(ITAC+irel,2)=DELTAZ(IT)
! 	   IF(NWR.EQ.1)THEN
! 	     WRITE(24,*)'TIEMPO',deltazf(ITAC+irel,1)
! 	     DO IK=1,NK
! 	      WRITE(24,*)Y(IK),XC(IK,IT)
! 	     ENDDO 
! 	    ENDIF
           GOTO 2006
          ENDIF
         ENDIF
!FIN CRITERIO DE CONVERGENCIA
         IF(IA.EQ.NITER) then
          WRITE(*,*) "NO CONVIRGIO",T(IT)
	  deltazf(ITAC+irel,2)=0.D0
	  return
         ENDIF
        ENDDO              !CICLO DE ITERACIONES         
 2006    continue		
	 irel=irel+1
      	ENDDO                !TERMINA INTERVALO TEMPORAL
CC
	ITAC=ITAC+250
	write(*,*)'tiempo acumulado: ',ITAC,deltazf(ITAC,1),
     ^                                 deltazf(ITAC,2)
      	IF(IRES.LE.NRES)GOTO 201
C
 1000  	continue
! 	IF(NWR.EQ.1) close(26)
       	Return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Funcion que calcula el coeficiente de difusion a timepos largos	C
C									C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function DLT(delz)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Este programa sirve para calcular el coeficiente de Difusion a tiemposC
C largos, a partir del archivo de salida de DeltaZ(t) de la teoria 	C
C autoconsistente.							C
C									C
C		Dl=1/(1+DeltaZ) 					C
C		con DeltaZ:=int_{0}^{infinity}DeltaZ(t)dt		C
C									C
C Febrero 2009								C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit double precision (a-h,o-z),integer*4(i-n)
	parameter(np=7000,nt=10**6)
	dimension dz(np),t(np),tv(nt),Ter(nt),delz(np,2)
c
	write(*,*)'calculando DL'
! 	lectura de la delta zeta estrella
	do i=1,7000
	 t(i)=delz(i,1)
	 dz(i)=delz(i,2)
	enddo
! mallas iniciales para el tiempo y w de t
	dt=1.d-3
	do i=0,nt
	 tv(i)=dble(i)*dt
	enddo
!  integral de DeltaZ(t) 
	do i=1,nt !ciclo temporal
	  Ter(i)=delta(tv(i),t,dz,np)
	enddo
	call SIMPSON(Ter,dt,nt,Dzstar)
	Dl=1.D0/(1.D0+Dzstar)
	write(*,*)'DL=',Dl,'integral=',Dzstar
	DLT=Dl
	Return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Subrutinas auxiliares							C
C									C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER function Lminima(Sk,N)

       Implicit none
       
       Integer N
       Integer h,m,j,kmin
       Double precision Sk(N,2)
       Double precision smin,smax

       smax=Sk(1,2)
               
       Do j=2, N
          if (Sk(j,2).GT.smax)then
            
             smax=Sk(j,2)
             m= j
          endif

       enddo
       
       smin= Sk(m,2)

       Do h=m , N

          if (Sk(h,2).LT.smin)then

             smin= Sk(h,2)
             kmin= h

          endif

       enddo

       Lminima=kmin
       Return

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Integer function Lmaximo(Sk,NK)
       Implicit none
       Integer NK
       Integer j,h
       Double precision Sk(NK,2)
       Double precision smax
     	smax=Sk(1,2)
       Do h=2,NK
	if (Sk(h,2).GT.smax)then
	  smax=Sk(h,2)
	  j=h
	endif
       enddo
       Lmaximo=j
       Return

       End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double precision function xlambda(x,xminima)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Funcion interpoladora lambda(k),se usa en la cerradura para		C 
C Cs(k,z) en la Teoria Autoconsistente. Ver detalles en la tesis	C
C de Laura Yeomans Reyna.						C
C Parametros: *x: valor de k en que queremos evaluar la funcion		C
C	      *xminima: valor de k en donde esta el primer minimo	C
C			de S(k)						C
C									C
C Agosto de 2006							C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       Implicit none

       Double precision x,xminima

       xlambda= 1.D0/(1.D0 + (x/xminima)**2)
       Return

      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Double precision function delta(tv,t,dz,np)
	implicit double precision (a-h,o-z),integer*4(i-n)
	dimension w(np),dz(np),t(np)
	do i=1,np
	 if(tv.eq.t(i))then
	  delta=dz(i)
	  return
	 endif
	 if(tv.lt.t(1))then
	  delta=dz(1)
	  return
	 endif
	 if(tv.lt.t(i))then
	  delta=(dz(i)+dz(i-1))/2.d0
	  return
	 endif
	enddo
	return
      end
