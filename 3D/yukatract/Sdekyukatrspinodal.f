      PROGRAM yukawaatrinter
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Este programa sirve para probar al subrutina de S(k) para el poten-C
C cial de yukawa atractivo.					     C
C								     C
C Enero 2007							     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none !toda variable se tiene que declarar
	Integer Nk,i
	parameter (Nk=2**12)
	Double precision fv,Ka,z
	Double precision dk,Rk(Nk),Sk(Nk)
	double precision fvg,fvl,Skg(Nk,2),Skl(Nk,2),Skef(Nk,2)
	double precision xi,pi,sp
	double precision Porl,Porg
	pi=4.D0*datan(1.D0)
CC
C Parametros
CC
	fv=0.5D0
	Ka=18.D0
	z=20.D0
	dk=0.01D0
CC
! 	open(90,file='spinodal.dat')
! 	do i=1,100
! 	 read(90,*)Ts,Ks,fvs
! 	 if((Ks-Ka).le.0.1D0)then
! 	  
! 	enddo		
CC
CC
C Ahora llamamos a la subrutina de S(k)
CC
	call msa(fv,Ka,z,dk,Nk,Rk,Sk,xi,sp)
CC
C salida
CC
	if(sp.gt.0.d0)then
	 open(50,file='Skyukatrz20fv0,5K18.dat')
 	 do i=1,Nk
	  write(50,*)2.d0*Rk(i),Sk(i)
	 enddo
	 close(50)
CC
C Calculamos extremos para interpolar
CC
	else
         call spinodal(fv,Ka,z,fvg,fvl)
CC
C calculamos porcentages
CC
! 	Porl=(fv-fvg)/(fvl-fvg)
! 	Porg=1.D0-Porl
CC
C Regla 1
CC
! 	 fvl=fvl+0.05*fvl
! 	 fvg=fvg-0.05*fvg
! 	 Porl=(fvl/fv)*(fv-fvg)/(fvl-fvg)
! 	 Porg=(fvg/fv)*(fvl-fv)/(fvl-fvg)
CC
C Regla 2
CC
! 	 Porl=(fv-fvg)/(fvl-fvg)
! 	 Porg=(fvl-fv)/(fvl-fvg)
! 	 write(*,*)Porl,'+',Porg,'=',Porl+Porg
CC
C calculamos los factores de estructura para interpolar
CC
	 call msa(fvg,Ka,z,dk,Nk,Rk,Sk,xi,sp)
	 do i=1,Nk
	  Skg(i,1)=2.D0*Rk(i)
	  Skg(i,2)=Sk(i)
! 	  write(80,*)Skg(i,1),Skg(i,2)
	 enddo
	 call msa(fvl,Ka,z,dk,Nk,Rk,Sk,xi,sp)
	 do i=1,Nk
	  Skl(i,1)=2.D0*Rk(i)
	  Skl(i,2)=Sk(i)
! 	  write(90,*)Skl(i,1),Skl(i,2)
CC
C Factor de estructura efectivo dentro de la espinodal
CC
	  Skef(i,1)=Skl(i,1)
CC
C La siguiente ecuacion vale para la regla 1 y 2
CC
! 	  Skef(i,2)= Porl*Skl(i,2)+Porg*Skg(i,2)
CC
C Regla 3
CC
	  Skef(i,2)=(fvl-fvg)*Skl(i,2)*Skg(i,2)/((fv-fvg)*Skg(i,2)
     &              +(fvl-fv)*Skl(i,2))
	 enddo

	 open(60,file='Skyukatrz20fv0,5K18guev.dat')
	 do i=1,Nk
	  write(60,*)Skef(i,1),Skef(i,2)
	 enddo
	 close(60)
	endif
	stop
      END
CC

      subroutine msa(vf,XDK,b,DQ,N,Q,SQ,chi,aa)
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
CCC      
      subroutine spinodal(fvo,Ka,z,fvg,fvl)
	implicit none
	Double precision fvo,fv,Ka,z,fvg,fvl
	Double precision Rk(2**12),Sk(2**12),xi
	Double precision spin,sping,spingn,spinl,spinln
	Double precision fvgn,fvln
	fv=fvo
	call msa(fv,Ka,z,0.001D0,2**12,Rk,Sk,xi,spin)
CC 
C Determinamos la fraccion de volumen de gas
CC
	fvg=0.03D0
	call msa(fvg,Ka,z,0.001D0,2**12,Rk,Sk,xi,sping)
	if((sping*spin).lt.0.D0)then
 50	 continue
	 fvgn=(fv+fvg)/2.D0
	 call msa(fvgn,Ka,z,0.001D0,2**12,Rk,Sk,xi,spingn)
CC
C si ya encontramos la espinodal en la parte gas
CC
	 if(dabs(spingn).le.1.D-6)then 
	  fvg=fvgn
	  write(*,*)'phi de gas',fvg,spingn
	  goto 100
	 endif
CC
C si no se ha encontrado bisecciones sucesivas
CC
	 if((sping*spingn).lt.0.D0)then
	  fv=fvgn
	  spin=spingn
	  goto 50
	 endif
	 if((spin*spingn).lt.0.D0)then
	  fvg=fvgn
	  sping=spingn
	 goto 50
	 endif
CC
 100	 continue
	endif
CC
C inicializamos para encontrar fraccion de volumen de liquido
CC	
	fv=fvo
	call msa(fv,Ka,z,0.001D0,2**12,Rk,Sk,xi,spin)
CC 
C Determinamos la fraccion de volumen de gas
CC
	fvl=0.6D0
	call msa(fvl,Ka,z,0.001D0,2**12,Rk,Sk,xi,spinl)
	if((spinl*spin).lt.0.D0)then
 60	 continue
	 fvln=(fv+fvl)/2.D0
	 call msa(fvln,Ka,z,0.001D0,2**12,Rk,Sk,xi,spinln)
CC
C si ya encontramos la espinodal en la parte gas
CC
	 if(dabs(spinln).le.1.D-6)then 
	  fvl=fvln
	  write(*,*)'phi de liquido',fvl,spinln
	  goto 200
	 endif
CC
C si no se ha encontrado bisecciones sucesivas
CC
	 if((spinl*spinln).lt.0.D0)then
	  fv=fvln
	  spin=spinln
	  goto 60
	 endif
	 if((spin*spinln).lt.0.D0)then
	  fvl=fvln
	  spinl=spinln
	 goto 60
	 endif
CC
 200	 continue
	endif
CC
C una vez que tenemos los dos nos vamos
CC
	return
      end