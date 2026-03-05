      PROGRAM DECIMATIONH
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

       IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
!tamaño de la malla en K     
       PARAMETER (NK=2001)
       
!NUMERO DE PUNTOS EN TIEMPOS CORTOS           
       PARAMETER (NT2=10)
             
!tamaño de la malla temporal      
       PARAMETER (MT=500)   
      
!PARAMETRO QUE ME MULTIPLICA EL kmin         
       PARAMETER (RN=1.) 
      
!NUMERO DE REESCALAMIENTOS      
C      PARAMETER (NRES=2)  
      
!NUMERO DE ITERACIONES           
       PARAMETER (NITER=50000)

!MALLA K, MALLA T, FACTOR DE ESTRUCTURA            
       DIMENSION Y(NK),T(0:MT+1),S(NK) 
      
!DELTAZ EN CADA REESCALAMIENTO      
       DIMENSION DELTAZ(MT),DELTAZOLD(MT)
      
!FUNCIONES EN CADA REESCALAMIENTO      
       DIMENSION FS(NK,MT),XC(NK,MT)
      
!MEMORIA SELF PARA CADA REESCALAMIENTO      
       DIMENSION CS(NK,MT)
      
!COEFICIENTES SEXP E INTERPOLADORA      
       DIMENSION AS(NK),BS(NK),AT(NK),BT(NK),XLAM(NK)
      
!MEMORIA COLECTIVA PARA CADA REESCALAMIENTO      
       DIMENSION CC(NK,MT)
      
!SEXP EN CADA REESCALAMIENTO      
       DIMENSION CSSEXP(NK,MT),CCSEXP(NK,MT)
      
!INTEGRALES DE COLA EN CADA REESCALAMIENTO      
       DIMENSION DFS(NK,MT),DXC(NK,MT)
       DIMENSION DCS(NK,MT),DCC(NK,MT)
      
!SUMAS EN CADA REESCALAMIENTO      
       DIMENSION A(MT),ASS(MT)
      
!G VINEYARD, GG CANTIDADES ESTATICAS      
       DIMENSION G(NK,MT),GG(NK)
      
!SALIDAS      
       DIMENSION XCOUT(NK,MT),FSOUT(NK,MT),DELTAZOUT(MT)
      
C auxiliares
       CHARACTER*2 DEC
       DIMENSION Sdek(NK,2),Gder(NK,2),GG0(NK),FH(NK)
!       OPEN (12,FILE='grsspy.dat')
       OPEN (13,FILE='facdesLJMfv5146nu18.dat')
       OPEN (14,FILE='selfLJMfv5146nu18.dat')	
       OPEN (15,FILE='colLJMfv5146nu18.dat')
!       OPEN (50,FILE='hdk.dat')
       OPEN (25,FILE='fktmaxLJMfv5146nu18.dat')
       OPEN (29,FILE='fktminLJMfv5146nu18.dat')
       OPEN (26,FILE='deltazLJMfv5146nu18.dat')
       OPEN (27,FILE='fktLJMfv5146nu18.dat')
       OPEN (28,FILE='sueltosLJMfv5146nu18.dat')
       OPEN (30,FILE='fsktLJMfv5146nu18.dat')
CC Mensaje inicial
       write(*,*)' '
       write(*,*)'Este programa usa el metodo de Fuchs para resolver'
       write(*,*)'el esquema autoconsistente'
       write(*,*)' '
       PI=4.*ATAN(1.0)
       IRES=0
C     DATOS DE ENTRADA:
       WRITE(*,*)'¿QUE LAMBDA QUIERES?(1:=1,2:LORENZIANA)'
       READ(*,*)LAMBEL
       WRITE(*,*)'¿QUE VINEYARD QUIERES?(1:MULTIPLICATIVA,2:ADITIVA)'
       READ(*,*)LVINE
       WRITE(*,*)'¿FRACCION DE VOLUMEN?'
       READ(*,*)PHI
       WRITE(*,*)'¿CUANTOS REESCALAMIENTOS?'
       READ(*,*)NRES
!        WRITE(*,*)'¿QUIERES CHECAR SOLO LA ESTATICA?(s o n)'
!        READ(*,*)DEC
	DEC='n'
C       PHI=0.55D0
CC
C Cambiando a unidades de esfera dura
CC
!         RHO=PHI*6.D0*(0.934313D0)**3/PI
! 	write(*,*)RHO
! 	PHI=PI*RHO/6.D0
CC
C     ESCALA DE DT
       DT=1.D-5
       TOL=0.00001
C      dk=100.D0/dble(NK)
C      Y(1)=0.001d0
       WRITE(*,*)'Leyendo insumos estaticos'
       DO IK=1,NK
!	READ(50,*)Y(IK),FH(IK)
	FH(IK)=1.D0
! 	READ(12,*)Gder(IK,1),Gder(IK,2)
CC 
C S(k) en unidades de esfera dura
CC
 	READ(13,*)XX,Sdek(IK,2)
	Sdek(IK,1)=(1.D0)*XX
!	Y(IK)=Sdek(IK,1)
	READ(14,*)XX,AS(IK),BS(IK)
	READ(15,*)XX,AT(IK),BT(IK)
	
! 	write(8,*)Sdek(IK,1),Sdek(IK,2)	
       ENDDO
!       CLOSE(12)
       CLOSE(13)
       CLOSE(14)
       CLOSE(15)
       write(*,*)'ya esta S(k), g(r) y H(k)'
       WRITE(*,*)'Calculando insumos estaticos.Tenga paciencia'
CCCCCCCCCCC FACTOR DE ESTRUCTURA Y CANTIDAES ESTATICASCCCCCCCCCCCCCCC
       write(*,*)'ya esta S(k), g(r) y Y(r)',PHI
CC minimo de S(k) para lambda y maximo para graficar
       IKMI=Lminima(Sdek,NK)
       YMIN=Sdek(IKMI,1)
       write(*,*)'el minimo de S(k)',IKMI,YMIN
       write(28,*)'el minimo de S(k)',IKMI,Sdek(IKMI,1),Sdek(IKMI,2)
       IKM=Lmaximo(Sdek,NK)
       write(*,*)'el maximo de S(k)',IKM,Sdek(IKM,1),Sdek(IKM,2)
       write(28,*)'el maximo de S(k)',IKM,Sdek(IKM,1),Sdek(IKM,2)

       DO 1 IK=1,NK
C      READ(12,*)Y(IK),S(IK)
	Y(IK)=Sdek(IK,1)
        S(IK)=Sdek(IK,2)
CC
C Valor inicial de F(k,t) y Correlador
CC
	WRITE(27,*)Y(IK),S(IK),1.D0
	
C	FH(IK)=S(IK)+1.d0
C     FUNCION LAMDA
        IF(LAMBEL.EQ.1)THEN
         XLAM(IK)=1.D0
        ELSE
         XLAM(IK)=1.0/(1.0+(Y(IK)/(RN*YMIN))**2)
        ENDIF
C     FUNCION G(k)
C      GG(IK)=(Y(IK)**4*(S(IK)-1.0)**2)/S(IK)
        GG(IK)=(Y(IK)**4*(S(IK)-1.D0)**2)/S(IK)
        GG0(IK)=GG(IK)/(36.D0*PI*PHI)
C     COEFICIENTES DE LA MEMORIA EXPONENCIAL SELF (AS(K) Y BS(K)):
C      WRITE(*,*)AS(IK),BS(IK)
C     COEFICIENTES DE LA MEMORIA EXPONENCIAL COLECTIVA AT(K) Y BT(K):
C      READ(11,*)XX,AT(IK),BT(IK)
1      CONTINUE 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     DIFERENCIAL DE K
        DY=Y(2)-Y(1)
C     calculo de deltaz(t=0)
        CALL SIMPSON(GG0,DY,NK,DEZ0)
        WRITE(*,*)'DELTAZ EN t IGUAL A CERO ES:',DEZ0
        WRITE(28,*)'DELTAZ EN t IGUAL A CERO ES:',DEZ0
CC
C ESCRIBIMOS EL VALOR INICIAL DEL CORRELADOR EN EL MAXIMO Y MINIMO
CC
	WRITE(25,*)0.D0,1.D0
	WRITE(29,*)0.D0,1.D0
C     MALLA DEL EJE TEMPORAL:
        T(0)=0.0
        DO IT=1,MT
         T(IT)=T(IT-1)+DT
 !     READ(13,*)T(IT)!,DELTAZ(IT)
         GRANDO = 0.
         DO IK=1,NK
          FS(IK,IT)=EXP(-Y(IK)**2*T(IT))
C         XC(IK,IT)=EXP(-Y(IK)**2*T(IT)/S(IK))
	  XC(IK,IT)=EXP(-Y(IK)**2*FH(IK)*T(IT)/S(IK))
          CCSEXP(IK,IT)=AT(IK)*EXP(-BT(IK)*T(IT))                    
          CSSEXP(IK,IT)=AS(IK)*EXP(-BS(IK)*T(IT)) 
	  CS(IK,IT)=CSSEXP(IK,IT)
	  CC(IK,IT)=CCSEXP(IK,IT)
          DFS(IK,IT)=FS(IK,IT)
          DXC(IK,IT)=XC(IK,IT)
          DCS(IK,IT)=CS(IK,IT)
          DCC(IK,IT)=CC(IK,IT)
          GRANDO = GRANDO + GG(IK)*FS(IK,IT)*XC(IK,IT)
         ENDDO
          DELTAZ(IT)=GRANDO*DY/(36.D0*PI*PHI)
	  if(IT.le.NT2/2)then
	   write(26,*)T(IT),DELTAZ(IT)
	   WRITE(25,*)T(IT),Y(IKM)**2*T(IT),XC(IKM,IT)
	   WRITE(29,*)T(IT),Y(IKMI)**2*T(IT),XC(IKMI,IT)
	  endif
        ENDDO
	if(DEC.eq.'s')goto 1000
C  AQUI EMPIEZA LO SABROSO>
	write(*,*)' '
	write(*,*)'Por fin termina estatica, empieza lo SABROSO'
	H=DT
        DO IT=NT2/2+1,MT,1  !COMIENZA INTERVALO TEMPORAL
         !write(*,*)"tiempo",it
         DELTAZ(IT)=DELTAZ(IT-1) !IMPUT INICIAL DE DELTAZ
         DO IA=1,NITER         !CICLO DE ITERACIONES        
          DELTAZOLD(IT)=DELTAZ(IT)
          !write(*,*)"iteracion",IA
          GRANDO=0.
          DO IK=1,NK         !CICLO SOBRE VECTORES DE ONDA
           CS(IK,IT)=CSSEXP(IK,IT)
     ^               +(DELTAZ(IT)-CSSEXP(IK,IT))*XLAM(IK)
c              convol=grando2(ik)+CS(IK,IT)exp(-BT(ik)*T(it))*H     !ojo con los reescal
CCCCCCC   vineyard CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	   IF(LVINE.EQ.1)THEN
            CONVOL=0.
            DO ITT=1,IT
             CONVOL=CONVOL+DELTAZ(ITT)*EXP(BT(IK)*(T(ITT)-T(IT)))
            ENDDO
            CC(IK,IT)=CS(IK,IT)+
     ^      (1.-XLAM(IK))*(AS(IK)*CCSEXP(IK,IT)/AT(IK)-CSSEXP(IK,IT))+
     ^      (BS(IK)-BT(IK))*XLAM(IK)*CONVOL*DT
            CC(IK,IT)=AT(ik)*CC(IK,IT)/AS(IK)
	   ELSE
	     CC(IK,IT)=CS(IK,IT)+CCSEXP(IK,IT)-CSSEXP(IK,IT)
	   ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           SUMAS=0.
           SUMAC=0.
           IF(MOD(IT,2).EQ.0.) THEN
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
C              DEN2=(1./H+(Y(IK)**2)/S(IK)+DCC(IK,1))
	    DEN2=(1./H+(FH(IK)*Y(IK)**2)/S(IK)+DCC(IK,1))
	    XC(IK,IT)=(CC(IK,IT)*(1.-DXC(IK,1))+A(IT))/DEN2
           ELSE
           DO IS=2,IT/2    !CICLO DE SUMATORIA
            SUMAS=SUMAS+DFS(IK,IS)*(CS(IK,IT-IS)-CS(IK,IT-IS+1))
     ^            +DCS(IK,IS)*(FS(IK,IT-IS)-FS(IK,IT-IS+1))
            SUMAC=SUMAC+DXC(IK,IS)*(CC(IK,IT-IS)-CC(IK,IT-IS+1))
     ^            +DCC(IK,IS)*(XC(IK,IT-IS)-XC(IK,IT-IS+1))
           ENDDO           !TERMINA CICLO DE SUMATORIA
           SUMAS=SUMAS+DFS(IK,IT-IT/2)*(CS(IK,IT-IT/2)-CS(IK,IT/2))/2.
     ^           +DCS(IK,IT-IT/2)*(FS(IK,IT-IT/2)-FS(IK,IT/2))/2.
           SUMAC=SUMAC+DXC(IK,IT-IT/2)*(CC(IK,IT-IT/2)-CC(IK,IT/2))/2.
     ^          +DCC(IK,IT-IT/2)*(XC(IK,IT-IT/2)-XC(IK,IT/2))/2.
      PROM=(CS(IK,IT-IT/2)+CS(IK,IT/2))*(FS(IK,IT-IT/2)+FS(IK,IT/2))/4.
           ASS(IT)=-PROM+DFS(IK,1)*CS(IK,IT-1)
     ^              +(1./H + DCS(IK,1))*FS(IK,IT-1)+SUMAS
           DEN1=(1./H+Y(IK)**2+DCS(IK,1))
      	   FS(IK,IT)=(CS(IK,IT)*(1.-DFS(IK,1))+ASS(IT))/DEN1
      PROM=(CC(IK,IT-IT/2)+CC(IK,IT/2))*(XC(IK,IT-IT/2)+XC(IK,IT/2))/4.
           A(IT)=-PROM+DXC(IK,1)*CC(IK,IT-1)
     ^            +(1./H + DCC(IK,1))*XC(IK,IT-1)+SUMAC
C              DEN2=(1./H+(Y(IK)**2)/S(IK)+DCC(IK,1))
	   DEN2=(1./H+(FH(IK)*Y(IK)**2)/S(IK)+DCC(IK,1))
	   XC(IK,IT)=(CC(IK,IT)*(1.-DXC(IK,1))+A(IT))/DEN2
          ENDIF
          DFS(IK,IT)=FS(IK,IT)
          DXC(IK,IT)=XC(IK,IT)
          DCS(IK,IT)=CS(IK,IT)
          DCC(IK,IT)=CC(IK,IT)
	  GRANDO=GRANDO+GG(IK)*FS(IK,IT)*XC(IK,IT)
         ENDDO             !CICLO SOBRE VECTORES DE ONDA
	 DELTAZ(IT)=GRANDO*DY/(36.D0*PI*PHI)
         GRANDO=0.
!CRITERIO DE CONVERGENCIA
	 IF((DELTAZ(IT)-DELTAZOLD(IT)).GE.0.D0)THEN
	  CRITERIO=(DELTAZ(IT)-DELTAZOLD(IT))/DELTAZ(IT)
	  IF(CRITERIO.LT.TOL)THEN
	   if(mod(it,100).eq.0.)then
	    WRITE(*,*)"ad",t(it)
	    do IK=2,NK
	     WRITE(27,*)Y(IK),sngl(XC(IK,IT)*S(IK)),sngl(XC(IK,IT))
	     WRITE(30,*)Y(IK),sngl(FS(IK,IT))
	    enddo
	   endif
	   write(26,*)T(IT),DELTAZ(IT)
	   WRITE(25,*)T(IT),Y(IKM)**2*T(IT),XC(IKM,IT)
	   WRITE(29,*)T(IT),Y(IKMI)**2*T(IT),XC(IKMI,IT)
           GOTO 1001
	  ENDIF
         ENDIF
	 IF((DELTAZOLD(IT)-DELTAZ(IT)).GT.0.D0)THEN
	  CRITERIO=(DELTAZOLD(IT)-DELTAZ(IT))/DELTAZOLD(IT)
	  IF(CRITERIO.LT.TOL)THEN
	   if(mod(it,100).eq.0.)then
	    WRITE(*,*)"ad",t(it)
	    do IK=2,NK
	     WRITE(27,*)Y(IK),sngl(XC(IK,IT)*S(IK)),sngl(XC(IK,IT))
	     WRITE(30,*)Y(IK),sngl(FS(IK,IT))
	    enddo
	   endif
	   write(26,*)T(IT),DELTAZ(IT)
	   WRITE(25,*)T(IT),Y(IKM)**2*T(IT),XC(IKM,IT)
	   WRITE(29,*)T(IT),Y(IKMI)**2*T(IT),XC(IKMI,IT)
           GOTO 1001
	  ENDIF
         ENDIF
!FIN CRITERIO DE CONVERGENCIA
         IF(IA.EQ.NITER) then
          WRITE(*,*) "NO CONVIRGIO",T(IT)
	  stop
         endif
         ENDDO              !CICLO DE ITERACIONES         
 1001    continue		
        ENDDO                !TERMINA INTERVALO TEMPORAL

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
COMIENZAN LOS BRINCOS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 201    IRES=IRES+1
        DT=DT+DT
        H=DT
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
           CSSEXP(IK,IT)=AS(IK)*EXP(-BS(IK)*T(IT))
           CCSEXP(IK,IT)=AT(IK)*EXP(-BT(IK)*T(IT))
           CS(IK,IT)=CSSEXP(IK,IT)
     ^                  +(DELTAZ(IT)-CSSEXP(IK,IT))*XLAM(IK)
c              if(IT.EQ.184)WRITE(*,*)'sexp',CSSEXP(IK,IT),CCSEXP(IK,IT)
c              if(IT.EQ.184)WRITE(*,*)'self',CS(IK,IT),DELTAZ(IT)
c              convol=grando2(ik)+CS(IK,IT)exp(-BT(ik)*T(it))*H     !ojo con los reescal
CCCCCCC   vineyard CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	   IF(LVINE.EQ.1)THEN
             CONVOL=0.
            DO ITT=1,IT
             CONVOL=CONVOL+DELTAZ(ITT)*EXP(BT(IK)*(T(ITT)-T(IT)))
            ENDDO
            CC(IK,IT)=CS(IK,IT)+
     ^      (1.-XLAM(IK))*(AS(IK)*CCSEXP(IK,IT)/AT(IK)-CSSEXP(IK,IT))+
     ^      (BS(IK)-BT(IK))*XLAM(IK)*CONVOL*DT
	    CC(IK,IT)=AT(ik)*CC(IK,IT)/AS(IK)
	   ELSE
            CC(IK,IT)=CS(IK,IT)+CCSEXP(IK,IT)-CSSEXP(IK,IT)
	   ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
C              DEN2=(1./H+(Y(IK)**2)/S(IK)+DCC(IK,1))
	   DEN2=(1./H+(FH(IK)*Y(IK)**2)/S(IK)+DCC(IK,1))
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
      PROM=(CS(IK,IT-IT/2)+CS(IK,IT/2))*(FS(IK,IT-IT/2)+FS(IK,IT/2))/4.

          ASS(IT)=-PROM+DFS(IK,1)*CS(IK,IT-1)
     ^             +(1./H + DCS(IK,1))*FS(IK,IT-1)+SUMAS
          DEN1=(1./H+Y(IK)**2+DCS(IK,1))
      	  FS(IK,IT)=(CS(IK,IT)*(1.-DFS(IK,1))+ASS(IT))/DEN1
      PROM=(CC(IK,IT-IT/2)+CC(IK,IT/2))*(XC(IK,IT-IT/2)+XC(IK,IT/2))/4.
          A(IT)=-PROM+DXC(IK,1)*CC(IK,IT-1)
     ^            +(1./H + DCC(IK,1))*XC(IK,IT-1)+SUMAC
C         DEN2=(1./H+(Y(IK)**2)/S(IK)+DCC(IK,1))
	  DEN2=(1./H+(FH(IK)*Y(IK)**2)/S(IK)+DCC(IK,1))
	  XC(IK,IT)=(CC(IK,IT)*(1.-DXC(IK,1))+A(IT))/DEN2
          ENDIF
c              if(IT.EQ.184)WRITE(*,*)'sum',SUMAS,SUMAC,ASS(IT),A(IT)
          DFS(IK,IT)=FS(IK,IT)
          DXC(IK,IT)=XC(IK,IT)
          DCS(IK,IT)=CS(IK,IT)
          DCC(IK,IT)=CC(IK,IT)
	  GRANDO=GRANDO+GG(IK)*FS(IK,IT)*XC(IK,IT)
c              if(IT.EQ.184)WRITE(*,*)'grando',GRANDO
         ENDDO             !CICLO SOBRE VECTORES DE ONDA
	 DELTAZ(IT)=GRANDO*DY/(36.D0*PI*PHI)
         GRANDO=0.
!CRITERIO DE CONVERGENCIA
	 IF((DELTAZ(IT)-DELTAZOLD(IT)).GE.0.D0)THEN
	  CRITERIO=(DELTAZ(IT)-DELTAZOLD(IT))/DELTAZ(IT)
	  IF(CRITERIO.LT.TOL)THEN
	   if(mod(it,100).eq.0.)then
	    WRITE(*,*)"ad",t(it)
	    do IK=2,NK
	     WRITE(27,*)Y(IK),sngl(XC(IK,IT)*S(IK)),sngl(XC(IK,IT))
	     WRITE(30,*)Y(IK),sngl(FS(IK,IT))
	    enddo
	   endif
	   write(26,*)T(IT),DELTAZ(IT)
	   WRITE(25,*)T(IT),Y(IKM)**2*T(IT),XC(IKM,IT)
	   WRITE(29,*)T(IT),Y(IKMI)**2*T(IT),XC(IKMI,IT)
           GOTO 2006
	  ENDIF
         ENDIF
	 IF((DELTAZOLD(IT)-DELTAZ(IT)).GT.0.D0)THEN
	  CRITERIO=(DELTAZOLD(IT)-DELTAZ(IT))/DELTAZOLD(IT)
	  IF(CRITERIO.LT.TOL)THEN
	   if(mod(it,100).eq.0.)then
	    WRITE(*,*)"ad",t(it)
	    do IK=2,NK
	     WRITE(27,*)Y(IK),sngl(XC(IK,IT)*S(IK)),sngl(XC(IK,IT))
	     WRITE(30,*)Y(IK),sngl(FS(IK,IT))
	    enddo
	   endif
	   write(26,*)T(IT),DELTAZ(IT)
	   WRITE(25,*)T(IT),Y(IKM)**2*T(IT),XC(IKM,IT)
	   WRITE(29,*)T(IT),Y(IKMI)**2*T(IT),XC(IKMI,IT)
           GOTO 2006
	  ENDIF
         ENDIF
!FIN CRITERIO DE CONVERGENCIA
         IF(IA.EQ.NITER) then
          WRITE(*,*) "NO CONVIRGIO",T(IT)
	  Stop
         endif
        ENDDO              !CICLO DE ITERACIONES         
 2006  continue		
       ENDDO                !TERMINA INTERVALO TEMPORAL

       write(*,*)ires,t(mt)
       write(28,*)'valor final de deltaz:',T(MT),DELTAZ(MT)
       IF(IRES.LE.NRES)GOTO 201
!        DO IT=1,MT,10
!         DO IK=1,NK
!          WRITE(27,*)Y(IK),sngl(XC(IK,IT)*S(IK)),sngl(XC(IK,IT))
!          !WRITE(14,*)sngl(Y(IK)),sngl(XCOUT(IK,IT+1)*S(IK))
!         ENDDO
!        ENDDO

      
 1000  close(25)
       close(26)
       close(27)
       close(29)
       STOP
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    C
C SUBRUTINAS PARA EL FACTOR DE ESTRUCUTURA                          C
C								    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    C
C FIN DE SUBRUTINAS PARA EL FACTOR DE ESTRUCTURA. INICIAN SUBRUTI-  C
C NAS AUXILIARES 						    C
C								    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C							             C
C La siguiente funcion sirve para calcular el primer minimo de S(k)  C
C para la funcion interpoladora.                                     C
C								     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER function Lminima(Sk,N)

       Implicit none
       
       Integer N
C       Double precision emp

       Integer h,m,j,kmin
       Double precision Sk(N,2)
C       Double precision y,dy
       Double precision smin,smax

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

       Do h=m , N

          if (Sk(h,2).LT.smin)then

             smin= Sk(h,2)
             kmin= h

          endif

       enddo

       Lminima=kmin
       Return

      End

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								      C
C Para observar el plateu localizamos el primer pico del correlador.  C
C la siguiente funcion hace esa chamba.				      C
C								      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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

