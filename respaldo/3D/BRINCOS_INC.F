      IRES=0
 201  IRES=IRES+1
c
c
C
C   REDEFINO FUNCIONES ("decimation")
C
	IF(IRES.LT.10.or.MOD(IRES,1).eq.0)THEN
      DT=DT+DT
	H=DT
      DO IT=1,NT/2
      T(IT)=T(2*IT)
      DELTAZ(1,1,IT)=DELTAZ(1,1,2*IT)
      DELTAZ(2,2,IT)=DELTAZ(2,2,2*IT)
!       T(IT)=T(2*IT)
       DO IK=1,NK
       FS(1,1,IK,IT)=FS(1,1,IK,2*IT)
      FS(2,2,IK,IT)=FS(2,2,IK,2*IT)
       FC(1,1,IK,IT)=FC(1,1,IK,2*IT)
      FC(1,2,IK,IT)=FC(1,2,IK,2*IT)
      FC(2,1,IK,IT)=FC(2,1,IK,2*IT)
      FC(2,2,IK,IT)=FC(2,2,IK,2*IT)
       CS(1,1,IK,IT)=CS(1,1,IK,2*IT)
      CS(2,2,IK,IT)=CS(2,2,IK,2*IT)
       CC(1,1,IK,IT)=CC(1,1,IK,2*IT)
      CC(1,2,IK,IT)=CC(1,2,IK,2*IT)
      CC(2,1,IK,IT)=CC(2,1,IK,2*IT)
      CC(2,2,IK,IT)=CC(2,2,IK,2*IT)
	ENDDO
      ENDDO
	ELSE
	DT=DT
	H=DT
	DO IT=1,NT/2
      T(IT)=T(NT/2+IT)!2*IT)
      DELTAZ(1,1,IT)=DELTAZ(1,1,NT/2+IT)!2*IT)
      DELTAZ(2,2,IT)=DELTAZ(2,2,NT/2+IT)!2*IT)
!       T(IT)=T(2*IT)
       DO IK=1,NK
       FS(1,1,IK,IT)=FS(1,1,IK,NT/2+IT)!2*IT)
      FS(2,2,IK,IT)=FS(2,2,IK,NT/2+IT)!2*IT)
       FC(1,1,IK,IT)=FC(1,1,IK,NT/2+IT)!2*IT)
      FC(1,2,IK,IT)=FC(1,2,IK,NT/2+IT)!2*IT)
      FC(2,1,IK,IT)=FC(2,1,IK,NT/2+IT)!2*IT)
      FC(2,2,IK,IT)=FC(2,2,IK,NT/2+IT)!2*IT)
       CS(1,1,IK,IT)=CS(1,1,IK,NT/2+IT)!2*IT)
      CS(2,2,IK,IT)=CS(2,2,IK,NT/2+IT)!2*IT)
       CC(1,1,IK,IT)=CC(1,1,IK,NT/2+IT)!2*IT)
      CC(1,2,IK,IT)=CC(1,2,IK,NT/2+IT)!2*IT)
      CC(2,1,IK,IT)=CC(2,1,IK,NT/2+IT)!2*IT)
      CC(2,2,IK,IT)=CC(2,2,IK,NT/2+IT)!2*IT)
	ENDDO
      ENDDO
	ENDIF
C
C
Ci    CALCULO DE CAPA [VER NOTAS]
      DO IK=1,NK !CICLO SOBRE K
      DO I=1,2 !RENGLONES
      CAPAS(I,I,IK)=ID(I,I)/H+AK(IK)**2*DCERO(I,I)+CS(I,I,IK,1)
      DO J=1,2 !COLUMNAS
      CAPA(I,J,IK)=ID(I,J)/H+AK(IK)**2*DCEROxSINV(I,J,IK)+
     ^            CC(I,J,IK,1)*ID(I,J)
      ENDDO !COLUNAS
      ENDDO !RENGLONES
C
      detCAPA=CAPA(1,1,IK)*CAPA(2,2,IK)-CAPA(1,2,IK)*CAPA(2,1,IK)
      CAPAINV(1,1,IK)=CAPA(2,2,IK)/detCAPA
      CAPAINV(2,2,IK)=CAPA(1,1,IK)/detCAPA
      CAPAINV(1,2,IK)=-CAPA(1,2,IK)/detCAPA
      CAPAINV(2,1,IK)=-CAPA(2,1,IK)/detCAPA
      ENDDO!CIERRE CICLO K
C
Cf    TERMINA CALCULO DE CAPA
Ci    DE UNA VEZ CALCULO S(k)-dF_1(K) Y 1/H+dC_1(k)
      DO IK=1,NK
      DO I=1,2
      SmDFS1(I,I,IK)=ID(I,I)-FS(I,I,IK,1)
      UEHMDCS1(I,I,IK)=ID(I,I)/H+CS(I,I,IK,1)
      DO J=1,2
      SmDFC1(I,J,IK)=S(I,J,IK)-FC(I,J,IK,1)
      UEHMDCC1(I,J,IK)=ID(I,J)/H+CC(I,J,IK,1)
      UEHMDCC1(I,J,IK)=UEHMDCC1(I,J,IK)*ID(I,J)
      ENDDO
      ENDDO
      ENDDO
Cf    TERMINA EL CALCULO DE S(k)-dF_1(K)
C
C
      DO IT=NT2/2+1,NT  !COMIENZA INTERVALO TEMPORAL
C
      T(IT)=T(IT-1)+DT
      DELTAZ(1,2,IT)=0.
      DELTAZ(2,1,IT)=0.
      DELTAZ(1,1,IT)=DELTAZ(1,1,IT-1) !IMPUT INICIAL DE DELTAZ
      DELTAZ(2,2,IT)=DELTAZ(2,2,IT-1)
C
      DO IA=1,NITER         !CICLO DE ITERACIONES        
C
      DELTAZOLD(1,1,IT)=DELTAZ(1,1,IT)
      DELTAZOLD(2,2,IT)=DELTAZ(2,2,IT)
C
      DO IK=1,NK         !CICLO SOBRE VECTORES DE ONDA
C
CALCULO DE DUNCIONES DE MEMORIA
      DO I=1,2
      DO J=1,2
      CCSEXP(I,J,IK,IT)=AC(I,J,IK)*exp(zc1(IK)*T(IT))
     ^                 +BC(I,J,IK)*exp(zc2(IK)*T(IT))
      CSSEXP(I,J,IK,IT)=ID(I,J)*AS(I,IK)*EXP(-BS(I,IK)*T(IT))
C
      CS(I,J,IK,IT)=CSSEXP(I,J,IK,IT)
     ^         +(DELTAZ(I,J,IT)-CSSEXP(I,J,IK,IT))*XLAM(I,J,IK)
! ! ! CCCCCCCCC NUEVA VINEYAR QUE SUSTITUYE LAS LINEAS ANTERIORRES------------
! ! ! 	CS(I,J,IK,IT)=ID(I,J)*DELTAZ(I,J,IT)*XLAM(I,J,IK)
! ! ! 	CC(I,J,IK,IT)=CS(I,J,IK,IT)
! ! ! CCCCCCCCC ---------------------------------------------------------------
CCCCC VINEYARD ADITIVA
      CC(I,J,IK,IT)=CS(I,J,IK,IT)*ID(I,J)+CCSEXP(I,J,IK,IT)
     ^                                   -CSSEXP(I,J,IK,IT)*ID(I,J)
c      vineyard multiplicativa
C            DO JT=1,IT
C            CONVC1(JT)=DELTAZ(I,J,JT)*EXP(ZC1(IK)*(T(IT)-T(JT)))
C            CONVC2(JT)=DELTAZ(I,J,JT)*EXP(ZC2(IK)*(T(IT)-T(JT)))
C            ENDDO
C            CALL TRAPE(CONVC1,DT,IT,RES1P)
C            CALL TRAPE(CONVC2,DT,IT,RES2P)
C      CC(I,J,IK,IT)=(AC(I,J,IK)+BC(I,J,IK))*XLAM(I,J,IK)*DELTAZ(I,J,IT)
C     &             /AS(I,IK)+(1.-XLAM(I,J,IK))*CCSEXP(I,J,IK,IT)+
C     &XLAM(I,J,IK)*AC(I,J,IK)*(BS(I,IK)+ZC1(IK))*RES1P/AS(I,IK)+
C     &XLAM(I,J,IK)*BC(I,J,IK)*(BS(I,IK)+ZC2(IK))*RES2P/AS(I,IK)
c     termina vineyard multiplicativa 
C
      ENDDO
      ENDDO
C
Ci    CALCULO ALGEBRAICO DE LA FC(k,t) Y LA FS(k,t)
Ci    CALCULOS AUXILIARES PARA FC Y FS
      DELTAZ(1,2,IT)=0.
      DELTAZ(2,1,IT)=0.
      UEHMDCC1(1,2,IK)=0.
      UEHMDCC1(2,1,IK)=0.
      DO I=1,2
      DO J=1,2
      !PARTE COLECTIVA
      AUXC1(I,J)=CAPAINV(I,1,IK)*CC(1,J,IK,IT)
     ^          +CAPAINV(I,2,IK)*CC(2,J,IK,IT)
      AUXC3(I,J)=CAPAINV(I,1,IK)*CC(1,J,IK,IT/2)
     ^          +CAPAINV(I,2,IK)*CC(2,J,IK,IT/2)
      AUXC5(I,J)=CAPAINV(I,1,IK)*CC(1,J,IK,IT-1)
     ^          +CAPAINV(I,2,IK)*CC(2,J,IK,IT-1)
      AUXC7(I,J)=CAPAINV(I,1,IK)*UEHMDCC1(1,J,IK)
     ^          +CAPAINV(I,2,IK)*UEHMDCC1(2,J,IK)
      ENDDO
      ENDDO
C
      DO I=1,2
      DO J=1,2
      !PARTE COLECTIVA
      AUXC2(I,J)=AUXC1(I,1)*SMDFC1(1,J,IK)
     ^          +AUXC1(I,2)*SMDFC1(2,J,IK)
      AUXC4(I,J)=AUXC3(I,1)*FC(1,J,IK,IT/2)
     ^          +AUXC3(I,2)*FC(2,J,IK,IT/2) 
      AUXC6(I,J)=AUXC5(I,1)*FC(1,J,IK,1)
     ^          +AUXC5(I,2)*FC(2,J,IK,1)
      AUXC8(I,J)=AUXC7(I,1)*FC(1,J,IK,IT-1)
     ^          +AUXC7(I,2)*FC(2,J,IK,IT-1)
      ENDDO
      ENDDO
      if((ia.eq.1).and.(it.eq.nt2/2+1))then
      endif
Ci    SUMA QUE APARECE EN LA F Y FS
      SUMAC(1,1)=0.
      SUMAC(1,2)=0.
      SUMAC(2,1)=0.
      SUMAC(2,2)=0.
      SUMAS(1,1)=0.
      SUMAS(1,2)=0.
      SUMAS(2,1)=0.
      SUMAS(2,2)=0.
      DO II=2,IT/2 !CICLO DE LA SUMA
      DO I=1,2
      DIFCS(I,I)=CS(I,I,IK,IT-II)-CS(I,I,IK,IT-II+1)
      DIFFS(I,I)=FS(I,I,IK,IT-II)-FS(I,I,IK,IT-II+1)
      DO J=1,2
      DIFCC(I,J)=CC(I,J,IK,IT-II)-CC(I,J,IK,IT-II+1)
      DIFCC(I,J)=DIFCC(I,J)*ID(I,J)
      DIFFC(I,J)=FC(I,J,IK,IT-II)-FC(I,J,IK,IT-II+1)
      ENDDO
      ENDDO
C
      DO I=1,2
      SUMS1(I,I)=DIFCS(I,I)*FS(I,I,IK,II)
      SUMS2(I,I)=CS(I,I,IK,II)*DIFFS(I,I)
      DO J=1,2
      SUMC1(I,J)=DIFCC(I,1)*FC(1,J,IK,II)+DIFCC(I,2)*FC(2,J,IK,II)
      SUMC2(I,J)=CC(I,1,IK,II)*DIFFC(1,J)+CC(I,2,IK,II)*DIFFC(2,J)
      ENDDO
      ENDDO
C
      DO I=1,2
      SUMAS(I,I)=SUMAS(I,I)+SUMS1(I,I)+SUMS2(I,I)
      DO J=1,2
      SUMAC(I,J)=SUMAC(I,J)+SUMC1(I,J)+SUMC2(I,J)
      ENDDO
      ENDDO
      ENDDO !CIERRE DEL CICLO SUMA
C
      DO I=1,2
      DO J=1,2
      SUMAF(I,J)=CAPAINV(I,1,IK)*SUMAC(1,J)+
     ^           CAPAINV(I,2,IK)*SUMAC(2,J)
      ENDDO
      ENDDO
Cf    TERMINA CALCULO DE LA SUMA
Cf    TERMINAN CALCULOS AUXILIARES PARA FC Y FS
C
      DO I=1,2
C
      FS(I,I,IK,IT)=CS(I,I,IK,IT)*(1.-FS(I,I,IK,1))-CS(I,I,IK,IT/2)
     &              *FS(I,I,IK,IT/2)+CS(I,I,IK,IT-1)*FS(I,I,IK,1)
     &              +(1./H+CS(I,I,IK,1))*FS(I,I,IK,IT-1)+SUMAS(I,I)
      FS(I,I,IK,IT)=FS(I,I,IK,IT)/CAPAS(I,I,IK)
      DO J=1,2
      FC(I,J,IK,IT)=AUXC2(I,J)-AUXC4(I,J)+AUXC6(I,J)+AUXC8(I,J)
     ^             +SUMAF(I,J)
      ENDDO
      ENDDO
C
      DO I=1,2
      DO J=1,2
      XI(I,J,IK,IT)=FC(I,1,IK,IT)*SINV(1,J,IK)+
     &             +FC(I,2,IK,IT)*SINV(2,J,IK)
      ENDDO
      ENDDO
c      if(mod(it,20).eq.0)WRITE(90,*)AK(IK),Fs(1,1,IK,IT),Fs(2,2,IK,IT)
c      if(mod(it,20).eq.0)WRITE(91,*)AK(IK),FC(1,2,IK,IT),FC(2,1,IK,IT)
c      if(mod(it,20).eq.0)WRITE(92,*)AK(IK),FC(1,1,IK,IT),FC(2,2,IK,IT)
      ENDDO !CIERRE DEL CICLO DE K's
C
Ci    CALCULO DE LA DELTAZ(t)
      DO IK=1,NK
	sum1=0.
	sum2=0.
	sum1=rho(1)*c_11(ik)*h_11(ik)*XI(1,1,IK,IT)+
     &	Sqrt(rho(1)*rho(2))*c_11(ik)*h_12(ik)*XI(1,2,IK,IT)+
     &	Sqrt(rho(1)*rho(2))*c_12(ik)*h_11(ik)*XI(2,1,IK,IT)+
     &	rho(2)*c_12(ik)*h_12(ik)*XI(2,2,IK,IT)
	sum2=rho(1)*c_12(ik)*h_12(ik)*XI(1,1,IK,IT)+
     &	Sqrt(rho(1)*rho(2))*c_12(ik)*h_22(ik)*XI(1,2,IK,IT)+
     &	Sqrt(rho(1)*rho(2))*c_22(ik)*h_12(ik)*XI(2,1,IK,IT)+
     &	rho(2)*c_22(ik)*h_22(ik)*XI(2,2,IK,IT)
! ! !       sum1=0.
! ! !       sum2=0.
! ! !       sum1=ache(1,1,IK)*(1.-SINV(1,1,IK))*XI(1,1,IK,IT)+
! ! !      &  (rho(2)/rho(1))*ache(1,2,IK)*(1.-SINV(1,1,IK))*XI(1,2,IK,IT)-
! ! !      &  sqrt(rho(1)/rho(2))*ache(1,1,IK)*SINV(1,2,IK)*XI(2,1,IK,IT)-
! ! !      &  sqrt(rho(2)/rho(1))*ache(1,2,IK)*SINV(1,2,IK)*XI(2,2,IK,IT)
! ! ! C
! ! !       sum2=ache(2,2,IK)*(1.-SINV(2,2,IK))*XI(2,2,IK,IT)+
! ! !      &  (rho(1)/rho(2))*ache(2,1,IK)*(1.-SINV(2,2,IK))*XI(2,1,IK,IT)-
! ! !      &  sqrt(rho(2)/rho(1))*ache(2,2,IK)*SINV(2,1,IK)*XI(1,2,IK,IT)-
! ! !      &  sqrt(rho(1)/rho(2))*ache(2,1,IK)*SINV(2,1,IK)*XI(1,1,IK,IT)  

      Zin1(ik)=(AK(IK)**4)*sum1*Fs(1,1,IK,IT)
      Zin2(ik)=(AK(IK)**4)*sum2*Fs(2,2,IK,IT)
C
      ENDDO
       call TRAPE(Zin1,dk,Nk,RES1)
       call TRAPE(Zin2,dk,Nk,RES2)
      DELTAZ(1,2,IT)=0.
      DELTAZ(2,1,IT)=0.
      DELTAZ(1,1,IT)=DS(1)*RES1/(6.*pi*pi)
      DELTAZ(2,2,IT)=DS(2)*RES2/(6.*pi**2)
Cf    TERMINA DELTAZ(t)
C
CCC-CRITERIOS DE CONVERGENCIA
      IF(POROSA.EQ.'NO')THEN
      INCLUDE 'CRITERIO_AMBAS_ESPECIES_BRINCOS.F'
      ENDIF
      IF(POROSA.EQ.'SI')THEN
      INCLUDE 'CRITERIO_POROSA_BRINCOS.F'
      ENDIF
! CCC-FIN DE CRITERIOS DE CONVERGENCIA
      ENDDO !CIERRE DEL CICLO DE ITERACIONES
2001  continue
c      write(73,*)T(it),deltaz(1,1,it),deltaz(2,2,it)
      ENDDO !CIERRE DEL CICLO TEMPORAL
Cf    TERMINA CALCULO ALGEBRAICO DE LA FC(k,t) Y LA FS(k,t)
C
C
       write(*,*)ires,t(nt)
C
CCC-MECANISMO PARA DETERMINAR LA TRANSICION
! 	IF(DELTAZ(1,1,IT).LT.TOL) GOTO 1
C
CCC-TERMINA MECANISMO PARA DETERMINAR LA TRANSICION
C
      IF(IRES.LE.NRES)GOTO 201