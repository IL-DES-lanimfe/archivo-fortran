      program SHS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      do q=0.001d0,100.d0,0.01d0
      sq=sqSHS(q,200.d0,0.4d0)
      write(1,*)q,sq
      enddo


      stop
      end
      DOUBLE PRECISION FUNCTION sqSHS(Q,TAU,PHI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      sqSHS=1.d0/
     -   ((1.d0 - (12.d0*phi)/q**3*
     -         (2.d0*(0.5d0*(1.d0 + 2.d0*phi - 
     -                 6.d0/phi*
     -                   ((tau + phi/(1.d0 - phi)) - 
     -                    ((tau + phi/(1.d0 - phi))**2 - 
     -                    (phi*(2.d0 + phi))/
     -                    (6.*(1.d0 - phi)**2))**0.5d0)*phi*
     -                  (1.d0 - phi)))/(1.d0 - phi)**2*
     -            (q*dCos(q) - dSin(q)) + 
     -           (-3.d0*phi + 
     -               6.d0/phi*
     -                 ((tau + phi/(1.d0 - phi)) - 
     -                   ((tau + phi/(1.d0 - phi))**2 - 
     -                    (phi*(2.d0 + phi))/
     -                    (6.d0*(1.d0 - phi)**2))**0.5d0)*phi*
     -                (1.d0 - phi))/((1.d0 - phi)**2*2.d0)*q*
     -            (dCos(q) - 1.d0) + 
     -           (6.d0/phi*
     -                ((tau + phi/(1.d0 - phi)) - 
     -                  ((tau + phi/(1.d0 - phi))**2 - 
     -                    (phi*(2.d0 + phi))/
     -                    (6.d0*(1.d0 - phi)**2))**0.5d0)*q**2)/
     -             12.*dSin(q)))**2 + 
     -     ((12.d0*phi)/q**3*
     -        (2.d0*(0.5d0*(1.d0 + 2.d0*phi - 
     -                6.d0/phi*
     -                  ((tau + phi/(1.d0 - phi)) - 
     -                    ((tau + phi/(1.d0 - phi))**2 - 
     -                    (phi*(2.d0 + phi))/
     -                    (6.d0*(1.d0 - phi)**2))**0.5d0)*phi*
     -                 (1.d0 - phi)))/(1.d0 - phi)**2*
     -           (q*dSin(q) + dCos(q) - 1.d0 - q**2/2.d0) + 
     -          (-3.d0*phi + 
     -              6.d0/phi*
     -                ((tau + phi/(1.d0 - phi)) - 
     -                  ((tau + phi/(1.d0 - phi))**2 - 
     -                    (phi*(2.d0 + phi))/
     -                    (6.d0*(1.d0 - phi)**2))**0.5d0)*phi*
     -               (1.d0 - phi))/((1.d0 - phi)**2*2.d0)*q*
     -           (dSin(q) - q) + 
     -          (6.d0/phi*((tau + phi/(1.d0 - phi)) - 
     -                 ((tau + phi/(1.d0 - phi))**2 - 
     -                    (phi*(2.d0 + phi))/
     -                    (6.d0*(1.d0 - phi)**2))**0.5d0)*q**2)/
     -            12.d0*(1.d0 - dCos(q))))**2)

      RETURN
      END