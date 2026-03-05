      Double precision function Fmaximo(Sk,NK,MT,IT)

       Implicit none
       
       Integer NK,MT,IT

       Integer j,m,h
       Double precision Sk(NK,MT)
       Double precision smin, smax


       smin=Sk(1,IT)
               
       Do j=2, NK/19
          if (Sk(j,IT).LT.smin)then
            
             smin=Sk(j,IT)
	     m=j

          endif
       enddo

	smax=Sk(m,IT)

       Do h=m,NK/9
	if (Sk(h,IT).GT.smax)then

	  smax=Sk(h,IT)

	endif
       enddo
       Fmaximo=smax
       Return

       End