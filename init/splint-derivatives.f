      SUBROUTINE SPLINT_DERIVATIVES(XA,YA,Y2A,N,X,Y,YP,YPP)
      implicit real*8 (a-h,o-z)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) Print *, 'Bad XA input.'
      A=(XA(KHI)-X)/H
      DA = -1/H
      B=(X-XA(KLO))/H
      DB = -1/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      YP=(-YA(KLO)/H)+(YA(KHI)/H)+((H**2/6)*((((1/H)-((3*A**2)/H))*
     *      Y2A(KLO))+(((-1/H)+((3*B**2)/H))*Y2A(KHI))))
      YPP=(((Y2A(KLO)*A))+((Y2A(KHI)*B)))
      RETURN
      END
