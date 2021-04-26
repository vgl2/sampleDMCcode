c
c  Reduction of a real symmetric matrix, A, to tridagonal
c  form by the Householder method. This is followed by the
c  evaluation of the eigenvalues and eigenvectors.
c
c    A - The (N,N) matrix to be diagonalized; physical dimension (NP,NP).
c    D - Array of dimension N containing eigenvalues on output; PD: (NP,NP).
c    
c    A is returned as the matrix of eigenvalues
c    E appears to be a scratch matrix used by the subroutine
      SUBROUTINE HOUSE(A,N,NP,D,E)
      implicit real*8(a-h,o-z)
      DIMENSION A(NP,NP),D(NP),E(NP)
      IF(N.GT.1)THEN
        DO 18 I=N,2,-1  
          L=I-1
          H=0.
          SCALE=0.
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+dABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-SIGN(dSQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=0.
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=0.
      E(1)=0.
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.0.)THEN
          DO 21 J=1,L
            G=0.
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1.
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=0.
            A(J,I)=0.
22        CONTINUE
        ENDIF
23    CONTINUE
c
c  Now diagonalize the triadiagonal matrix produced above.
c
      IF (N.GT.1) THEN
        DO 41 I=2,N
          E(I-1)=E(I)
41      CONTINUE
        E(N)=0.
        DO 45 L=1,N
          ITER=0
1         DO 42 M=L,N-1
            DD=dABS(D(M))+dABS(D(M+1))
            IF (dABS(E(M))+DD.EQ.DD) GO TO 2
42        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30)PAUSE 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.*E(L))
            R=SQRT(G**2+1.)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.
            C=1.
            P=0.
            DO 44 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(dABS(F).GE.dABS(G))THEN
                C=G/F
                R=dSQRT(C**2+1.)
                E(I+1)=F*R
                S=1./R
                C=C*S
              ELSE
                S=F/G
                R=dSQRT(S**2+1.)
                E(I+1)=G*R
                C=1./R  
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 43 K=1,N
                F=A(K,I+1)
                A(K,I+1)=S*A(K,I)+C*F
                A(K,I)=C*A(K,I)-S*F
43            CONTINUE
44          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.
            GO TO 1
          ENDIF
45      CONTINUE
      ENDIF
c
c  Now sort the eigenvalues and eigenvectors increasing order.
c
      DO 33 I=1,N-1
        K=I
        P=D(I)
        DO 31 J=I+1,N
          IF(D(J).LE.P)THEN
            K=J
            P=D(J)
          ENDIF
31      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 32 J=1,N
            P=A(J,I)
            A(J,I)=A(J,K)
            A(J,K)=P
32        CONTINUE
        ENDIF
33    CONTINUE
      RETURN
      END
