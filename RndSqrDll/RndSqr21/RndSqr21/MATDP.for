      !A(N,N)- ¬ходна€ матрица, размерностью (N,N)
      !E(N,N)- ќбращенна€ матрица, размерностью (N,N)
      subroutine matr(A,E,N)
      implicit complex*16(c,w,z), real*8(a-b,d-h,o-v,x-y)
      DIMENSION A(N,N),E(N,N)
25    FORMAT(5X,'—бой при обращении матрицы')
      DO 11 I=1,N
      DO 11 J=1,N
      IF(I.EQ.J)GO TO 10
      E(I,J)=0.
      GO TO 11
10    E(I,J)=1.
11    CONTINUE
      DO 27 I=1,N
      J=I
      L=I
      LD=0
12    IF(ABS(A(L,J)))16,13,16
13    LD=1
      IF(L-N)15,14,14
14    WRITE(3,25)
      GOTO 27
15    L=L+1
      GO TO 12
16    IF(LD-1)26,17,26
17    DO 18 K=1,N
      T=A(L,K)
      A(L,K)=A(I,K)
      A(I,K)=T
      T=E(L,K)
      E(L,K)=E(I,K)
18    E(I,K)=T
26    BM=1./A(I,I)
      DO 19 K=1,N
19    E(I,K)=E(I,K)*BM
      I1=I+1
      DO 20 KI=I1,N
      K=N+I1-KI
20    A(I,K)=A(I,K)*BM
      DO 23 L=1,N
      IF(L.EQ.I)GOTO 23
      BM=A(L,I)
      DO 21 K=1,N
21    E(L,K)=E(L,K)-E(I,K)*BM
      I1=I+1
      DO 22 KI=I1,N
      K=N+I1-KI
22    A(L,K)=A(L,K)-A(I,K)*BM
23    CONTINUE
27    CONTINUE
      RETURN
      END
