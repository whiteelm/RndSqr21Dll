      SUBROUTINE NROOT (M,A,B,XL,X)                                     NROO 370
      DIMENSION A(1),B(1),XL(1),X(1)                                    NROO 380
      DOUBLE PRECISION A,B,XL,X,SUMV                                    NROO 460
      K=1                                                               NROO 610
      DO 100 J=2,M                                                      NROO 620
      L=M*(J-1)                                                         NROO 630
      DO 100 I=1,J                                                      NROO 640
      L=L+1                                                             NROO 650
      K=K+1                                                             NROO 660
  100 B(K)=B(L)                                                         NROO 670
      MV=0                                                              NROO 710
      CALL EIGEN (B,X,M,MV)                                             NROO 720
      L=0                                                               NROO 770
      DO 110 J=1,M                                                      NROO 780
      L=L+J                                                             NROO 790
  110 XL(J)=1.0/ DSQRT( DABS(B(L)))                                     NROO 800
      K=0                                                               NROO 810
      DO 115 J=1,M                                                      NROO 820
      DO 115 I=1,M                                                      NROO 830
      K=K+1                                                             NROO 840
  115 B(K)=X(K)*XL(J)                                                   NROO 850
      DO 120 I=1,M                                                      NROO 890
      N2=0                                                              NROO 900
      DO 120 J=1,M                                                      NROO 910
      N1=M*(I-1)                                                        NROO 920
      L=M*(J-1)+I                                                       NROO 930
      X(L)=0.0                                                          NROO 940
      DO 120 K=1,M                                                      NROO 950
      N1=N1+1                                                           NROO 960
      N2=N2+1                                                           NROO 970
  120 X(L)=X(L)+B(N1)*A(N2)                                             NROO 980
      L=0                                                               NROO 990
      DO 130 J=1,M                                                      NROO1000
      DO 130 I=1,J                                                      NROO1010
      N1=I-M                                                            NROO1020
      N2=M*(J-1)                                                        NROO1030
      L=L+1                                                             NROO1040
      A(L)=0.0                                                          NROO1050
      DO 130 K=1,M                                                      NROO1060
      N1=N1+M                                                           NROO1070
      N2=N2+1                                                           NROO1080
  130 A(L)=A(L)+X(N1)*B(N2)                                             NROO1090
      CALL EIGEN (A,X,M,MV)                                             NROO1130
      L=0                                                               NROO1140
      DO 140 I=1,M                                                      NROO1150
      L=L+I                                                             NROO1160
  140 XL(I)=A(L)                                                        NROO1170
      DO 150 I=1,M                                                      NROO1210
      N2=0                                                              NROO1220
      DO 150 J=1,M                                                      NROO1230
      N1=I-M                                                            NROO1240
      L=M*(J-1)+I                                                       NROO1250
      A(L)=0.0                                                          NROO1260
      DO 150 K=1,M                                                      NROO1270
      N1=N1+M                                                           NROO1280
      N2=N2+1                                                           NROO1290
  150 A(L)=A(L)+B(N1)*X(N2)                                             NROO1300
      L=0                                                               NROO1310
      K=0                                                               NROO1320
      DO 180 J=1,M                                                      NROO1330
      SUMV=0.0                                                          NROO1340
      DO 170 I=1,M                                                      NROO1350
      L=L+1                                                             NROO1360
  170 SUMV=SUMV+A(L)*A(L)                                               NROO1370
  175 SUMV= DSQRT(SUMV)                                                 NROO1380
      DO 180 I=1,M                                                      NROO1390
      K=K+1                                                             NROO1400
  180 X(K)=A(K)/SUMV                                                    NROO1410
      RETURN                                                            NROO1420
      END                                                               NROO1430
