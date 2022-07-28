        SUBROUTINE GHIONE(X,C,N,M)
*   N - количество линий передачи 1...19                    21-06-99
* X(2N+2) - вектор координат краевых точек проводников
*   М - количество узлов квадратурной формулы Гаусса M=1..1000
*
*   Gl    1     2     ...    N        Gr (экран)
* =====--===--=====--====--======-----==========
*      1 2  3 4   5       2N   2N+1  2N+2
*
*   C(N,N) - матрица емкостей многопроводных линий передачи
*            нормированная для верхней полуплоскости
*
*-----------------------------------------------------------------------
      real*8 C(1), xi(40000),sigma(40000),X(40),    ! (2N+2)*M = 40*1000 = 40000
     *       F(361),G1(361),U1(361),signF,signG,sumF,sumG,
     *       pi,co,pf1,pf2,pf3,pf4,d,xx
       pi=3.14159265359
* ==================== Начало циклов ===================================
             do 1 j=1,N
        signF=(-1)**(N+1-j)*pi/M
        signG=(-1)**(N+2-j)*pi/M
             do 1 i=1,N
        sumF=0.
        sumG=0.
             do 2 k=1,M
        co=cos((2.*k-1.)*pi/(2.*M))
        jk=(k-1)*N+j
        xi(jk)   =(x(2*j+1)+ x(2*j))/2.+(x(2*j+1)- x(2*j))*co/2.
        sigma(jk)=(x(2*j)+ x(2*j-1))/2.+(x(2*j)- x(2*j-1))*co/2.
          pf1=1.
        do 11 mm=1,2*j-1
11      pf1=pf1*sqrt(xi(jk)-x(mm))
          pf2=1.
        do 12 mm=2*j+2,2*N+2
12      pf2=pf2*sqrt(x(mm)-xi(jk))
          pf3=1.
        if(j.ne.1)then
        do 13 mm=1,2*j-2
13      pf3=pf3*sqrt(sigma(jk)-x(mm))
        endif
          pf4=1.
        do 14 mm=2*j+1,2*N+2
14      pf4=pf4*sqrt(x(mm)-sigma(jk))
             sumF=sumF+   xi(jk)**(i-1)/(pf1*pf2)
2            sumG=sumG+sigma(jk)**(i-1)/(pf3*pf4)
        ij=(j-1)*N+i
        F(ij)=signF*sumF
1       G1(ij)=signG*sumG

* ====================== Конец циклов ==================================
        CALL DMINV(G1,N,D)
          do3 i=1,N
          do3 j=1,N
        xx=0.
        ij=(j-1)*N+i
            U1(ij)=0.
            if(i.eq.j)U1(ij)=1.
            if(i+1.eq.j)U1(ij)=-1.
          do4 k=1,N
        ik=(k-1)*N+i
        kj=(j-1)*N+k
4       xx=xx+G1(ik)*F(kj)
3       C(ij)=xx

          do7 i=1,N              ! Уточнениние 21-06-99
          do7 j=1,N
        ij=(j-1)*N+i
7       F(ij)=C(ij)

          do5 i=1,N
          do5 j=1,N
        xx=0.
        ij=(j-1)*N+i
          do6 k=1,N
        ik=(k-1)*N+i
        kj=(j-1)*N+k
6       xx=xx+U1(ik)*F(kj)
5       C(ij)=-xx
        return
        end