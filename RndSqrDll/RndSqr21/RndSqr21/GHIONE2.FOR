        SUBROUTINE GHIONE2(X,C,N,M) 
! Расчет многопроводных копланарных полосок (МКПП) 	       31.03.2021						27-03-2021
! Multiconductor Coplanar Stripline (MCPS)	
* N - количество линий передачи N = 1..19                       
* X(2N+2) - вектор координат краевых точек проводников  4..40
* М - количество узлов квадратурной формулы Гаусса M=1..1000
*
*       1     2    ...      N      N+1(экран)
* -----===---===--====----======----========-------
*      1 2   3 4        2N-1 2N   2N+1  2N+2
*
*   C(N,N) - матрица емкостей многопроводных линий передачи
*            нормированная для верхней полуплоскости 
*
*-----------------------------------------------------------------------
      real*8 X(1), C(1), xi,sigma,          ! 19 x 19 = 361    Nmax = 19
     *       F1(361),G(361),L1(361),sumF,sumG,
     *       pi,co,pf1,pf2,pf3,pf4,d,xx
       pi=3.14159265359
* ==================== Начало циклов ===================================
      do 1 j=1,N
      do 1 i=1,N
        sumF=0.d0
        sumG=0.d0

            do 2 k=1,M
             co = dcos((2.*k-1.)*pi/(2.*M))
             xi    = (x(2*j+1)+ x(2*j))/2. + (x(2*j+1)- x(2*j))*co/2.
             sigma = (x(2*j)+ x(2*j-1))/2. + (x(2*j)- x(2*j-1))*co/2.
         pf1=1.
         do 11 mm=1,2*j-1
11       pf1=pf1*dsqrt(xi-x(mm))
             pf2=1.
             do 12 mm=2*j+2,2*N+2
12           pf2=pf2*dsqrt(x(mm)-xi)
         pf3=1.
         if(j.ne.1)then
         do 13 mm=1,2*j-2
13       pf3=pf3*dsqrt(sigma-x(mm))
         endif
              pf4=1.
              do 14 mm=2*j+1,2*N+2
14            pf4=pf4*dsqrt(x(mm)-sigma)

         sumF=sumF+   xi**(i-1)/(pf1*pf2)
2        sumG=sumG+sigma**(i-1)/(pf3*pf4)
        ij=(j-1)*N+i
        F1(ij) =((-1)**(N+1-j)*pi/M)*sumF
1       G(ij)=((-1)**(N+2-j)*pi/M)*sumG

* ====================== Конец циклов ==================================
        CALL DMINV(F1,N,D)
          do3 i=1,N
          do3 j=1,N
        xx=0.
        ij=(j-1)*N+i
            L1(ij)=0.
            if(i.eq.j)L1(ij)=1.
            if(i.eq.j+1)L1(ij)=-1.
          do4 k=1,N
        ik=(k-1)*N+i
        kj=(j-1)*N+k
4       xx=xx+F1(ik)*G(kj)
3       C(ij)=xx

          do7 i=1,N              
          do7 j=1,N
        ij=(j-1)*N+i
7       F1(ij)=C(ij)

        do5 i=1,N
        do5 j=1,N
            xx=0.
            ij=(j-1)*N+i
         do6 k=1,N
             ik=(k-1)*N+i
             kj=(j-1)*N+k
6       xx=xx+L1(ik)*F1(kj)
5       C(ij)=-xx
        return
        end