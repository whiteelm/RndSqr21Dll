    subroutine IRndSqr21(a, b, h1, h2, s, s0, d1, d2, er, dC, dL, Um, em, dZ0)
    !dec$ attributes stdcall, dllexport :: IRndSqr21
    !DEC$ ATTRIBUTES VALUE :: a, b, h1, h2, s, s0, d1, d2, er
    !DEC$ ATTRIBUTES REFERENCE :: dC, dL, Um, em, dZ0
    implicit complex*16(c,w,z), real*8(a-b,d-h,o-v,x-y)
    dimension dC(2,2), dL(2,2), Um(2,2), em(2), dZ0(2,2), dCC(2,2), dLL(2,2)
    call RndSqr21(a, b, h1, h2, s, s0, d1, d2, er, dC, dL)
    n=2
    dCC=dC 
    dLL=dL
    call dminv(dLL,n,ad)
    call nroot(n,dCC,11.127*dLL,em,Um)
    call impedance(n,dC,Um,em,dZ0)
    dC=transpose(dC);
    dL=transpose(dL);
    Um=transpose(Um);
    dZ0=transpose(dZ0);
    end