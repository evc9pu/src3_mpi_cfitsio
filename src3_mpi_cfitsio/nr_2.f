      FUNCTION erfc(x)
      real(8) erfc,x
CU    USES gammp,gammq
      real(8) gammp,gammq
      if(x.lt.0.d0)then
        erfc=1.d0+gammp(.5d0,x**2)
      else
        erfc=gammq(.5d0,x**2)
      end if
      return
      END
      
      FUNCTION gammln(xx)
      real(8) gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
      
      FUNCTION gammp(a,x)
      real(8) a,gammp,x
CU    USES gcf,gser
      real(8) gammcf,gamser,gln
      if(x.lt.0.d0.or.a.le.0.d0) then
         print*,'bad arguments in gammp'
         stop
      end if
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.d0-gammcf
      end if
      return
      END
      
      FUNCTION gammq(a,x)
      real(8) a,gammq,x
CU    USES gcf,gser
      real(8) gammcf,gamser,gln
      if(x.lt.0.d0.or.a.le.0.d0) then
         print*,'bad arguments in gammq'
         stop
      end if
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammq=1.d0-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      end if
      return
      END
      
      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      real(8) a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
CU    USES gammln
      INTEGER i
      real(8) an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.d0-a
      c=1.d0/FPMIN
      d=1.d0/b
      h=d
      do 11 i=1,ITMAX
        an=-dble(i)*(dble(i)-a)
        b=b+2.d0
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=d*c
        h=h*del
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      print*,'a too large, ITMAX too small in gcf'
      stop
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END
      
      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      real(8) a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)
CU    USES gammln
      INTEGER n
      real(8) ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.d0)then
        if(x.lt.0.d0) then
           print*, 'x < 0 in gser'
           stop
        end if
        gamser=0.d0
        return
      end if
      ap=a
      sum=1.d0/a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.d0
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      print*, 'a too large, ITMAX too small in gser'
      stop
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
      
      SUBROUTINE locate(xx,n,x,j)
c     from numerical recipes
c     searches an ordered table, using bisection
      INTEGER j,n
      real(8) x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
!         print*,jl,ju
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        end if
      goto 10
      end if
      j=jl
      if(j==0) then
        if(x==xx(1)) j=1
      end if
      return
      END

      SUBROUTINE locate_int(xx,n,x,j)
c     from numerical recipes
c     searches an ordered table, using bisection
      INTEGER j,n
      INTEGER x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
!         print*,jl,ju
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        end if
      goto 10
      end if
      j=jl
      if(j==0) then
        if(x==xx(1)) j=1
      end if
      return
      END
      
      SUBROUTINE locate2(xx,nn,n,x,j)
c     from numerical recipes
c     searches an ordered table, using bisection

c     baw 8/1/2000, allow dimensioned array size to be bigger than size used.
c     dimensioned array size is nn; usable is n.

      INTEGER j,n,nn
      real(8) x,xx(nn)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        end if
      goto 10
      end if
      j=jl
      if(j==0) then
        if(x==xx(1)) j=1
      end if
      return
      END
          
      SUBROUTINE qsimp(func,a,b,s,converg)
      INTEGER JMAX
      real(8) a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-5, JMAX=20)
C     USES trapzd
      INTEGER j
      real(8) os,ost,st
      LOGICAL converg
      converg=.false.
      ost=-1.d30
      os= -1.d30
      do 11 j=1,JMAX
        call trapzd(func,a,b,st,j)
        s=(4.d0*st-ost)/3.d0
        if (abs(s-os).lt.EPS*abs(os)) return
        os=s
        ost=st
11    continue
      converg=.true.
c      print*,'too many steps in qsimp'
c      pause 'too many steps in qsimp'
      END

      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      real(8) a,b,s,func
      EXTERNAL func
      INTEGER it,j
      real(8) del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=dble(it)
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      end if
      return
      END
      
      SUBROUTINE qsimp1(func,a,b,s,converg)
      INTEGER JMAX
      real(8) a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-5, JMAX=20)
C     USES trapzd
      INTEGER j
      real(8) os,ost,st
      LOGICAL converg
      converg=.false.
      ost=-1.d30
      os= -1.d30
      do 11 j=1,JMAX
        call trapzd1(func,a,b,st,j)
        s=(4.d0*st-ost)/3.d0
        if (abs(s-os).lt.EPS*abs(os)) return
        os=s
        ost=st
11    continue
      converg=.true.
c      print*,'too many steps in qsimp'
c      pause 'too many steps in qsimp'
      END

      SUBROUTINE trapzd1(func,a,b,s,n)
      INTEGER n
      real(8) a,b,s,func,s1,xexp
      EXTERNAL func
      INTEGER it,j
      real(8) del,sum,tnm,x,dx,x1
      
      xexp=2.5d0
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=dble(it)
        del=(b-a)/tnm**xexp
        s1=0.d0
        do 11 j=1,it
           if (func(a).gt.func(b)) then
              x=a+dx*j**xexp
              x1=a+dx*(j-1)**xexp
           else
              x=b-dx*j**xexp
              x1=b-dx*(j-1)**xexp
           endif
           s1=s1+func((x+x1)/2.d0)*abs(x-x1)
 11     continue
        s=s1
      end if
      return
      END
      

