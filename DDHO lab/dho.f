      subroutine dho
      implicit none
      include 'dho.inc'
      mass=1
      springk=1
      fkinetic=0.01
      bdrag=0
      delta0=0
      omega0=sqrt(springk/mass)
      c0=fkinetic/springk
      gamma=bdrag/(2*mass)
      print *, 'dho package initialized with',
     +         ' mass=',mass,
     +         ' springk=',springk,
     +         ' fkinetic=',fkinetic,
     +         ' bdrag=',bdrag
      period=twopi/omega0
      ncoeffs=0
      end

      real function x(t)
      implicit none
      real t
      include 'dho.inc'
      integer nper,n
      real frac
      real tsub
      frac=t/period-delta0/twopi
      if (frac.lt.0) then
        print *, 'cannot sample x before the start of the first cycle'
        stop
      endif
      nper=frac
      frac=frac-nper
      nper=nper+1
      if (ncoeffs.eq.0) then
        a(1)=1
        b(1)=a(1)-2*c0
        ncoeffs=1
      endif
      if (nper.gt.999999) then
        print *, 'cannot sample past period 999999'
        stop
      endif
      do n=ncoeffs+1,nper
        a(n)=max(b(n-1)-2*c0,0.)
        b(n)=max(a(n)-2*c0,0.)
      enddo
      ncoeffs=nper
      if (frac.lt.0.5) then
        if (a(nper).gt.c0) then
          x=a(nper)*cos(frac*twopi)+c0
        else
          x=0
        endif
      else
        if (b(nper).gt.c0) then
          x=b(nper)*cos(frac*twopi)-c0
        else
          x=0
        endif
      endif
      end
