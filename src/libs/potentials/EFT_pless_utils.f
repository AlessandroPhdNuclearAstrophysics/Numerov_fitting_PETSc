      subroutine lagm(g,j0,jf,r0,h,rin,in0,nin,gin,i0)
      implicit real*8(a-h,o-z)
c       la funcion g es tabulada con paso h a partir de g(j0)=g(r0)
c       hasta g(jf) con jf-j0 mayor de 4. gin es la funcion inter-
c       polada en los nin puntos de abcisa rin, el primer valor es
c       rin(in0) al que corresponde gin(in0). rin(in0), ...,
c       rin(in0+nin-1) debe estar comprendida entre r0 y r0+h*(jf-j0).
c
        dimension g(1),rin(1),gin(1)
c
        n=jf-j0+1
        if(n-6) 1,10,10
 1      write(6,1000)
        stop
 10     r=r0+h*(n-1)
        if(rin(in0).lt.r0.or.rin(in0+nin-1).gt.r) then
           write(*,*)rin(in0),r0,rin(in0+nin-1),r
           goto 100
        endif
        do 80 i=1,nin
        amr=(rin(in0+i-1)-r0)/h
        mr=amr+1
        mr=amin0(mr,n-3)
        if(mr-3) 20,20,30
 20     p=amr-2.d0
        g0=g(j0)
        l=j0+1
        goto 50
 30     l=mr+j0-2
        p=amr-(l-j0+1)
        g0=g(l-1)
 50     p3=p-3.d0
        p2=p-2.d0
        p1=p-1.d0
        p4=p+1.d0
        p5=p+2.d0
        p23=p3*p5
        p12=p4*p2
        pp1=p*p1
 80     gin(i0+i-1)=(pp1*(0.1d0*p12*(p5*g(l+4)-p3*g0)+0.5d0*p23*
     >              (p2*g(l)-p4*g(l+3)))+p12*p23*(-p1*g(l+1)+p*
     >               g(l+2)))/12.d0
        return
 100    write(6,1100)
        stop
 1000   format(1x,'   n demasiado pequeo')
 1100   format(1x,'   punto fuera del rango')
        end
c
c
         double precision function b5(jb,m1,m2,m3,h1,h2,h3,a,ias)
         implicit real*8(a-h,o-z)
c       jb=1,2,3 para integrar en 1,2,3 bloques a partir del punto
c       de abcisas ias.
c       m1,m2,m3 numero de pasos en cada bloque, h1,h2,h3 paso corres
c       pondiente. hi=hi/22.5
c
        dimension a(1)
c
        b5=0.d0
        n=m1
        h5=h1
        j0=ias
        ir=jb-2
 1      if(n-4) 44,2,2
 2      hl5=h5/32.d0
        nb=n/4
        n1=n-4*nb+1
        s2=0.d0
        s3=s2
        s4=s2
        l=j0
        do 4 i=1,nb
        s2=s2+a(l+2)
        s3=s3+a(l+1)+a(l+3)
        s4=s4+a(l)+a(l+4)
 4      l=l+4
        goto (8,12,16,20),n1
 8      ax=0.d0
        goto 40
 12     ax=hl5*(-19.d0*a(l-3)+106.d0*a(l-2)-264.d0*a(l-1)+646.d0*a(l)+
     >          251.d0*a(l+1))
        goto 40
 16     ax=hl5*(-8.d0*a(l-2)+32.d0*a(l-1)+192.d0*a(l)+992.d0*a(l+1)+
     >           232.d0*a(l+2))
        goto 40
 20     ax=hl5*(-27.d0*a(l-1)+378.d0*a(l)+648.d0*a(l+1)+918.d0*a(l+2)+
     >           243.d0*a(l+3))
 40     b5=b5+ax+h5*(7.d0*s4+32.d0*s3+12.d0*s2)
        if(ir) 60,50,55
 50     j0=m1+ias
        n=m2
        h5=h2
        ir=-1
        goto 1
 55     j0=m2+m1+ias
        n=m3
        h5=h3
        ir=0
        goto 1
 60     return
c
 44     write(4,200)
 200    format(//1x,'  numero de intervalos menor que 4'/)
        stop
        end




c----------------------------------------------------------------------------------
      real*8 function dgamma( x )
c----------------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      save
c
      pi=4.d0*datan( 1.d0 )
      z=x
      n=idint( x )
c
      if( dabs( x-dble( n ) ) .lt. 1.d-8 ) then
      f=1.d0
c
 10   if( z .lt. 2.d0 ) go to 20
      z=z-1.d0
      f=f*z
      go to 10
c
 20   dgamma=f
      else
c
      f=1.d0
c
      do 30 i=1,n
      f=f*0.5d0*(2*i-1.d0)
 30   continue
c
      dgamma=f*dsqrt(pi)
      endif
c
      return
      end



c_______________________________________________________________________
      subroutine gauss( ng , ainf , asup , xg , wg )
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      dimension xg(ng),wg(ng)
c
      data eps / 1.d-12 /
c
      apl=.5d0*( asup+ainf )
      amn=.5d0*( asup-ainf )
c
      c=.5d0*dacos(-1.d0)/dble( 2*ng+1 )
c
      do 30 i=1,( ng+1 )/2
      x=dcos( c*dble( 4*i-1 ) )
c
10    p0=1.d0
      p1=0.d0
c
      do 20 j=1,ng
      c1=2.d0-1.d0/dble(j)
      c2=1.d0-1.d0/dble(j)
      p2=p1
      p1=p0
      p0=x*c1*p1-c2*p2
20    continue
c
      dp0=dble( ng )*( x*p0-p1 )/( x**2-1.d0 )
      x0=x
      x=x0-p0/dp0
c
      xchk=x-x0
      if( xchk .lt. 0.d0 ) xchk=-xchk
      if( xchk .gt. eps ) go to 10
c
      xg(i     )=apl-amn*x
      xg(ng+1-i)=apl+amn*x
      wg(i     )=2.d0*amn/( (1.d0-x**2)*dp0**2 )
      wg(ng+1-i)=2.d0*amn/( (1.d0-x**2)*dp0**2 )
30    continue
c
      return
      end


c
      subroutine ng(vec,h,n)
      implicit real*8(a-h,o-z)
c
c     accelerazione del processo di convergenza
c
      dimension vec(4,200),dff(3,200),a(200),b(200),c(6)
      n1=n-1
      pi4=4.d0*3.1415925635898d0
      h1=h/22.5d0
      do 10 i=1,n
      dn=vec(4,i)-vec(3,i)
      dn1=vec(3,i)-vec(2,i)
      dn2=vec(2,i)-vec(1,i)
      dff(1,i)=dn-dn1
      dff(2,i)=dn-dn2
      dff(3,i)=dn
 10   continue
c
      do 20 j=1,3
      do 30 i=1,n
      a(i)=dff(1,i)*dff(j,i)
      if(j.gt.1) b(i)=dff(2,i)*dff(j,i)
 30   continue
      c(j)=pi4*b5(1,n1,0,0,h1,0.d0,0.d0,a,1)
      if(j.gt.1) c(j+2)=pi4*b5(1,n1,0,0,h1,0.d0,0.d0,b,1)
 20   continue
c
      det=c(1)*c(4)-c(2)**2
      if(dabs(det).lt.(1.d-30)) go to 200
      c1=(c(3)*c(4)-c(2)*c(5))/det
      c2=(c(1)*c(5)-c(2)*c(3))/det
      do 100 i=1,n
      x4=vec(4,i)
      x3=vec(3,i)
      x2=vec(2,i)
      vec(4,i)=(1.d0-c1-c2)*x4+c1*x3+c2*x2
 100  continue
c
 200  return
      end
c
c

c
c_______________________________________________________________________
      real*8 function xg( l1 , m1 , l2 , m2 , l , m )
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
c
c     clebsch-gordan coefficients
c
c     n.b. l1, m1, l2, m2, l, m are doubled
c
c
      xg=0.d0
c
      if( l1 .lt. 0 .or. l2 .lt. 0 .or. l .lt. 0 ) return
      if( m1+m2 .ne. m ) return
      if( l .lt. iabs(l1-l2) .or. l .gt. (l1+l2) ) return
c
      if( iabs(m1) .gt. l1 ) return
      if( iabs(m2) .gt. l2 ) return
      if( iabs(m ) .gt. l  ) return
c
      ll=( l1-l2+m )/2
      sgn=dble( (-1)**ll )
c
      xg=sgn*dsqrt( dble(l+1) )*threej( l1 , l2 , l , m1 , m2 ,-m )
c
      return
      end
c
c
c
c_______________________________________________________________________
      real*8 function threej( j1 , j2 , j , m1 , m2 , m )
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      integer z,zmin,zmax
      save
c
c
c     three-j symbol--j's and m's ====> 2j's and 2m's
c
c     it must be ( j1+j2+j )/2+1 < nf
c
c
      parameter( nf=78 )
c
      common/facc/xfac(0:nf)
c
      cc=0.d0
c
      if( m1+m2 .ne.-m .or. iabs(m1) .gt. iabs(j1)
     x                 .or. iabs(m2) .gt. iabs(j2)
     x                 .or. iabs(m ) .gt. iabs(j )
     x                 .or. j .gt. j1+j2
     x                 .or. j .lt. iabs(j1-j2) ) go to 20
c
      zmin=0
      if( j-j2+m1 .lt. 0 ) zmin=-j+j2-m1
      if( j-j1-m2+zmin .lt. 0 ) zmin=-j+j1+m2
c
      zmax=j1+j2-j
      if( j2+m2-zmax .lt. 0 ) zmax=j2+m2
      if( j1-m1-zmax .lt. 0 ) zmax=j1-m1
c
      do 10 z=zmin,zmax,2
      ja=z/2
      jb=( j1+j2-j-z )/2
      jc=( j1-m1-z )/2
      jd=( j2+m2-z )/2
      je=( j-j2+m1+z )/2
      jf=( j-j1-m2+z )/2
c
      sgn=dble( (-1)**(z/2) )
      cc=cc+sgn/( xfac(ja)*xfac(jb)*xfac(jc)*
     x            xfac(jd)*xfac(je)*xfac(jf) )
10    continue
c
      ja=( j1+j2-j )/2
      jb=( j1-j2+j )/2
      jc=(-j1+j2+j )/2
c
      jd=( j1+m1 )/2
      je=( j1-m1 )/2
      jf=( j2+m2 )/2
      jg=( j2-m2 )/2
c
      jh=( j-m )/2
      ji=( j+m )/2
c
      jj=( j1+j2+j+2 )/2
c
      if( jj .gt. nf ) then
      write(6,*) 'error: max factorial exceeded'
      stop
      endif
c
      cc=dsqrt( dble(j+1)*xfac(ja)*xfac(jb)*
     x                    xfac(jc)*xfac(jd) )*cc
c
      cc=dsqrt( xfac(je)*xfac(jf)*xfac(jg)*
     x          xfac(jh)*xfac(ji)/xfac(jj) )*cc
c
20    threej= dble( (-1)**( ( j1-j2-m )/2 ) )
     x       *cc/dsqrt( 10.d0*dble( j+1 ) )
c
      return
      end
c
c
