c
      subroutine av14pw(lemp,emcoul,l,s,j,t1z,t2z,r,vpw)
      implicit real*8(a-h,o-z)
      integer t1z, t2z, s
      double precision vpw(2,2)
      logical emcoul
c
      ipte=1
      ipls=1
      ipqq=1
      ipbb=1
      ipsb=1
      ur=1.d0/r
c
      vpw=0
c
      if ((t1z+t2z).lt.0) inn=1
      if ((t1z+t2z).eq.0) inn=2
      if ((t1z+t2z).gt.0) inn=3
c
      call potl(r,ur,ipte,ipls,ipqq,ipbb,ipsb,emcoul,lemp,j,
     >          vpp1,vpp2,vpp3,vpp4,vpp5,
     >          vpn1,vpn2,vpn3,vpn4,vpn5,
     >          vnn1,vnn2,vnn3,vnn4,vnn5)
c
      if(inn.eq.1) then
         vx1=vpp1
         vx2=vpp2
         vx3=vpp3
         vx4=vpp4
         vx5=vpp5
      else if(inn.eq.2) then
         vx1=vpn1
         vx2=vpn2
         vx3=vpn3
         vx4=vpn4
         vx5=vpn5
      else if(inn.eq.3) then
         vx1=vnn1
         vx2=vnn2
         vx3=vnn3
         vx4=vnn4
         vx5=vnn5
      end if
c
      if (j.eq.l.and.s.eq.0) then
        vpw(1,1)=vx1
        return
      else if (j.eq.l.and.s.eq.1) then
        vpw(1,1)=vx2
        return
      else 
        if (j.eq.0) then
          vpw(1,1)=vx4
          return
        endif
        vpw(1,1)=vx3
        vpw(1,2)=vx5
        vpw(2,1)=vx5
        vpw(2,2)=vx4
        return
      endif
      stop "ERROR in av14pw"
      return
      end
c
c
c
       subroutine pot(r,ur,ipls,ipqq,ipbb,
     >                potc,pott,pots,potm,ptec,ptet,plsc,plst,
     >                pqqc,pqqt,pqqs,pqqm,pbbc,pbbt,
     >                psb1,psb2,psb3,psb4)
       implicit real*8(a-h,o-z)
c
c      v= potc+pots sigmaij+pott tauij+potm sigmaij tauij +
c         (ptec + ptet tauij) sij + (plsc + plst tauij) lsij +
c         (pqqc+pqqs sigmaij+pqqt tauij+pqqm sigmaij tauij)lij**2 +
c         (pbbc + pbbt tauij) (lsij)**2 +
c         (psb1 + psb2 sigmaij + psb3 sij ) tij + 
c          psb4 (tauiz+taujz)
c
c      ipls  = 1 to include the ls potential
c      ipqq  = 1 to include the l**2 potential
c      ipbb  = 1 to include the (ls)**2 potential
c
      parameter (pi=3.141592653589793d0)
      parameter (dmu=0.7d0,hh=10.463d0,dmu2=.699488167d0)
      parameter (udmu=1.d0/dmu,udmu2=1.d0/dmu2)
      parameter (u2=0.5d0,u4=0.25d0,u16=0.0625d0,u3=1.d0/3.d0)
c
c
      potc=0.d0
      pots=0.d0
      pott=0.d0
      potm=0.d0
      ptec=0.d0
      ptet=0.d0
      plsc=0.d0
      plst=0.d0
      pqqc=0.d0
      pqqt=0.d0
      pqqs=0.d0
      pqqm=0.d0
      pbbc=0.d0
      pbbt=0.d0
      psb1=0.d0
      psb2=0.d0
      psb3=0.d0
      psb4=0.d0
      if(r.gt.90.d0) return

      x=r*dmu2
      ux=ur*udmu2
      x2=2.d0*r*r
      if(x2.gt.40.d0) cutoff=1.d0
      if(x2.le.40.d0) cutoff=1.d0-dexp(-x2)
      yx=dexp(-x)*ux
      yp=yx*cutoff
      tp=(1.d0+3.d0*ux+3.d0*ux*ux)*yp*cutoff
      tp2=tp*tp
      ww=1.d0/(1.d0+dexp((r-0.5d0)*5.d0))
      potc=-4.801125d0*tp2 +2061.5625d0*ww
      pott= 0.798925d0*tp2 - 477.3125d0*ww
      pots= 1.189325d0*tp2 - 502.3125d0*ww
      potm= 0.182875d0*tp2 +  97.0625d0*ww +3.72681d0*yp
      ptec=-0.1575d0  *tp2 + 108.75d0  *ww
      ptet=-0.7525d0  *tp2 + 297.25d0  *ww +3.72681d0*tp
      if(ipls.eq.0) return
      plsc= 0.5625d0  *tp2 - 719.75d0  *ww
      plst= 0.0475d0  *tp2 - 159.25d0  *ww
      if(ipqq.eq.0) return
      pqqc= 0.070625d0*tp2 +   8.625d0 *ww
      pqqt=-0.148125d0*tp2 +   5.625d0 *ww
      pqqs=-0.040625d0*tp2 +  17.375d0 *ww
      pqqm=-0.001875d0*tp2 -  33.625d0 *ww
      if(ipbb.eq.0) return
      pbbc=-0.5425d0  *tp2 + 391.0d0   *ww
      pbbt= 0.0025d0  *tp2 + 145.0d0   *ww
c
      return
      end
c
c
c
       subroutine potl(r,ur,ipte,ipls,ipqq,ipbb,ipsb,emcoul,lemp,j,
     >                 vpp1,vpp2,vpp3,vpp4,vpp5,
     >                 vpn1,vpn2,vpn3,vpn4,vpn5,
     >                 vnn1,vnn2,vnn3,vnn4,vnn5)
       implicit real*8(a-h,o-z)
c------------------------------------------------------------------------
c      nuclear potentials calculated in the point r (ur=1/r)
c      between states with total angular momentum j;
c      vpp=proton-proton potentials = <T=1 Tz=1|V|T=1 Tz=1>
c      vpn=proton-neutron potentials = <T Tz=0|V|T Tz=0>
c      vnn=neutron-neutron potentials = <T=1 Tz=-1|V|T=1 Tz=-1>
c
c      n.b. <T Tz|V|T' Tz'>=0 unless T=T' and Tz=Tz'
c
c      for a given (2S+1)L_J state, T is fixed by the antisymmetry
c      condition (-)**(L+S+T)=-1; clearly for bra,ket states with T=0
c      vpp=vnn=0 and vpn= <T=0 Tz=0|V|T=0 Tz=0>; for bra,ket states
c      of T=1 vpn= <T=1 Tz=0|V|T=1 Tz=0>;
c
c      vNN1=V(1J_J)
c      vNN2=V(3J_J)
c      vNN3=V(3[J-1]_J)
c      vNN4=V(3[J+1]_J)
c      vNN5=V(3[J-1]_J - 3[J+1]_J)
c
c      ipte  = 1 to include the tensor potential
c      ipls  = 1 to include the ls potential
c      ipqq  = 1 to include the l**2 potential
c      ipbb  = 1 to include the (ls)**2 potential
c      ipsb  = 1 to include the isospin symmetry breaking terms
c      emcoul= .TRUE. to include the coulomb potential
c      lemp  = 1-14 to include the e.m. potential for the av18
c------------------------------------------------------------------------
c
      data icont/0/
c
c------------------------------------------------------------------------
c
c    coefficienti per potenziale coulombiano
      parameter (eqmf=197.327053d0/137.035989d0)
c--------------------------------------------------------------
c
      logical emcoul
      dimension vem(14),lts0(0:1),lts1(0:1),ltin(0:1)
c
      data zero/0.d0/,one/1.d0/,u2/2.d0/,u3/3.d0/,u4/4.d0/,u6/6.d0/
c
      data lts0/1,0/     !lts0(mod(L))=value of T: case S=0
      data lts1/0,1/     !lts1(mod(L))=value of T: case S=1
      data ltin/1,0/     !ltin(T=0)=1 and ltin(T=1)=0
c
      jp=mod(j,2)
      uj=one/dfloat(2*j+1)
      xjj=j*(j+1)
c
      vpp1=zero
      vpp2=zero
      vpp3=zero
      vpp4=zero
      vpp5=zero
      vpn1=zero
      vpn2=zero
      vpn3=zero
      vpn4=zero
      vpn5=zero
      vnn1=zero
      vnn2=zero
      vnn3=zero
      vnn4=zero
      vnn5=zero
c
      call pot(r,ur,ipls,ipqq,ipbb,
     >           potc,pott,pots,potm,ptec,ptet,plsc,plst,
     >           pqqc,pqqt,pqqs,pqqm,pbbc,pbbt,
     >           psb1,psb2,psb3,psb4)
      if(ipte.eq.0) then
        ptec=zero
        ptet=zero
      end if
      if(ipls.eq.0) then
        plsc=zero
        plst=zero
      end if
      if(ipqq.eq.0) then
        pqqc=zero
        pqqt=zero
        pqqs=zero
        pqqm=zero
      end if
      if(ipbb.eq.0) then
        pbbc=zero
        pbbt=zero
      end if
      if(ipsb.eq.0) then
        psb1=zero
        psb2=zero
        psb3=zero
        psb4=zero
      end if
      icont=icont+1
      if(icont.eq.1) then
        write(6,*) ipte,ipls,ipqq,ipbb,ipsb
        write(6,*) potc,pott,pots,potm
        write(6,*) ptec,ptet,plsc,plst
        write(6,*) pqqc,pqqt,pqqs,pqqm
        write(6,*) pbbc,pbbt
        write(6,*) psb1,psb2,psb3,psb4
      end if
c
      vcoul=zero
      if(emcoul) vcoul=eqmf*ur        
      if(lemp.gt.0) call empot(lemp,r,vem)
c
c=======case 1 si=sj=0 li=lj=j
c
      lti=lts0(jp)
      p1=potc-u3*pots+xjj*(pqqc-u3*pqqs)
      p2=pott-u3*potm+xjj*(pqqt-u3*pqqm)
      p3=psb1-u3*psb2
      p4=psb4
c
      if(lti.eq.1) then
        vempp=vcoul+vem(1)+vem(2)+vem(3)+vem(4)-u3*vem(6)
        vempn=vem(5)-u3*vem(8)
        vemnn=-u3*vem(7)
        vpp1=p1+p2+u2*(p3+p4)+vempp
        vpn1=p1+p2-u4*p3+vempn
        vnn1=p1+p2+u2*(p3-p4)+vemnn
      else
        vempn=vem(5)-u3*vem(8)
        vpn1=p1-u3*p2+vempn
        vnn1=0.d0
        vpp1=0.d0
      end if
c
c=======case 2 si=sj=1 li=lj=j
c
      lti=lts1(jp)
      if(j.gt.0) then
        p1=potc+pots+u2*ptec-plsc+xjj*(pqqc+pqqs)+pbbc
        p2=pott+potm+u2*ptet-plst+xjj*(pqqt+pqqm)+pbbt
        p3=psb1+psb2+u2*psb3
        p4=psb4
c
        if(lti.eq.1) then
          vempp=vcoul+vem(1)+vem(2)+vem(3)+vem(4)+vem(6)
     >           +u2*vem(9)-vem(12)
          vempn=vem(5)+vem(8)+u2*vem(11)-vem(14)
          vemnn=vem(7)+u2*vem(10)-vem(13)
          vpp2=p1+p2+u2*(p3+p4)+vempp
          vpn2=p1+p2-u4*p3+vempn
          vnn2=p1+p2+u2*(p3-p4)+vemnn
        else
          vempn=vem(5)+vem(8)+u2*vem(11)-vem(14)
          vpp2=0.d0
          vpn2=p1-u3*p2+vempn
          vnn2=0.d0
        end if
      end if
c
c=======case 3 si=sj=1 li=lj=j-1
c
      lti=ltin(lti)
c
      if(j.gt.0) then
        xte=-u2*(j-1)*uj
        xls=j-1
        xll=j*(j-1)
        xls2=xls*xls
        p1=potc+pots+xte*ptec+xls*plsc+xll*(pqqc+pqqs)+xls2*pbbc
        p2=pott+potm+xte*ptet+xls*plst+xll*(pqqt+pqqm)+xls2*pbbt
        p3=psb1+psb2+xte*psb3
        p4=psb4
c
        if(lti.eq.1) then
          vempp=vcoul+vem(1)+vem(2)+vem(3)+vem(4)+vem(6)
     >           +xte*vem(9)+xls*vem(12)
          vempn=vem(5)+vem(8)+xte*vem(11)+xls*vem(14)
          vemnn=vem(7)+xte*vem(10)+xls*vem(13)
          vpp3=p1+p2+u2*(p3+p4)+vempp
          vpn3=p1+p2-u4*p3+vempn
          vnn3=p1+p2+u2*(p3-p4)+vemnn
        else
          vempn=vem(5)+vem(8)+xte*vem(11)+xls*vem(14)
          vpn3=p1-u3*p2+vempn
          vnn3=0.d0
          vpp3=0.d0
        end if
      end if
c
c=======case 4 si=sj=1 li=lj=j+1
c
      xte=-u2*(j+2)*uj
      xls=-j-2
      xll=(j+1)*(j+2)
      xls2=xls*xls
      p1=potc+pots+xte*ptec+xls*plsc+xll*(pqqc+pqqs)+xls2*pbbc
      p2=pott+potm+xte*ptet+xls*plst+xll*(pqqt+pqqm)+xls2*pbbt
      p3=psb1+psb2+xte*psb3
      p4=psb4
c
      if(lti.eq.1) then
        vempp=vcoul+vem(1)+vem(2)+vem(3)+vem(4)+vem(6)
     >         +xte*vem(9)+xls*vem(12)
        vempn=vem(5)+vem(8)+xte*vem(11)+xls*vem(14)
        vemnn=vem(7)+xte*vem(10)+xls*vem(13)
        vpp4=p1+p2+u2*(p3+p4)+vempp
        vpn4=p1+p2-u4*p3+vempn
        vnn4=p1+p2+u2*(p3-p4)+vemnn
      else
        vempn=vem(5)+vem(8)+xte*vem(11)+xls*vem(14)
        vpn4=p1-u3*p2+vempn
        vnn4=0.d0
        vpp4=0.d0
      end if
c
c=======case 5 si=sj=1 li=j-1 lj=j+1
c
      if(j.gt.0) then
        xte=u6*dsqrt(xjj)*uj
        p1=xte*ptec
        p2=xte*ptet
        p3=xte*psb3
        p4=zero
c
        if(lti.eq.1) then
          vempp=xte*vem(9)
          vempn=xte*vem(11)
          vemnn=xte*vem(10)
          vpp5=p1+p2+u2*(p3+p4)+vempp
          vpn5=p1+p2-u4*p3+vempn
          vnn5=p1+p2+u2*(p3-p4)+vemnn
        else
          vempn=xte*vem(11)
          vpn5=p1-u3*p2+vempn
          vnn5=0.d0
          vpp5=0.d0
        end if
      end if
c
      return
      end

