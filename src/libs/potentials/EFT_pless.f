      subroutine setpot(ilb,itz)
      implicit real*8(a-h,o-z)
      save
c
      include 'params'
c
      parameter( j0=0 , jx=3 , nf=78 , ngx=101, mpx=41 )
c
      common/facc/xfac(0:nf)
      common/gyx/ygx(ngx),wgy(ngx)
      common/eftctc_st/vc (0:1,0:1),vt  (0:1,0:1),
     x                 vb (0:1,0:1),vq  (0:1,0:1),
     x                 vp (0:1,0:1),vpd (0:1,0:1),
     x                 vtp(0:1,0:1),vtpd(0:1,0:1),vbb(0:1,0:1),
     x                 vcv(0:1,0:1),vtv (0:1,0:1),vbv(0:1,0:1),
     x                 vct(0:1,0:1),vtt (0:1,0:1),vbt(0:1,0:1)
      common/eftpot_rv/vcn (nr+1,0:1,0:1),vtn  (nr+1,0:1,0:1),
     x                 vbn (nr+1,0:1,0:1),vqn  (nr+1,0:1,0:1),
     x                 vpn (nr+1,0:1,0:1),vpdn (nr+1,0:1,0:1),
     x                 vtpn(nr+1,0:1,0:1),vtpdn(nr+1,0:1,0:1),
     x                 vbbn(nr+1,0:1,0:1),
     x                 vcvn(nr+1,0:1,0:1),vtvn (nr+1,0:1,0:1),
     x                 vbvn(nr+1,0:1,0:1),
     x                 vctn(nr+1,0:1,0:1),vttn (nr+1,0:1,0:1),
     x                 vbtn(nr+1,0:1,0:1)
c
      call set_eft_ctc( 0 , itz, ilb )
c
      do is=0,1
      do it=0,1
c
      do i=1,nr
      rr=hr*dble(i)
      call eft_ctc( itz, rr )
c
      vcn  (i+1,is,it)=vc  (is,it)
      vtn  (i+1,is,it)=vt  (is,it)
      vbn  (i+1,is,it)=vb  (is,it)
      vqn  (i+1,is,it)=vq  (is,it)
      vpn  (i+1,is,it)=vp  (is,it)
      vpdn (i+1,is,it)=vpd (is,it)
      vtpn (i+1,is,it)=vtp (is,it)
      vtpdn(i+1,is,it)=vtpd(is,it)
      vbbn (i+1,is,it)=vbb (is,it)
      vcvn (i+1,is,it)=vcv (is,it)
      vtvn (i+1,is,it)=vtv (is,it)
      vbvn (i+1,is,it)=vbv (is,it)
      vctn (i+1,is,it)=vct (is,it)
      vttn (i+1,is,it)=vtt (is,it)
      vbtn (i+1,is,it)=vbt (is,it)
      enddo
      vcn  (1,is,it)=3.d0*(vcn  (2,is,it)-vcn  (3,is,it))+vcn  (4,is,it)
      vtn  (1,is,it)=3.d0*(vtn  (2,is,it)-vtn  (3,is,it))+vtn  (4,is,it)
      vbn  (1,is,it)=3.d0*(vbn  (2,is,it)-vbn  (3,is,it))+vbn  (4,is,it)
      vqn  (1,is,it)=3.d0*(vqn  (2,is,it)-vqn  (3,is,it))+vqn  (4,is,it)
      vpn  (1,is,it)=3.d0*(vpn  (2,is,it)-vpn  (3,is,it))+vpn  (4,is,it)
      vpdn (1,is,it)=3.d0*(vpdn (2,is,it)-vpdn (3,is,it))+vpdn (4,is,it)
      vtpn (1,is,it)=3.d0*(vtpn (2,is,it)-vtpn (3,is,it))+vtpn (4,is,it)
      vtpdn(1,is,it)=3.d0*(vtpdn(2,is,it)-vtpdn(3,is,it))+vtpdn(4,is,it)
      vbbn (1,is,it)=3.d0*(vbbn (2,is,it)-vbbn (3,is,it))+vbbn (4,is,it)
      vcvn (1,is,it)=3.d0*(vcvn (2,is,it)-vcvn (3,is,it))+vcvn (4,is,it)
      vtvn (1,is,it)=3.d0*(vtvn (2,is,it)-vtvn (3,is,it))+vtvn (4,is,it)
      vbvn (1,is,it)=3.d0*(vbvn (2,is,it)-vbvn (3,is,it))+vbvn (4,is,it)
      vctn (1,is,it)=3.d0*(vctn (2,is,it)-vctn (3,is,it))+vctn (4,is,it)
      vttn (1,is,it)=3.d0*(vttn (2,is,it)-vttn (3,is,it))+vttn (4,is,it)
      vbtn (1,is,it)=3.d0*(vbtn (2,is,it)-vbtn (3,is,it))+vbtn (4,is,it)
c
      enddo
      enddo
 1000 continue
      return
      end


      SUBROUTINE EFT_PLESS_PW(LEMP, ILB, L, S, J, T1Z, T2Z, R, VPW)
      IMPLICIT NONE 
      INTEGER, INTENT(IN) :: LEMP, ILB
      INTEGER, INTENT(IN) :: L, S, J
      INTEGER, INTENT(IN) :: T1Z, T2Z
      DOUBLE PRECISION, INTENT(IN)  :: R
      DOUBLE PRECISION, INTENT(OUT) :: VPW(2,2)
c
      INTEGER :: TZ
      INTEGER :: INN, JBN, IPSB
      LOGICAL :: HEFORM, SING, TRIP, COUP, ENDEP
      CHARACTER*4 :: LABEL
      DOUBLE PRECISION :: V, XX
      LOGICAL, SAVE :: FIRST_CALL = .TRUE.
c      
      COMMON /CPOTR/   V(6), XX
      COMMON /CSTATE/ JBN,HEFORM,SING,TRIP,COUP,ENDEP,LABEL
      COMMON /CNN/ INN

      TZ = T1Z + T2Z
      IF (FIRST_CALL) THEN 
            CALL setpot(ILB, TZ)
            FIRST_CALL = .FALSE.
      ENDIF

      IPSB = 1
      INN = -TZ/2 + 2
      JBN = J
      XX = R

      CALL rpotr_eft(LEMP, IPSB)

      VPW = 0.D0
      IF (J.EQ.0) THEN
            IF (L.EQ.0.AND.S.EQ.0) VPW(1,1) = V(1)
            IF (L.EQ.1.AND.S.EQ.1) VPW(1,1) = V(3)
            RETURN
      ELSE IF (MOD(J,2).EQ.1) THEN
            IF (L.EQ.J.AND.S.EQ.0) THEN
                  IF (TZ.EQ.0) VPW(1,1) = V(1)
                  RETURN
            ENDIF
            IF (L.EQ.J.AND.S.EQ.1) THEN
                  VPW(1,1) = V(2)
                  RETURN
            ENDIF
            IF (TZ.NE.1) RETURN
            VPW(1,1) = V(4)
            VPW(1,2) = V(5)
            VPW(2,1) = V(6)
            VPW(2,2) = V(3)
            RETURN
      ELSE IF (MOD(J,2).EQ.0) THEN
            IF (L.EQ.J.AND.S.EQ.0) THEN
                  VPW(1,1) = V(1)
                  RETURN
            ENDIF
            IF (L.EQ.J.AND.S.EQ.1) THEN
                  IF (TZ.EQ.0) VPW(1,1) = V(2)
                  RETURN
            ENDIF
            VPW(1,1) = V(4)
            VPW(1,2) = V(5)
            VPW(2,1) = V(6)
            VPW(2,2) = V(3)
            RETURN
      ELSE
            STOP "ERROR IN EFT_PLESS_PW"
      ENDIF
      RETURN
      
      END SUBROUTINE EFT_PLESS_PW


      subroutine rpotr_eft(lemp,ipsb)
      implicit real*8(a-h,o-z)
      character*4 label
      logical heform,sing,trip,coup,endep
c
      common /cpotr/   v(6),xx
      common /cstate/ jbn,heform,sing,trip,coup,endep,label
      common /cnn/ inn
c      common /cpotr/inn,jbn,xx,v(6)
c
      r=xx
      ur=1.d0/r
c
      j=jbn
c
      do i=1,6
         v(i)=0.d0
      enddo
c
      call potl_eft(r,ur,lemp,ipsb,j,
     >              vpp1,vpp2,vpp3,vpp4,vpp5,
     >              vpn1,vpn2,vpn3,vpn4,vpn5,
     >              vnn1,vnn2,vnn3,vnn4,vnn5)
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
      if(j.eq.0) then
         v(1)=vx1
         v(3)=vx4
      else
         if(mod(j,2).eq.1) then
            v(2)=vx2
            if(inn.eq.2) then
               v(1)=vx1
               v(3)=vx4
               v(4)=vx3
               v(5)=vx5
               v(6)=vx5
            end if
         else if(mod(j,2).eq.0) then
            v(1)=vx1
            v(3)=vx4
            v(4)=vx3
            v(5)=vx5
            v(6)=vx5
            if(inn.eq.2) then
               v(2)=vx2
            end if
         end if
      end if
c
      return
      end




c-----------------------------------------------------------------------
      subroutine potl_eft(r,ur,lemp,ipsb,j,
     >                    vpp1,vpp2,vpp3,vpp4,vpp5,
     >                    vpn1,vpn2,vpn3,vpn4,vpn5,
     >                    vnn1,vnn2,vnn3,vnn4,vnn5)
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
c     if ipsb=0 no operators 12-19
c     the eft_r potential is used, and ilb select the version
c     ! ONLY ILB=4 IS READY
c      lemp  = 1-14 to include the e.m. potential for the av18
c------------------------------------------------------------------------
c
      data icont/0/
c
c------------------------------------------------------------------------
c
c    coefficienti per potenziale coulombiano
c      parameter (eqmf=197.327053d0/137.03599d0)  !as in *empot*
      parameter (eqmf=197.327053d0/137.035989d0) !new *empot*
      include 'params'
c--------------------------------------------------------------
c
      dimension vem(14),lts0(0:1),lts1(0:1),ltin(0:1)
c
      data zero/0.d0/,one/1.d0/,u2/2.d0/,u3/3.d0/,u4/4.d0/,u6/6.d0/
c
      data lts0/1,0/     !lts0(mod(L))=value of T: case S=0
      data lts1/0,1/     !lts1(mod(L))=value of T: case S=1
      data ltin/1,0/     !ltin(T=0)=1 and ltin(T=1)=0
c
      common/eftpot_rv/vc (nr+1,0:1,0:1),vt  (nr+1,0:1,0:1),
     x                 vb (nr+1,0:1,0:1),vq  (nr+1,0:1,0:1),
     x                 vp (nr+1,0:1,0:1),vpd (nr+1,0:1,0:1),
     x                 vtp(nr+1,0:1,0:1),vtpd(nr+1,0:1,0:1),
     x                 vbb(nr+1,0:1,0:1),
     x                 vcv(nr+1,0:1,0:1),vtv (nr+1,0:1,0:1),
     x                 vbv(nr+1,0:1,0:1),
     x                 vct(nr+1,0:1,0:1),vtt (nr+1,0:1,0:1),
     x                 vbt(nr+1,0:1,0:1)
c
      dimension fin(nr+1),
     x          vcr (0:1,0:1),vtr  (0:1,0:1),
     x          vbr (0:1,0:1),vqr  (0:1,0:1),
     x          vpr (0:1,0:1),vpdr (0:1,0:1),
     x          vtpr(0:1,0:1),vtpdr(0:1,0:1),vbbr(0:1,0:1),
     x          vcvr(0:1,0:1),vtvr (0:1,0:1),vbvr(0:1,0:1),
     x          vctr(0:1,0:1),vttr (0:1,0:1),vbtr(0:1,0:1)
c
      jp=mod(j,2)
      uj=one/dfloat(2*j+1)
      xjj=j*(j+1)
c
c     interpolation
c
      do is=0,1
      do it=0,1
         if(r.gt.hr*nr)then
         vcr  (is,it)=zero
         vtr  (is,it)=zero
         vbr  (is,it)=zero
         vqr  (is,it)=zero
         vpr  (is,it)=zero
         vpdr (is,it)=zero
         vtpr (is,it)=zero
         vtpdr(is,it)=zero
         vbbr (is,it)=zero
         vcvr (is,it)=zero
         vtvr (is,it)=zero
         vbvr (is,it)=zero
         vctr (is,it)=zero
         vttr (is,it)=zero
         vbtr (is,it)=zero
c
         else
c
         do i=1,nr+1
            fin(i)=vc(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vcr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vt(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vtr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vb(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vbr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vq(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vqr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vp(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vpr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vpd(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vpdr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vtp(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vtpr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vtpd(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vtpdr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vbb(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vbbr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vcv(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vcvr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vtv(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vtvr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vbv(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vbvr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vct(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vctr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vtt(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vttr(is,it)=fout
c
         do i=1,nr+1
            fin(i)=vbt(i,is,it)
         enddo
         call lagm(fin,1,nr+1,0.d0,hr,r,1,1,fout,1)
         vbtr(is,it)=fout
c
      endif
c
      enddo
      enddo
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
      if(ipsb.eq.0)then
         pc0=zero
         pcp=zero
         pcm=zero
      endif
      vcoul=zero
      if(lemp.eq.0) vcoul=eqmf*ur
      if(lemp.gt.0) call em_pot(lemp,r,vem)
c
c=======case 1 si=sj=0 li=lj=j
c
        lti=lts0(jp)
        p0=vcr(0,0)+xjj*vqr(0,0)
        p1=vcr(0,1)+xjj*vqr(0,1)
c
        if(lti.eq.1) then
          if(ipsb.ne.0)pcd=vcvr(0,1)+vctr(0,1)
c
          vempp=vcoul+vem(1)+vem(2)+vem(3)+vem(4)-u3*vem(6)
          vempn=vem(5)-u3*vem(8)
          vemnn=-u3*vem(7)
          vpp1=p1+pcd+vempp
          vpn1=p1+pcd+vempn
          vnn1=p1+pcd+vemnn
        else
          vempn=vem(5)-u3*vem(8)
          vpn1=p0+vempn
          vnn1=0.d0
          vpp1=0.d0
        end if
c
c=======case 2 si=sj=1 li=lj=j
c
        lti=lts1(jp)
        if(j.gt.0) then
          p0=vcr(1,0)+u2*vtr(1,0)-vbr(1,0)+xjj*vqr(1,0)+vbbr(1,0)
          p1=vcr(1,1)+u2*vtr(1,1)-vbr(1,1)+xjj*vqr(1,1)+vbbr(1,1)
c
          if(lti.eq.1) then
             if(ipsb.ne.0)
     x         pcd=     vcvr(1,1)+vctr(1,1)
     x            +u2* (vtvr(1,1)+vttr(1,1))
     x            +xjj*(vbvr(1,1)+vbtr(1,1))
c
            vempp=vcoul+vem(1)+vem(2)+vem(3)+vem(4)+vem(6)
     >           +u2*vem(9)-vem(12)
            vempn=vem(5)+vem(8)+u2*vem(11)-vem(14)
            vemnn=vem(7)+u2*vem(10)-vem(13)
            vpp2=p1+pcd+vempp
            vpn2=p1+pcd+vempn
            vnn2=p1+pcd+vemnn
          else
            vempn=vem(5)+vem(8)+u2*vem(11)-vem(14)
            vpp2=0.d0
            vpn2=p0+vempn
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
          p0=vcr(1,0)
     x      +xte*vtr(1,0)+xls*vbr(1,0)+xll*vqr(1,0)+xls2*vbbr(1,0)
          p1=vcr(1,1)
     x      +xte*vtr(1,1)+xls*vbr(1,1)+xll*vqr(1,1)+xls2*vbbr(1,1)
c
          if(lti.eq.1) then
             if(ipsb.ne.0)
     x        pcd=      vcvr(1,1)+vctr(1,1)
     x            +xte*(vtvr(1,1)+vttr(1,1))
     x            +xls*(vbvr(1,1)+vbtr(1,1))
c
            vempp=vcoul+vem(1)+vem(2)+vem(3)+vem(4)+vem(6)
     >           +xte*vem(9)+xls*vem(12)
            vempn=vem(5)+vem(8)+xte*vem(11)+xls*vem(14)
            vemnn=vem(7)+xte*vem(10)+xls*vem(13)
            vpp3=p1+pcd+vempp
            vpn3=p1+pcd+vempn
            vnn3=p1+pcd+vemnn
          else
            vempn=vem(5)+vem(8)+xte*vem(11)+xls*vem(14)
            vpn3=p0+vempn
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
        p0=vcr(1,0)+xte*vtr(1,0)+xls *vbr (1,0)
     x             +xll*vqr(1,0)+xls2*vbbr(1,0)
        p1=vcr(1,1)+xte*vtr(1,1)+xls *vbr (1,1)
     x             +xll*vqr(1,1)+xls2*vbbr(1,1)
c
        if(lti.eq.1) then
           if(ipsb.ne.0)
     x      pcd=     vcvr(1,1)+vctr(1,1)
     x         +xte*(vtvr(1,1)+vttr(1,1))
     x         +xls*(vbvr(1,1)+vbtr(1,1))
c
          vempp=vcoul+vem(1)+vem(2)+vem(3)+vem(4)+vem(6)
     >         +xte*vem(9)+xls*vem(12)
          vempn=vem(5)+vem(8)+xte*vem(11)+xls*vem(14)
          vemnn=vem(7)+xte*vem(10)+xls*vem(13)
          vpp4=p1+pcd+vempp
          vpn4=p1+pcd+vempn
          vnn4=p1+pcd+vemnn
        else
          vempn=vem(5)+vem(8)+xte*vem(11)+xls*vem(14)
          vpn4=p0+vempn
          vnn4=0.d0
          vpp4=0.d0
        end if
c
c=======case 5 si=sj=1 li=j-1 lj=j+1
c
        if(j.gt.0) then
          xte=u6*dsqrt(xjj)*uj
          p0=xte*vtr(1,0)
          p1=xte*vtr(1,1)
c
          if(lti.eq.1) then
           if(ipsb.ne.0)
     x      pcd=xte*(vtvr(1,1)+vttr(1,1))
c
            vempp=xte*vem(9)
            vempn=xte*vem(11)
            vemnn=xte*vem(10)
            vpp5=p1+pcd+vempp
            vpn5=p1+pcd+vempn
            vnn5=p1+pcd+vemnn
          else
            vempn=xte*vem(11)
            vpn5=p0+vempn
            vnn5=0.d0
            vpp5=0.d0
          end if
        end if
c
c=====================================================================
c
      return
      end











c_______________________________________________________________________
      subroutine em_pot( lemp , rr , vvem )
c_______________________________________________________________________
c
c lemp: 0=Coulomb w/ff (pp) only in MeV units
c       1=full electromagnetic potential in MeV units
c
c order of operators in vvem(l)
c l:     1=Coulomb w/ff (pp)   2=DF    (pp)          3=OO      (pp)
c        4=VP    (pp)                                5=Coulomb (np)
c        6=s1.s2 (pp)          7=s1.s2 (nn)          8=s1.s2   (np)
c        9=S12   (pp)         10=S12   (nn)         11=S12     (np)
c       12=L.S   (pp)         13=L.S   (nn)         14=L.S     (np)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      real*8 me,mp,mn,mr,mup,mun
      save
c
      parameter ( alpha=1.d0/137.03599d0 , pi=dacos(-1.d0 ) )
c
      dimension vvem(14)
c
      data hbarc / 197.32697d0 /
      data me / 0.510999d0 / , mp / 938.2720d0 / , mup / 2.79285d0 / ,
     x                         mn / 939.5653d0 / , mun /-1.91304d0 / ,
     x     bb / 4.27d0 / , beta / 0.0189d0 /
c
      br=bb*rr
      fcoul=1.d0-( 1.d0+11.d0*br   /16.d0
     x                 + 3.d0*br**2/16.d0
     x                 +      br**3/48.d0 )*dexp(-br )
      vvem(1)=alpha*hbarc*fcoul/rr
c-----------------------ATTENTION--------------------------------------
c     vvem(1)=vvem(1)-alpha*hbarc/rr ! vvem(1) defined as the correction to point Coulomb
c-----------------------ATTENTION--------------------------------------
c
      if( lemp .eq. 1 ) then
c
      mr=mp*mn/( mp+mn )
      ha=alpha*hbarc**3
c
      ftensor=1.d0-( 1.d0+br+br**2/  2.d0
     x                      +br**3/  6.d0
     x                      +br**4/ 24.d0
     x                      +br**5/144.d0 )*dexp(-br )
      fspinor=1.d0-( 1.d0+br+     br**2/ 2.d0
     x                      +7.d0*br**3/48.d0
     x                      +     br**4/48.d0 )*dexp(-br )
      fdelta=bb**3*( 1.d0+br+br**2/3.d0 )*dexp(-br )/16.d0
      fnp   =bb**2*( 15.d0*br   +15.d0*br**2
     x              + 6.d0*br**3+      br**4 )*dexp(-br )/384.d0
      fivp=ffvp( 2.d0*me*rr/hbarc )
c-------------------------------------------
c  central pp cpmponents follow
c
      vvem(2)=-ha*fdelta/( 4.d0*mp**2 )
      vvem(3)=-vvem(1)**2/mp
      vvem(4)=2.d0*alpha*vvem(1)*fivp/( 3.d0*pi )
c-------------------------------------------
c  central np component follows
c
      vvem(5)=hbarc*alpha*beta*fnp/rr
c-------------------------------------------
c  s1.s2 components for pp, nn, and np follow
c
      vvem(6)=-ha*mup*mup*fdelta/( 6.d0*mp*mp )
      vvem(7)=-ha*mun*mun*fdelta/( 6.d0*mn*mn )
      vvem(8)=-ha*mup*mun*fdelta/( 6.d0*mn*mp )
c-------------------------------------------
c  tensor components for pp, nn, and np follow
c
      vvem( 9)=-ha*mup*mup*ftensor/( 4.d0*mp*mp*rr**3 )
      vvem(10)=-ha*mun*mun*ftensor/( 4.d0*mn*mn*rr**3 )
      vvem(11)=-ha*mup*mun*ftensor/( 4.d0*mp*mn*rr**3 )
c-------------------------------------------
c  spin-orit components for pp, nn, and np follow
c
      vvem(12)=-ha*( 4.d0*mup-1.d0 )*fspinor/( 2.d0*mp*mp*rr**3 )
      vvem(13)=0.d0
      vvem(14)=-ha*       mun       *fspinor/( 2.d0*mn*mr*rr**3 )
c
      elseif( lemp .eq. 2 ) then         ! optionally the 1/r**2 part is omitted
c
      vvem(3)=0
      vvem(4)=0
c
      else
c
      do 10 i=2,14
      vvem(i)=0.d0
 10   continue
c
      endif
c
      return
      end
c
c
c
c_______________________________________________________________________
      real*8 function ffvp( zz )
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      parameter ( nx=201 , zinf=0.01d0 )
c
      common/gp/gg(nx),wgg(nx)
c
      dimension xx(nx),wx(nx)
c
      data gamma / 0.577215664901533d0 / ,
     x        pi / 3.141592653589793d0 /
c
      ffvp=0.d0
c
      if( zz .le. zinf ) then
      ffvp=-gamma-5.d0/6.d0-dlog( 0.50d0*zz )+3.d0*pi*zz/8.d0
c
      else
c
      alf=1.0d0/zz
      bet=0.5d0*pi*alf
c
      do 10 i=1,nx
      xxi=0.5d0*pi*gg(i)
c
      xx(i)=1.d0+alf*dtan( xxi )
      wx(i)=     bet*wgg(i)/dcos( xxi )**2
 10   continue
c
      do 20 i=1,nx
      x1=xx(i)
      x2=xx(i)*xx(i)
c
      fx=dexp(-zz*( x1-1.d0 ) )*dsqrt( x2-1.d0 )*( 1.d0+0.5d0/x2 )/x2
c
      ffvp=ffvp+wx(i)*fx
 20   continue
c
      ffvp=dexp(-zz )*ffvp
      endif
c
      return
      end











c------------------------------------------------------------------------
      subroutine set_eft_ctc( iprt , itz , ictc )
c------------------------------------------------------------------------
c
c EFT contact interactions at LO, NLO and N3LO with four possible
c combinations of cutoffs in the T=0 and T=1 channels
c
c IPRT = 0 or 1 not to or to print
c
c ITZ = 1 for pp; itz= 0 for np; itz=-1 for nn
c
c ICTC selects the interaction with the following values for the flag;
c for example, the 1.7/1.5 N3LO interaction corresponds to ICTC = 9
c
c-----------------------------------------------------------
c R0/R1    1.7/1.5   1.9/2.0   2.1/2.5   2.3/3.0  optimized
c-----------------------------------------------------------
c    LO       1        2         3         4       5
c
c   NLO       6        7         8         9      10
c
c  N3LO      11       12        13        14      15
c-----------------------------------------------------------
c
c   The subroutine SET_EFT_CTC is called once to set up the correct
c   set of low-energy constants. The subroutine EFT_CTC generates the
c   (strong-interaction) potential at the given inter-nucleon separation
c   rr (in fm).  The components of the potential in MeV are stored in
c   common block /EFTCTC_ST/ with the convention
c
c                   s=0 or 1 and t=0 or 1
c
c   vc  (s,t) ---> central term
c
c   vt  (s,t) ---> tensor term (is zero if s=0 )
c
c   vb  (s,t) ---> spin-orbit term (is zero if s=0 )
c
c   vq  (s,t) ---> L^2-dependent term
c
c   vp  (s,t) ---> p^2-dependent central term (is set to zero here)
c
c   vpd (s,t) ---> derivative p^2-dependent central term (is set to zero here)
c
c   vtp (s,t) ---> p^2-dependent term term (is set to zero here)
c
c   vtpd(s,t) ---> derivative p^2-dependent tensor term (is set to zero here)
c
c   vbb (s,t) ---> quadratic-spin-orbit term (is zero if s=0 )
c
c   vcv (s,t) ---> central csb term (is zero if t=0)
c
c   vtv (s,t) ---> tensor csb term (is zero if t=0)
c
c   vbv (s,t) ---> spin-orbit csb term (is zero if t=0)
c
c   vct (s,t) ---> central cd  term (is zero if t=0)
c
c   vtt (s,t) ---> tensor cd term (is set if t=0)
c
c   vbt (s,t) ---> tensor cd term (is set if t=0)
c
c--------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      character*5 clecs
      save
c
      parameter ( mpx=36 , nx=201 )
c
      common/gp/gg(nx),wgg(nx)
      common/lecs/rcf0,rcf1,c0c0,c0c1,c2c00,c2c10,c2c01,c2c11,
     x                                c4c00,c4c10,c4c01,c4c11,
     x                c2t0,c2t1,c4t0,c4t1,c2b,c4b0,c4b1,c4bb,
     x                c4q0,c4q1,c4p0,c4p1,c4tp0,c4tp1,
     x                c0cv,c2c0v,c2c1v,c2tv,c2bv,
     x                c0ct,c2c0t,c2c1t,c2tt,c2bt
      common/lecs_st/cc0cx  (0:1,0:1),cc2cx  (0:1,0:1),cc4cx (0:1,0:1),
     x               ct2tx  (0:1,0:1),ct4tx  (0:1,0:1),cb2bx (0:1,0:1),
     x               cb4bx  (0:1,0:1),cq4qx  (0:1,0:1),cp4px (0:1,0:1),
     x               cbb4bbx(0:1,0:1),ctp4tpx(0:1,0:1),cc0cvx(0:1,0:1),
     x               cc2cvx (0:1,0:1),ct2tvx (0:1,0:1),cb2bvx(0:1,0:1),
     x               cc0ctx (0:1,0:1),cc2ctx (0:1,0:1),ct2ttx(0:1,0:1),
     x               cb2btx (0:1,0:1)
c
      dimension clecs(mpx),dd(mpx)
      equivalence (dd( 1),rcf0 ),(dd( 2),rcf1 ),
     x            (dd( 3),c0c0 ),(dd( 4),c0c1 ),
     x            (dd( 5),c2c00),(dd( 6),c2c10),
     x            (dd( 7),c2c01),(dd( 8),c2c11),
     x            (dd( 9),c4c00),(dd(10),c4c10),
     x            (dd(11),c4c01),(dd(12),c4c11),
     x            (dd(13),c2t0 ),(dd(14),c2t1 ),
     x            (dd(15),c4t0 ),(dd(16),c4t1 ),
     x            (dd(17),c2b  ),(dd(18),c4b0 ),
     x            (dd(19),c4b1 ),(dd(20),c4bb ),
     x            (dd(21),c4q0 ),(dd(22),c4q1 ),
     x            (dd(23),c4p0 ),(dd(24),c4p1 ),
     x            (dd(25),c4tp0),(dd(26),c4tp1),
     x            (dd(27),c0cv ),(dd(28),c2c0v),
     x            (dd(29),c2c1v),(dd(30),c2tv ),
     x            (dd(31),c2bv ),(dd(32),c0ct ),
     x            (dd(33),c2c0t),(dd(34),c2c1t),
     x            (dd(35),c2tt ),(dd(36),c2bt )
c
      data clecs /'rcf0 ','rcf1 ',
     x            'c0c0 ','c0c1 ',
     x            'c2c00','c2c10','c2c01','c2c11',
     x            'c4c00','c4c10','c4c01','c4c11',
     x            'c2t0 ','c2t1 ',
     x            'c4t0 ','c4t1 ',
     x            'c2b  ','c4b0 ','c4b1 ','c4bb ',
     x            'c4q0 ','c4q1 ',
     x            'c4p0 ','c4p1 ','c4tp0','c4tp1',
     x            'c0cv ','c2c0v','c2c1v','c2tv ','c2bv ',
     x            'c0ct ','c2c0t','c2c1t','c2tt ','c2bt '/
      data hbarc / 197.32697d0 /
c
      call gauss( nx , 0.d0 , 1.d0 , gg , wgg )
c
      do i=1,mpx
      dd(i)=0.d0
      enddo
c
      if( ictc .gt. 22 ) then
      write(6,5000) ictc
      return
      endif
c
      if( ictc .eq. 1 ) then
c
c lecs for LO R0/R1 = 1.7/1.5 follow
c
      dd(1)=1.7d0
      dd(2)=1.5d0
c
      dd(3)=-.438524414d+01 ! c0c0
      dd(4)=-.800783936d+01 ! c0c1
c
      elseif( ictc .eq. 2 ) then
c
c lecs for LO R0/R1 = 1.9/2.0 follow
c
      dd(1)=1.9d0
      dd(2)=2.0d0
c
      dd(3)=-.572220536d+01 ! c0c0
      dd(4)=-.934392090d+01 ! c0c1
c
      elseif( ictc .eq. 3 ) then
c
c lecs for LO R0/R1 = 2.1/2.5 follow
c
      dd(1)=2.1d0
      dd(2)=2.5d0
c
      dd(3)=-.700250932d+01 ! c0c0
      dd(4)=-.107734100d+02 ! c0c1
c
      elseif( ictc .eq. 4 ) then
c
c lecs for LO R0/R1 = 2.3/3.0 follow
c
      dd(1)=2.3d0
      dd(2)=3.0d0
c
      dd(3)=-.822926713d+01 ! c0c0
      dd(4)=-.122993164d+02 ! c0c1
c
      elseif( ictc .eq. 5 ) then
c
c lecs for optimized LO with R0/R1 obtained by fitting effective range expansion
c
      dd(1)=0.154592984d+01
      dd(2)=0.183039397d+01
c
      dd(3)=-.527518671d+01 ! c0c0
      dd(4)=-.704040080d+01 ! c0c1
c
      elseif( ictc .eq. 6 ) then
c
c lecs for NLO R0/R1 = 1.7/1.5 follow
c
      dd(1)=1.7d0
      dd(2)=1.5d0
c
      dd( 3)=-.511051122d+01 ! c0c0
      dd( 4)=-.422732988d+01 ! c0c1
      dd( 5)=-.828234407d+01 ! c2c00
      dd( 6)=-.237961671d+01 ! c2c10
      dd( 7)=0.106696423d+01 ! c2c01
      dd( 8)=-.646100475d+00 ! c2c11
      dd(13)=0.258995384d+01 ! c2t0
      dd(14)=-.797433238d+00 ! c2t1
      dd(17)=-.155550814d+01 ! c2b

      dd(32)=0.190747072d-01 ! c0ct
c
      elseif( ictc .eq. 7 ) then
c
c lecs for NLO R0/R1 = 1.9/2.0 follow
c
      dd(1)=1.9d0
      dd(2)=2.0d0
c
      dd( 3)=-.515205193d+01 ! c0c0
      dd( 4)=-.486213195d+01 ! c0c1
      dd( 5)=-.856027246d+01 ! c2c00
      dd( 6)=-.414687057d+01 ! c2c10
      dd( 7)=-.997027624d+00 ! c2c01
      dd( 8)=-.583022650d+00 ! c2c11
      dd(13)=0.299664427d+01 ! c2t0
      dd(14)=-.100664567d+01 ! c2t1
      dd(17)=-.138788868d+01 ! c2b
c
      dd(32)=0.242061782d-01 ! c0ct
c
      elseif( ictc .eq. 8 ) then
c
c lecs for NLO R0/R1 = 2.1/2.5 follow
c
      dd(1)=2.1d0
      dd(2)=2.5d0
c
      dd( 3)=-.524089036d+01 ! c0c0
      dd( 4)=-.147490885d+01 ! c0c1
      dd( 5)=-.168474897d+02 ! c2c00
      dd( 6)=-.993610707d+01 ! c2c10
      dd( 7)=-.514045901d+01 ! c2c01
      dd( 8)=-.428220919d+00 ! c2c11
      dd(13)=0.608990650d+01 ! c2t0
      dd(14)=-.111195915d+01 ! c2t1
      dd(17)=-.150745124d+01 ! c2b
c
      dd(32)=0.343911021d-01 ! c0ct
c
      elseif( ictc .eq. 9 ) then
c
c lecs for NLO R0/R1 = 2.3/3.0 follow
c
      dd(1)=2.3d0
      dd(2)=3.0d0
c
      dd( 3)=-.534645335d+01 ! c0c0
      dd( 4)=-.442765927d+01 ! c0c1
      dd( 5)=-.556714095d+00 ! c2c00
      dd( 6)=-.120137656d+02 ! c2c10
      dd( 7)=-.120651283d+02 ! c2c01
      dd( 8)=-.504387462d+00 ! c2c11
      dd(13)=0.380297934d+01 ! c2t0
      dd(14)=-.122136856d+01 ! c2t1
      dd(17)=-.153475063d+01 ! c2b
c
      dd(32)=0.488093390d-01 ! c0ct

      elseif( ictc .eq. 10 ) then
c
c lecs for optimized NLO with R0/R1 from optimized LO
c
      dd(1)=0.154592984d+01
      dd(2)=0.183039397d+01
c
      dd( 3)=-.514608170d+01 ! c0c0
      dd( 4)=-.564430900d+01 ! c0c1
      dd( 5)=-.838477662d+01 ! c2c00
      dd( 6)=-.389762585d+00 ! c2c10
      dd( 7)=-.744446852d-01 ! c2c01
      dd( 8)=-.582484599d+00 ! c2c11
      dd(13)=0.156269012d+01 ! c2t0
      dd(14)=-.924209886d+00 ! c2t1
      dd(17)=-.136793827d+01 ! c2b
c
      dd(32)=0.219960910d-01 ! c0ct
c
      elseif( ictc .eq. 11 ) then
c
c lecs for N3LO R0/R1 = 1.7/1.5 follow
c
      dd(1)=1.7d0
      dd(2)=1.5d0
c
      dd( 3)=-.511424764d+01 ! c0c0
      dd( 4)=-.425743601d+01 ! c0c1
      dd( 5)=-.830665210d+01 ! c2c00
      dd( 6)=-.239136315d+01 ! c2c10
      dd( 7)=0.109812283d+01 ! c2c01
      dd( 8)=-.710230431d+00 ! c2c11
      dd( 9)=-.181868453d-01 ! c4c00
      dd(10)=0.221901443d-01 ! c4c10
      dd(11)=-.153123742d-01 ! c4c01
      dd(12)=0.617231418d-02 ! c4c11
      dd(13)=0.260926378d+01 ! c2t0
      dd(14)=-.788483914d+00 ! c2t1
      dd(15)=-.125947898d-01 ! c4t0
      dd(16)=0.398120601d-01 ! c4t1
      dd(17)=-.148260980d+01 ! c2b
      dd(18)=-.190612708d-01 ! c4b0
      dd(19)=0.258095198d-01 ! c4b1
      dd(20)=0.111357163d+00 ! c4bb
      dd(21)=0.199158467d+00 ! c4q0
      dd(22)=-.229235566d-01 ! c4q1

c
      dd(32)=0.616726547d-02 ! c0ct
      dd(33)=0.195747500d-01 ! c2c0t
      dd(34)=-.126831643d-01 ! c2c1t
      dd(35)=-.261236310d-01 ! c2tt
      dd(36)=0.156812161d-02 ! c2bt
c
      elseif( ictc .eq. 12 ) then
c
c lecs for N3LO R0/R1 = 1.9/2.0 follow
c
      dd(1)=1.9d0
      dd(2)=2.0d0
c
      dd( 3)=-.508230349d+01 ! c0c0
      dd( 4)=-.473602278d+01 ! c0c1
      dd( 5)=-.860411467d+01 ! c2c00
      dd( 6)=-.423065440d+01 ! c2c10
      dd( 7)=-.906565836d+00 ! c2c01
      dd( 8)=-.625041622d+00 ! c2c11
      dd( 9)=-.990927063d-01 ! c4c00
      dd(10)=-.215651571d-01 ! c4c10
      dd(11)=-.201602182d+00 ! c4c01
      dd(12)=0.326112893d+00 ! c4c11
      dd(13)=0.298626502d+01 ! c2t0
      dd(14)=-.997041036d+00 ! c2t1
      dd(15)=0.162098906d-01 ! c4t0
      dd(16)=0.567306211d-01 ! c4t1
      dd(17)=-.121337175d+01 ! c2b
      dd(18)=0.769565092d-01 ! c4b0
      dd(19)=-.939920749d-01 ! c4b1
      dd(20)=0.112830624d+00 ! c4bb
      dd(21)=0.392181450d+00 ! c4q0
      dd(22)=0.716506159d-02 ! c4q1
c
      dd(32)=-.221853840d-01 ! c0ct
      dd(33)=0.913635316d-01 ! c2c0t
      dd(34)=-.190925764d-01 ! c2c1t
      dd(35)=-.156366250d-01 ! c2tt
      dd(36)=0.583713002d-01 ! c2bt
c
      elseif( ictc .eq. 13 ) then
c
c lecs for N3LO R0/R1 = 2.1/2.5 follow
c
      dd(1)=2.1d0
      dd(2)=2.5d0
c
      dd( 3)=-.503452047d+01 ! c0c0
      dd( 4)=-.822959678d+00 ! c0c1
      dd( 5)=-.169235588d+02 ! c2c00
      dd( 6)=-.964326070d+01 ! c2c10
      dd( 7)=-.437509314d+01 ! c2c01
      dd( 8)=-.352440064d+01 ! c2c11
      dd( 9)=-.967846481d+00 ! c4c00
      dd(10)=-.425444862d+00 ! c4c10
      dd(11)=-.172710802d+01 ! c4c01
      dd(12)=-.987551702d+00 ! c4c11
      dd(13)=0.508667430d+01 ! c2t0
      dd(14)=-.591092927d+00 ! c2t1
      dd(15)=0.949731880d+00 ! c4t0
      dd(16)=-.165746282d+01 ! c4t1
      dd(17)=-.149303317d+01 ! c2b
      dd(18)=-.154418998d+01 ! c4b0
      dd(19)=-.224728260d+00 ! c4b1
      dd(20)=-.145021865d+00 ! c4bb
      dd(21)=0.980175820d+00 ! c4q0
      dd(22)=0.993835808d+00 ! c4q1
c
      dd(32)=-.406049402d-01 ! c0ct
      dd(33)=0.198860745d+00 ! c2c0t
      dd(34)=0.117430362d+01 ! c2c1t
      dd(35)=-.372938280d-01 ! c2tt
      dd(36)=0.137474019d+00 ! c2bt
c
      elseif( ictc .eq. 14 ) then
c
c lecs for N3LO R0/R1 = 2.3/3.0 follow
c
      dd(1)=2.3d0
      dd(2)=3.0d0
c
      dd( 3)=-.503178655d+01 ! c0c0
      dd( 4)=-.376510267d+01 ! c0c1
      dd( 5)=-.138856138d+02 ! c2c00
      dd( 6)=-.115155939d+02 ! c2c10
      dd( 7)=-.927114901d+01 ! c2c01
      dd( 8)=-.534391071d+01 ! c2c11
      dd( 9)=-.768444931d+01 ! c4c00
      dd(10)=-.816662965d+00 ! c4c10
      dd(11)=-.700061242d+01 ! c4c01
      dd(12)=-.178012661d+02 ! c4c11
      dd(13)=0.243913791d+01 ! c2t0
      dd(14)=0.417191713d+00 ! c2t1
      dd(15)=0.140894086d+01 ! c4t0
      dd(16)=-.425038125d+01 ! c4t1
      dd(17)=-.175175554d+01 ! c2b
      dd(18)=-.832846394d+00 ! c4b0
      dd(19)=-.625130230d+01 ! c4b1
      dd(20)=0.564521782d+01 ! c4bb
      dd(21)=-.118230699d+02 ! c4q0
      dd(22)=0.346894956d+01 ! c4q1
c
      dd(32)=-.127554450d+00 ! c0ct
      dd(33)=0.642256720d+00 ! c2c0t
      dd(34)=0.293333551d+01 ! c2c1t
      dd(35)=-.299306179d-01 ! c2tt
      dd(36)=-.218783861d+00 ! c2bt
c
      elseif( ictc .eq. 15 ) then
c
c lecs for optimized N3LO with R0/R1 from optimized LO
c
      dd(1)=0.154592984d+01
      dd(2)=0.183039397d+01
c
      dd( 3)=-.512268575d+01 ! c0c0
      dd( 4)=-.569749961d+01 ! c0c1
      dd( 5)=-.839788290d+01 ! c2c00
      dd( 6)=-.373568342d+00 ! c2c10
      dd( 7)=-.680778104d-01 ! c2c01
      dd( 8)=-.614743769d+00 ! c2c11
      dd( 9)=0.137148891d-01 ! c4c00
      dd(10)=-.286999005d-01 ! c4c10
      dd(11)=-.378665625d-01 ! c4c01
      dd(12)=0.875705263d-01 ! c4c11
      dd(13)=0.161976638d+01 ! c2t0
      dd(14)=-.972924744d+00 ! c2t1
      dd(15)=-.491377039d-01 ! c4t0
      dd(16)=0.268116821d-01 ! c4t1
      dd(17)=-.134307736d+01 ! c2b
      dd(18)=0.218026841d-01 ! c4b0
      dd(19)=-.369118294d-01 ! c4b1
      dd(20)=0.226506657d-01 ! c4bb
      dd(21)=-.624395862d-02 ! c4q0
      dd(22)=0.312122677d-01 ! c4q1
c
      dd(32)=0.713292586d-02 ! c0ct
      dd(33)=0.264716400d-01 ! c2c0t
      dd(34)=-.239979852d-01 ! c2c1t
      dd(35)=0.374105167d-03 ! c2tt
      dd(36)=0.298742271d-01 ! c2bt
c
      elseif( ictc .eq. 16 ) then
c
c lecs for optimized NLO by fitting phase shifts with R0/R1 from optimized LO
c
      dd(1)=0.154592984d+01
      dd(2)=0.183039397d+01
c
      dd( 3)=-.493643068d+01 ! c0c0
      dd( 4)=-.656370228d+01 ! c0c1
      dd( 5)=-.820312500d+01 ! c2c00
      dd( 6)=-.608265239d-01 ! c2c10
      dd( 7)=-.455950114d+00 ! c2c01
      dd( 8)=0.508124677d-01 ! c2c11
      dd(13)=0.996063497d+00 ! c2t0
      dd(14)=-.809327695d+00 ! c2t1
      dd(17)=-.155823261d+01 ! c2b
c
      dd(32)=0.206041226d-01 ! c0ct
c
      elseif( ictc .eq. 17 ) then
c
c optimized NLO obtained by fitting R0/R1 + LECs
c up to 50 MeV (chi2=3.02 on 1652 data)
c
      dd(1)=0.162418606d+01
      dd(2)=0.170728665d+01
c
      dd( 3)=-.503069500d+01 ! c0c0
      dd( 4)=-.594950343d+01 ! c0c1
      dd( 5)=-.841654055d+01 ! c2c00
      dd( 6)=-.639748984d+00 ! c2c10
      dd( 7)=0.257616266d+00 ! c2c01
      dd( 8)=-.365824374d-01 ! c2c11
      dd(13)=0.174319470d+01 ! c2t0
      dd(14)=-.812753053d+00 ! c2t1
      dd(17)=-.157526292d+01 ! c2b
c
      dd(32)=0.191454558d-01 ! c0ct
c
      elseif( ictc .eq. 18 ) then
c
c optimized NLO obtained by fitting R0/R1 + LECs
c up to 25 MeV (chi2=1.29 on 1102 data)
c
      dd(1)=0.155209236d+01
      dd(2)=0.183468010d+01
c
      dd( 3)=-.514472335d+01 ! c0c0
      dd( 4)=-.567791015d+01 ! c0c1
      dd( 5)=-.840904110d+01 ! c2c00
      dd( 6)=-.425830770d+00 ! c2c10
      dd( 7)=-.976609508d-01 ! c2c01
      dd( 8)=-.505843262d+00 ! c2c11
      dd(13)=0.155337935d+01 ! c2t0
      dd(14)=-.926615204d+00 ! c2t1
      dd(17)=-.144029429d+01 ! c2b
c
      dd(32)=0.220327143d-01 ! c0ct
c
      elseif( ictc .eq. 19 ) then
c
c optimized NLO obtained by fitting R0/R1 + LECs
c up to 15 MeV (chi2=1.92 on 1102 data)

      dd(1)=0.112551429d+01
      dd(2)=0.124357315d+01
c
      dd( 3)=-.228869385d+02 ! c0c0
      dd( 4)=-.558807983d+01 ! c0c1
      dd( 5)=-.375752246d+01 ! c2c00
      dd( 6)=0.138687186d+01 ! c2c10
      dd( 7)=-.214137298d+01 ! c2c01
      dd( 8)=-.939496128d+00 ! c2c11
      dd(13)=0.733571886d+00 ! c2t0
      dd(14)=-.673165532d+00 ! c2t1
      dd(17)=-.134221054d+01 ! c2b
c
      dd(32)=0.127009276d-01 ! c0ct
c
      elseif ( ictc .eq. 20 ) then
c
c lecs arxiv 2408.02480 pionless LO full
c
      dd(1)=2.0d0
      dd(2)=2.0d0
c      
      dd(3)=-.571710d1 ! c0c0
      dd(4)=-.100468d2 ! c0c1
c
      elseif ( ictc .eq. 21 ) then 
c
c lecs arxiv 2408.02480 pionless NLO full
c
      dd(1)=2.0d0
      dd(2)=2.0d0
c
      dd( 3)=-.52794d1 ! c0c0
      dd( 4)=-.18177d1 ! c0c1
      dd( 5)= .14973d1 ! c2c00
      dd( 6)=-.90723d1 ! c2c10
      dd( 7)=-.76604d0 ! c2c01
      dd( 8)=-.37466d1 ! c2c11
      dd(13)= .33174d1 ! c2t0
      dd(14)= .29972d1 ! c2t1
      dd(17)=-.23540d0 ! c2b
c
      dd(27)= .01500d0 ! c0cv
c
      dd(32)= .02000d0 ! c0ct
c      
      elseif ( ictc .eq. 22 ) then
c
c lecs arxiv 2408.02480 pionless NLO no P-waves
c
      dd(1)=2.0d0
      dd(2)=2.0d0
c
      dd( 3)=-.53502d1 ! c0c0
      dd( 4)=-.35152d1 ! c0c1
      dd( 6)=-.78473d1 ! c2c10
      dd( 7)=-.70560d0 ! c2c01
      dd(13)= .18046d1 ! c2t0
c
      dd(27)= .01570d0 ! c0cv
c
      dd(32)= .01990d0 ! c0ct
c      
      endif
      

c
      if( iprt .eq. 1 ) then
      write(6,5010) ictc
      write(6,5020) ( dd(i),clecs(i), i=1,mpx )
      write(6,5030)
      endif
c
c   S=0 and T=0 channel LECs follow
c----------ATTENTION------------------
c LO only acts in ST=01 or 10
c
c     cc0cx(0,0)=c0c0
c-------------------------------------
      cc0cx(0,0)=0.d0
c-------------------------------------
      cc2cx(0,0)=c2c00
      cc4cx(0,0)=c4c00
c
      ct2tx(0,0)=0.d0
      ct4tx(0,0)=0.d0
c
      cb2bx(0,0)=0.d0
      cb4bx(0,0)=0.d0
c
      cq4qx(0,0)=c4q0
      cp4px(0,0)=c4p0
c
      cbb4bbx(0,0)=0.d0
      ctp4tpx(0,0)=0.d0
c
      cc0cvx(0,0)=0.d0
      cc2cvx(0,0)=0.d0
      ct2tvx(0,0)=0.d0
      cb2bvx(0,0)=0.d0
c
      cc0ctx(0,0)=0.d0
      cc2ctx(0,0)=0.d0
      ct2ttx(0,0)=0.d0
      cb2btx(0,0)=0.d0
c
c   S=0 and T=1 channel LECs follow
c
      cc0cx(0,1)=c0c0
      cc2cx(0,1)=c2c01
      cc4cx(0,1)=c4c01
c
      ct2tx(0,1)=0.d0
      ct4tx(0,1)=0.d0
c
      cb2bx(0,1)=0.d0
      cb4bx(0,1)=0.d0
c
      cq4qx(0,1)=c4q0
      cp4px(0,1)=c4p0
c
      cbb4bbx(0,1)=0.d0
      ctp4tpx(0,1)=0.d0
c
      cc0cvx(0,1)=dble( 2*itz )*c0cv
      cc2cvx(0,1)=dble( 2*itz )*c2c0v
      ct2tvx(0,1)=0.d0
      cb2bvx(0,1)=0.d0
c
      cc0ctx(0,1)=dble(-4+6*iabs( itz ) )*c0ct
      cc2ctx(0,1)=dble(-4+6*iabs( itz ) )*c2c0t
      ct2ttx(0,1)=0.d0
      cb2btx(0,1)=0.d0
c
c   S=1 and T=0 channel LECs follow
c
      cc0cx(1,0)=c0c1
      cc2cx(1,0)=c2c10
      cc4cx(1,0)=c4c10
c
      ct2tx(1,0)=c2t0
      ct4tx(1,0)=c4t0
c
      cb2bx(1,0)=c2b
      cb4bx(1,0)=c4b0
c
      cq4qx(1,0)=c4q1
      cp4px(1,0)=c4p1
c
      cbb4bbx(1,0)=c4bb
      ctp4tpx(1,0)=c4tp0
c
      cc0cvx(1,0)=0.d0
      cc2cvx(1,0)=0.d0
      ct2tvx(1,0)=0.d0
      cb2bvx(1,0)=0.d0
c
      cc0ctx(1,0)=0.d0
      cc2ctx(1,0)=0.d0
      ct2ttx(1,0)=0.d0
      cb2btx(1,0)=0.d0
c
c   S=1 and T=1 channel LECs follow
c
c----------ATTENTION------------------
c LO only acts in ST=01 or 10
c
c     cc0cx(1,1)=c0c1
c-------------------------------------
      cc0cx(1,1)=0.d0
c-------------------------------------
      cc2cx(1,1)=c2c11
      cc4cx(1,1)=c4c11
c
      ct2tx(1,1)=c2t1
      ct4tx(1,1)=c4t1
c
      cb2bx(1,1)=c2b
      cb4bx(1,1)=c4b1
c
      cq4qx(1,1)=c4q1
      cp4px(1,1)=c4p1
c
      cbb4bbx(1,1)=c4bb
      ctp4tpx(1,1)=c4tp1
c
      cc0cvx(1,1)=dble( 2*itz )*c0cv
      cc2cvx(1,1)=dble( 2*itz )*c2c1v
      ct2tvx(1,1)=dble( 2*itz )*c2tv
      cb2bvx(1,1)=dble( 2*itz )*c2bv
c
      cc0ctx(1,1)=dble(-4+6*iabs( itz ) )*c0ct
      cc2ctx(1,1)=dble(-4+6*iabs( itz ) )*c2c1t
      ct2ttx(1,1)=dble(-4+6*iabs( itz ) )*c2tt
      cb2btx(1,1)=dble(-4+6*iabs( itz ) )*c2bt
c
      do it=0,1
      do is=0,1
      cc0cx(is,it)=hbarc*cc0cx(is,it)
      cc2cx(is,it)=hbarc*cc2cx(is,it)
      cc4cx(is,it)=hbarc*cc4cx(is,it)
c
      ct2tx(is,it)=hbarc*ct2tx(is,it)
      ct4tx(is,it)=hbarc*ct4tx(is,it)
c
      cb2bx(is,it)=hbarc*cb2bx(is,it)
      cb4bx(is,it)=hbarc*cb4bx(is,it)
c
      cq4qx(is,it)=hbarc*cq4qx(is,it)
      cp4px(is,it)=hbarc*cp4px(is,it)
c
      cbb4bbx(is,it)=hbarc*cbb4bbx(is,it)
      ctp4tpx(is,it)=hbarc*ctp4tpx(is,it)
c
      cc0cvx(is,it)=hbarc*cc0cvx(is,it)
      cc2cvx(is,it)=hbarc*cc2cvx(is,it)
      ct2tvx(is,it)=hbarc*ct2tvx(is,it)
      cb2bvx(is,it)=hbarc*cb2bvx(is,it)
c
      cc0ctx(is,it)=hbarc*cc0ctx(is,it)
      cc2ctx(is,it)=hbarc*cc2ctx(is,it)
      ct2ttx(is,it)=hbarc*ct2ttx(is,it)
      cb2btx(is,it)=hbarc*cb2btx(is,it)
      enddo
      enddo
c
      return
 5000 format(//,'  wrong value for ictc=',i5)
 5010 format(//,10x,'ictc=',i3,//)
 5020 format(5x,e16.9,2x,a5)
 5030 format(//)
      end
c
c
c
c_______________________________________________________________________
      subroutine eft_ctc( itz , rr  )
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      parameter( ncf=2 , rmax=20.d0 )
c
      common/lecs_st/cc0cx  (0:1,0:1),cc2cx  (0:1,0:1),cc4cx (0:1,0:1),
     x               ct2tx  (0:1,0:1),ct4tx  (0:1,0:1),cb2bx (0:1,0:1),
     x               cb4bx  (0:1,0:1),cq4qx  (0:1,0:1),cp4px (0:1,0:1),
     x               cbb4bbx(0:1,0:1),ctp4tpx(0:1,0:1),cc0cvx(0:1,0:1),
     x               cc2cvx (0:1,0:1),ct2tvx (0:1,0:1),cb2bvx(0:1,0:1),
     x               cc0ctx (0:1,0:1),cc2ctx (0:1,0:1),ct2ttx(0:1,0:1),
     x               cb2btx (0:1,0:1)
      common/eftctc_st/vc (0:1,0:1),vt  (0:1,0:1),
     x                 vb (0:1,0:1),vq  (0:1,0:1),
     x                 vp (0:1,0:1),vpd (0:1,0:1),
     x                 vtp(0:1,0:1),vtpd(0:1,0:1),vbb(0:1,0:1),
     x                 vcv(0:1,0:1),vtv (0:1,0:1),vbv(0:1,0:1),
     x                 vct(0:1,0:1),vtt (0:1,0:1),vbt(0:1,0:1)
c
      if( rr .gt. rmax ) then
c
      do it=0,1
      do is=0,1
      vc(is,it)=0.d0
      vt(is,it)=0.d0
      vb(is,it)=0.d0
      vq(is,it)=0.d0
c
      vp  (is,it)=0.d0
      vpd (is,it)=0.d0
      vtp (is,it)=0.d0
      vtpd(is,it)=0.d0
      vbb (is,it)=0.d0
c
      vcv(is,it)=0.d0
      vtv(is,it)=0.d0
      vbv(is,it)=0.d0
c
      vct(is,it)=0.d0
      vtt(is,it)=0.d0
      vbt(is,it)=0.d0
      enddo
      enddo
c
      elseif( rr .le. rmax ) then
c
      do it=0,1
      do is=0,1
c
      call ctc( it , ncf , rr , fc0  , fc2  , ft2   , fb2 ,
     x                          fc4  , ft4  , fb4 ,
     x                          fbb4 , fq4  ,
     x                          fp4  , fp4d  , ftp4 , ftp4d )
c
      vc(is,it)=cc0cx(is,it)*fc0+cc2cx(is,it)*fc2+cc4cx(is,it)*fc4
      vt(is,it)=                 ct2tx(is,it)*ft2+ct4tx(is,it)*ft4
      vb(is,it)=                 cb2bx(is,it)*fb2+cb4bx(is,it)*fb4
      vq(is,it)=                                  cq4qx(is,it)*fq4
c
      vp  (is,it)=cp4px  (is,it)*fp4
      vpd (is,it)=cp4px  (is,it)*fp4d
      vtp (is,it)=ctp4tpx(is,it)*ftp4
      vtpd(is,it)=ctp4tpx(is,it)*ftp4d
      vbb (is,it)=cbb4bbx(is,it)*fbb4
c
c  isospin-symmetry violating terms porportional to ( tau_{1z}+tau_{2z} )
c
      vcv(is,it)=cc0cvx(is,it)*fc0+cc2cvx(is,it)*fc2
      vtv(is,it)=                  ct2tvx(is,it)*ft2
      vbv(is,it)=                  cb2bvx(is,it)*fb2
c
c  isospin-symmetry violating terms porportional to  T_{12}
c
      vct(is,it)=cc0ctx(is,it)*fc0+cc2ctx(is,it)*fc2
      vtt(is,it)=                  ct2ttx(is,it)*ft2
      vbt(is,it)=                  cb2btx(is,it)*fb2
      enddo
      enddo
c
      endif
c
      return
      end
c
c
c
c_______________________________________________________________________
      subroutine ctc( it , n , rr ,
     x                fc0  , fc2   , ft2 , fb2  ,
     x                fc4  , ft4   , fb4 , fbb4 , fq4 ,
     x                fp4  , fp4d  ,
     x                ftp4 , ftp4d )
c_______________________________________________________________________
c
c Contact radial functions at LO, NLO, N3LO:
c
c n = exponent for r-space cutoff
c
c rr = inter-nucleon separation in fm units
c
c fc0   = v^c(r) at LO in fm^{-3} units
c fc2   = v^c(r) at NLO in fm^{-5} units
c ft2   = v^t(r) at NLO in fm^{-5} units
c fb2   = v^b(r) at NLO in fm^{-5} units
c fc4   = v^c(r) at N3LO in fm^{-7} units
c ft4   = v^t(r) at N3LO in fm^{-7} units
c fb4   = v^b(r) at N3LO in fm^{-7} units
c fbb4  = v^bb(r) at N3LO in fm^{-7} units
c fq4   = v^{q}(r) an N3LO in fm^{-7} units
c fp4   = v^{p}(r) an N3LO in fm^{-5} units
c fp4d  = first derivative of v^{p}(r) an N3LO in fm^{-6} units
c ftp4  = v^{tp}(r) an N3LO in fm^{-5} units
c ftp4d = first derivative of v^{tp}(r) an N3LO in fm^{-6} units
c
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      common/lecs/rcf0,rcf1,c0c0,c0c1,c2c00,c2c10,c2c01,c2c11,
     x                                c4c00,c4c10,c4c01,c4c11,
     x                c2t0,c2t1,c4t0,c4t1,c2b,c4b0,c4b1,c4bb,
     x                c4q0,c4q1,c4p0,c4p1,c4tp0,c4tp1,
     x                c0cv,c2c0v,c2c1v,c2tv,c2bv,
     x                c0ct,c2c0t,c2c1t,c2tt,c2bt
c
      data pi / 3.141592653589793d0 /
c
      if( it .eq. 0 ) rcfx=rcf0
      if( it .eq. 1 ) rcfx=rcf1
c
      rcfn=rcfx**n
c
      n1=n-1
      n2=n-2
      n3=n-3
      n4=n-4
c
      rrn0=rr**n
      rrn1=rr**n1
      rrn2=rr**n2
      rrn3=rr**n3
      rrn4=rr**n4
c
      rr2=rr*rr
      rr3=rr*rr2
      rr4=rr*rr3
c
      rcff=dble( n )/rcfn
c
      ann=0.25d0*n/( pi*dgamma( 3.d0/dble( n ) )*rcfx**3 )
      wcf=ann*dexp(-rrn0/rcfn )
c
      wcfd1=-rcff*wcf*rrn1
c
      wcfd2=-rcff*(            wcfd1*rrn1
     x             +dble( n1 )*wcf  *rrn2 )
c
      wcfd3=-rcff*(                    wcfd2*rrn1
     x             +2.d0*dble( n1    )*wcfd1*rrn2
     x             +     dble( n1*n2 )*wcf  *rrn3 )
c
      wcfd4=-rcff*(                       wcfd3*rrn1
     x             +3.d0*dble( n1       )*wcfd2*rrn2
     x             +3.d0*dble( n1*n2    )*wcfd1*rrn3
     x             +     dble( n1*n2*n3 )*wcf  *rrn4 )
c
      fc0=wcf
c
      fc2=-( wcfd2+2.d0*wcfd1/rr )
      ft2=-( wcfd2-     wcfd1/rr )
      fb2=-             wcfd1/rr
c
      fc4=  ( wcfd4+4.d0*wcfd3/rr                               )
      ft4=  ( wcfd4+     wcfd3/rr-6.d0*wcfd2/rr2+6.d0*wcfd1/rr3 )
      fb4=  (            wcfd3/rr+2.d0*wcfd2/rr2-2.d0*wcfd1/rr3 )
      fbb4=-(                          wcfd2/rr2-     wcfd1/rr3 )
      fq4 =-(                          wcfd2/rr2-     wcfd1/rr3 )
c
      fp4 =fc2
      ftp4=ft2
c
      fp4d =-wcfd3-2.d0*wcfd2/rr+2.d0*wcfd1/rr2
      ftp4d=-wcfd3+     wcfd2/rr-     wcfd1/rr2
c
      return
      end
