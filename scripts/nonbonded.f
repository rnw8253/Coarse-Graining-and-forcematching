CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCC Main Program CCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM Force_Match
        include 'omp_lib.h'
        integer nmax
        parameter(nmax=100000)
C paramtop values
        integer natom,ntyp,nbonh,nbona,ntheth,ntheta
        integer nphih,nphia,nnb,nres
        double precision charge(nmax),mass(nmax)
        integer ityp(nmax),nbparm(nmax),nrpnt(nmax)
        double precision kbnd(nmax),rbnd(nmax),kang(nmax),tang(nmax)
        double precision kdih(nmax),pdih(nmax),ndih(nmax)
        double precision scee(nmax),scnb(nmax)
        double precision alj(nmax),blj(nmax)
        integer bhlist(3,nmax),balist(3,nmax)
        integer ahlist(4,nmax),aalist(4,nmax)
        integer dhlist(5,nmax),dalist(5,nmax)
        integer excl(nmax),nexc(nmax)
        double precision kumb1,kumb2,rumb
        integer umblist(nmax,2)
C solvent values
        integer ma(55),inext,inextp
        integer npts,ncpu,resin(nmax),isel
        double precision w(nmax),gr2(2,nmax)
        double precision alp(nmax),g0(nmax),x0(nmax)
        double precision As(nmax),Bs(nmax),vtot(nmax)
C
        double precision x(3,nmax),v(3,nmax),f(3,nmax)
        double precision lbox(3),hlbox(3),dt,hdt,T,ran3,try,pnu
C
        character*64 rst,frc,xyz,vel,ene,com,fname,temp
        character*64 frst,prmtop,outfile
C     
        integer i,j,igo,jgo,kgo,tid,ip,ivel,step,iex,jex
        double precision com_r(3),r2,mtot,com_dist
        double precision hist_min,hist_max,bin_size
        double precision force_hist(7,7,580), F_AB(7,7)
        integer hist_counter(7,7,580),nsolute,cutoff
        integer dist_bin,counter(7,7)
        integer nspnt(8),num_bins,m,k
        double precision bin_dist
C     
        common /param/ natom,ntyp,nbonh,nbona,ntheth,ntheta,nphih,nphia,
     &       nnb,nres,charge,mass,ityp,nbparm,nrpnt,kbnd,rbnd,kang,tang,
     &       kdih,pdih,ndih,scee,scnb,alj,blj,bhlist,balist,ahlist,
     &       aalist,dhlist,dalist,excl,kumb1,kumb2,umblist,rumb,nspnt,
     &       nexc
C
        common /solvent/ gr2,alp,x0,g0,As,Bs,vtot,w,npts,ncpu
C
        common /random/ ma,inext,inextp


C create selections
        natom = 1080

C loop over files
        write(prmtop,888)
888     format('chx60_vac.no_coul.prmtop')  
        call setparam(prmtop)
        write(fname,999)
        open(20,file=fname,status='old')
        read (20,*)
        write(outfile,997)
        open(25,FILE=outfile, STATUS='unknown')
        do step=1,5000
           call readcrd(x, lbox, hlbox, natom)           
           call getfrc(x,lbox,hlbox,f)
           write(25,995) "Step:", step
           do i=1,natom
              write(25,996) f(1,i), f(2,i), f(3,i)
           enddo                     
           write(6,998)  step
        enddo
        close(20)

999     format('chx60_vac.no_coul.run0.crd7')
998     format(1(i5.5))
995     format(a,i5.5)
997     format('chx60_vac.no_coul.run0.vdw.dat')
996     format(f12.6,1x,f12.6,1x,f12.6)

        stop
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readcrd(x,lbox,hlbox,natom)
        integer nmax
        parameter(nmax=100000)

        integer i,j,l,k
        double precision x(3,nmax), lbox(3), hlbox(3)
        integer natom
        double precision Xo(10)


      i = 1
      j = 1
      do while(i.le.natom)
         read(20,991) Xo(1),Xo(2),Xo(3),Xo(4),Xo(5),Xo(6),Xo(7),Xo(8),
     &        Xo(9),Xo(10)
         do k=1,10
            x(j,i) = Xo(k)
            j = j + 1
            if(j.eq.4) then 
               j = 1
               i = i + 1
            endif
         enddo
      enddo
      read(20,992) lbox(1),lbox(2),lbox(3)
      do l=1,3
         hlbox(l) = lbox(l)/2.0
      enddo

991   format(10(f8.3))
992   format(3(f8.3))

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getfrc(x,lbox,hlbox,f)
        integer nmax
        parameter(nmax=100000)
C paramtop values
        integer natom,nbonh,nbona,ntheth,ntheta
        integer nphih,nphia,nnb,nres,ntyp
        double precision charge(nmax),mass(nmax)
        integer ityp(nmax),nbparm(nmax),nrpnt(nmax)
        double precision kbnd(nmax),rbnd(nmax),kang(nmax),tang(nmax)
        double precision kdih(nmax),pdih(nmax),ndih(nmax)
        double precision scee(nmax),scnb(nmax)
        double precision alj(nmax),blj(nmax)
        integer bhlist(3,nmax),balist(3,nmax)
        integer ahlist(4,nmax),aalist(4,nmax)
        integer dhlist(5,nmax),dalist(5,nmax)
        integer excl(nmax),nexc(nmax)
        double precision kumb1,kumb2,rumb
        integer umblist(nmax,2)
        integer resin(nmax)
        double precision a,b,c11,c12,c13,c22,c23,c33,f1,f14e,f14v
        double precision f4,fang,fbnd,fdih,t1,t2,t3,t4,t5,t6,th
        double precision d1(3),d2(3),d3(3)
        
        double precision x(3,nmax)
        double precision lbox(3),hlbox(3)
C
        integer it,jt,nlj,i,j,k,iex,jex
        integer ksel,lsel,nsel
        double precision f(3,nmax)
        integer i1,i2,i3,i4,i5,nspnt(7)
        double precision flj,fc
        double precision r2,r6,r(3)
C
        double precision dx
        double precision gnow,fs,try,prob
C
        common /param/ natom,ntyp,nbonh,nbona,ntheth,ntheta,nphih,nphia,
     &       nnb,nres,charge,mass,ityp,nbparm,nrpnt,kbnd,rbnd,kang,tang,
     &       kdih,pdih,ndih,scee,scnb,alj,blj,bhlist,balist,ahlist,
     &       aalist,dhlist,dalist,excl,kumb1,kumb2,umblist,rumb,nspnt,
     &       nexc

        do i=1,natom
           do j=1,3
              f(j,i)=0.d0
           enddo
        enddo

C  Non-bonded
        iex=1
        if (excl(1).eq.0) then
          jex=1
        else
          jex=0
        endif
        do i=1,natom
          it=ityp(i)
          do j=i+1,natom
            if (excl(iex).ne.j) then
              jt=ityp(j)
              nlj=ntyp*(it-1)+jt
              nlj=nbparm(nlj)
              r2=0.d0
              do k=1,3
                dx=x(k,i)-x(k,j)
                if (dx.gt.hlbox(k)) then
                  dx=dx-lbox(k)
                endif
                if (dx.lt.-hlbox(k)) then
                  dx=dx+lbox(k)
                endif
                r(k)=dx
                r2=r2+dx*dx
              enddo
              r6=r2**(-3)
              flj=r6*(12.d0*alj(nlj)*r6-6.d0*blj(nlj))/r2
              fc=charge(i)*charge(j)/dsqrt(r2)/r2
              fc = 0
              if (r2.gt.144.d0) then
                 flj=0
                 fc=0
              endif
              do k=1,3
                f(k,i)=f(k,i)+(flj+fc)*r(k)
                f(k,j)=f(k,j)-(flj+fc)*r(k)
              enddo
            else
              iex=iex+1
            endif
          enddo
          if (excl(iex).eq.0) then
            if (jex.eq.0) then
              jex=1
            else
              iex=iex+1
              if (excl(iex).eq.0) then
                jex=1
              else
                jex=0
              endif
           endif
          endif
       enddo

C     Bonds      
C        do i=1,nbona
C           i1=balist(1,i)
C           i2=balist(2,i)
C           it=balist(3,i)
C           write(6,*) i1,i2,it
C           r2=0.d0
C           do k=1,3
C              r(k)=x(k,i1)-x(k,i2)
C              r2=r2+r(k)*r(k)
C           enddo
C           fbnd=2.d0*kbnd(it)*(rbnd(it)/dsqrt(r2)-1.d0)
CC       if(dabs(fbnd).gt.20.d0) write(6,*) r2,i1,i2,fbnd,kbnd(it),rbnd(it)
C           do k=1,3
C              f(k,i1)=f(k,i1)+fbnd*r(k)
C              f(k,i2)=f(k,i2)-fbnd*r(k)
C           enddo
C        enddo
CC        write(6,*) natom,nbona
C        do i=1,nbonh
C           i1=bhlist(1,i)
C           i2=bhlist(2,i)
C           it=bhlist(3,i)
CC           write(6,*) i1,i2,it
C           r2=0.d0
C           do k=1,3
C              r(k)=x(k,i1)-x(k,i2)
C              r2=r2+r(k)**2
C           enddo
C           fbnd=2.d0*kbnd(it)*(rbnd(it)/dsqrt(r2)-1.d0)
C           do k=1,3
C              f(k,i1)=f(k,i1)+fbnd*r(k)
C              f(k,i2)=f(k,i2)-fbnd*r(k)
C           enddo
C        enddo
C
C
CC  Angles
C        do i=1,ntheta
C          i1=aalist(1,i)
C          i2=aalist(2,i)
C          i3=aalist(3,i)
C          it=aalist(4,i)
C          c11=0.d0
C          c22=0.d0
C          c12=0.d0
C          do k=1,3
C             d1(k)=x(k,i1)-x(k,i2)
C             d2(k)=x(k,i2)-x(k,i3)
C             c11=c11+d1(k)*d1(k)
C             c22=c22+d2(k)*d2(k)
C             c12=c12+d1(k)*d2(k)
C          enddo
C          b=-c12/dsqrt(c11*c22)
C          if(b.ge.1.d0) then
C             th=1.d-16
C          elseif(b.le.-1.d0) then
C             th=3.1415926535d0
C          else
C             th=dacos(b)
C          endif
C          fang=2.d0*kang(it)*(th-tang(it))/dsqrt(c11*c22-c12*c12)
C          do k=1,3
C             f(k,i1)=f(k,i1)+fang*(c12/c11*d1(k)-d2(k))
C             f(k,i2)=f(k,i2)+fang*
C     &            ((1.d0+c12/c22)*d2(k)-(1.d0+c12/c11)*d1(k))
C             f(k,i3)=f(k,i3)+fang*(d1(k)-c12/c22*d2(k))
C             enddo
C        enddo
C        do i=1,ntheth
C          i1=ahlist(1,i)
C          i2=ahlist(2,i)
C          i3=ahlist(3,i)
C          it=ahlist(4,i)
C          c11=0.d0
C          c22=0.d0
C          c12=0.d0
C          do k=1,3
C             d1(k)=x(k,i1)-x(k,i2)
C             d2(k)=x(k,i2)-x(k,i3)
C             c11=c11+d1(k)*d1(k)
C             c22=c22+d2(k)*d2(k)
C             c12=c12+d1(k)*d2(k)
C          enddo
C          b=-c12/dsqrt(c11*c22)
C          if(b.ge.1.d0) then
C             th=1.d-16
C          elseif(b.le.-1.d0) then
C             th=3.1415926535d0
C          else
C             th=dacos(b)
C          endif
C          fang=2.d0*kang(it)*(th-tang(it))/dsqrt(c11*c22-c12*c12)
C          do k=1,3
C             f(k,i1)=f(k,i1)+fang*(c12/c11*d1(k)-d2(k))
C             f(k,i2)=f(k,i2)+fang*
C     &            ((1.d0+c12/c22)*d2(k)-(1.d0+c12/c11)*d1(k))
C             f(k,i3)=f(k,i3)+fang*(d1(k)-c12/c22*d2(k))
C          enddo
C        enddo
C
C  Dihedrals
C        do i=1,nphih
C           i1=dhlist(1,i)
C           i2=dhlist(2,i)
C           i3=dhlist(3,i)
C           i4=dhlist(4,i)
C           i5=dhlist(5,i)
C           if(i4.lt.0) i4=-i4
C           if(i3.lt.0) then
C              i3=-i3
C           else
C              r2=0.d0
C              do k=1,3
C                 r(k)=x(k,i1)-x(k,i4)
C                 r2=r2+r(k)*r(k)
C              enddo
C              r6=r2**(-3)
C              it=ityp(i1)
C              jt=ityp(i4)
C              nlj=ntyp*(it-1)+jt
C              nlj=nbparm(nlj)
C              f14e=charge(i1)*charge(i4)/r2/dsqrt(r2)/scee(i5)
C              f14v=r6*(12.d0*alj(nlj)*r6-6.d0*blj(nlj))/scnb(i5)/r2
C              do k=1,3
C                 f(k,i1)=f(k,i1)+(f14e+f14v)*r(k)
C                 f(k,i4)=f(k,i4)-(f14e+f14v)*r(k)
C              enddo
C           endif                 
CC
C           if(i2.gt.0) then
C              c11=0.d0
C              c12=0.d0
C              c13=0.d0
C              c22=0.d0
C              c23=0.d0
C              c33=0.d0
C              do k=1,3
C                 d1(k)=x(k,i1)-x(k,i2)
C                 d2(k)=x(k,i2)-x(k,i3)
C                 d3(k)=x(k,i3)-x(k,i4)
C                 c11=c11+d1(k)*d1(k)
C                 c12=c12+d1(k)*d2(k)
C                 c13=c13+d1(k)*d3(k)
C                 c22=c22+d2(k)*d2(k)
C                 c23=c23+d2(k)*d3(k)
C                 c33=c33+d3(k)*d3(k)
C              enddo
CC     
C              t1=c13*c22-c12*c23
C              t2=c11*c23-c12*c13
C              t3=c12*c12-c11*c22
C              t4=c22*c33-c23*c23
C              t5=c13*c23-c12*c33
C              t6=-t1
CC     
C              b=dsqrt(-t3*t4)
CC     
C              a=t6/b
C              if(a.le.-1.d0) then
C                 fdih=0.d0
C              elseif(a.ge.1.d0) then
C                 fdih=0.d0
C              else
C                 th=dacos(a)
C                 fdih=ndih(i5)*kdih(i5)*dsin(ndih(i5)*th-pdih(i5))
C     &                /dsin(th)*c22/b
C              endif
C           else
C              i2=-i2
C              if(a.le.-1.d0) then
C                 fdih=0.d0
C              elseif(a.ge.1.d0) then
C                 fdih=0.d0
C              else
C                 th=dacos(a)
C                 fdih=ndih(i5)*kdih(i5)*dsin(ndih(i5)*th-pdih(i5))
C     &                /dsin(th)*c22/b
C              endif
C           endif
C           do k=1,3
C              f1=fdih*(t1*d1(k)+t2*d2(k)+t3*d3(k))/t3
C              f4=-fdih*(t4*d1(k)+t5*d2(k)+t6*d3(k))/t4
C              f(k,i1)=f(k,i1)+f1
C              f(k,i2)=f(k,i2)-(1.d0+c12/c22)*f1+c23/c22*f4
C              f(k,i3)=f(k,i3)+c12/c22*f1-(1.d0+c23/c22)*f4
C              f(k,i4)=f(k,i4)+f4
C           enddo
C        enddo        
C        do i=1,nphia
C           i1=dalist(1,i)
C           i2=dalist(2,i)
C           i3=dalist(3,i)
C           i4=dalist(4,i)
C           i5=dalist(5,i)
C           if(i4.lt.0) i4=-i4
C           if(i3.lt.0) then
C              i3=-i3
C           else
C              r2=0.d0
C              do k=1,3
C                 r(k)=x(k,i1)-x(k,i4)
C                 r2=r2+r(k)*r(k)
C              enddo
C              r6=r2**(-3)
C              it=ityp(i1)
C              jt=ityp(i4)
C              nlj=ntyp*(it-1)+jt
C              nlj=nbparm(nlj)
C              f14e=charge(i1)*charge(i4)/r2/dsqrt(r2)/scee(i5)
C              f14v=r6*(12.d0*alj(nlj)*r6-6.d0*blj(nlj))/scnb(i5)/r2
C              do k=1,3
C                 f(k,i1)=f(k,i1)+(f14e+f14v)*r(k)
C                 f(k,i4)=f(k,i4)-(f14e+f14v)*r(k)
C              enddo
C           endif
C     
C           if(i2.gt.0) then
C              c11=0.d0
C              c12=0.d0
C              c13=0.d0
C              c22=0.d0
C              c23=0.d0
C              c33=0.d0
C              do k=1,3
C                 d1(k)=x(k,i1)-x(k,i2)
C                 d2(k)=x(k,i2)-x(k,i3)
C                 d3(k)=x(k,i3)-x(k,i4)
C                 c11=c11+d1(k)*d1(k)
C                 c12=c12+d1(k)*d2(k)
C                 c13=c13+d1(k)*d3(k)
C                 c22=c22+d2(k)*d2(k)
C                 c23=c23+d2(k)*d3(k)
C                 c33=c33+d3(k)*d3(k)
C              enddo
CC     
C              t1=c13*c22-c12*c23
C              t2=c11*c23-c12*c13
C              t3=c12*c12-c11*c22
C              t4=c22*c33-c23*c23
C              t5=c13*c23-c12*c33
C              t6=-t1
CC     
C              b=dsqrt(-t3*t4)
CC     
C              a=t6/b
C              if(a.le.-1.d0) then
C                 fdih=0.d0
C              elseif(a.ge.1.d0) then
C                 fdih=0.d0
C              else
C                 th=dacos(a)
C                 fdih=ndih(i5)*kdih(i5)*dsin(ndih(i5)*th-pdih(i5))
C     &                /dsin(th)*c22/b
C              endif
C           else
C              i2=-i2
C              if (a.le.-1.d0) then
C                 fdih=0.d0
C              elseif (a.ge.1.d0) then
C                 fdih=0.d0
C              else
C                 th=dacos(a)
C                 fdih=ndih(i5)*kdih(i5)*dsin(ndih(i5)*th-pdih(i5))
C     &                /dsin(th)*c22/b
C              endif
C           endif
C           do k=1,3
C              f1=fdih*(t1*d1(k)+t2*d2(k)+t3*d3(k))/t3
C              f4=-fdih*(t4*d1(k)+t5*d2(k)+t6*d3(k))/t4
C              f(k,i1)=f(k,i1)+f1
C              f(k,i2)=f(k,i2)-(1.d0+c12/c22)*f1+c23/c22*f4
C              f(k,i3)=f(k,i3)+c12/c22*f1-(1.d0+c23/c22)*f4
C              f(k,i4)=f(k,i4)+f4
C           enddo
C        enddo

        return
        end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine setparam(prmtop)
        integer nmax
        parameter(nmax=100000)
        integer natom,ntyp,nbonh,nbona,ntheth,ntheta,nphih,nphia,natyp
        integer nnb,nres,numbnd,numang,nptra
        double precision charge(nmax),mass(nmax)
        integer ityp(nmax),nbparm(nmax),nrpnt(nmax)
        double precision kbnd(nmax),rbnd(nmax),kang(nmax),tang(nmax)
        double precision kdih(nmax),pdih(nmax),ndih(nmax)
        double precision scee(nmax),scnb(nmax)
        double precision alj(nmax),blj(nmax)
        integer bhlist(3,nmax),balist(3,nmax)
        integer ahlist(4,nmax),aalist(4,nmax)
        integer dhlist(5,nmax),dalist(5,nmax)
        integer excl(nmax),nexc(nmax)
        double precision kumb1,kumb2,rumb
        integer umblist(nmax,2),nspnt(7)
        integer i,j,inext,i1,i2,i3,i4
        character*80 str
        character*64 prmtop
C
        common /param/ natom,ntyp,nbonh,nbona,ntheth,ntheta,nphih,nphia,
     &       nnb,nres,charge,mass,ityp,nbparm,nrpnt,kbnd,rbnd,kang,tang,
     &       kdih,pdih,ndih,scee,scnb,alj,blj,bhlist,balist,ahlist,
     &       aalist,dhlist,dalist,excl,kumb1,kumb2,umblist,rumb,nspnt,
     &       nexc
C
        open(20,FILE=prmtop,STATUS='old')
C  numbers
        do i=1,6
          read(20,*)
        enddo
        read(20,999) natom,ntyp,nbonh,nbona,ntheth,ntheta,nphih,nphia
        read(20,998) nnb,nres,numbnd,numang,nptra,natyp
C  charge
        inext=7+(natom-1)/20
        do i=1,inext
          read(20,*)
        enddo
        inext=(natom-1)/5+1
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) charge(j),charge(j+1),charge(j+2),
     &         charge(j+3),charge(j+4)
        enddo
C  mass
        inext=5+(natom-1)/10
        do i=1,inext
          read(20,*)
        enddo
        inext=(natom-1)/5+1
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) mass(j),mass(j+1),mass(j+2),
     &         mass(j+3),mass(j+4)
        enddo
C  ityp
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=(natom-1)/10+1
        do i=1,inext
          j=(i-1)*10+1
          read(20,996) ityp(j),ityp(j+1),ityp(j+2),ityp(j+3),ityp(j+4),
     &         ityp(j+5),ityp(j+6),ityp(j+7),ityp(j+8),ityp(j+9)
        enddo
C  n-excluded
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=(natom-1)/10+1
        do i=1,inext
          j=(i-1)*10+1
          read(20,996) nexc(j),nexc(j+1),nexc(j+2),nexc(j+3),nexc(j+4),
     &         nexc(j+5),nexc(j+6),nexc(j+7),nexc(j+8),nexc(j+9)
        enddo
C non-bonded parm
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=(ntyp*ntyp-1)/10+1
        do i=1,inext
          j=(i-1)*10+1
          read(20,996) nbparm(j),nbparm(j+1),nbparm(j+2),nbparm(j+3),
     &         nbparm(j+4),nbparm(j+5),nbparm(j+6),nbparm(j+7),
     &         nbparm(j+8),nbparm(j+9)
        enddo
C  residue pointers
        inext=5+(nres-1)/20
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nres-1)/10
        do i=1,inext
          j=(i-1)*10+1
          read(20,996) nrpnt(j),nrpnt(j+1),nrpnt(j+2),nrpnt(j+3),
     &         nrpnt(j+4),nrpnt(j+5),nrpnt(j+6),nrpnt(j+7),
     &         nrpnt(j+8),nrpnt(j+9)
        enddo
        nrpnt(nres+1)=natom+1
C  bond force
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(numbnd-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) kbnd(j),kbnd(j+1),kbnd(j+2),kbnd(j+3),kbnd(j+4)
        enddo
C  bond distance
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(numbnd-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) rbnd(j),rbnd(j+1),rbnd(j+2),rbnd(j+3),rbnd(j+4)
        enddo
C  angle force
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(numang-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) kang(j),kang(j+1),kang(j+2),kang(j+3),kang(j+4)
        enddo
C  angle values
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(numang-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) tang(j),tang(j+1),tang(j+2),tang(j+3),tang(j+4)
        enddo
C  dihedral force
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nptra-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) kdih(j),kdih(j+1),kdih(j+2),kdih(j+3),kdih(j+4)
        enddo
C  dihedral preriodicity
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nptra-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) ndih(j),ndih(j+1),ndih(j+2),ndih(j+3),ndih(j+4)
        enddo
C  dihedral phase
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nptra-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) pdih(j),pdih(j+1),pdih(j+2),pdih(j+3),pdih(j+4)
        enddo
C  SCEE scale factor
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nptra-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) scee(j),scee(j+1),scee(j+2),scee(j+3),scee(j+4)
        enddo
C  SCNB scale factor
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nptra-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) scnb(j),scnb(j+1),scnb(j+2),scnb(j+3),scnb(j+4)
        enddo
C  Lennard-Jones A
        inext=5+(natyp-1)/5
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(ntyp*(ntyp+1)/2-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) alj(j),alj(j+1),alj(j+2),alj(j+3),alj(j+4)
        enddo
C  Lennard-Jones B
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(ntyp*(ntyp+1)/2-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) blj(j),blj(j+1),blj(j+2),blj(j+3),blj(j+4)
        enddo
C  Bond list w/ H
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(3*nbonh-1)/10
        do i=1,inext
          j=(i-1)*10+3
          i1=mod((j+1),3)
          if (i1.eq.0) i1=3
          i2=i1+1
          if (i2.eq.4) i2=1
          i3=i2+1
          if (i3.eq.4) i3=1
          read(20,996) bhlist(i1,j/3),bhlist(i2,(j+1)/3),
     &         bhlist(i3,(j+2)/3),bhlist(i1,(j+3)/3),bhlist(i2,(j+4)/3),
     &         bhlist(i3,(j+5)/3),bhlist(i1,(j+6)/3),bhlist(i2,(j+7)/3),
     &         bhlist(i3,(j+8)/3),bhlist(i1,(j+9)/3)
        enddo
        do i=1,nbonh
          bhlist(1,i)=bhlist(1,i)/3+1
          bhlist(2,i)=bhlist(2,i)/3+1
        enddo
C  Bond list w/o H
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(3*nbona-1)/10
        do i=1,inext
          j=(i-1)*10+3
          i1=mod((j+1),3)
          if (i1.eq.0) i1=3
          i2=i1+1
          if (i2.eq.4) i2=1
          i3=i2+1
          if (i3.eq.4) i3=1
          read(20,996) balist(i1,j/3),balist(i2,(j+1)/3),
     &         balist(i3,(j+2)/3),balist(i1,(j+3)/3),balist(i2,(j+4)/3),
     &         balist(i3,(j+5)/3),balist(i1,(j+6)/3),balist(i2,(j+7)/3),
     &         balist(i3,(j+8)/3),balist(i1,(j+9)/3)
        enddo
        do i=1,nbona
          balist(1,i)=balist(1,i)/3+1
          balist(2,i)=balist(2,i)/3+1
        enddo
C  Angle list w/ H
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(4*ntheth-1)/10
        do i=1,inext
          j=(i-1)*10+4
          i1=mod((j+1),4)
          if (i1.eq.0) i1=4
          i2=i1+1
          if (i2.eq.5) i2=1
          i3=i2+1
          if (i3.eq.5) i3=1
          i4=i3+1
          if (i4.eq.5) i4=1
          read(20,996) ahlist(i1,j/4),ahlist(i2,(j+1)/4),
     &         ahlist(i3,(j+2)/4),ahlist(i4,(j+3)/4),ahlist(i1,(j+4)/4),
     &         ahlist(i2,(j+5)/4),ahlist(i3,(j+6)/4),ahlist(i4,(j+7)/4),
     &         ahlist(i1,(j+8)/4),ahlist(i2,(j+9)/4)
        enddo
        do i=1,ntheth
          ahlist(1,i)=ahlist(1,i)/3+1
          ahlist(2,i)=ahlist(2,i)/3+1
          ahlist(3,i)=ahlist(3,i)/3+1
        enddo
C  Angle list w/o H
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(4*ntheta-1)/10
        do i=1,inext
          j=(i-1)*10+4
          i1=mod((j+1),4)
          if (i1.eq.0) i1=4
          i2=i1+1
          if (i2.eq.5) i2=1
          i3=i2+1
          if (i3.eq.5) i3=1
          i4=i3+1
          if (i4.eq.5) i4=1
          read(20,996) aalist(i1,j/4),aalist(i2,(j+1)/4),
     &         aalist(i3,(j+2)/4),aalist(i4,(j+3)/4),aalist(i1,(j+4)/4),
     &         aalist(i2,(j+5)/4),aalist(i3,(j+6)/4),aalist(i4,(j+7)/4),
     &         aalist(i1,(j+8)/4),aalist(i2,(j+9)/4)
        enddo
        do i=1,ntheta
          aalist(1,i)=aalist(1,i)/3+1
          aalist(2,i)=aalist(2,i)/3+1
          aalist(3,i)=aalist(3,i)/3+1
        enddo
C  Dihedral list w/ H
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(5*nphih-1)/10
        do i=1,inext
          j=(i-1)*2+1
          read(20,996) dhlist(1,j),dhlist(2,j),
     &         dhlist(3,j),dhlist(4,j),dhlist(5,j),
     &         dhlist(1,j+1),dhlist(2,j+1),dhlist(3,j+1),
     &         dhlist(4,j+1),dhlist(5,j+1)
        enddo
        do i=1,nphih
          dhlist(1,i)=dhlist(1,i)/3+1
          dhlist(2,i)=dhlist(2,i)/3+1
          if(dhlist(3,i).lt.0) then
            dhlist(3,i)=dhlist(3,i)/3-1
          else
            dhlist(3,i)=dhlist(3,i)/3+1
          endif
          if(dhlist(4,i).lt.0) then
            dhlist(4,i)=dhlist(4,i)/3-1
          else
            dhlist(4,i)=dhlist(4,i)/3+1
          endif
        enddo
        do i=2,nphih
          j=i-1
          if(iabs(dhlist(1,i)).eq.iabs(dhlist(1,j))) then
            if(iabs(dhlist(2,i)).eq.iabs(dhlist(2,j))) then
              if(iabs(dhlist(3,i)).eq.iabs(dhlist(3,j))) then
                if(iabs(dhlist(4,i)).eq.iabs(dhlist(4,j))) then
                  dhlist(2,i)=-dhlist(2,i)
                endif
              endif
            endif
          endif
        enddo
C  Dihedral list w/o H
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(5*nphia-1)/10
        do i=1,inext
          j=(i-1)*2+1
          read(20,996) dalist(1,j),dalist(2,j),
     &         dalist(3,j),dalist(4,j),dalist(5,j),
     &         dalist(1,j+1),dalist(2,j+1),dalist(3,j+1),
     &         dalist(4,j+1),dalist(5,j+1)
        enddo
        do i=1,nphia
          dalist(1,i)=dalist(1,i)/3+1
          dalist(2,i)=dalist(2,i)/3+1
          if(dalist(3,i).lt.0) then
            dalist(3,i)=dalist(3,i)/3-1
          else
            dalist(3,i)=dalist(3,i)/3+1
          endif
          if(dalist(4,i).lt.0) then
            dalist(4,i)=dalist(4,i)/3-1
          else
            dalist(4,i)=dalist(4,i)/3+1
          endif
        enddo
        do i=2,nphia
          j=i-1
          if(iabs(dalist(1,i)).eq.iabs(dalist(1,j))) then
            if(iabs(dalist(2,i)).eq.iabs(dalist(2,j))) then
              if(iabs(dalist(3,i)).eq.iabs(dalist(3,j))) then
                if(iabs(dalist(4,i)).eq.iabs(dalist(4,j))) then
                  dalist(2,i)=-dalist(2,i)
                endif
              endif
            endif
          endif
        enddo
C  Excluded atom list
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nnb-1)/10
        do i=1,inext
          j=(i-1)*10+1
          read(20,996) excl(j),excl(j+1),excl(j+2),excl(j+3),excl(j+4),
     &         excl(j+5),excl(j+6),excl(j+7),excl(j+8),excl(j+9)
        enddo
C
999     format(8(i8))
998     format(2(i8),3(8X),4(i8))
997     format(5(e16.8))
996     format(10(i8))
        close(20)
C
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
