      program bdyoutcompare
      print *,'link file 1 to fort.91'
      print *,'link file 2 to fort.92'
      print *,'enter ix,jx,kx'
      read (*,*) ix,jx,kx
      print *,'enter inhyd,imoist,iice,igraupel'
      read (*,*) inhyd,imoist,iice,igraupel
      print *,'enter tolerance
      call rdbc(ix,jx,kx,91,92,inhyd,imoist,iice,igraupel)

      subroutine rdbc(ix,jx,kx,iunit,junit,inhyd,imoist,iice,igraupel)

      parameter(nd=5)

      dimension tresc(ix,jx),tresc2(ix,jx)
      real eb3d(ix,kx,nd),wb3d(ix,kx,nd),nb3d(jx,kx,nd),
     $     sb3d(jx,kx,nd),et3d(ix,kx,nd),wt3d(ix,kx,nd),
     $     nt3d(jx,kx,nd),st3d(jx,kx,nd)
      real peb(ix,nd),pwb(ix,nd),pnb(jx,nd),
     $     psb(jx,nd),pet(ix,nd),pwt(ix,nd),
     $     pnt(jx,nd),pst(jx,nd)
      real web(ix,kx+1,nd),wwb(ix,kx+1,nd),
     $     wnb(jx,kx+1,nd),wsb(jx,kx+1,nd),
     $     wet(ix,kx+1,nd),wwt(ix,kx+1,nd),
     $     wnt(jx,kx+1,nd),wst(jx,kx+1,nd)
      real eb3d2(ix,kx,nd),wb3d2(ix,kx,nd),nb3d2(jx,kx,nd),
     $     sb3d2(jx,kx,nd),et3d2(ix,kx,nd),wt3d2(ix,kx,nd),
     $     nt3d2(jx,kx,nd),st3d2(jx,kx,nd)
      real peb2(ix,nd),pwb2(ix,nd),pnb2(jx,nd),
     $     psb2(jx,nd),pet2(ix,nd),pwt2(ix,nd),
     $     pnt2(jx,nd),pst2(jx,nd)
      real web2(ix,kx+1,nd),wwb2(ix,kx+1,nd),
     $     wnb2(jx,kx+1,nd),wsb2(jx,kx+1,nd),
     $     wet2(ix,kx+1,nd),wwt2(ix,kx+1,nd),
     $     wnt2(jx,kx+1,nd),wst2(jx,kx+1,nd)

      read(iunit,err=888) tresc
      read(junit,err=888) tresc2
      print*,'read mm5v1 boundary file output'
      print*,'user has inhyd,imoist,iice,igraupel=',
     $     inhyd,imoist,iice,igraupel
 10   continue
      read(iunit,end=999) tbegin,tend,m1,m2
      print *,'read fort.91 between',tbegin,' and',
     $     tend,'minutes.  dates ',m1,' and',m2
      read(junit,end=999) tbegin,tend,m1,m2
      print *,'read fort.92 between',tbegin,' and',
     $     tend,'minutes.  dates ',m1,' and',m2
      
c
c     p* and tendency
c
      read(iunit) peb,pwb,pnb,psb
      read(iunit) pet,pwt,pnt,pst
      read(junit) peb2,pwb2,pnb2,psb2
      read(junit) pet2,pwt2,pnt2,pst2
c
c     p*u,p*v,p*t,p*q and tendencies
c
      do 700 ivariabs=1,4
c u,v,t,q at boundaries and tendencies at boundaries
         read(iunit) eb3d,wb3d,nb3d,sb3d
         read(iunit) et3d,wt3d,nt3d,st3d
         read(junit) eb3d2,wb3d2,nb3d2,sb3d2
         read(junit) et3d2,wt3d2,nt3d2,st3d2
 700  continue
      print*,'read p*u,p*v,p*t,p*q and tendencies'
c
c     p*w,p*pp and tendencies
c
      if(inhyd.eq.1) then
         read(iunit) web,wwb,wnb,wsb         
         read(iunit) wet,wwt,wnt,wst
         read(iunit) eb3d,wb3d,nb3d,sb3d
         read(iunit) et3d,wt3d,nt3d,st3d
         print*,'read p*w,p*pp,and tendencies from 91'
         read(iunit) web2,wwb2,wnb2,wsb2
         read(iunit) wet2,wwt2,wnt2,wst2
         read(iunit) eb3d2,wb3d2,nb3d2,sb3d2
         read(iunit) et3d2,wt3d2,nt3d2,st3d2
         print*,'read p*w,p*pp,and tendencies from 92'
      endif
c
c     p*clw,p*rnw,and tendencies
c
      if(imoist.eq.2) then
         read(iunit) eb3d,wb3d,nb3d,sb3d
         read(iunit) et3d,wt3d,nt3d,st3d
         read(iunit) eb3d,wb3d,nb3d,sb3d
         read(iunit) et3d,wt3d,nt3d,st3d
         print*,'read p*clw,p*rnw,and tendencies from 91'
         read(iunit) eb3d2,wb3d2,nb3d2,sb3d2
         read(iunit) et3d2,wt3d2,nt3d2,st3d2
         read(iunit) eb3d2,wb3d2,nb3d2,sb3d2
         read(iunit) et3d2,wt3d2,nt3d2,st3d2
         print*,'read p*clw,p*rnw,and tendencies from 92'
      endif
c
c     p*icew,p*snoww,and tendencies
c
      if (iice.eq.1) then
         read(iunit) eb3d,wb3d,nb3d,sb3d
         read(iunit) et3d,wt3d,nt3d,st3d
         read(iunit) eb3d,wb3d,nb3d,sb3d
         read(iunit) et3d,wt3d,nt3d,st3d
         print*,'read p*icew,p*snoww,and tendencies from 91'
         read(iunit) eb3d2,wb3d2,nb3d2,sb3d2
         read(iunit) et3d2,wt3d2,nt3d2,st3d2
         read(iunit) eb3d2,wb3d2,nb3d2,sb3d2
         read(iunit) et3d2,wt3d2,nt3d2,st3d2
         print*,'read p*icew,p*snoww,and tendencies from 92'
      endif
c
c     p*graw,p*npart,and tendencies
c
      if(igraupel.eq.1) then
         read(iunit) eb3d,wb3d,nb3d,sb3d
         read(iunit) et3d,wt3d,nt3d,st3d
         read(iunit) eb3d,wb3d,nb3d,sb3d
         read(iunit) et3d,wt3d,nt3d,st3d
         print*,'read p*graw,p*npart,and tendencies from 91'
         read(iunit) eb3d2,wb3d2,nb3d2,sb3d2
         read(iunit) et3d2,wt3d2,nt3d2,st3d2
         read(iunit) eb3d2,wb3d2,nb3d2,sb3d2
         read(iunit) et3d2,wt3d2,nt3d2,st3d2
         print*,'read p*graw,p*npart,and tendencies from 92'
      endif
      goto 1
c
 888  print*,'input data is not in mm5v1 format'
      stop
 999  continue
      return
      end
