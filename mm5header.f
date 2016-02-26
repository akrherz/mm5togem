      program mm5togem
      include 'GEMINC:GEMPRM.PRM'
      parameter(ixl=120,jyl=130,kzl=40,maxtimes=20)
      parameter(kzp=kzl,num_out_flds=40)

      character*80 mrfc(1000,20),mifc(1000,20),
     $     filein,outnam, innam,tempfile,chdum
      character*72 oldproj,cproj
      character*40 tmpa
      character*20 gdattm(2)
      character*12 parm
      character dummy

      real mrf(1000,20),rnvblk(LLNNAV),anlblk(LLNANL)
c  Raw data 3-d fields:
      real uu(ixl,jyl,kzl),      vv(ixl,jyl,kzl),   tt(ixl,jyl,kzl)
      real qv(ixl,jyl,kzl),      qc(ixl,jyl,kzl),   qr(ixl,jyl,kzl)
      real cice(ixl,jyl,kzl),    sice(ixl,jyl,kzl), ww(ixl,jyl,kzl+1)
      real pp(ixl,jyl,kzl),      rh(ixl,jyl,kzl),   rtten(ixl,jyl,kzl)
c  Raw precip 2-d fields
      real re(ixl,jyl), rc(ixl,jyl)
      real tot_prec(ixl,jyl), prec_3hr(ixl,jyl)
     $     ,re_3hr(ixl,jyl),rc_3hr(ixl,jyl)
     $     ,last_output_hr_tot_prec(ixl,jyl)
     $     ,last_output_hr_exp_prec(ixl,jyl)
     $     ,last_output_hr_con_prec(ixl,jyl)
c  Raw data 2-d fields:
      real pstx(ixl,jyl),        tg(ixl,jyl)
      real ter(ixl,jyl),         xmap(ixl,jyl),     dotmap(ixl,jyl)
      real dotcor(ixl,jyl),      irst(ixl,jyl),     xlat(ixl,jyl)
      real xluse(ixl,jyl),       xsnowc(ixl,jyl),   xlon(ixl,jyl)
c  Derived fields:
      real prs(ixl,jyl,kzl),     pstd(ixl,jyl),     dew(ixl,jyl,kzl)
      real sfp(ixl,jyl),         slp(ixl,jyl),      terd(ixl,jyl)
      real dotvar(ixl,jyl),      dat(ixl*jyl),      dluse(ixl,jyl)
      real wwh(ixl,jyl,kzl),     tmp2d(ixl,jyl),    prsd(ixl,jyl,kzl)
      real pdata(ixl,jyl,kzp),   ttp(ixl,jyl,kzp),  prsp(ixl,jyl,kzp)
      real tmp3dp(ixl,jyl,kzp),  tmp3ds(ixl,jyl,kzl),hlhsl(ixl,jyl,2)
c  Miscellaneous 2-d fields:
      real tvirt(ixl,jyl)
c  Miscellaneous 1-d fields:
      real plevs(kzp),maxp(kzl),minp(kzl),maxpd(kzl),minpd(kzl)
      real sih(kzl),             sif(kzl)
c     optional fields                    
      REAL senhtflx(ixl,jyl), lathtflx(ixl,jyl), 
     &     lonwvrad(ixl,jyl),
     &     shtwvrad(ixl,jyl), pblht(ixl,jyl), ustar(ixl,jyl),
     &     pblregime(ixl,jyl), moisture(ixl,jyl)
      character i_output(num_out_flds)*1,sig_fld(num_out_flds)*4
     $     ,prs_fld(num_out_flds)*4,fld_name(num_out_flds)*70


      integer mif(1000,20),ilvl(2)
      integer aa,bb,cc,dd
      integer cmn(maxtimes),iptimes(20),output_hr,last_output_hr
      integer dom_special, special_hour
c  Frequencies at which precip may be written out:
c  Note:  0 writes out precip since last output time as P00M if output
c         is more often than 1-hourly.  99 is skipped
      data navsz,ianlsz,ihdrsz/LLNNAV,LLNANL,2/
      data iptimes/0,1,99,3,99,6,99,12,24,99,99,99,99,99,99,99,99,
     $     99,99,99/
      data ighttype/1/ !flag determining how geopotential height is calculated
c                      ighttype=1 means calculate height directly on pressure
c                               levels.
c                      ighttype=2 means calculate height on sigma levels and
c                               interpolate to pressure levels.  Extrapolation
c                               below surface uses standard lapse rate.
c       I Suggest using ighttype=1 with many pressure levels, and few sigma
c       levels, or 2 with many sigma levels and few pressure levels.  The 
c       Difference is small at least in case of 36 sigma levels and 19 
c       pressure levels.
      logical respond,exists,write,rewrite,ifprs,ifsig,ifchoose
      logical ifuu,ifvv,iftt,ifqv,ifqc,ifqr,ifcice,ifsice,ifww,ifpp
     $     ,ifrtten
      logical ifcon,ifres,quiet,angflg

c     additions for gempak library calls
      parameter (nhdrsz=2)
      integer ighdr(nhdrsz)
      real gbnds(4),ebnds(4)
      integer iebnds(4),ihdrsz
      real angle1,angle2,angle3 

      eps = 0.622
c
c   Check for command line arguements
c
      itapfrq = 0
      maxgrd  = 0
      igemfrq = 0
      innam   = ''
      outnam  = ''
      quiet   = .false.
      i = 1
c
c     get information from input file, rather than hard-wire these
c
      filein = 'input_file'
      open (unit=10,file=filein
     $     ,form='formatted',status='unknown',err=999)
c     read info about how different types of observations are stored
c     in the ftype variable in the program
      read(10,1005) innam
 1005 format(/,a80)
      read(10,1005) outnam
      read (10,*) dummy
      read(10,*) tempfile,dom_special,special_hour
      read(10,1006) ioverwrite,issig,isprs,iprecise
 1006 format(/,i1,//,i1,//,i1,//,i1)
      read(10,1007) issfp,num_pres_levs
      read(10,*) (plevs(kdum),kdum=1,num_pres_levs)
      ifchoose=.false.
      if (issig.eq.1) then
         ifsig = .true.
      else
         ifsig = .false.
      endif
      if (isprs.eq.1) then
         ifprs = .true.
      else
         ifprs = .false.
      endif
      read(10,*) dummy
      read(10,*) dummy
      read(10,*) dummy
      read(10,*) dummy
      read(10,*) dummy
      read(10,*) num_3d_outputs
      if (num_3d_outputs.gt.num_out_flds) then
         print*,'too many outputs'
         goto 999
      endif
      do iline = 1,num_3d_outputs
         read(10,1008) i_output(iline),sig_fld(iline),prs_fld(iline)
     $        ,fld_name(iline)
      enddo
      read(10,*) dummy
      read(10,*) dummy
      read(10,*) dummy
      read(10,*) dummy
      read(10,*) dummy
      read(10,*) num_2d_outputs
      if ((num_3d_outputs+num_2d_outputs).gt.num_out_flds) then
         print*,'too many 2d outputs'
         goto 999
      endif
      do iline = num_3d_outputs+1,num_3d_outputs+num_2d_outputs
         read(10,1008) i_output(iline),sig_fld(iline),prs_fld(iline)
     $        ,fld_name(iline)
      enddo
 1007 format(/,/,i1,x,i3)
 1008 format(a1,x,a4,x,a4,x,a)
      close (10)

      do while (i.le.iargc())
        call getarg(i,tmpa)
        if(tmpa(1:2).eq.'-i')then       !override output frequency from header
          call getarg(i+1,tmpa)
          read(tmpa,fmt="(i7)") itapfrq
          i=i+1
        elseif(tmpa(1:2).eq.'-g')then   !specify maximum number of grids
          call getarg(i+1,tmpa)
          read(tmpa,fmt="(i5)") maxgrd
          if(maxgrd.gt.MMHDRS)maxgrd=MMHDRS
          i=i+1
        elseif(tmpa(1:2).eq.'-f')then   !specify mmout file on command line
          call getarg(i+1,innam)
          i=i+1
        elseif(tmpa(1:2).eq.'-o')then   !specify gempak file on command line
          call getarg(i+1,outnam)
          i=i+1
        elseif(tmpa(1:2).eq.'-q')then   !specify quiet mode
          quiet=.true.
        elseif(tmpa(1:2).eq.'-t')then   !specify time interval to write to
          call getarg(i+1,tmpa)         !gempak file.
          read(tmpa,fmt="(i5)") igemfrq  
          i=i+1
        else
          print*,'Argument error.  Available options are:',
     &'-i <mm5 output freq in minutes>  To override tapfrq in header.',
     &'-g <max # of grids in file>  Create gempak file with number',
     &'                             of grids less than maximum.',
     &'-f <mm5 output file>         User specified file for input',
     &'-o <gempak output file>      User specified file for output',
     &'-q                           Quiet mode. Just write out',
     &'                             all times without prompting user.',
     &'-t (interval(min.)for output> Specify time interval to write',
     &'                                      to gempak file.',
     &'conv can also be run with no commmand line arguements'
          stop
        endif
        i=i+1
      enddo
      open(unit=22,file=tempfile,form='unformatted')
c
c Initialize all logical variables defining output fields to true
c
      ifuu=.false.
      ifvv=.false.
      iftt=.false.
      ifqv=.false.
      ifqc=.false.
      ifqr=.false.
      ifqcice=.false.
      ifww=.false.
      ifrtten=.false.
      ifcon=.false.
      ifres=.false.
c

      print *, '          MM5 HEADER'
      print *, ' '
      print *, 'This Program prints out raw MM5 Version 1 and 2'
      print *, ' '
 
10    continue
c
c   Read header:
c
      open(unit=20,name=innam,status='old',form='unformatted')
      read (20,END=999) mif,mrf,mifc,mrfc
      iprog=mif(1,1)
      if(iprog.ne.1.and.iprog.ne.2.and.iprog.ne.3)then
c       deleted stuff to put in input_file
      elseif(iprog.eq.1)then
        ifprs=.false.
        ifsig=.false.
      elseif(iprog.eq.2.or.iprog.eq.3)then
        ifprs=.true.
        ifsig=.false.
      endif
c     from /usr/local/gempak/gempak5.4/source/gemlib/ip/ipinit.f
c     initialize GEMPAK user interface
      call ip_init(respond,iret)
      if(iret.ne.0)stop 1
c     from  /usr/local/gempak/gempak5.4/source/gemlib/fl/flinqr.f
c     gempak determination of file existence
      call fl_inqr(outnam,exists,dummy,iret)
      if(iret.ne.0)stop 2
      if(exists .and. ioverwrite.eq.0) goto 999
      if(exists .and. ioverwrite.eq.1) then
         write(chdum,846) outnam(1:72)
 846     format('/bin/rm ',a72)
         call system(chdum)
         exists=.false.
      endif
c
c   Hardwire some constants:
c
      ilvl(2) = -1
      write   = .true.
      rewrite = .true.
      ipktyp  = mdgdec
c      ipktyp  = MDGNON          ! no grid packing
c     iprecise   = 2             ! read in from input_file now
      itime   = 0
      iout    = 0
      idot    = 0
      icross  = 1
c
c   Process some header information:
c
      ix = mif(104,1)
      jy = mif(105,1)
      kz = nint(mrf(101,iprog))
      if(iprog.eq.2.or.iprog.eq.3)kz=mif(101,iprog)
      if(itapfrq.ne.0)then
        print*,'Output every',itapfrq,' minutes according to user.'
      else
        itapfrq= nint(mrf(302,6))
        if(iprog.eq.2.or.iprog.eq.3)itapfrq=mif(5,iprog)*60
      endif
      ptop=mrf(1,2)
      print *, 'I = ',ix,' J = ',jy,' K = ', kz
      if(ix.gt.ixl.or.jy.gt.jyl.or.kz.gt.kzl) then
         print *, 'Sorry, Model Arrays too big!'
         print *, 'Need to recompile program'
         goto 999
      endif
c
c   Find half sigma levels:
c
      do i = 1,kz
         sih(i) = mrf(101+i,iprog)
      enddo
c      
c   Find full sigma levels:             ! calculated, but never used
c
      sif(1)=0.
      sif(kz+1)=1000.
      do i = 2,kz
        sif(i)=(2*nint(mrf(101+(i-1),iprog)*10000.0)) - sif(i-1) 
      enddo
c
c  Find the highest pressure level to use
c
      kp=num_pres_levs
      do k=num_pres_levs,1,-1
        if(ptop.gt.plevs(k))kp=k-1
      enddo
c
c     Calculate map area           ! simpler than Ken's old version using
c      refrat=mrf(1,1)/mrf(101,1)
c      yllc=mif(106,1)+(0)/refrat
c      yurc=mif(106,1)+(ix-1.)/refrat
c      xllc=mif(107,1)+(0)/refrat
c      xurc=mif(107,1)+(jy-1.)/refrat
c
      yllc=mrf(102,1)
      xllc=mrf(103,1)
      yurc=mrf(104,1)
      xurc=mrf(105,1)
      call xyll(xllc,yllc,xlatll,xlonll,mrf(4,1),mrf(1,1),
     +          mrf(2,1),mrf(3,1),mif(2,1),mif(3,1))
      call xyll(xurc,yurc,xlatur,xlonur,mrf(4,1),mrf(1,1),
     +          mrf(2,1),mrf(3,1),mif(2,1),mif(3,1))
      print *,'your corners are: ',xlatll,xlonll
      print *,'and: ',xlatur,xlonur
      print *,'continue?'
      read (*,*) chdum
      if (chdum(1:1).ne.'y') then
         stop
      endif
c
c   Create new gempak grid file if necessary
c
      if(exists) then  !check dimensions of old gempak file
        call gd_opnf(outnam,write,igdfln,navsz,rnvblk,ianlsz,anlblk,
     $               ihdrsz,maxgrd,iret)

        call gr_rnav(rnvblk,oldproj,ioldkx,ioldky,iret)
        if(iret.ne.0)then
          print*,'gr_rnav: iret=',iret
          goto 999
        endif
        if(ioldkx.ne.jy.or.ioldky.ne.ix) then
          print*,'Model array size does not match grid file.'
          print*,'Quit, or retry(q/r)?'
          read(5, '(A)') tmpa
          if(tmpa(1:1).eq.'r') then
             close(20)
             goto 10
          endif
          if(tmpa(1:1).eq.'q')goto 999
        endif
      else                          !must create new gempak file
c
c#######################################################################
c
c     Build navigation block
c
c#######################################################################
c
C************************************************************************
C* GR_MNAV								*
C*									*
C* This subroutine makes a navigation block for a grid file.  The 	*
C* projection may be any simple, full or graph projection.  If 		*
C* ANGFLG is set, the projection must be a full map projection.  	*
C* Otherwise, a simple map projection will be defined.			*
C*									*
C* GR_MNAV  ( PROJ, KX, KY, RLAT1, RLON1, RLAT2, RLON2, ANGL1, 		*
C*            ANGL2, ANGL3, ANGFLG, RNVBLK, IRET )			*
C*									*
C* Input parameters:							*
C*	PROJ		CHAR*		Projection name			*
C*	KX		INTEGER		Number of x grid points		*
C*	KY		INTEGER		Number of y grid points 	*
C*	RLAT1		REAL		Lower left latitude/x		*
C*	RLON1		REAL		Lower left longitude/y		*
C*	RLAT2		REAL		Upper right latitude/x		*
C*	RLON2		REAL		Upper right longitude/y		*
C*	ANGL1		REAL		Projection angle 1		*
C*	ANGL2		REAL		Projection angle 2		*
C*	ANGL3		REAL		Projection angle 3		*
C*	ANGFLG		LOGICAL		Full projection flag		*
C*									*
C* Output parameters:							*
C*	RNVBLK (LLNNAV)	REAL		Navigation block		*
C*	IRET		INTEGER		Return code			*
C*					  0 = normal return		*
C************************************************************************
         if(mif(4,1).eq.1)then  !lambert conic conformal projection
            cproj = 'LCC'       !projection type (lambert conformal=LCC)
            angflg=.true.       !full map projection
            angle1=  mrf(6,1)   !standard latitude 1 for GEMPAK (true lat 2 MM5)
            angle2=  mrf(3,1)   !central longitude for GEMPAK (pole longitude in MM5)
            angle3=  mrf(5,1)   !standard latitude 2 for GEMPAK (true lat 1 MM5)
         elseif(mif(4,1).eq.2)then !polar stereographic projection
            cproj = 'STR'       !projection type (polar stereographic=STR)
            angflg=.true.       !full map projection
            angle1= mrf(7,1)    !latitude of point of tangency
            angle2= mrf(3,1)    !longitude of point of tangency
            angle3= 0.          !not used
         elseif(mif(4,1).eq.3)then !Mercator projection
            cproj = 'MER'       !projection type (Mercator=STR)
            angflg=.true.       !full map projection
            angle1= mrf(7,1)    !latitude of point of tangency
            angle2= mrf(3,1)    !longitude of point of tangency
            angle1= 0.          !not used
            angle3= 0.          !not used
         else
            print*,'This map projection is not supported'
            goto 999
         endif 
         CALL GR_MNAV  ( cproj, jy, ix, xlatll, xlonll, xlatur, xlonur,
     +        angle1, angle2, angle3, angflg,
     +        rnvblk, iret )
         IF( iret.ne.0 ) THEN
            write(6,820) iret
 820        FORMAT('  Error building GEMPAK navigation block:',i6)
            STOP
         END IF
         print *,'built GEMPAK navigation block ok'
c
c#######################################################################
c
c     Build analysis block
c
c#######################################################################
c
C************************************************************************
C* GR_MBAN								*
C*									*
C* This subroutine makes a Barnes analysis block.  The analysis block	*
C* generated is LLNANL words long.  All the bounds must be entered in 	*
C* the order:  lower left latitude; lower left longitude; upper		*
C* right latitude; upper right longitude.				*
C*									*
C* GR_MBAN  ( DELTAN, DELTAX, DELTAY, GBNDS, EBNDS, DBNDS, ANLBLK,	*
C*            IRET )							*
C*									*
C* Input parameters:							*
C*	DELTAN		REAL		Station spacing			*
C*	DELTAX		REAL		Grid spacing in x dir		*
C*	DELTAY		REAL		Grid spacing in y dir		*
C*	GBNDS (4)	REAL		Grid area bounds		*
C*	EBNDS (4)	REAL		Extended area bounds		*
C*	DBNDS (4)	REAL		Data area bounds		*
C*									*
C* Output parameters:							*
C*	ANLBLK (LLNANL)	REAL		Analysis block			*
C*	IRET		INTEGER		Return code			*
C*					  0 = normal return		*
C**									*
C* Log:									*
C* M. desJardins/GSFC	 2/85						*
C* M. desJardins/GSFC	 9/88	GEMPAK4					*
C* G. Huffman/GSC	 4/89	Document removing GAMMA			*
C* K. Brill/NMC		02/92	Use LLNANL				*
C************************************************************************
         deltan=2.*mrf(1,1)     !station spacing(deg) for barnes analysis
     &        /111.19842        !roughly twice the grid spacing
         deltax= RMISSD
         deltay= RMISSD
         gbnds(1)=xlatll
         gbnds(2)=xlonll
         gbnds(3)=xlatur
         gbnds(4)=xlonur
         ebnds(1)= xlatll - 1.
         ebnds(2)= xlonll
         ebnds(3)= xlatur
         ebnds(4)= xlonur + 1.

         CALL GR_MBAN  ( deltan, deltax, deltay,
     +        gbnds, ebnds, ebnds, anlblk, iret )
         IF( iret.ne.0 ) THEN
            write(6,825) iret
 825        FORMAT('  Error GEMPAK analysis block:',i6)
            STOP
         END IF
         print *,'built GEMPAK analysis block ok'
         if(maxgrd.eq.0)then
            maxgrd=MMHDRS
         endif
c
c#######################################################################
c
c     Create the grid file
c
c#######################################################################
c
C************************************************************************
C* GD_CREF								*
C*									*
C* This subroutine creates a new GEMPAK5 grid file.  If MAXGRD is zero	*
C* or negative, it will default to 400.  IHDRSZ is the length of	*
C* the grid header which will be stored with every grid.  This		*
C* header is intended to save offsets from a base grid, but is not	*
C* currently used.  IHDRSZ should usually be set to 2.			*
C*									*
C* GD_CREF  ( FILNAM, NAVSZ, RNVBLK, IANLSZ, ANLBLK, IHDRSZ,		*
C*            MAXGRD, IGDFLN, IRET )					*
C*									*
C* Input parameters:							*
C*	FILNAM		CHAR*		File name			*
C*	NAVSZ		INTEGER		Navigation blck length (LLNNAV)	*
C*	RNVBLK (NAVSZ)	REAL		Navigation block		*
C*	IANLSZ		INTEGER		Analysis block length (LLNANL)	*
C*	ANLBLK (IANLSZ)	REAL		Analysis block			*
C*	IHDRSZ		INTEGER		Grid header length		*
C*	MAXGRD		INTEGER		Max number of grids in file	*
C*									*
C* Output parameters:							*
C*	IGDFLN		INTEGER		Grid file number		*
C*	IRET		INTEGER		Return code			*
C*					  0 = normal return		*
C*					 -1 = file cannot be created	*
C*					-13 = grid header too long	*
C**									*
         ighdr(1)=0             ! used later in gd_wdat
         ighdr(2)=0             ! used later in gd_wdat
         call gd_cref(outnam,navsz,rnvblk,ianlsz,anlblk,ihdrsz,maxgrd,
     &        igdfln,iret)
         if(iret.ne.0)then
            print*,'gd_cref:  problem creating output gempak file:'
     $           ,outnam
            goto 999
         endif
      endif                     ! end of creating new gempak file or reading old one

30    if(itime.ne.0)then
        read (20,END=999) mif,mrf,mifc,mrfc
      endif
      itime = itime + 1
      num3d=mif(201,iprog)
      num2d=mif(202,iprog)
c
c   Read data
c
c  3-D fields:
      do ifld=205,204+num3d
      if(mifc(ifld,iprog)(1:8).eq.'U       ')then       ! U COMPONENT OF WIND
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((( uu(i,j,k),i=1,ix),j=1,jy),k=1,kz)
        if(itime.eq.1)then
          if(ifchoose)then
            print*,'Would you like this field written out(y/n)?'
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')ifuu=.true.
          else
            ifuu=.true.
          endif
        endif
      elseif(mifc(ifld,iprog)(1:8).eq.'V       ')then   ! V COMPONENT OF WIND
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((( vv(i,j,k),i=1,ix),j=1,jy),k=1,kz)
        if(itime.eq.1)then
          if(ifchoose)then
            print*,'Would you like this field written out(y/n)?'
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')ifvv=.true.
          else
            ifvv=.true.
          endif
        endif
      elseif(mifc(ifld,iprog)(1:8).eq.'T       ')then   ! TEMPERATURE
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((( tt(i,j,k),i=1,ix),j=1,jy),k=1,kz)
        if(itime.eq.1)then
          if(ifchoose)then
            print*,'Would you like this field written out(y/n)?'
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')iftt=.true.
          else
            iftt=.true.
          endif
        endif
      elseif(mifc(ifld,iprog)(1:8).eq.'Q       '.or.
     &       mifc(ifld,iprog)(1:8).eq.'RH      ')then   ! QV or RH
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((( qv(i,j,k),i=1,ix),j=1,jy),k=1,kz)
        if(itime.eq.1)then
          if(ifchoose)then
            if(iprog.eq.2.or.iprog.eq.3)then

               print*,
     $           'Would you like RELATIVE HUMIDITY written out(y/n)?'
            else
              print*,'Would you like DEWPOINT written out(y/n)?'
            endif
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')ifqv=.true.
          else
            ifqv=.true.
          endif
        endif
      elseif(mifc(ifld,iprog)(1:8).eq.'CLW     ')then   ! QV CLOUD
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((( qc(i,j,k),i=1,ix),j=1,jy),k=1,kz)
        if(itime.eq.1)then
          if(ifchoose)then
            print*,'Would you like this field written out(y/n)?'
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')ifqc=.true.
          else
            ifqc=.true.
          endif
        endif
      elseif(mifc(ifld,iprog)(1:8).eq.'RNW     ')then   ! QV RAIN
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((( qr(i,j,k),i=1,ix),j=1,jy),k=1,kz)
        if(itime.eq.1)then
          if(ifchoose)then
            print*,'Would you like this field written out(y/n)?'
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')ifqr=.true.
          else
            ifqr=.true.
          endif
        endif
      elseif(mifc(ifld,iprog)(1:8).eq.'ICE     ')then   ! CLOUD ICE
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((( cice(i,j,k),i=1,ix),j=1,jy),k=1,kz)
        if(itime.eq.1)then
          if(ifchoose)then
            print*,'Would you like this field written out(y/n)?'
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')ifcice=.true.
          else
            ifcice=.true.
          endif
        endif
      elseif(mifc(ifld,iprog)(1:8).eq.'SNOW    ')then   ! SNOW ICE
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((( sice(i,j,k),i=1,ix),j=1,jy),k=1,kz)
        if(itime.eq.1)then
          if(ifchoose)then
            print*,'Would you like this field written out(y/n)?'
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')ifsice=.true.
          else
            ifsice=.true.
          endif
        endif
      elseif(mifc(ifld,iprog)(1:8).eq.'GRAW    ')then   ! GRAUPEL
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        print*,'(this field will not be written to gempak file)'
        read (20) ((( tmp3ds(i,j,k),i=1,ix),j=1,jy),k=1,kz)
      elseif(mifc(ifld,iprog)(1:8).eq.'NICEPART')then   !  Number concentration
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        print*,'(this field will not be written to gempak file)'
        read (20) ((( tmp3ds(i,j,k),i=1,ix),j=1,jy),k=1,kz)
      elseif(mifc(ifld,iprog)(1:8).eq.'W       ')then   ! VERTICAL WIND COMP.
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((( ww(i,j,k),i=1,ix),j=1,jy),k=1,kz+1)
        if(itime.eq.1)then
          if(ifchoose)then
            print*,'Would you like this field written out(y/n)?'
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')ifww=.true.
          else
            ifww=.true.
          endif
        endif
      elseif(mifc(ifld,iprog)(1:8).eq.'PP      '.or.    !PRESSURE PERTURBATION
     &       mifc(ifld,iprog)(1:8).eq.'H       ')then   !or GEOPOTENTIAL HEIGHT
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((( pp(i,j,k),i=1,ix),j=1,jy),k=1,kz)
        if(itime.eq.1)then
          if(ifchoose)then
            if(iprog.eq.2.or.iprog.eq.3)then
              print*,
     $         'Would you like GEOPOTENTIAL HEIGHT written out(y/n)?'
            else
              print*,
     $         'Would you like PRESSURE and/or GEOPOTENTIAL HEIGHT',
     &         ' written out(y/n)?'
            endif
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')ifpp=.true.
          else
            ifpp=.true.
          endif
        endif
      elseif(mifc(ifld,iprog)(1:8).eq.'RTTEN   ')then   ! Radiative Tendency
        ifrtten = .true.
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((( rtten(i,j,k),i=1,ix),j=1,jy),k=1,kz)
        if(itime.eq.1)then
          if(ifchoose)then
            print*,'Would you like this field written out(y/n)?'
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')ifrtten=.true.
          else
            ifrtten=.true.
          endif
        endif
      else                                              ! Unknown field
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        print*,'(this field will not be written to gempak file)'
        read (20) ((( tmp3ds(i,j,k),i=1,ix),j=1,jy),k=1,kz)
      endif
      enddo
C  2-D fields:
      do ifld=205+num3d,204+num3d+num2d
      if(mifc(ifld,iprog)(1:8).eq.'PSTARCRS')then       ! PSTAR
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((pstx(i,j),i=1,ix),j=1,jy)
      elseif(mifc(ifld,iprog)(1:8).eq.'GROUND T')then   ! GROUND TEMPERATURE
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((tg(i,j),i=1,ix),j=1,jy)
      elseif(mifc(ifld,iprog)(1:8).eq.'RAIN CON')then   ! ACCUMULATED CONV PPT
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((rc(i,j),i=1,ix),j=1,jy)
        if(itime.eq.1)then
          if(ifchoose)then
            print*,'Would you like this field written out(y/n)?'
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')ifcon=.true.
          else
            ifcon=.false.
          endif
        endif
      elseif(mifc(ifld,iprog)(1:8).eq.'RAIN NON')then   ! ACCUMULATED RES PPT
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((re(i,j),i=1,ix),j=1,jy)
        if(itime.eq.1)then
          if(ifchoose)then
c            print*,'Would you like this field written out(y/n)?'
            print*,'Would you like 3-hr precip written out(y/n)?'
            read(5, '(A)') tmpa
            if(tmpa(1:1).eq.'y')ifres=.true.
          else
            ifres=.true.
          endif
        endif
      elseif(mifc(ifld,iprog)(1:8).eq.'TERRAIN ')then   ! TERRAIN
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((ter(i,j),i=1,ix),j=1,jy)
      elseif(mifc(ifld,iprog)(1:8).eq.'MAPFACCR')then   ! MAP FACTOR CROSS
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((xmap(i,j),i=1,ix),j=1,jy)
      elseif(mifc(ifld,iprog)(1:8).eq.'MAPFACDT')then   ! MAP FACTOR DOT
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((dotmap(i,j),i=1,ix),j=1,jy)
      elseif(mifc(ifld,iprog)(1:8).eq.'CORIOLIS')then   ! CORIOLIS
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((dotcor(i,j),i=1,ix),j=1,jy)
      elseif(mifc(ifld,iprog)(1:8).eq.'RES TEMP')then   ! RESERVOIR TEMPERATURE
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((irst(i,j),i=1,ix),j=1,jy)
      elseif(mifc(ifld,iprog)(1:8).eq.'LATITCRS')then   ! LATITUDE
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((xlat(i,j),i=1,ix),j=1,jy)
      elseif(mifc(ifld,iprog)(1:8).eq.'LONGICRS')then   ! LONGITUDE
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((xlon(i,j),i=1,ix),j=1,jy)
      elseif(mifc(ifld,iprog)(1:8).eq.'LAND USE')then   ! LAND USE CATEGORY
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((xluse(i,j),i=1,ix),j=1,jy)
      elseif(mifc(ifld,iprog)(1:8).eq.'SNOWCOVR')then   ! SNOW COVER
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        read (20) ((xsnowc(i,j),i=1,ix),j=1,jy)
c
c     Ken's Additions for Boundary Layer Fields
c
C        14) SENSIBLE HEAT FLUX
      elseif(mifc(ifld,iprog)(1:18).eq.'SENSIBLE HEAT FLUX')then
        print*,'Reading ', mifc(ifld,iprog)(1:67)
        read (20) ((senhtflx(i,j),i=1,ix),j=1,jy)
C        15) LATENT HEAT FLUX
      elseif(mifc(ifld,iprog)(1:16).eq.'LATENT HEAT FLUX')then
        print*,'Reading ', mifc(ifld,iprog)(1:67)
        read (20) ((lathtflx(i,j),i=1,ix),j=1,jy)
C        16) LONGWAVE RADIATION
      elseif(mifc(ifld,iprog)(1:8).eq.'LONGWAVE')then
        print*,'Reading ', mifc(ifld,iprog)(1:67)
        READ (20) ((lonwvrad(i,j),i=1,ix),j=1,jy)
C        17) SHORTWAVE RADIATION
      elseif(mifc(ifld,iprog)(1:9).eq.'SHORTWAVE')then
        print*,'Reading ', mifc(ifld,iprog)(1:67)
        READ (20) ((shtwvrad(i,j),i=1,ix),j=1,jy)
C        18) PLANETARY BOUNDARY LAYER HEIGHT
      elseif(mifc(ifld,iprog)(1:9).eq.'PLANETARY')then
        print*,'Reading ', mifc(ifld,iprog)(1:67)
        READ (20) ((pblht(i,j),i=1,ix),j=1,jy)
C        19) FRICTION VELOCITY
      elseif(mifc(ifld,iprog)(1:8).eq.'FRICTION')then
        print*,'Reading ', mifc(ifld,iprog)(1:67)
        READ (20) ((ustar(i,j),i=1,ix),j=1,jy)
C        20) PBL REGIME
      elseif(mifc(ifld,iprog)(1:10).eq.'PBL REGIME')then
        print*,'Reading ', mifc(ifld,iprog)(1:67)
        READ (20) ((pblregime(i,j),i=1,ix),j=1,jy)
C        21) MOISTURE AVAILABILITY
      elseif(mifc(ifld,iprog)(1:21).eq.'MOISTURE AVAILABILITY')then
        print*,'Reading ', mifc(ifld,iprog)(1:67)
        READ (20) ((moisture(i,j),i=1,ix),j=1,jy)
c
c     end of Ken's additions to MM5 output files
c
      elseif(mifc(ifld,iprog)(1:8).eq.'PSEALVLD')then   !Sea-level pres
        print*,'Reading ', mifc(ifld,iprog)(27:67)      !from datagrid
        read (20) ((slp(i,j),i=1,ix),j=1,jy)
      else                                              ! Unknown field
        print*,'Reading ', mifc(ifld,iprog)(27:67)
        print*,'(this field will not be written to gempak file)'
        read (20) ((tmp2d(i,j),i=1,ix),j=1,jy)
      endif                     ! enf of scanning the field names
      enddo                     ! end of looping through the 2-d fields

      if (itime.eq.1)then       ! adjust end points of land use
        do aa = 1,ix
        do bb = 1,jy
           if(aa.eq.ix) then
              xluse(aa,bb) = xluse(aa-1,bb)
           endif
           if(bb.eq.jy) then
              xluse(aa,bb) = xluse(aa,bb-1)
           endif
        enddo
        enddo
        call xtodot2(ter,terd,ixl,jyl,ix,jy)
        call xtodot2(xluse,dluse,ixl,jyl,ix,jy)
      endif
c
c     ?? I don't really understand yet what cmn(itime) is doing
c
      if(iprog.eq.5)then
        cmn(itime) = nint(mrf(1,iprog))
      else
        cmn(itime) = (itime-1)*itapfrq
      endif
      print*,' '
c
c  Calculate date/time info and possibly prompt user
c
      if(iprog.eq.5)then        !interp output file / mminput file
         if(mod(nint(mrf(1,5)),60).ne.0)then
            print *, 'Date =',mif(1,iprog), ' Forecast Time =',
     &           int(mrf(1,iprog)/60.),' hours and',nint(mrf(1,iprog)
     &           - int(mrf(1,iprog)/60.)*60),' minutes'
         else
            print *, 'Date =',mif(1,iprog), ' Forecast Time =',
     &           nint(mrf(1,iprog)/60.0),' hours'
         endif
         iyr = mif(4,5)/100
         ihr = mod(mif(4,5),100)*100
         iftim=nint(mrf(1,5))
      elseif(iprog.eq.6)then    !output from mm5 - - more time info
         idate1=mif(13,6)+100*mif(14,6)+mif(15,6)*10000+mif(16,6)*1000000
         idate2=mif(11,6)+100*mif(12,6)
         if(mod(itapfrq,60).ne.0)then
            write(6,fmt='(a6,1x,i8.8,a1,i4.4,1x,a15,1x,i4,1x,a9,1x,i2
     $           ,1x,a7)')
     &           'Date =',idate1,'/',idate2,'Forecast Time =',
     &           int(mrf(1,iprog)/60.),'hours and',nint(mrf(1,iprog)
     &           - int(mrf(1,iprog)/60.)*60),'minutes'
         else
            write(6,fmt='(a6,1x,i8.8,a1,i4.4,1x,a15,1x,i4,1x,a5)')
     &           'Date =',idate1,'/',idate2,'Forecast Time =',
     &           nint(mrf(1,iprog)/60.0),'hours'
         endif
         iyr = mif(1,5)/100
         ihr = mod(mif(1,5),100)*100
         iftim=(itime-1)*itapfrq + int(mrf(1,6))
      elseif(iprog.eq.1)then    !output from program TERRAIN
         iyr=010101
         ihr=0
         ifhr=0
         ifmn=0
         write(gdattm(1),fmt="(i6.6,a1,i4.4,a1,i3.3,i2.2,'   ')")
     &        iyr,'/',ihr,'F',ifhr,ifmn
         goto 321               !can skip all the following calculations for terrain
      elseif(iprog.eq.2.or.iprog.eq.3)then
         print*,'Date=',mif(1,iprog)
         iyr = mif(4,iprog)/100
         ihr = mod(mif(4,iprog),100)*100
         iftim=(itime-1)*itapfrq
      elseif(iprog.eq.4.or.iprog.gt.6.or.iprog.lt.1)then
         print*,'This output is not from interp, mm5, or terrain, ',
     &        'datagrid, or rawins . . .stopping'
         goto 999
      endif

      ifhr  = nint(float(iftim)/60.)
c     new addition 8/13/97 to fix for Brad
c      ifmn  = iftim - ifhr*60
c      write(gdattm(1),fmt="(i6.6,a1,i4.4,a1,i3.3,i2.2,'   ')")
c     &     iyr,'/',ihr,'F',ifhr,ifmn
      write(gdattm(1),fmt="(i6.6,a1,i4.4,a1,i3.3,'   ')")
     &     iyr,'/',ihr,'F',ifhr


      tmpa='n'
      if(quiet)then             !don't prompt user
         tmpa='n'
      elseif(igemfrq.ne.0)then  !also don't prompt user
         if((iprog.eq.6.and.mod(itapfrq*(itime-1),igemfrq).eq.0).or.
     &        (iprog.eq.5.and.mod(nint(mrf(1,5)),   igemfrq).eq.0).or.
     &        (iprog.eq.2.and.mod(itapfrq*(itime-1),igemfrq).eq.0).or.
     &        (iprog.eq.3.and.mod(itapfrq*(itime-1),igemfrq).eq.0))then
            tmpa='n'
         else
            tmpa='y'
         endif
      endif

      if(tmpa(1:1).eq.'y') then
         goto 30
      endif
      if(tmpa(1:1).eq.'q') goto 999
      iout=iout+1
      if(iprog.eq.2.or.iprog.eq.3)goto 678
CCCCC
CC calculate tot_prec and convert units
CCCCC
      do i=1,ix-1
        do j=1,jy-1
           tot_prec(i,j) = (10.0*rc(i,j))+(10.0*re(i,j))
           re(i,j) = 10.0*re(i,j)
           rc(i,j) = 10.0*rc(i,j)
        enddo
      enddo
c
c  calculate pstar on dot points (pstd)
c
      do i=1,ix
      do j=1,jy
        aa = min(i,ix-1)
        bb = min(j,jy-1)
        cc = max(i-1,1)
        dd = max(j-1,1)
        pstd(i,j)= 0.25*(pstx(aa,bb)+pstx(cc,bb)
     &       +pstx(aa,dd)+pstx(cc,dd))
      enddo
      enddo
c
c Decouple Data
c
      do i=1,ix
      do j=1,jy
      do k=1,kz
        uu(i,j,k) = uu(i,j,k)/pstd(i,j)
        vv(i,j,k) = vv(i,j,k)/pstd(i,j)
      enddo
      enddo
      enddo
 
      do i=1,ix-1
      do j=1,jy-1
      do k=1,kz
        tt(i,j,k) = tt(i,j,k)/pstx(i,j)
        qv(i,j,k) = qv(i,j,k)/pstx(i,j)
        qc(i,j,k) = qc(i,j,k)/pstx(i,j)
        qr(i,j,k) = qr(i,j,k)/pstx(i,j)
        cice(i,j,k) = cice(i,j,k)/pstx(i,j)
        sice(i,j,k) = sice(i,j,k)/pstx(i,j)
        pp(i,j,k) = 0.01*pp(i,j,k)/pstx(i,j)
      enddo
      enddo
      enddo
 
      if((ifprs.or.ifsig).and.ifww)then
      do i=1,ix-1
      do j=1,jy-1
      do k=1,kz+1
        ww(i,j,k) = ww(i,j,k)/pstx(i,j)
      enddo
      enddo
      enddo
      endif
 
      do i=1,ix-1
      do j=1,jy-1
        pstx(i,j)=pstx(i,j)*10.0
      enddo
      enddo

      do i=1,ix
      do j=1,jy
        pstd(i,j)=pstd(i,j)*10.0
      enddo
      enddo
c
c  Calculate total pressure on sigma levels
c
      do k=1,kz
      maxp(k)=0.
      minp(k)=10000.
      do i=1,ix-1
      do j=1,jy-1
        prs(i,j,k)=pstx(i,j)*sih(k)+mrf(1,2)+pp(i,j,k)
        if(prs(i,j,k).gt.maxp(k))maxp(k)=prs(i,j,k)
        if(prs(i,j,k).lt.minp(k))minp(k)=prs(i,j,k)
      enddo
      enddo
      enddo
c
c  Find pressure on dot point domain for use in interpolation to pressure levs
c
      if(ifprs)then
      do k=1,kz
        maxpd(k)=0.
        minpd(k)=10000.
        call xtodot(prs(1,1,k),dotvar,ixl,jyl,ix,jy)
        do j=1,jy
        do i=1,ix
          prsd(i,j,k)=dotvar(i,j)
          if(prsd(i,j,k).gt.maxpd(k))maxpd(k)=prsd(i,j,k)
          if(prsd(i,j,k).lt.minpd(k))minpd(k)=prsd(i,j,k)
        enddo
        enddo
      enddo
      endif
c
c   Calculate Surface and sea-level Pressure
c
      do i=1,ix-1
      do j=1,jy-1
        sfp(i,j)=pstx(i,j)+mrf(1,2)+pp(i,j,kz)
        wmixrat=qv(i,j,kz)
        tvll=tt(i,j,kz) * (1.+0.608*wmixrat )
        beta=287.04/9.81*
     &       log((pstx(i,j)+mrf(1,2)+pp(i,j,kz))/
     &       (sih(kz)*pstx(i,j)+mrf(1,2)+pp(i,j,kz)))
        gpll=ter(i,j)+beta*tvll/(1.-.5*beta*0.0065)
        tvavg = tvll + 0.0065 * (gpll - .5 * ter(i,j))
        slp(i,j)=(pstx(i,j)+mrf(1,2)+pp(i,j,kz))*
     &           exp(9.81*ter(i,j)/287.04/tvavg)
      enddo
      enddo
c
c Calculate dewpoint
c
      do k=1,kz
      do j=1,jy-1
      do i=1,ix-1
        q=max(qv(i,j,k),1.e-15)
        t=tt(i,j,k)
        p=pstx(i,j)*sih(k)+mrf(1,2)+pp(i,j,k)
        e=q*p/(0.622+(1.-0.622)*q)
        dew(i,j,k)=5418.12/(19.84659-alog(e/6.112))
c     and the relative humidity
        tcels=t-273.15
        esat=6.112*exp((17.67*tcels)/(tcels+243.5))
        qsat=eps*esat/p 
        rh(i,j,k)=q/qsat*100
      enddo
      enddo
      enddo
      print*,'t',t,' p ',p,' q ',q,' qsat ',qsat
      print*,'eps ',eps,' esat ',esat
      PRINT*,'rh(5,5,5) = ',rh(5,5,5) 

CCCCCC
C Calculate Rain totals
      output_hr = NINT(mrf(1,6)/60.0)
      do j=1, jy-1
         do i=1, ix-1
            prec_3hr(i,j)=0.0
            re_3hr(i,j)=0.0
            rc_3hr(i,j)=0.0
         enddo
      enddo
c
c     unless this is the initialization time, compute the precip
c
      if (output_hr.eq.0 .or. (output_hr.eq.special_hour .and.
     $     mif(101,1).eq.dom_special) ) then
         print 1027,mif(101,1),output_hr,dom_special,special_hour
 1027    format('skipping precip calc.s for domain',i2
     $        ,' output hour = ',i3,'. dom_special = ',i2
     $        ' and special_hour = ',i3)
         continue
      else
         print*,'reading the last output from tempfile'
         read(22) last_output_hr
         if ((output_hr-last_output_hr) .ne. 3) then
            print*,'missing the last output file precip'
         else
            read(22) last_output_hr_tot_prec
            read(22) last_output_hr_exp_prec
            read(22) last_output_hr_con_prec
            do j=1,jy-1            
               do i=1,ix-1
                  prec_3hr(i,j)=tot_prec(i,j)-
     &                 last_output_hr_tot_prec(i,j)
                  re_3hr(i,j)=re(i,j)-
     &                 last_output_hr_exp_prec(i,j)
                  rc_3hr(i,j)=rc(i,j)-
     &                 last_output_hr_con_prec(i,j)
               enddo
            enddo
         endif
      endif
c   write the current total precip grid to the tempfile      
      rewind (22)
      write(22) output_hr
      write(22) tot_prec
      write(22) re
      write(22) rc

c
c  Calculate the height above ground level of lowest half sigma level
c
c     if(.not.ifsig.and.iout.eq.1)then
      if(.not.ifsig.and.itime.eq.1)then
        call ghtcalc(tt(1,1,kz-1),prs(1,1,kz-1),qv(1,1,kz-1),sfp,ter
     $        ,hlhsl,ix,jy,2,ixl,jyl,2,2)
        do j=1,jy-icross
        do i=1,ix-icross
          hlhsl(i,j,1)=hlhsl(i,j,2)-ter(i,j)
          hlhsltot=hlhsltot+hlhsl(i,j,2)-ter(i,j)
        enddo
        enddo
        zagllhsl=hlhsltot/((ix-1)*(jy-1))
        if(zagllhsl.ge.40.)then
          izagllhsl=10*nint(zagllhsl/10.)
        elseif(zagllhsl.ge.5.)then
          izagllhsl=5*nint(zagllhsl/5.)
        else
          izagllhsl=nint(zagllhsl)
        endif
      endif
c
c  Interpolate vertical velocity to half sigma levels
c
      if((ifprs.or.ifsig).and.ifww)then
      do k = 1,kz
      do j = 1,jy-icross
      do i = 1,ix-icross
        wwh(i,j,k)=(ww(i,j,k)+ww(i,j,k+1))/2.
      enddo
      enddo
      enddo
      endif
c
c   Dump out arrays for gempak
c
 678  continue
c  
c  3-d fields first:
c
      if(iprog.eq.2.or.iprog.eq.3)then  !datagrid or rawins output
        parm='UWND        '
        ivcor=1
        print*,'Writing U on pressure levels - ',uu(ix/2,jy/2,1),'M/S'
        do kk = 1,kz
          ilvl(1)=mif(101+kk,iprog)
          do ii=1,ix
          do jj=1,jy
            jjii=jj+((ii-1)*jy)
            dat(jjii)=uu(ii,jj,kk)
          enddo
          enddo
c          call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
c     $                 parm,rewrite,ipktyp,iprecise,iret)
C* GD_WPGD  ( IGDFLN, GRID, IGX, IGY, IGHDR, GDATTM, LEVEL, IVCORD,	*
C*            PARM, REWRIT, IPKTYP, NBITS, IRET )			*
C*									*
C* Input parameters:							*
C*	IGDFLN		INTEGER		Grid file number		*
C*	GRID (IGX,IGY)	REAL		Grid data			*
C*	IGX		INTEGER		Number of horizontal points	*
C*	IGY		INTEGER		Number of vertical points 	*
C*	IGHDR (IHDRSZ)	INTEGER		Grid header			*
C*	GDATTM (2)	CHAR*20		GEMPAK times			*
C*	LEVEL  (2)	INTEGER		Vertical levels			*
C*	IVCORD		INTEGER		Vertical coordinate		*
C*				  	   0 = NONE			*
C*				  	   1 = PRES			*
C*					   2 = THTA			*
C*					   3 = HGHT			*
C*	PARM		CHAR*12		Parameter name			*
C*	REWRIT		LOGICAL		Flag to replace existing grid	*
C*	IPKTYP		INTEGER		Packing type			*
C*	NBITS		INTEGER		Number of bits / precision	*
C*									*
C* Output parameters:							*
C*	IRET		INTEGER		Return code			*
C*					  0 = normal return		*
C*					 -4 = file not open		*
C*					 -5 = no write access		*
C*					 -6 = read/ write error		*
C*					 -9 = invalid grid size		*
C*					-10 = grid already exists	*
C*					-11 = grid file is full		*
          call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
          if(iret.ne.0) then
            print *,'gd_wpgd: iret=',iret,' Stopping . . .'
            goto 999
          endif
        enddo
        parm='VWND        '
        ivcor=1
        print*,'Writing V on pressure levels - ',vv(ix/2,jy/2,1),'M/S'
        do kk = 1,kz
          ilvl(1)=mif(101+kk,iprog)
          do ii=1,ix
          do jj=1,jy
            jjii=jj+((ii-1)*jy)
            dat(jjii)=vv(ii,jj,kk)
          enddo
          enddo
          call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
          if(iret.ne.0) then
            print *,'gd_wpgd: iret=',iret,' Stopping . . .'
            goto 999
          endif
        enddo

        parm='TMPK        '
        ivcor=1
        print*,'Writing T on pressure levels -', tt(ix/2,jy/2,1)
        do kk = 1,kz
          call xtodot(tt(1,1,kk),dotvar,ixl,jyl,ix,jy)
          ilvl(1)=mif(101+kk,iprog)
          do ii=1,ix
          do jj=1,jy
            jjii=jj+((ii-1)*jy)
            dat(jjii)=dotvar(ii,jj)
          enddo
          enddo
          call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
          if(iret.ne.0) then
            print*,'gd_wpgd: iret=',iret,' Stopping . . .'
            goto 999
          endif
        enddo

        parm='RELH        '
        ivcor=1
        print*,'Writing RH on pressure levels -', qv(ix/2,jy/2,1)
        do kk = 1,kz
          call xtodot(qv(1,1,kk),dotvar,ixl,jyl,ix,jy)
          ilvl(1)=mif(101+kk,iprog)
          do ii=1,ix
          do jj=1,jy
            jjii=jj+((ii-1)*jy)
            dat(jjii)=dotvar(ii,jj)
          enddo
          enddo
          call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
          if(iret.ne.0) then
            print*,'gd_wpgd: iret=',iret,' Stopping . . .'
            goto 999
          endif
        enddo

        parm='HGHT        '
        ivcor=1
        print*,'Writing geo. hght. on pressure levels - '
     $       ,pp(ix/2,jy/2,1),'M/S'
        do kk = 1,kz
          ilvl(1)=mif(101+kk,iprog)
          do ii=1,ix
          do jj=1,jy
            jjii=jj+((ii-1)*jy)
            dat(jjii)=pp(ii,jj,kk)
          enddo
          enddo
          call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
          if(iret.ne.0) then
            print *,'gd_wpgd: iret=',iret,' Stopping . . .'
            goto 999
          endif
        enddo

      goto 321 ! plot out surface fields, same as for mmout and interp
      endif    !end of stuff for datagrid output

      do 421 inum = 1,num_3d_outputs
c      if(ifuu)then 
         if(sig_fld(inum).eq.'UWND'.or.sig_fld(inum).eq.'UREL' .or.
     $        prs_fld(inum).eq.'UWND'.or.prs_fld(inum).eq.'UREL')then 
            if(ifprs)then
               parm=prs_fld(inum)//'        '
               ivcor=1
               call pintp (uu,pdata,prsd,tvirt,plevs,maxpd,minpd,ix,jy
     $              ,kz,ixl,jyl,kzl,kp,kzp,idot,0)
               print*,'Writing U on pressure levels - '
     $              ,pdata(ix/2,jy/2,1),'M/S'
               do kk = 1,kp
                  ilvl(1)=int(plevs(kk))
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=pdata(ii,jj,kk)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print *,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(ifsig)then
               print*,'Writing U on sigma levels -',uu(ix/2,jy/2,kz)
     $              ,'M/S'
               parm=sig_fld(inum)//'        '
               ivcor=4
               do kk = 1,kz
                  ilvl(1) = nint(mrf(101+kk,iprog)*10000.0)
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=uu(ii,jj,kk)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(.not.ifsig.and.iprog.ne.2)then
               print*,'Writing U on ',izagllhsl,' meter level -'
     $              ,uu(ix/2,jy/2,kz),'M/S'
               parm='UWND        '
               ivcor=1514227532
               kk = kz
               ilvl(1) = izagllhsl
               do ii=1,ix
                  do jj=1,jy
                     jjii=jj+((ii-1)*jy)
                     dat(jjii)=uu(ii,jj,kk)
                  enddo
               enddo
               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $              parm,rewrite,ipktyp,iprecise,iret)
               if(iret.ne.0) then
                  print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                  goto 999
               endif
            endif
         endif


         if(sig_fld(inum).eq.'VWND'.or.sig_fld(inum).eq.'VREL' .or.
     $        prs_fld(inum).eq.'VWND'.or.prs_fld(inum).eq.'VREL')then 
c      if(ifvv)then
            if(ifprs)then
               parm=prs_fld(inum)//'        '
               ivcor=1
               call pintp (vv,pdata,prsd,tvirt,plevs,maxpd,minpd
     $              ,ix,jy,kz,ixl,jyl,kzl,kp,kzp,idot,0)
               print*,'Writing V on pressure levels -'
     $              ,pdata(ix/2,jy/2,1),'M/S'
               do kk = 1,kp
                  ilvl(1)=int(plevs(kk))
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=pdata(ii,jj,kk)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(ifsig)then
               print *, 'Writing V on sigma levels - '
     $              ,vv(ix/2,jy/2,kz),'M/S'
               parm=sig_fld(inum)//'        '
               ivcor=4
               do kk = 1,kz
                  ilvl(1) = nint(mrf(101+kk,iprog)*10000.0)
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=vv(ii,jj,kk)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(.not.ifsig)then
               print*,'Writing V on ',izagllhsl,' meter level -'
     $              ,vv(ix/2,jy/2,kz),'M/S'
               parm='VWND        '
               ivcor=1514227532
               kk = kz
               ilvl(1) = izagllhsl
               do ii=1,ix
                  do jj=1,jy
                     jjii=jj+((ii-1)*jy)
                     dat(jjii)=vv(ii,jj,kk)
                  enddo
               enddo
               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $              parm,rewrite,ipktyp,iprecise,iret)
               if(iret.ne.0) then
                  print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                  goto 999
               endif
            endif
         endif


c      if(ifww)then
         if(sig_fld(inum).eq.'WWND'.or.prs_fld(inum).eq.'WWND')then 
            if(ifprs)then
               parm=prs_fld(inum)//'        '
               ivcor=1
               call pintp (wwh,pdata,prsd,tvirt,plevs,maxpd,minpd
     $              ,ix,jy,kz,ixl,jyl,kzl,kp,kzp,icross,0)
               print*,'Writing W on pressure levels -'
     $              ,pdata(ix/2,jy/2,1),'M/S'
               do kk = 1,kp
                  call xtodot(pdata(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1)=int(plevs(kk))
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=dotvar(ii,jj)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(ifsig)then
               print *, 'Writing W on sigma levels -'
     $              ,wwh(ix/2,jy/2,kz),'M/S'
               parm=sig_fld(inum)//'        '
               ivcor=4
               do kk = 1,kz
                  call xtodot(wwh(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1) = nint(mrf(101+kk,iprog)*10000.0)
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=dotvar(ii,jj)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
     $                 ivcor,parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo 
            endif
         endif

 
c      if(iftt)then
         if(sig_fld(inum).eq.'TMPK'.or.prs_fld(inum).eq.'TMPK')then 
            if(ifprs)then
               parm=prs_fld(inum)//'        '
               ivcor=1
               call pintp (tt,ttp,prs,tvirt,plevs,maxp,minp
     $              ,ix,jy,kz,ixl,jyl,kzl,kp,kzp,icross,3)
               print*,'Writing T on pressure levels -', ttp(ix/2,jy/2,1)
               do kk = 1,kp
                  call xtodot(ttp(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1)=int(plevs(kk))
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=dotvar(ii,jj)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(ifsig)then
               print *, 'Writing T on sigma levels - ', tt(ix/2,jy/2,kz)
               parm=sig_fld(inum)//'        '
               ivcor=4
               do kk = 1,kz
                  call xtodot(tt(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1) = nint(mrf(101+kk,iprog)*10000.0)
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=dotvar(ii,jj)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(.not.ifsig)then
               print*,'Writing T on ',izagllhsl,' meter level -'
     $              ,tt(ix/2,jy/2,kz)
               parm='TMPK        '
               ivcor=1514227532
               kk = kz
               call xtodot(tt(1,1,kk),dotvar,ixl,jyl,ix,jy)
               ilvl(1) = izagllhsl
               do ii=1,ix
                  do jj=1,jy
                     jjii=jj+((ii-1)*jy)
                     dat(jjii)=dotvar(ii,jj)
                  enddo
               enddo
               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $              parm,rewrite,ipktyp,iprecise,iret)
               if(iret.ne.0) then
                  print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                  goto 999
               endif
            endif
         endif


         if(sig_fld(inum).eq.'DWPK'.or.prs_fld(inum).eq.'DWPK')then 
c     if(ifqv)then
            if(ifprs)then
               parm=prs_fld(inum)//'        '
               ivcor=1
               call pintp (dew,pdata,prs,tvirt,plevs,maxp,minp
     $              ,ix,jy,kz,ixl,jyl,kzl,kp,kzp,icross,3)
               print*,'Writing Dewpoint on pressure levels -'
     $              , pdata(ix/2,jy/2,1)
               do kk = 1,kp
                  call xtodot(pdata(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1)=int(plevs(kk))
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=dotvar(ii,jj)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(ifsig)then
               print *, 'Writing Dewpoint on sigma levels - '
     $              , dew(ix/2,jy/2,kz)
               parm=sig_fld(inum)//'        '
               ivcor=4
               do kk = 1,kz
                  call xtodot(dew(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1) = nint(mrf(101+kk,iprog)*10000.0)
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=dotvar(ii,jj)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(.not.ifsig)then
               print*,'Writing Dewpoint on ',izagllhsl,' meter level -',
     &              dew(ix/2,jy/2,kz),'M/S'
               parm='DWPK        '
               ivcor=1514227532
               kk = kz
               call xtodot(dew(1,1,kk),dotvar,ixl,jyl,ix,jy)
               ilvl(1) = izagllhsl
               do ii=1,ix
                  do jj=1,jy
                     jjii=jj+((ii-1)*jy)
                     dat(jjii)=dotvar(ii,jj)
                  enddo
               enddo
               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $              parm,rewrite,ipktyp,iprecise,iret)
               if(iret.ne.0) then
                  print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                  goto 999
               endif
            endif
         endif


         if(sig_fld(inum).eq.'PRES')then 
c     if(ifpp)then
            if(ifsig)then
               print *, 'Writing Pressure on sigma levels - '
     $              , prs(ix/2,jy/2,kz)
               parm=sig_fld(inum)//'        '
               ivcor=4
               do kk = 1,kz
                  call xtodot(prs(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1) = nint(mrf(101+kk,iprog)*10000.0)
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=dotvar(ii,jj)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
         endif
         if(prs_fld(inum).eq.'HGHT')then 
            if(ifprs)then
               print*, 'Calculating Geopotential Height . . .'
               if(ighttype.eq.2)then
                  call ghtcalc(tt,prs,qv,sfp,ter,tmp3ds
     $                 ,ix,jy,kz,ixl,jyl,kzl,2)
c     Find virtual temperature on lowest sigma level:
                  do j=1,jy-1
                     do i=1,ix-1
                        tvirt(i,j)=tt(i,j,kz)*(1.+qv(i,j,kz)/.62197)
     $                       /(1.+qv(i,j,kz))
                     enddo
                  enddo
                  call pintp(tmp3ds,pdata,prs,tvirt,plevs,maxp,minp
     $                 ,ix,jy,kz,ixl,jyl,kzl,kp,kzp,icross,4)
               elseif(ighttype.eq.1)then
                  do k=1,kp
                     do j=1,jy-1
                        do i=1,ix-1
                           prsp(i,j,k)=plevs(k)
                        enddo
                     enddo
                  enddo
c     :skip following loop to do interpolation of q rather than log q
                  do k=1,kz
                     do j=1,jy-1
                        do i=1,ix-1
                           qv(i,j,k)=log(max(qv(i,j,k),1.e-15))
                        enddo
                     enddo
                  enddo
                  call pintp(qv,tmp3dp,prs,tvirt,plevs,maxp,minp
     $                 ,ix,jy,kz,ixl,jyl,kzl,kp,kzp,icross,1)
c     :skip following loop to do interpolation of q rather than log q
                  do k=1,kp
                     do j=1,jy-1
                        do i=1,ix-1
                           tmp3dp(i,j,k)=exp(tmp3dp(i,j,k))
                        enddo
                     enddo
                  enddo
                  call ghtcalc(ttp,prsp,tmp3dp,sfp,ter,pdata
     $                 ,ix,jy,kp,ixl,jyl,kzp,1)
               endif
               print *, 'Writing Geopotential Height - '
     $              , pdata(ix/2,jy/2,1)
               parm='HGHT        '
               ivcor=1
               do kk = 1,kp
                  call xtodot(pdata(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1) = int(plevs(kk))
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=dotvar(ii,jj)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif               !done with Geopotential height
         endif

         if(sig_fld(inum).eq.'QCLD'.or.prs_fld(inum).eq.'QCLD')then 
c     if(ifqc)then
            if(ifprs)then
               parm=prs_fld(inum)//'        '
               ivcor=1
               call pintp (qc,pdata,prs,tvirt,plevs,maxp,minp
     $              ,ix,jy,kz,ixl,jyl,kzl,kp,kzp,icross,0)
               print*,'Writing cloud water mixing ratio(g/kg) '
     $              ,'on pressure levels -',
     &              pdata(ix/2,jy/2,1)*1000.
               do kk = 1,kp
                  call xtodot2(pdata(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1)=int(plevs(kk))
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=1000.*dotvar(ii,jj) !also change to g/kg
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(ifsig)then
               print*,'Writing cloud water mixing ratio(g/kg)'
     $              ,' on sigma levels -',
     &              qc(ix/2,jy/2,kz)*1000.
               parm=sig_fld(inum)//'        '
               ivcor=4
               do kk = 1,kz
                  call xtodot2(qc(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1) = nint(mrf(101+kk,iprog)*10000.0)
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=1000.*dotvar(ii,jj) !also change to g/kg
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
         endif

         if(sig_fld(inum).eq.'MIXR'.or.prs_fld(inum).eq.'MIXR')then 
            if(ifprs)then
               parm=prs_fld(inum)//'        '
               ivcor=1
               call pintp (qv,pdata,prs,tvirt,plevs,maxp,minp
     $              ,ix,jy,kz,ixl,jyl,kzl,kp,kzp,icross,0)
               print*,'Writing total mixing ratio(g/kg) '
     $              ,'on pressure levels -',
     &              pdata(ix/2,jy/2,1)*1000.
               do kk = 1,kp
                  call xtodot2(pdata(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1)=int(plevs(kk))
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=1000.*dotvar(ii,jj) !also change to g/kg
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(ifsig)then
               print*,'Writing total mixing ratio(g/kg)'
     $              ,' on sigma levels -',
     &              qv(ix/2,jy/2,kz)*1000.
               parm=sig_fld(inum)//'        '
               ivcor=4
               do kk = 1,kz
                  call xtodot2(qv(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1) = nint(mrf(101+kk,iprog)*10000.0)
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=1000.*dotvar(ii,jj) !also change to g/kg
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
         endif

         if(sig_fld(inum).eq.'QRAI'.or.prs_fld(inum).eq.'QRAI')then 
c     if(ifqr)then
            if(ifprs)then
               parm=prs_fld(inum)//'        '
               ivcor=1
               call pintp (qr,pdata,prs,tvirt,plevs,maxp,minp
     $              ,ix,jy,kz,ixl,jyl,kzl,kp,kzp,icross,0)
               print*,'Writing rain water mixing ratio(g/kg)'
     $              ,' on pressure levels -',
     &              pdata(ix/2,jy/2,1)*1000.
               do kk = 1,kp
                  call xtodot2(pdata(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1)=int(plevs(kk))
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=1000.*dotvar(ii,jj) !also change to g/kg
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(ifsig)then
               print*,'Writing rain water mixing ratio(g/kg)'
     $              ,' on sigma levels -',
     &              qr(ix/2,jy/2,kz)*1000.
               parm=sig_fld(inum)//'        '
               ivcor=4
               do kk = 1,kz
                  call xtodot2(qr(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1) = nint(mrf(101+kk,iprog)*10000.0)
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=1000.*dotvar(ii,jj)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(.not.ifsig)then
               print*,'Writing Rainwater mixing ratio(g/kg) on ',izagllhsl,
     &              ' meter level -',qr(ix/2,jy/2,kz)*1000.
               parm='QRAI        '
               ivcor=1514227532
               kk = kz
               call xtodot2(qr(1,1,kk),dotvar,ixl,jyl,ix,jy)
               ilvl(1) = izagllhsl
               do ii=1,ix
                  do jj=1,jy
                     jjii=jj+((ii-1)*jy)
                     dat(jjii)=1000.*dotvar(ii,jj)
                  enddo
               enddo
               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $              parm,rewrite,ipktyp,iprecise,iret)
               if(iret.ne.0) then
                  print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                  goto 999
               endif
            endif
         endif

         if(sig_fld(inum).eq.'QICE'.or.prs_fld(inum).eq.'QICE')then 
c     if(ifcice)then
            print *, 'Writing ice water mixing ratio(g/kg) - ', 
     &           cice(ix/2,jy/2,kz)*1000.
            if(ifprs)then
               parm=prs_fld(inum)//'        '
               ivcor=1
               call pintp (cice,pdata,prs,tvirt,plevs,maxp,minp
     $              ,ix,jy,kz,ixl,jyl,kzl,kp,kzp,icross,0)
               do kk = 1,kp
                  call xtodot2(pdata(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1)=int(plevs(kk))
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=1000.*dotvar(ii,jj) !also change to g/kg
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(ifsig)then
               parm=sig_fld(inum)//'        '
               ivcor=4
               do kk = 1,kz
                  call xtodot2(cice(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1) = nint(mrf(101+kk,iprog)*10000.0)
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=1000.*dotvar(ii,jj) !also change to g/kg
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
         endif
         
         if(sig_fld(inum).eq.'QSNO'.or.prs_fld(inum).eq.'QSNO')then 
            if(ifsice)then
               print *, 'Writing snow water mixing ratio(g/kg) - ', 
     &              sice(ix/2,jy/2,kz)*1000.
               if(ifprs)then
                  parm=prs_fld(inum)//'        '
                  ivcor=1
                  call pintp (sice,pdata,prs,tvirt,plevs,maxp,minp
     $                 ,ix,jy,kz,ixl,jyl,kzl,kp,kzp,icross,0)
                  do kk = 1,kp
                     call xtodot2(pdata(1,1,kk),dotvar,ixl,jyl,ix,jy)
                     ilvl(1)=int(plevs(kk))
                     do ii=1,ix
                        do jj=1,jy
                           jjii=jj+((ii-1)*jy)
                           dat(jjii)=1000.*dotvar(ii,jj) !also change to g/kg
                        enddo
                     enddo
                     call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl
     $                    ,ivcor,parm,rewrite,ipktyp,iprecise,iret)
                     if(iret.ne.0) then
                        print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                        goto 999
                     endif
                  enddo
               endif
               if(ifsig)then
                  parm=sig_fld(inum)//'        '
                  ivcor=4
                  do kk = 1,kz
                     call xtodot2(sice(1,1,kk),dotvar,ixl,jyl,ix,jy)
                     ilvl(1) = nint(mrf(101+kk,iprog)*10000.0)
                     do ii=1,ix
                        do jj=1,jy
                           jjii=jj+((ii-1)*jy)
                           dat(jjii)=1000.*dotvar(ii,jj)
                        enddo
                     enddo
                     call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                    parm,rewrite,ipktyp,iprecise,iret)
                     if(iret.ne.0) then
                        print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                        goto 999
                     endif
                  enddo
               endif
            endif
         endif

         if(sig_fld(inum).eq.'RTEN'.or.prs_fld(inum).eq.'RTEN')then 
c     if(ifrtten)then
            print *, 'Writing Radiative Tendency(K/sec) - '
     $           , rtten(ix/2,jy/2,kz)
            if(ifprs)then
               parm=prs_fld(inum)//'        '
               ivcor=1
               call pintp (rtten,pdata,prs,tvirt,plevs,maxp,minp
     $              ,ix,jy,kz,ixl,jyl,kzl,kp,kzp,icross,0)
               do kk = 1,kp
                  call xtodot2(pdata(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1)=int(plevs(kk))
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=dotvar(ii,jj)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
            if(ifsig)then
               parm=sig_fld(inum)//'        '
               ivcor=4
               do kk = 1,kz
                  call xtodot2(rtten(1,1,kk),dotvar,ixl,jyl,ix,jy)
                  ilvl(1) = nint(mrf(101+kk,iprog)*10000.0)
                  do ii=1,ix
                     do jj=1,jy
                        jjii=jj+((ii-1)*jy)
                        dat(jjii)=dotvar(ii,jj)
                     enddo
                  enddo
                  call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,ivcor,
     $                 parm,rewrite,ipktyp,iprecise,iret)
                  if(iret.ne.0) then
                     print*,'gd_wpgd: iret=',iret,' Stopping . . .'
                     goto 999
                  endif
               enddo
            endif
         endif
 421  continue                  ! end of loop for 3-d output fields
c
c  Now do surface fields
c
 321  ivcor=0
      ilvl(1) = 0
      
      if(iprog.eq.2.or.iprog.eq.3)then !datagrid/rawins fields here
         print *, 'Writing sea level Pressure'
         parm='PMSL        '
         call xtodot(slp,dotvar,ixl,jyl,ix,jy)
         do ii=1,ix
            do jj=1,jy
               jjii=jj+((ii-1)*jy)
               dat(jjii)=dotvar(ii,jj)
            enddo
         enddo
         call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
     &        ivcor,parm,rewrite,ipktyp,iprecise,iret)
c     
c     add sea-surface temperature,and snow-cover here in future.
c     
      endif                     !end of datagrid 2-d fields
      
      do 521 inum = num_3d_outputs+1,num_3d_outputs+num_2d_outputs
         if(iprog.gt.4)then     !mm5 or interp output
            if(sig_fld(inum).eq.'PRES'.or.prs_fld(inum).eq.'PRES')then 
               print *, 'Writing surface Pressure'
               parm='PRES        '
               call xtodot(sfp,dotvar,ixl,jyl,ix,jy)
               do ii=1,ix
                  do jj=1,jy
                     jjii=jj+((ii-1)*jy)
                     dat(jjii)=dotvar(ii,jj)
                  enddo
               enddo
               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
     $              ivcor,parm,rewrite,ipktyp,iprecise,iret)
            endif
            
            if(sig_fld(inum).eq.'PMSL'.or.prs_fld(inum).eq.'PMSL')then 
               print *, 'Writing sea level Pressure'
               parm='PMSL        '
               call xtodot(slp,dotvar,ixl,jyl,ix,jy)
               do ii=1,ix
                  do jj=1,jy
                     jjii=jj+((ii-1)*jy)
                     dat(jjii)=dotvar(ii,jj)
                  enddo
               enddo
               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
     $              ivcor,parm,rewrite,ipktyp,iprecise,iret)
            endif

            if(sig_fld(inum).eq.'TMPK'.or.prs_fld(inum).eq.'TMPK'.or.
     $           sig_fld(inum).eq.'TGRD'.or.prs_fld(inum).eq.'TGRD')
     $           then 
               print *, 'Writing Ground Temp as though at sigma = 1- '
     $              , tg(ix/2,jy/2)
c     parm='TGRD        '
               parm='TMPK        '
               call xtodot(tg,dotvar,ixl,jyl,ix,jy)
               do ii=1,ix
                  do jj=1,jy
                     jjii=jj+((ii-1)*jy)
                     dat(jjii)=dotvar(ii,jj)
                  enddo
               enddo
c               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,10000,
c     $              4,parm,rewrite,ipktyp,iprecise,iret)
               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
     $              ivcor,parm,rewrite,ipktyp,iprecise,iret)
            endif

c     print *, 'Writing lowest half sig level temp as sfc temp'
c     parm='TMPK        '
c     call xtodot(tt(1,1,kz),dotvar,ixl,jyl,ix,jy)
c     do ii=1,ix
c     do jj=1,jy
c     jjii=jj+((ii-1)*jy)
c     dat(jjii)=dotvar(ii,jj)
c     enddo
c     enddo
c     call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
c     $             ivcor,parm,rewrite,ipktyp,iprecise,iret)

c     print *, 'Writing lowest half sig level dwpk as sfc dwpk'
c     parm='DWPK        '
c     call xtodot(dew(1,1,kz),dotvar,ixl,jyl,ix,jy)
c     do ii=1,ix
c     do jj=1,jy
c     jjii=jj+((ii-1)*jy)
c     dat(jjii)=dotvar(ii,jj)
c     enddo
c     enddo
c     call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
c     $             ivcor,parm,rewrite,ipktyp,iprecise,iret)

            if(sig_fld(inum).eq.'P03M'.or.prs_fld(inum).eq.'P03M')then 
               parm='P03M        '
               call xtodot2(prec_3hr(1,1),dotvar,ixl,jyl,ix,jy)
               print*,'Writing accumulated 3hr Ppt. -'
     $              ,dotvar(ix/2,jy/2),'mm'
               do ii=1,ix
                  do jj=1,jy
                     jjii=jj+((ii-1)*jy)
                     dat(jjii)=dotvar(ii,jj)
                  enddo
               enddo
               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
     $              ivcor,parm,rewrite,ipktyp,iprecise,iret)
            endif

            if(sig_fld(inum).eq.'E03M'.or.prs_fld(inum).eq.'E03M')then 
               parm='E03M        '
               call xtodot2(re_3hr(1,1),dotvar,ixl,jyl,ix,jy)
               print*,'Writing 3hr explicit Ppt. -'
     $              ,dotvar(ix/2,jy/2),'mm'
               do ii=1,ix
                  do jj=1,jy
                     jjii=jj+((ii-1)*jy)
                     dat(jjii)=dotvar(ii,jj)
                  enddo
               enddo
               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
     $              ivcor,parm,rewrite,ipktyp,iprecise,iret)
            endif

            if(sig_fld(inum).eq.'C03M'.or.prs_fld(inum).eq.'C03M')then 
               parm='C03M        '
               call xtodot2(rc_3hr(1,1),dotvar,ixl,jyl,ix,jy)
               print*,'Writing 3hr convective Ppt. -'
     $              ,dotvar(ix/2,jy/2),'mm'
               do ii=1,ix
                  do jj=1,jy
                     jjii=jj+((ii-1)*jy)
                     dat(jjii)=dotvar(ii,jj)
                  enddo
               enddo
               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
     $              ivcor,parm,rewrite,ipktyp,iprecise,iret)
            endif

            if(sig_fld(inum).eq.'HGHT'.or.prs_fld(inum).eq.'HGHT')then 
               print *, 'Writing Terrain - ', terd(ix/2,jy/2)
               parm='HGHT        '
               do ii=1,ix
                  do jj=1,jy
                     jjii=jj+((ii-1)*jy)
                     dat(jjii)=terd(ii,jj)
                  enddo
               enddo
               call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
     $              ivcor,parm,rewrite,ipktyp,iprecise,iret)
            endif

c     print*,'Writing water use -', dluse(ix/2,jy/2)
c     parm='WATR        '
c     do ii=1,ix
c     do jj=1,jy
c     jjii=jj+((ii-1)*jy)
c     dat(jjii)=dluse(ii,jj)
c     enddo
c     enddo
c     call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
c     $             ivcor,parm,rewrite,ipktyp,iprecise,iret)

c     -----------------------------------------------------
c     This section is unused in U.W. implementation
C     if (cmn(itime).gt.0) then
C     if(ifcon)then
C     parm='CXXM        '
C     call xtodot2(rc(1,1),dotvar,ixl,jyl,ix,jy)
C     print*,'Writing accumulated convective Ppt. -'
C     $       ,dotvar(ix/2,jy/2),'mm'
C     do ii=1,ix
C     do jj=1,jy
C     jjii=jj+((ii-1)*jy)
C     dat(jjii)=10.*dotvar(ii,jj)
C     enddo
C     enddo
C     call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
C     $               ivcor,parm,rewrite,ipktyp,iprecise,iret)
C     parm='EXXM        '
C     call xtodot2(re(1,1),dotvar,ixl,jyl,ix,jy)
C     print*,'Writing accumulated explicit Ppt. -'
C     $       ,dotvar(ix/2,jy/2),'mm'
C     do ii=1,ix
C     do jj=1,jy
C     jjii=jj+((ii-1)*jy)
C     dat(jjii)=10.*dotvar(ii,jj)
C     enddo
C     enddo
C     call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
C     $               ivcor,parm,rewrite,ipktyp,iprecise,iret)
C     endif
C     
C     if(ifcon.or.ifres)then
C     parm='PXXM        '
C     do jj=1,jy
C     do ii=1,ix
C     tmp2d(ii,jj)=re(ii,jj)+rc(ii,jj)
C     enddo
C     enddo
C     call xtodot2(tmp2d,dotvar,ixl,jyl,ix,jy)
C     print*,'Writing accumulated total Ppt. -'
C     $       ,dotvar(ix/2,jy/2),'mm'
C     do ii=1,ix
C     do jj=1,jy
C     jjii=jj+((ii-1)*jy)
C     dat(jjii)=10.*dotvar(ii,jj)
C     enddo
C     enddo
C     call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
C     $               ivcor,parm,rewrite,ipktyp,iprecise,iret)
C     
C     do iptime=1,20     !loop through up to 20 precip intervals
C     iph=iptimes(iptime)
C     ipm=iph*60
C     ipflag=0
C     if(ipm.ne.0)then
C     if(mod(ipm,itapfrq).eq.0.and.mod(cmn(itime),ipm).eq.0.and.
C     &        (itime-ipm/itapfrq).ge.1) ipflag=1
C     else
C     if(itime.gt.1.and.ipm.eq.0.and.mod(itapfrq,60).ne.0) ipflag=1
C     endif
C     if(ipflag.eq.1)then
C     
C     c      write x-hourly precip:
C     c
C     c   Total precip
C     c
C     write(parm,fmt="('P',i2.2,'M        ')")iph
C     do jj=1,jy
C     do ii=1,ix
C     c original              itime2=itime-ipm/itapfrq
C     itime2=maxtimes
C     if(ipm.eq.0)itime2=itime-1
C     tmp2d(ii,jj)= re(ii,jj)+ rc(ii,jj)
C     &             -re(ii,jj,itime2)-rc(ii,jj,itime2)
C     enddo
C     enddo
C     call xtodot2(tmp2d,dotvar,ixl,jyl,ix,jy)
C     print*,'Writing',ipm/60,' hour total precip(',
C     &             parm(1:4),') -',dotvar(ix/2,jy/2),'mm'
C     do ii=1,ix
C     do jj=1,jy
C     jjii=jj+((ii-1)*jy)
C     dat(jjii)=10.*dotvar(ii,jj)
C     enddo
C     enddo
C     call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
C     $                   ivcor,parm,rewrite,ipktyp,iprecise,iret)
C     c
C     c    Convective precip
C     c
C     if(ifcon)then
C     write(parm,fmt="('C',i2.2,'M        ')")iph
C     do jj=1,jy
C     do ii=1,ix
C     itime2=itime-ipm/itapfrq
C     if(ipm.eq.0)itime2=itime-1
C     tmp2d(ii,jj)= rc(ii,jj) - rc(ii,jj,itime2)
C     enddo
C     enddo
C     call xtodot2(tmp2d,dotvar,ixl,jyl,ix,jy)
C     print*,'Writing ',ipm/60,' hour convective precip(',
C     &             parm(1:4),') -',dotvar(ix/2,jy/2),'mm'
C     do ii=1,ix
C     do jj=1,jy
C     jjii=jj+((ii-1)*jy)
C     dat(jjii)=10.*dotvar(ii,jj)
C     enddo
C     enddo
C     call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
C     $                   ivcor,parm,rewrite,ipktyp,iprecise,iret)
C     c
C     c    Explicit precip
C     c
C     write(parm,fmt="('E',i2.2,'M        ')")iph
C     do jj=1,jy
C     do ii=1,ix
C     itime2=itime-ipm/itapfrq
C     if(ipm.eq.0)itime2=itime-1
C     tmp2d(ii,jj)= re(ii,jj) - re(ii,jj,itime2)
C     enddo
C     enddo
C     call xtodot2(tmp2d,dotvar,ixl,jyl,ix,jy)
C     print*,'Writing ',ipm/60,' hour explicit precip(',
C     &             parm(1:4),') -',dotvar(ix/2,jy/2),'mm'
C     do ii=1,ix
C     do jj=1,jy
C     jjii=jj+((ii-1)*jy)
C     dat(jjii)=10.*dotvar(ii,jj)
C     enddo
C     enddo
C     call gd_wpgd(igdfln,dat,jy,ix,ighdr,gdattm,ilvl,
C     $                   ivcor,parm,rewrite,ipktyp,iprecise,iret)
C     endif
C     endif
C     enddo
C     endif
c     endif
c     -----------------------------------------------------
         endif
 521  continue

c
c  Return for next time
      goto 30
999   print *,'End'
      call gd_clos(igdfln,iret)
c     call gendp(1,iret)
c     if  ( iret .ne. 0 )  call er_wmsg ( 'GEMPLT', iret, ' ', ier )
      call ip_exit(iret)
      if(iret.ne.0)print*,'ip_exit:  iret=',iret
      stop 'Conversion finished'
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ghtcalc(temp,pres,q,psfc,hsfc,hght,ix,jx,kx,
     &                   ixl,jxl,kxl,istart)
C This subroutine builds a geopotential height field.  Above model
C surface elevation, heights are integrated up from surface elevation.
C Below model surface elevation, heights are integrated down from surface.
C Input:
C    temp: temperature array
C    pres: pressure levels
C    q   : mixing ratio array
C    psfc: surface pressure array
C    hsfc: surface elevation array
C    ix,jx,kx: dimensions of arrays actually containing data
C    ixl,jxl,kxl: fortran dimensions of arrays
C    istart: determines the vertical ordering of arrays.
C                 istart=2 if input levels start at highest and go to lowest
C                 istart=1 if input levels start at lowest and go to highest
C                 (MM5 sigma coordinate uses 2, conv pressure coord uses 1
C Output:
C    hght: geopotential height array
C Assumptions:
C    At least two levels completely above ground
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real temp(ixl,jxl,kxl),pres(ixl,jxl,kxl),q(ixl,jxl,kxl)
      real psfc(ixl,jxl)    ,hsfc(ixl,jxl),    hght(ixl,jxl,kxl)
      logical havelower,haveupper
      eps = .622
      R=287.04
      g=9.80616
      icross=1   !change to zero if using dot point data
      if(istart.eq.1)then
        k1=   1
        k2=  kx
        k3=   1
        k4=kx-2
        k5=   1
        k6=  -1
        j1=   1
        j2=   2
      elseif(istart.eq.2)then
        k1=  kx
        k2=   1
        k3=  -1
        k4=   3
        k5=  kx
        k6=   1
        j1=  -1
        j2=  -2
      endif
 
      do 10 j=1,jx-icross
         do 20 i=1,ix-icross
            haveupper=.false.
            havelower=.false.
c
c  Deal with above ground points first
c
            do k=k1,k2,k3
               if(pres(i,j,k).le.psfc(i,j))then !above ground - - build ght up
                  if(.not.havelower)then !extrapolate virtual temperature to surface
                     tvl= temp(i,j,k+j1)*(1.+q(i,j,k+j1)/eps)
     $                    /(1.+q(i,j,k+j1)) +
     &                    (log(psfc(i,j)) - log(pres(i,j,k+j1))) *
     &                    (temp(i,j,k)*(1.+q(i,j,k)/eps)/
     $                    (1.+q(i,j,k)) -
     &                    temp(i,j,k+j1)*(1.+q(i,j,k+j1)/eps)
     $                    /(1.+q(i,j,k+j1)))/
     &                    (log(pres(i,j,k)) - log(pres(i,j,k+j1)))
                     pl=psfc(i,j)
                     hl=hsfc(i,j)
                  else
                     tvl=tvh
                     pl=ph
                     hl=hght(i,j,k-j1)
                  endif
                  tvh=temp(i,j,k)*(1.+q(i,j,k)/eps)/(1.+q(i,j,k))
                  ph =pres(i,j,k)
                  tv=.5*(tvh+tvl)
                  dlnp=log(ph) - log(pl)
                  dz=(-R*tv/g)*dlnp
                  hght(i,j,k)=hl+dz
                  havelower=.true.
               endif
            enddo
c
c  Now deal with below ground points
c
            do k=k4,k5,k6
               if(pres(i,j,k).gt.psfc(i,j))then !below ground - - build ght down
                  if(.not.haveupper)then !extrapolate virtual temperature to surface
                     tvh= temp(i,j,k+j2)*(1.+q(i,j,k+j2)/eps)
     $                    /(1.+q(i,j,k+j2)) +
     &                    (log(psfc(i,j)) -  log(pres(i,j,k+j2))) *
     &                    (temp(i,j,k+j1)*(1.+q(i,j,k+j1)/eps)/
     $                    (1.+q(i,j,k+j1)) - temp(i,j,k+j2)*
     $                    (1.+q(i,j,k+j2)/eps)/(1.+q(i,j,k+j2)))/
     &                    (log(pres(i,j,k+j1)) -  log(pres(i,j,k+j2)))
                     ph=psfc(i,j)
                     hh=hsfc(i,j)
                  else
                     tvh=tvl
                     ph=pl
                     hh=hght(i,j,k+j1)
                  endif
                  tvl=temp(i,j,k)*(1.+q(i,j,k)/eps)/(1.+q(i,j,k))
                  pl =pres(i,j,k)
                  tv=.5*(tvh+tvl)
                  dlnp=log(ph) - log(pl)
                  dz=(-R*tv/g)*dlnp
                  hght(i,j,k)=hh-dz
                  haveupper=.true.
               endif
            enddo
 20      continue
 10   continue
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine pintp (array,parray,pres,tvirt,plevs,maxp,minp,ix,jx,
     &                   kxs,ixl,jxl,kxl,kxp,kxpl,icross,iextrap)
C interpolation algorithm from sigma to pressure levels.
C Below the lowest level (i.e. below ground) pressure level data is
C extrapolated using 1 of 6 options.  Otherwise, lowest sigma 
C level value used below ground.
C Above ptop, values are extrapolated since the highest half sigma level
C is slightly below ptop.
c  array:       Array on sigma levels to be interpolated(lowest sigma value 
c                    has third index of 1)
c  parray:      The array after interpolation
c  pres:        Array of pressure on sigma levels
c  tvirt:       Virtual temperature array on lowest sigma level.  Only used
c               if iextrap = 4
c  plevs:       Pressure levels to interpolate to (lowest pres. at position 1)
c  ix,jx,kxs:   The dimensions of the part of sigma array containing data
c  ixl,jxl,kxl: The fortran dimensions of the sigma array
c  kxp:         The number of pressure levels to output
c  kxpl:        The fortran dimensions of the pressure level array
c  iextrap:     Flag to determine type of extrapolation below surface:
c               0: Use lowest sigma values below lowest sigma level
c               1: Use linear log(pressure) extrapolation based on model lapse
c               2: Use linear pressure extrapolation based on model lapse rate
c               3: extrapolation based on standard lapse rate ( 6.5 K/km )
c               4: Extrapolate geopotential heights below ground level using
c                standard environmental lapse rate. Sigma level temperature
c               grid passed in is used. Equation: z2=z1+(T1/gamma)*[1-(P2/P1)**
c                                                     (gamma*R/g)]
c                  Derived from dz=-(RT/Pg)dP  and T=T1(P/P1)**(gamma*R/g)
c               5: set to zero below lowest sigma level
c               6: set to missing below lowest sigma level
c               Most other variables should use 0 since things like negative
c               or over 100% Relative Humidity could occur otherwise.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       include 'GEMINC:GEMPRM.PRM'
       real array(ixl,jxl,kxl),parray(ixl,jxl,kxpl),pres(ixl,jxl,kxl)
       real plevs(kxpl),maxp(kxl),tvirt(ixl,jxl),minp(kxl),alr
       integer etop
       alr=.0065                !degrees K / Meter.  U.S. average lapse rate
       etop=10                  !top level to use for determining lapse rates in
c     extrapolations(1 is the highest sigma level)
       ebot=kxs-1               !bottom level to use for determining lapse rates in
c     extrapolations(kxs is the lowest sigma level)
       if(etop.gt.ebot-1)etop=ebot-1

       do 10 kp=1,kxp
          do 20 ks=1,kxs-1
             if(plevs(kp).ge.minp(ks).or.plevs(kp).le.maxp(ks))then
c
c Case 1: pressure level partly or completely between two sigma levels
c
                do j=1,jx-icross
                   do i=1,ix-icross
                      if(plevs(kp).ge.pres(i,j,ks).and.
     $                     plevs(kp).le.pres(i,j,ks+1))then
c
c Case a: two sigma levels surround pressure level at a point - - interpolate
c
                         parray(i,j,kp)= (array(i,j,ks))     +
     &                        (log(plevs(kp)) - log(pres(i,j,ks))) *
     &                        ((array(i,j,ks+1) - array(i,j,ks))     /
     &                        (log(pres(i,j,ks+1))-log(pres(i,j,ks))))
                      elseif(plevs(kp).gt.pres(i,j,kxs))then
c
c Case b: pressure level below lowest sigma level at a point - - extrapolate
c
                         if(iextrap.eq.1)then
                            parray(i,j,kp)= (array(i,j,etop))    +
     &                           (log(plevs(kp)) -log(pres(i,j,etop))) *
     &                           ((array(i,j,ebot) -array(i,j,etop))   /
     &                           (log(pres(i,j,ebot))-
     $                           log(pres(i,j,etop))))
                         elseif(iextrap.eq.2)then
                            parray(i,j,kp)= (array(i,j,etop))        +
     &                           ((plevs(kp)) -(pres(i,j,etop))) *
     &                           ((array(i,j,ebot) -array(i,j,etop)) /
     &                           ((pres(i,j,ebot)) -(pres(i,j,etop))))
                         elseif(iextrap.eq.3)then !use: t2=t1*(p2/p1)**(gamma*R/g)
                            parray(i,j,kp)= (array(i,j,kxs))   *
     &                           (plevs(kp)/pres(i,j,kxs))**
     &                           (alr*RDGAS/GRAVTY)      
                         elseif(iextrap.eq.4)then !extrapolate geopotential heights
                            parray(i,j,kp)= array(i,j,kxs)+
     $                           (tvirt(i,j)/alr)*
     &                           (1 - (plevs(kp)/pres(i,j,kxs))**
     &                           (alr*RDGAS/GRAVTY)                )
                         elseif(iextrap.eq.5)then
                            parray(i,j,kp)=0.
                         elseif(iextrap.eq.6)then
                            parray(i,j,kp)=-9999.
                         else   !iextrap=0 or input error - do not extrapolate
                            parray(i,j,kp)= array(i,j,kxs)
                         endif
                      elseif(plevs(kp).lt.pres(i,j,1))then
c
c Case c: pressure level above highest sigma level at a point - shouldn't be
c         here, but deal with it.
c
                         parray(i,j,kp)= (array(i,j,2))        +
     &                        (log(plevs(kp)) -log(pres(i,j,2))) *
     &                        ((array(i,j,1)  -array(i,j,2   ))  /
     &                        (log(pres(i,j,1))  -log(pres(i,j,2))))
                      endif
                   enddo
                enddo
             endif
 20       continue
          if(plevs(kp).gt.maxp(kxs))then
c
c Case 2: pressure level completely below lowest sigma level - extrapolate
c
             do j=1,jx-icross
                do i=1,ix-icross
                   if(iextrap.eq.1)then
                      parray(i,j,kp)= (array(i,j,etop))             +
     &                     (log(plevs(kp))  -log(pres(i,j,etop))) *
     &                     ((array(i,j,ebot)    -array(i,j,etop))  /
     &                     (log(pres(i,j,ebot))-log(pres(i,j,etop))))
                   elseif(iextrap.eq.2)then
                      parray(i,j,kp)= (array(i,j,etop))            +
     &                     ((plevs(kp))      -(pres(i,j,etop))) *
     &                     ((array(i,j,ebot) -array(i,j,etop))   /
     &                     ((pres(i,j,ebot))-(pres(i,j,etop))))
                   elseif(iextrap.eq.3)then
                      parray(i,j,kp)= (array(i,j,kxs))         *
     &                     (plevs(kp)/pres(i,j,kxs))**
     &                     (alr*RDGAS/GRAVTY)      
                   elseif(iextrap.eq.4)then !extrapolate geopotential heights
                      parray(i,j,kp)= array(i,j,kxs)+
     $                     (tvirt(i,j)/alr)*
     &                     (1 - (plevs(kp)/pres(i,j,kxs))**
     &                     (alr*RDGAS/GRAVTY)                    )
                   elseif(iextrap.eq.5)then
                      parray(i,j,kp)=0.
                   elseif(iextrap.eq.6)then
                      parray(i,j,kp)=-9999.
                   else         !iextrap=0 or an input error - do not extrapolate
                      parray(i,j,kp)= array(i,j,kxs)
                   endif
                enddo
             enddo
          elseif(plevs(kp).lt.minp(1))then
c     
c Case 3: pressure level completely above highest sigma level - shouldn't be
c           here, but deal with it.
c
             do j=1,jx-icross
                do i=1,ix-icross
                   parray(i,j,kp)= (array(i,j,2))                 +
     &                  (log(plevs(kp))    -log(pres(i,j,2))) *
     &                  ((array(i,j,1)      -array(i,j,2   ))  /
     &                  (log(pres(i,j,1))  -log(pres(i,j,2))))
                enddo
             enddo
          endif
 10    continue

       return
       end



      SUBROUTINE SMTHER(SLAB,slabnew,ix,IMX,jx,JMX,kx,KMX,NPASS,ICRSDOT)
C
C     SECTION  TOOLS
C     PURPOSE  SPATIALLY SMOOTH (USUALLY SLAB) TO REMOVE HIGH
C              FREQUENCY WAVES (2 delta x per pass?)
C     slab is 3-d array to smooth
c     imx is x dimension,ix is actual size in x
c     jmx is y dimension,jx is actual size in y
c     kmx is z dimension,kx is actual size in z
c     npass is number of passes
c     icrsdot=1 for cross points,=0 for dot points 
C
      DIMENSION SLABNEW(IMX,JMX)
      DIMENSION SLAB(IMX,JMX,KMX),XNU(2)
C
c      IE=IMX-1-ICRSDOT
      ie=ix-1-icrsdot
c      JE=JMX-1-ICRSDOT
      je=jx-1-icrsdot
      XNU(1)=0.50
      XNU(2)=-0.52
c      DO 100 K=1,KMX
      do 100 k=1,kx
      DO 100 LOOP=1,NPASS*2
c        if (loop.eq.1)print*,'smoothing level ',k
         N=2-MOD(LOOP,2)
C
C        ... FIRST SMOOTH IN THE IMX DIRECTION
C
         DO 20 I=2,IE
         DO 20 J=2,JE
            SLABNEW(I,J)=SLAB(I,J,K)+XNU(N)*
     *         ((SLAB(I,J+1,K)+SLAB(I,J-1,K))*0.5-SLAB(I,J,K))
20       CONTINUE
         DO 30 I=2,IE
         DO 30 J=2,JE
            SLAB(I,J,K)=SLABNEW(I,J)
30       CONTINUE
C
C        ... NOW SMOOTH IN THE JMX DIRECTION
C
         DO 40 J=2,JE
         DO 40 I=2,IE
            SLABNEW(I,J)=SLAB(I,J,K)+XNU(N)*
     *         ((SLAB(I+1,J,K)+SLAB(I-1,J,K))*0.5-SLAB(I,J,K))
40       CONTINUE
         DO 50 I=2,IE
         DO 50 J=2,JE
            SLAB(I,J,K)=SLABNEW(I,J)
50       CONTINUE
100   CONTINUE
      RETURN
      END



C                                                             C
C*************************************************************C
C                                                             C
C            SUBROUTINE XYLL(X,Y,XLAT,XLON)                   C
C                                                             C
C                                                             C
C   THIS SUBROUTINE TRANSFORMS X AND Y BACK TO LAT AND LON    C
C                                                             C
C                                                             C
C   JULY 28,1980           1ST VERSION BY BILL KUO            C
C                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE XYLL(X,Y,XLAT,XLON,XN,DS,
     &                XLATC,XLONC,IMAX,JMAX)
C   XN: CONE FACTOR =1.0 FOR POLAR STEREOGRAPHIC.
C                   =0.716 FOR LAMBERT CONFORMAL.
C   DS: GRID SPACING IN KM.
C   XLATC,XLONC: CENTRAL LATITUDE AND LONGITUDE.
C   IMAX,JMAX: GRID NUMBER IN Y & X DIRECTION
C   THE X AND Y ARE ON THE SAME COORDINATE AS THE GRID SYSTEM.
C   WHICH MEANS IF X=1.5 IT FALLS BETWEEN J=1 AND J=2
C
C   NOTE, FOR DOT POINT VARIABLES (U AND V) IMAX IS IMAX
C   (SUCH AS 41). FOR CROSS POINT VARIABLES (T, Q, P*) IMAX
C   IMAX - 1 (SUCH AS 40). SAME ARGUMENT HOLDS FOR JMAX.
C
      A=6370.
      CONV=45./ATAN(1.)
      PHI1=30.0
      PHIC=90.0-XLATC
      PHI1=PHI1/CONV
      PHIC=PHIC/CONV
      CNST=A*SIN(PHI1)/(XN*TAN(PHI1/2.)**XN)
      XXC=0.0
      YYC=-1.*CNST*(TAN(PHIC/2.)**XN)
      CENTRI=(FLOAT(IMAX)+1.)/2.
      CENTRJ=(FLOAT(JMAX)+1.)/2.
      XX=(X-CENTRJ)*DS+XXC
      YY=(Y-CENTRI)*DS+YYC
      R=SQRT(XX*XX+YY*YY)
      YXN=1./XN
      PHI=(R*XN/(SIN(PHI1)*A))**YXN
      PHI=PHI*TAN(PHI1/2.)
      PHI=2.0*ATAN(PHI)
      PHI=PHI*CONV
      XLAT=90.0-PHI
      S=ASIN(XX/R)
      S=S*CONV
      S=S/XN
      XLON=XLONC+S
      RETURN
      END

c                                                                       c
c***********************************************************************c
c                                                                       c
      subroutine llxy(xlat,xlon,x,y,xn,ds,xlatc,xlonc,imax,jmax)
c                                                                          
c     Subprogram to calculate x and y given latitude and longitude.        
c     Peter Howells, NCAR, 1984                                            
c                                                                          
      CONV=45./atan(1.)
      a = 6370.0                                                           
      phi1 = 30.0                                                          
      pole = 90.0                                                          
      if ( xlatc.lt.0.0 ) then                                             
         phi1 = -phi1                                                      
         pole = -pole                                                      
      end if                                                               
      phic = ( pole - xlatc )/conv                                         
      phi1 = phi1/conv                                                     
      xc = 0.0                                                             
      yc = -a/xn*sin(phi1)*(tan(phic/2.0)/tan(phi1/2.0))**xn               
      centri = float(imax + 1)/2.0                                         
      centrj = float(jmax + 1)/2.0                                         
c                                                                          
c     Calculate x,y coords. relative to pole                               
c                                                                          
      flp = xn*( xlon - xlonc )/conv
      psx = ( pole - xlat )/conv                      
      r = -a/xn*sin(phi1)*(tan(psx/2.0)/tan(phi1/2.0))**xn                 
      if ( xlatc.lt.0.0 ) then                                             
         xx = r*sin(flp)                                                   
         yy = r*cos(flp)                                                   
      else                                                                 
         xx = -r*sin(flp)                                                  
         yy = r*cos(flp)                                                   
      end if                                                               
c                                                                          
c  Transform (1,1) to the origin                                           
c                                                                          
      x = ( xx - xc )/ds + centrj
      y = ( yy - yc )/ds + centri
 100  continue
c
      return                                                               
      end                                                                  

      subroutine xtodot2(xslab,dotslab,ix,jy,maxix,maxjy)
c  This subroutine interpolates to dot points, but along edges of
c  dotpoint domain, only uses the same values as the crosspoint edge.
c  The result is that fields like precip can not obtain negative values
c  as a result of extrapolation
      dimension xslab(ix,jy),dotslab(ix,jy)
      integer aa,bb,cc,dd
      do i=1,maxix
      do j=1,maxjy
        aa = min(i,maxix-1)
        bb = min(j,maxjy-1)
        cc = max(i-1,1)
        dd = max(j-1,1)
        dotslab(i,j) = 0.25*(xslab(aa,bb)+xslab(cc,bb)
     &                 +xslab(aa,dd)+xslab(cc,dd))
      enddo
      enddo
      return
      end


C
c**********************************************************************c
c                                                                      c
      subroutine xtodot(xslab,dotslab,iy,jx,maxiy,maxjx)
c
      dimension xslab(iy,jx), dotslab(iy,jx)
c
c   Interpolate in the interior.
c
      do 100 j=2,maxjx-1
      do 100 i=2,maxiy-1
         dotslab(i,j)=.25*(xslab(i-1,j-1)+xslab(i,j-1)+xslab(i-1,j)+
     &      xslab(i,j))
  100    continue
c
c   Extrapolate out to top and bottom edges.
c
      do 200 j=2,maxjx-1
         dotslab(1,j)=(3.*(xslab(1,j-1)+xslab(1,j))-
     &      (xslab(2,j-1)+xslab(2,j)))/4.
         dotslab(maxiy,j)=(3.*(xslab(maxiy-1,j-1)+xslab(maxiy-1,j))-
     &      (xslab(maxiy-2,j-1)+xslab(maxiy-2,j)))/4.
  200 continue
c
c   Extrapolate out to left and right edges.
c
      do 300 i=2,maxiy-1
         dotslab(i,1)=(3.*(xslab(i-1,1)+xslab(i,1))-
     &      (xslab(i-1,2)+xslab(i,2)))/4.
         dotslab(i,maxjx)=(3.*(xslab(i-1,maxjx-1)+xslab(i,maxjx-1))-
     &      (xslab(i-1,maxjx-2)+xslab(i,maxjx-2)))/4.
  300 continue
c
c   Extrapolate out to corners.
c
      dotslab(1,1)=(3.*xslab(1,1)-xslab(2,2))/2.
      dotslab(maxiy,1)=(3.*xslab(maxiy-1,1)-xslab(maxiy-2,2))/2.
      dotslab(1,maxjx)=(3.*xslab(1,maxjx-1)-xslab(2,maxjx-2))/2.
      dotslab(maxiy,maxjx)=(3.*xslab(maxiy-1,maxjx-1)-
     &   xslab(maxiy-2,maxjx-2))/2.
      return
      end



