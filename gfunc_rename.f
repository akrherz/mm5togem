      program gfunc_rename
c     -----------------------------------------------------------------  
c     Written by David Ovens (ovens@atmos.washington.edu) 1/21/2000
c
c     -----------------------------------------------------------------  
c     Program Description
c     -----------------------------------------------------------------  
c     This program will read a gempak file and append a string to
c       all variables that it copies to a second gempak file.
c     f77 -o read_gempak_file.x -g read_gempak_file.f  \
c        /usr/local/ldm/NAWIPS-5.4/lib/sol/gemlib.a \
c        /usr/local/ldm/NAWIPS-5.4/lib/sol/appl.a \
c        /usr/local/ldm/NAWIPS-5.4/lib/sol/syslib.a \
c        /usr/local/ldm/NAWIPS-5.4/lib/sol/gemlib.a
c

      implicit none
c     GEMPRM.PRM types
      real RMISSD, RDIFFD, PI, HALFPI, TWOPI, PI4TH, DTR, RTD
     $     ,RADIUS, OMEGA, GRAVTY, RDGAS, RKAP, RKAPPA, AKAPPA
     $     ,GAMUSD, TMCK, RADCLM, RADSKY, RSZPTN
      integer IMISSD, MMKEY, MMHDRS, MMPRT, MMLIST, MMFREE
     $     , MMFILE, MBLKSZ, MCACHE, MMPARM, MMFHDR, MMSRCH
     $     , MMFLDP, MTVAX, MTSUN, MTIRIS, MTAPOL, MTIBM
     $     , MTIGPH, MTULTX, MTHP, MTALPH, MTMACH, MMRECL
     $     , MDREAL, MDINTG, MDCHAR, MDRPCK, MDGRID, MDGNON
     $     , MDGGRB, MDGNMC, MDGDIF, MDGDEC, MFSF, MFSN, MFGD
     $     , MFUNKN, MFAIRW, MFMETR, MFSHIP, MFBUOY, MFSYNP
     $     , MFRAOB, MFVAS, MFGRID, MFTEXT, LLMXLV, LLMXTM
     $     , LLMXGT, LLMXST, LLMXDT, LLMXPT, LLSTFL
     $     , LLMXGD, LLMDGG, MXLOOP, LLNNAV, LLNANL, LLSTHL
     $     , LLGDHD, LLOAGD, LLCLEV, LLAXIS, LLTMCX, IFINVD
     $     , IFAREA, IFGINI, IFNIDS, IFNOWR, JOFLST, JOFLDT
     $     , IGBSZM, IGBSIZ, IGTBSZ, IGDSZM, IGDSIZ, MCCYL
     $     , MPCEQU, MSCEQU, MCAZM, MPAEQU, MPASTR, MPAORT
     $     , MPALAM, MPAGNO, MSANOR, MSASOU, MPCMER, MPCMCD
     $     , MCCON, MPCNOR, MPCSOU, MCMER, MPTMER, MPUTM
     $     , MPOBLQ, MCGOES, MPMCI, MXCLNM, NDVCHR
c     GEMPAK 5.6 Additions
     $     , MXFLSZ, MXNMFL, MFCOUN, LLMXLN, LLMAXD
c     GEMPAK 5.6.A Additions
     $     , IFNCDF, IFNEXZ

c     grdcmn.cmn types
      integer kgrid, khdrln, kbfhdr, lnavbl, lanlbl, igdatm, ndattm
     $     , ktgrid, ksrtl, navsz, ianlsz, ihdrsz
      real savnav, savanl
c
      include 'GEMINC:GEMPRM.PRM'
      include 'grdcmn.cmn'
      
      integer dbug,level(2),IGHDR(2), maxtimes, ipktyp, iargc
     $     , iarglen, isuflen, ilensuffix, iret, istat, ilenout
     $     , ilenin, i, ihd, iflno, inav, ianl, maxg, ier, kx, ky
     $     , igdfil, igdfln, igdflo, ntimes, itime, fhour, igrid
     $     , ivcord, igx, igy, iplen
      parameter (maxtimes=240)    
      real grid(LLMXGD)
      real rnvblk(LLNNAV),anlblk(LLNANL)
      character arg*256,filein*256,fileout*256,suffix*8,ch256*256
      character*12 parm
      CHARACTER	cprj*24
      character*20 gdattm(2),timarr(maxtimes)
      character*20 dummy
      logical respond,write,exists
      data navsz,ianlsz,ihdrsz/LLNNAV,LLNANL,2/

      dbug = 1
      ipktyp  = mdgdec
      if (iargc().gt.0) then
         call getarg(1,arg)
      endif
      if (arg .eq. '-h' .or. arg .eq. '-help'.or.iargc().lt.2) then
         print *,
     $        'gfunc_rename  this program will copy all of the fields'
         print *,
     $        '   in a given GEMPAK input file into a given GEMPAK'
         print *,
     $        '   output file adding an optional suffix to their gfunc'
         print *,
     $        '   names.'
         print *,
     $        'If the output file does not exist, it will be created'
         print *,
     $        '   using the navigation and analysis blocks of the input'
         print *,
     $        '   file.'
         print *,
     $        'Input and output files must have the same navigation and'
         print *,
     $        '   analysis information.'
         print *,
     $        'Due to GEMPAK limitations, input and output files must'
         print *,
     $        '   be no more than 72 characters in length'
         print *,' '
         print *,
     $        'Usage: gfunc_rename infile outfile [suffix]'
         stop 1
      endif
      call trim(arg,iarglen)
      filein = arg(1:iarglen)
      call getarg(2,arg)
      call trim(arg,iarglen)
      fileout = arg(1:iarglen)
      if (iargc().lt.3) then
         suffix = ' '
         isuflen = 1
      else
         call getarg(3,arg)
         call trim(arg,iarglen)
         ilensuffix = len(suffix)
         if (iarglen .gt. ilensuffix) then
            print *,'suffix ',arg,' too long, truncating it to '
     $           ,arg(1:ilensuffix)
            suffix = arg(1:ilensuffix)
            isuflen = ilensuffix
         else
            suffix = arg(1:iarglen)
            isuflen = iarglen
         endif
      endif

c     initialize GEMPAK user interface
      call ip_init(respond,iret)
      if(iret.ne.0)stop 1

      write = .false.
      call ginitp(0,istat,iret)
c     from  /usr/local/gempak/gempak5.4/source/gemlib/fl/flinqr.f
c     gempak determination of file existence
      call fl_inqr(filein,exists,dummy,iret)
      if (exists .ne. .true.) then
         print *,'input file ',filein
         call bail_out('input file does not exist')
      endif
      call fl_inqr(fileout,exists,dummy,iret)
      if(iret.ne.0)stop 2
      if (exists .ne. .true.) then
         call trim(fileout,ilenout)
         call trim(filein,ilenin)
         print *,'creating output file: ', 
     $        fileout(1:ilenout),' from input file: ',filein(1:ilenin)
         if (ilenout .gt. 72) then
            print *,'WARNING: output filename may be'
     $           //' too long ',ilenout,' chars for dg_ofil'
            print *,'errors of [FL 2 or 1018] [DM -2] [GD -2] [DG -30]'
     $           //' may follow'
         endif
         print *,'using gg_init'
         CALL GG_INIT  ( 1, iret )
         IF  ( iret .ne. 0 )  THEN
            print *,iret
            call bail_out('gg_init failure')
         endif
         DO  i = 1, LLNANL
            anlblk (i) = 0.
         END DO
         ihd = 2
C*	      CASE 3: Get the navigation and analysis blocks from the
C*	      old file.
         print *,'using gd_opnf'
         CALL GD_OPNF  ( filein, .false., iflno, inav, rnvblk,
     +        ianl,  anlblk, ihd, maxg, iret )
         IF  ( iret .ne. 0 )  THEN
            print *,iret
            call bail_out('gd_opnf failure')
         endif
         print *,'using gd_clos,gr_rnav,gd_cref'
         CALL GD_CLOS  ( iflno, ier )
         CALL GR_RNAV  ( rnvblk, cprj, kx, ky, ier )
         CALL GD_CREF  ( fileout, LLNNAV, rnvblk, LLNANL,
     +        anlblk, ihd, MMHDRS, igdfil, iret )
         IF  ( iret .eq. 0 )  THEN
            CALL GD_CLOS ( igdfil, ier )
            print *,'succeeded in creating new file ',fileout(1:ilenout)
         else
            print *,iret
            call bail_out('gd_cref failure')
         END IF
c         open(unit=11,file='gdcfilinput',status='unknown')
c         write (11,*) 'gdoutf = ',fileout(1:70)
c         write (11,*) 'proj = '
c         write (11,*) 'grdarea = '
c         write (11,*) 'kxky = '
c         write (11,*) 'maxgrd = MMHDRS'
c         write (11,*) 'cpyfil = ',filein(1:70)
c         write (11,*) 'anlyss = '
c         write (11,*) 'run'
c         write (11,*) ''
c         write (11,*) 'exit'
c         close (11)
c         call system ('$GEMEXE/gdcfil < gdcfilinput')
      endif

      call dg_ofil(filein,fileout,.true.,igdfln,igdflo,iret)
      print *,'DG_OFIL: iret = ',iret,' igdfln,igdflo = ',
     $     igdfln,igdflo
      call gd_gnav(igdfln,rnvblk,navsz,iret)
      if(iret.ne.0)then
         print*,'gd_gnav: iret=',iret
         goto 999
      endif
      call gr_snav(navsz,rnvblk,iret)      
      if(iret.ne.0)then
         print*,'gd_snav: iret=',iret
         goto 999
      endif
      call gd_gtim  ( igdfln, maxtimes, timarr, ntimes, iret )
      if(iret.ne.0)then
         print*,'gd_gtim: iret=',iret
         goto 999
      endif
      do itime = 1,ntimes
c     create reverse time mapping array (0,3,6...45,48 --> 1,2,3,..)
         read(timarr(itime)(13:15),'(i3)') fhour
         print *,'Time #',itime,' = ',timarr(itime),' fhour ='
     $        ,fhour
      enddo

c     kgrid(igdfln), the number of grids in this file, is from the
c       grdcmn.cmn common block
      do 900 igrid = 1,kgrid(igdfln)
C************************************************************************
C* GD_GIDN								*
C*									*
C* This subroutine returns a grid identifier given the grid number.	*
C*									*
C* GD_GIDN  ( IGDFLN, IGNUM, GDATTM, LEVEL, IVCORD, PARM, IRET )	*
C*									*
C* Input parameters:							*
C*	IGDFLN		INTEGER		Grid file number		*
C*	IGNUM		INTEGER		Grid number			*
C*									*
C* Output parameters:							*
C*	GDATTM (2)	CHAR*20		GEMPAK times			*
C*	LEVEL  (2)	INTEGER 	Vertical levels			*
C*	IVCORD		INTEGER		Vertical coordinate		*
C*	PARM		CHAR*12		Parameter name			*
C*	IRET		INTEGER		Return code			*
C*					  0 = normal return		*
C*					 -4 = file not open		*
C*					 -6 = read/write error		*
C*					-12 = invalid grid number	*
         call gd_gidn (igdfln,igrid,gdattm,level,ivcord,parm,iret)
         if (dbug.gt.1) print*,
     $        '1igrid,parm,level,ivcord=',igrid,parm,level,ivcord
         if(iret.ne.0)then
            print*,'gd_gidn: iret=',iret
            goto 999
         endif
C************************************************************************
C* GD_GGRD								*
C*									*
C* This subroutine reads the requested grid from a grid file given	*
C* the grid number.							*
C*									*
C* GD_GGRD  ( IGDFLN, IGNUM, GDATTM, LEVEL, IVCORD, PARM, GRID, IGX,	*
C*            IGY, IGHDR, IRET )					*
C*									*
C* Input parameters:							*
C*	IGDFLN		INTEGER		Grid file number		*
C*	IGNUM		INTEGER		Grid number			*
C*									*
C* Output parameters:							*
C*	GDATTM (2)	CHAR*20		GEMPAK times			*
C*	LEVEL  (2)	INTEGER		Vertical levels			*
C*	IVCORD		INTEGER		Vertical coordinate		*
C*	PARM		CHAR*12		Parameter name			*
C*	GRID (IGX,IGY)	REAL		Grid data			*
C*	IGX		INTEGER		Number of horizontal points	*
C*	IGY		INTEGER		Number of vertical points 	*
C*	IGHDR (IHDRSZ)	INTEGER		Grid header			*
C*	IRET		INTEGER		Return code			*
C*					  0 = normal return		*
C*					 -4 = file not open		*
C*					 -6 = read error		*
C*					-12 = grid does not exist	*
         call gd_ggrd (igdfln,igrid,gdattm,level,ivcord,parm,grid
     $        ,igx,igy,ighdr,iret)
         if (dbug.gt.1) print*,
     $        '2igrid,parm,level,ivcord,kx,ky,iret=',igrid,parm,level
     $        ,ivcord,igx,igy,iret
         if(iret.ne.0)then
c
c     Some model fields are not output at all times and for all levels
c
            if (iret.eq.-6) then
               print *,'no ',parm,
     $              ' actually found at this level and time'
               goto 900
            else
               print*,'gd_ggrd: iret=',iret
               goto 999
            endif
         endif
         call trim(parm,iplen)
         parm = parm(1:iplen)//suffix(1:isuflen)
         print *,'writing new parm name = ',parm
     $        ,' grid(400) = ',grid(400)

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
         call gd_wpgd(igdflo,grid,igx,igy,ighdr,gdattm,level,ivcord,
     $                 parm,.true.,ipktyp,3,iret)
         if(iret.ne.0) then
            print *,'gd_wpgd: iret=',iret,' Stopping . . .'
            call bail_out('error in write_data')
         endif
 900  continue

 999  call gd_clos(igdfln,iret)
      call gd_clos(igdflo,iret)
      call gendp(1,iret)
      if  ( iret .ne. 0 )  call er_wmsg ( 'GEMPLT', iret, ' ', ier )
      call ip_exit(iret)
      if(iret.ne.0)print*,'ip_exit:  iret=',iret

      stop
      end

      subroutine bail_out(gerr)
      integer ier
      character*(*) gerr

      print *,'ERROR: bail_out ',gerr
      call gendp(1,ier)         ! used for GEMPLT
      if  ( ier .ne. 0 )  call er_wmsg ( 'GEMPLT', ier, ' ', ier )
      call ip_exit(ier)
      if(ier.ne.0)print*,'ip_exit:  ier=',ier
      stop
      end

      subroutine trim (word,ilength)
      character*(*) word
      integer ilength,il,ifirst,ilast
      
      il = len(word)
      ilast = 0
      ifirst = 0
      do ii = 1,il
         if (word(ii:ii) .eq. ' ') then
            if (word(ii-1:ii-1) .ne. ' ') ilast = ii - 1
         else if (ifirst .eq. 0) then
            ifirst = ii
         endif
      enddo
      ilength = (ilast - ifirst) + 1
      ii = 0
      do i = ifirst,ilast
         ii = ii + 1
         word(ii:ii) = word(i:i)
      enddo
      do ii = ilength+1,il
         word(ii:ii) = ' '
      enddo
      return
      end
