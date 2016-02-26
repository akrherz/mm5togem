      program read_gempak_file
c     -----------------------------------------------------------------  
c     NOTE:  THE SIMPLEST VERSION OF THIS PROGRAM IS CONTAINED IN
c            ~ovens/gempak/read_gempak.tar.gz 
c
c     Written by David Ovens (ovens@atmos.washington.edu) 6/4/97
c
c     -----------------------------------------------------------------  
c     Program Description
c     -----------------------------------------------------------------  
c     This program will read a gempak file and fill the 2-D grid()
c       array with the data for each successive grid in the file.  The
c       grid order is the same as that displayed by running gdinfo for
c       the file.
c     In addition to filling the grid() array, the analysis block array,
c       anlblk(), and the navigation block array, rnvblk(), are also filled.
c     An input file called "read_gem_input_file" is used to indicate
c       the name of the gempak file as well as other parameters that
c       might be used in the future to choose what fields and what
c       levels the user wants output in another file.  At present, 6/4/97,
c       the only item in this input file that is used is the name of the
c       gempak file (innam).
c     -----------------------------------------------------------------  
c     Dependent Files and Compiling
c     -----------------------------------------------------------------  
c     The makefile in this directory indicates how to compile this
c       code.  Essentially, you need to link to $LIBDIR (here it is
c       /usr/local/ldm/NAWIPS-5.4) with its three main
c       libraries in the following order (yes, gemlib.a is used twice):
c        $LIBDIR/lib/sol/gemlib.a
c        $LIBDIR/lib/sol/appl.a
c        $LIBDIR/lib/sol/syslib.a
c        $LIBDIR/lib/sol/gemlib.a
c       and the following include files should be linked to the directory
c       where you are compiling this program:
c          Local Link            Actual file      
c        GEMINC:ADBUFF.CMN -> $LIBDIR/gempak5.4/include/adbuff.cmn
c        GEMINC:ERMISS.FNC -> $LIBDIR/gempak5.4/include/ermiss.fnc
c        GEMINC:GEMPRM.PRM -> $LIBDIR/gempak5.4/include/gemprm.SunOS
c        GEMINC:GMBDTA.CMN -> $LIBDIR/gempak5.4/include/gmbdta.cmn
c        grdcmn.cmn        -> $LIBDIR/gempak5.4/source/gemlib/gd/grdcmn.cmn
c     Remember that the gempak file you wish to read is declared in the
c       "read_gem_input_file" file in this directory.
c     The final compile line is this on Solaris:
c     f77 -o read_gempak_file.x -g read_gempak_file.f  \
c        /usr/local/ldm/NAWIPS-5.4/lib/sol/gemlib.a \
c        /usr/local/ldm/NAWIPS-5.4/lib/sol/appl.a \
c        /usr/local/ldm/NAWIPS-5.4/lib/sol/syslib.a \
c        /usr/local/ldm/NAWIPS-5.4/lib/sol/gemlib.a
c

      include 'GEMINC:GEMPRM.PRM'
      include 'grdcmn.cmn'
      parameter(ixl=140,jyl=140,kzl=50,maxtimes=20)
      parameter(kzp=kzl,num_out_flds=40)

      integer mif(1000,20),levarr(2,kzl),level(2)
      real mrf(1000,20),rnvblk(LLNNAV),anlblk(LLNANL),grid(ixl,jyl)
      real angle1,angle2,angle3 
      real plevs(kzp),rlat(ixl,jyl),rlon(ixl,jyl)
      character*12 parm
      character*20 gdattm(2),firstm,lasttm,timarr(maxtimes)
      character*72 cproj
      character*80 filein,outnam, innam
      character i_output(num_out_flds)*1,sig_fld(num_out_flds)*4
     $     ,prs_fld(num_out_flds)*4,fld_name(num_out_flds)*70
      character dummy
      logical write,respond
      data navsz,ianlsz,ihdrsz/LLNNAV,LLNANL,2/

      filein = 'read_gem_input_file'
      open (unit=10,file=filein
     $     ,form='formatted',status='unknown',err=999)
c     read info about how different types of observations are stored
c     in the ftype variable in the program
      read(10,1005) innam
 1005 format(/,a80)
      read(10,1005) outnam
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

c     initialize GEMPAK user interface
      call ip_init(respond,iret)
      if(iret.ne.0)stop 1

      write = .false.
C************************************************************************
C* GINITP                                                               *
C*                                                                      *
C* This subroutine initializes GEMPLT for an application program.       *
C* The GPLT subprocess is begun if it has not been started by a         *
C* previously-executed program.  This subroutine must be the first      *
C* subroutine called by any program that uses GEMPLT.  MODE specifies   *
C* whether map or graph plots are to be made.                           *
C*                                                                      *
C* GINITP  ( MODE, ISTAT, IRET )                                        *
C*                                                                      *
C* Input parameters:                                                    *
C*      MODE            INTEGER         Plot mode                       *
C*                                         0 = no change                *
C*                                         1 = map coordinates          *
C*                                         2 = graph coordinates        *
C*                                                                      *
C* Output parameters:                                                   *
C*      ISTAT           INTEGER         Status code                     *
C*                                         0 = GPLT started             *
C*                                         1 = GPLT previously started  *
C*                                                                      *
C*      IRET            INTEGER         Return code                     *
C**                                                                     *
C* Log:                                                                 *
C* I. Graffman/RDS       4/85   GEMPLT Version 3.1                      *
C* M. desJardins/GSFC    5/88   Documentation                           *
C* I. Graffman/RDS       6/88   Clean up                                *
C* M. desJardins/NMC     7/91   UNIX version                            *
C* M. desJardins/NMC     1/92   Make UNIX & VMS calls identical         *
C* L. Williams/EAi       3/94   Removed blank comments from header      *
C************************************************************************
       call ginitp(0,istat,iret)
C************************************************************************
C* GD_OPNF								*
C*									*
C* This subroutine opens an existing GEMPAK grid file.  If the		*
C* file requires shared, write access, the subroutine GD_OPNR		*
C* should be used.							*
C*									*
C* GD_OPNF ( FILNAM, WRTFLG, IGDFLN, NAVSZ, RNVBLK, IANLSZ, ANLBLK, 	*
C*           IHDRSZ, MAXGRD, IRET )					*
C*									*
C* Input parameters:							*
C*	FILNAM		CHAR*		File name 			*
C*	WRTFLG		LOGICAL		Flag for write access		*
C*									*
C* Output parameters:							*
C*	IGDFLN		INTEGER		File number			*
C*	NAVSZ		INTEGER		Navigation block length		*
C*	RNVBLK (NAVSZ)	REAL		Navigation block		*
C*	IANLSZ		INTEGER		Analysis block length		*
C*	ANLBLK (IANLSZ)	REAL		Analysis block			*
C*	IHDRSZ		INTEGER		Grid header length		*
C*	MAXGRD		INTEGER		Maximum number of grids		*
C*	IRET		INTEGER		Return code			*
C*					  0 = normal return		*
C*					 -2 = file cannot be opened	*
C*					 -7 = not a GEMPAK5 grid file	*
C*					 -8 = nav cannot be read	*
C*					-13 = grid header too long	*
C*					-14 = file name is blank	*
      call gd_opnf(innam,write,igdfln,navsz,rnvblk,ianlsz,anlblk,
     $     ihdrsz,maxgrd,iret)
      if(iret.ne.0)then
         print*,'gd_opnf: iret=',iret
         goto 999
      endif

C************************************************************************
C* GR_RNAV								*
C*									*
C* This subroutine gets the projection and grid size from a grid	*
C* navigation block.							(
C*									*
C* GR_RNAV  ( RNVBLK, PROJ, KX, KY, IRET )				*
C*									*
C* Input parameters:							*
C*	RNVBLK (LLNNAV)	REAL		Navigation block		*
C*									*
C* Output parameters:							*
C*	PROJ		CHAR*		Projection name			*
C*	KX		INTEGER		Number of points in x dir	*
C*	KY		INTEGER		Number of points in y dir	*
C*	IRET		INTEGER		Return code			*
C*					  0 = normal return		*
C*					 -6 = invalid navigation	*
      call gr_rnav(rnvblk,cproj,kx,ky,iret)
      if(iret.ne.0)then
         print*,'gr_rnav: iret=',iret
         goto 999
      endif
c
c     This section fills header arrays for MM5
c
c     kx = rnvblk (5)  ! redundant
c     ky = rnvblk (6)  ! redundant
      rlat1 = rnvblk (7)  
      rlon1 = rnvblk (8)  
      rlat2 = rnvblk (9)  
      rlon2 = rnvblk (10) 
      angle1 = rnvblk (11) 
      angle2 = rnvblk (12) 
      angle3 = rnvblk (13) 
      if (cproj.eq.'LCC') then !projection type (lambert conformal=LCC)
         mif(4,1) = 1
         mrf(6,1) = angle1      !standard latitude 1 for GEMPAK (true lat 2 MM5)
         mrf(3,1) = angle2      !central longitude for GEMPAK (pole longitude in MM5)
         mrf(5,1) = angle3      !standard latitude 2 for GEMPAK (true lat 1 MM5)
      elseif(cproj.eq.'STR')then !polar stereographic projection
         mif(4,1) = 2
         mrf(7,1) = angle1      !latitude of point of tangency
         mrf(3,1) = angle2      !longitude of point of tangency
      elseif(cproj.eq.'MER')then !Mercator projection
         mif(4,1) = 3
         mrf(7,1) = angle1      !latitude of point of tangency
         mrf(3,1) = angle2      !longitude of point of tangency
      else
         print*,'This map projection is not supported'
         goto 999
      endif 
c     end of MM5 section
      
C************************************************************************
C* GR_SNAV                                                              *
C*                                                                      *
C* This subroutine sets up a grid coordinate system in GEMPLT.  The     *
C* navigation block should be sent as it was received from the grid     *
C* file open subroutine.  Note that the graphics projection and mode    *
C* must be defined before GR_SNAV is called.  This subroutine will fail *
C* if the grid mode is not the same as the current GEMPLT mode.         *
C*                                                                      *
C* GR_SNAV  ( NAVSZ, RNVBLK, IRET )                                     *
C*                                                                      *
C* Input parameters:                                                    *
C*      NAVSZ           INTEGER         Length of navigation block      *
C*      RNVBLK (NAVSZ)  REAL            Navigation block                *
C*                                                                      *
C* Output parameters:                                                   *
C*      IRET            INTEGER         Return code                     *
C*                                        0 = normal return             *
C*                                       -6 = invalid navigation type   *
C*                                       -7 = GEMPLT error              *
C**                                                                     *
C* Log:                                                                 *
C* M. desJardins/GSFC   12/84                                           *
C* M. desJardins/GSFC    8/88   Fixed for GEMPAK4                       *
C* K. Brill/NMC         01/92   Replace GERROR with ER_WMSG             *
C************************************************************************
       call gr_snav(navsz,rnvblk,iret)      
      
C************************************************************************
C* GD_GTIM								*
C*									*
C* This subroutine returns all the times present in a grid file.	*
C* Only the first times are returned.  They are sorted from earliest	*
C* to latest.								*
C*									*
C* GD_GTIM  ( IGDFLN, MAXTIM, TIMARR, NTIMES, IRET )			*
C*									*
C* Input parameters:							*
C*	IGDFLN		INTEGER		Grid file number		*
C*	MAXTIM		INTEGER		Maximum number of times		*
C*									*
C* Output parameters:							*
C*	TIMARR (NTIMES)	CHAR*		GEMPAK times			*
C*	NTIMES		INTEGER		Number of times			*
C*	IRET		INTEGER		Return code			*
C*					  0 = normal return		*
C*					 -4 = file not open		*
C*					 -6 = read/write error		*
      call gd_gtim  ( igdfln, maxtimes, timarr, ntimes, iret )
      if(iret.ne.0)then
         print*,'gd_gtim: iret=',iret
         goto 999
      endif
      do it = 1,ntimes
         print *,'Time #',it,' = ',timarr(it)
         gdattm(it) = timarr(it)
      enddo

c     kgrid(igdfln), the number of grids in this file, is from the
c       grdcmn.cmn common block

C************************************************************************
C* GR_LTLN                                                              *
C*                                                                      *
C* This subroutine computes the latitude and longitude at each grid     *
C* point.  The grid must be defined in GEMPLT before this subroutine    *
C* is called.                                                           *
C*                                                                      *
C* GR_LTLN  ( KX, KY, RLAT, RLON, IRET )                                *
C*                                                                      *
C* Input parameters:                                                    *
C*      KX              INTEGER         Number of points in x dir       *
C*      KY              INTEGER         Number of points in y dir       *
C*                                                                      *
C* Output parameters:                                                   *
C*      RLAT (KX,KY)    REAL            Latitudes in degrees            *
C*      RLON (KX,KY)    REAL            Longitudes in degrees           *
C*      IRET            INTEGER         Return code                     *
C*                                        0 = normal return             *
C*                                       -6 = grid projection error     *
C**                                                                     *
C* Log:                                                                 *
C* M. desJardins/GSFC   11/88   From DG_LTLN                            *
C* K. Brill/NMC         01/92   Replace GERROR with ER_WMSG             *
C************************************************************************
       call gr_ltln( kx, ky, rlat, rlon, iret )
       if(iret.ne.0)then
          print*,'gr_ltln: iret=',iret
          goto 999
       endif


      do igrid = 1,kgrid(igdfln)

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
         call gd_ggrd (igdfln,igrid,gdattm,level,ivcord,parm,grid,igx,
     $        igy,ighdr,iret)
         if(iret.ne.0)then
c
c     Some model fields are not output at all times and for all levels
c
            if (parm(1:4).eq.'QCLD' .and. iret.eq.-6) then
               print *,'no QCLD found at this level and time'
            else
               print*,'gd_ggrd: iret=',iret
               goto 999
            endif
         endif
         print *,'gdattm,level,ivcord,parm,igx,igy,grid(10,10)='
     $        ,gdattm(1),level(1),ivcord,parm,igx,igy,grid(10,10)


C************************************************************************
C* GD_GLEV								*
C*									*
C* This subroutine returns all the levels present in a grid file for	*
C* a given date and vertical coordinate.  The levels returned are	*
C* not sorted.								*
C*									*
C* GD_GLEV  ( IGDFLN, GDATTM, IVCORD, MAXLEV, LEVARR, NLEV, IRET )	*
C*									*
C* Input parameters:							*
C*	IGDFLN		INTEGER		Grid file number		*
C*	GDATTM (2)	CHAR*20		GEMPAK times			*
C*	IVCORD		INTEGER		Vertical coordinate		*
C*	MAXLEV		INTEGER		Maximum number of levels	*
C*									*
C* Output parameters:							*
C*	LEVARR (2,NLEV)	INTEGER		Levels found			*
C*	NLEV		INTEGER		Number of levels found		*
C*	IRET		INTEGER		Return code			*
C*					  0 = normal return		*
C*					 -4 = file not open		*
C*					 -6 = read/write error		*
c         call gd_glev (igdfln, gdattm, ivcord, kzl, levarr, nlev, iret )
c         if(iret.ne.0)then
c            print*,'gd_glev: iret=',iret
c            goto 999
c         endif
c         do it = 1,nlev
c            print *,'levarr ',it,' = ',levarr(1,it)
c         enddo
         
      enddo

C************************************************************************
C* GD_RDAT								*
C*									*
C* This subroutine reads the requested grid from a grid file.		*
C*									*
C* GD_RDAT  ( IGDFLN, GDATTM, LEVEL, IVCORD, PARM, GRID, IGX, IGY,	*
C*            IGHDR,  IRET )						*
C*									*
C* Input parameters:							*
C*	IGDFLN		INTEGER		Grid file number		*
C*	GDATTM (2)	CHAR*20		GEMPAK times			*
C*	LEVEL  (2)	INTEGER		Vertical levels			*
C*	IVCORD		INTEGER		Vertical coordinate		*
C*					  0 = NONE			*
C*					  1 = PRES			*
C*					  2 = THTA			*
C*					  3 = HGHT			*
C*					  4 = SGMA			*
C*	PARM		CHAR*12		Parameter name			*
C*									*
C* Output parameters:							*
C*	GRID (IGX,IGY)	REAL		Grid data			*
C*	IGX		INTEGER		Number of horizontal points	*
C*	IGY		INTEGER		Number of vertical points	*
C*	IGHDR (IHDRSZ)	INTEGER		Grid header			*
C*	IRET		INTEGER		Return code			*
C*					  0 = normal return		*
C*					 -4 = file not open		*
C*					 -6 = read/write error		*
C*					-12 = grid does not exist	*
c      call gd_rdat(igdfln, gdattm, level, ivcord, parm, grid, igx, igy,
c     $     ighdr,  iret )
      
999   print *,'End'
      call gd_clos(igdfln,iret)
c     call gendp(1,iret)   ! used for GEMPLT
c     if  ( iret .ne. 0 )  call er_wmsg ( 'GEMPLT', iret, ' ', ier )
      call ip_exit(iret)
      if(iret.ne.0)print*,'ip_exit:  iret=',iret
      stop 'File reading finished'
      end



