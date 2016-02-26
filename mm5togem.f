      program mm5togem
      implicit none

c     GEMPRM.PRM types
      real RMISSD, RDIFFD, PI, HALFPI, TWOPI, PI4TH, DTR, RTD
     $     ,RADIUS, OMEGA, GRAVTY, RDGAS, RKAP, RKAPPA, AKAPPA
     $     ,GAMUSD, TMCK, RADCLM, RADSKY, RSZPTN
      integer IMISSD, MMKEY, MMHDRS, MMPRT, MMLIST, MMFREE
     $     , MMFILE, MBLKSZ, MCACHE, MMPARM, MMFHDR, MMSRCH
     $     , MMFLDP, MTVAX, MTSUN, MTLNUX, MTIRIS, MTAPOL, MTIBM
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
      INCLUDE 'GEMINC:GEMPRM.PRM'
c
      integer igdatm,kbfhdr,kgrid,khdrln,ksrtl,ktgrid,lanlbl
     $     ,lnavbl,ndattm,ntimes,ndates
      real savanl,savnav,ftimes(60),lasttime
      INCLUDE 'grdcmn.cmn'
      character chdates(60)*24,lastdate*24
      character*80 mrfc(1000,20),mifc(1000,20),
     $     outnam, innam,tempfile,newtempfile,chdum
      character*128 filein
      character*72 oldproj,cproj
      character*40 tmpa
      character*20 gdattm(2)
      character*12 parm
      character yymmdd*6,chdate14*14
      character dummy
      character dataform*8

      integer ighdr(2),idomid

      real mrfv1(1000,20),rnvblk(LLNNAV),anlblk(LLNANL)
      integer bhi(50,20), flag
      real bhr(20,20)
      character*80 bhic(50,20),bhrc(20,20)
      integer :: ndim
      real :: time, sample
      integer, dimension(4) :: start_index, end_index
      character (len= 4) :: staggering
      character (len= 4) :: ordering
      character (len=24) :: start_date
      character (len=24) :: current_date
      character (len= 9) :: name
      character (len=25) :: units
      character (len=46) :: description

c  Raw data 3-d fields (ivert,ihorz,kz):
      real, allocatable, dimension(:,:,:) ::  uu, vv, tt, qv, qc,
     $     qr, cice, sice, pp, rh, rtten

c     3-d (ivert,ihorz,kz+1):
      real, allocatable, dimension(:,:,:) :: ww

c     3-d (ivert,ihorz,6)
      real, allocatable, dimension(:,:,:) :: soilt

c  Raw precip 2-d fields (ivert,ihorz):
      real, allocatable, dimension(:,:)   :: re, rc, prec,
     $     prec_3hr, re_3hr, rc_3hr, tot_prec,
     $     last_output_hr_tot_prec, last_output_hr_exp_prec,
     $     last_output_hr_con_prec

c  Raw data 2-d fields (ivert,ihorz):
      real, allocatable, dimension(:,:)   ::  pstx, tg, ter,
     $     xmap, dotmap, dotcor, sst, slpcrs, irst, xlat, xluse,
     $     xsnowc, xlon

c  Derived 3d fields (ivert,ihorz,kz):
      real, allocatable, dimension(:,:,:) ::  prs, dew, wwh,
     $     prsd, pdata, ttp, prsp, qvlog,
     $     tmp3dp, tmp3ds
      real, allocatable, dimension(:,:,:) :: hlhsl ! ivert*ihorz*2
c     derived 3-d fields (ihorz,ivert,kz):
      real, allocatable, dimension(:,:,:) :: unor,vnor

c  Derived 2d fields (ivert,ihorz):
      real, allocatable, dimension(:,:)   :: pstd, sfp, slp,
     $     dotvar, dat, dluse, tmp2d

c  Miscellaneous 2-d fields (ivert,ihorz):
      real, allocatable, dimension(:,:)   :: tvirt

c     Miscellaneous 1-d fields (kz):
      real, allocatable, dimension(:)     :: plevs, maxp, minp,
     $     maxpd, minpd, sih, modplevs

c     optional 2-d fields (ivert,ihorz):
      real, allocatable, dimension(:,:)   :: senhtflx, lathtflx,
     $     lonwvrad, shtwvrad, pblht, ustar,
     $     pblregime, moisture, u10, v10, t2

      character i_output(150)*1,sig_fld(150)*12
     $     ,prs_fld(150)*12,fld_name(150)*70

      integer mifv1(1000,20),ilvl(2)
      integer aa,bb,cc,dd
      integer iptimes(20),output_hr,last_output_hr,livert,lihorz
      integer dom_special, special_hour
c  Frequencies at which precip may be written out:
c  Note:  0 writes out precip since last output time as P00M if output
c         is more often than 1-hourly.  99 is skipped
      integer navsz,ianlsz,ihdrsz,ighttype,itapfrq,maxgrd,igemfrq,i,
     $     iargc,ioverwrite,iprecise,isprs,issig,issfp,num_pres_levs
     $     ,kdum,num_3d_outputs,iline,num_2d_outputs,iprog,ivert,ihorz
     $     ,kz,iexpanded,iret,istat,ipktyp,itime,iout,idot,icross,ivertl
     $     ,ihorzl,ntotlevels,kp,k,igdfln,ioldkx,ioldky,iend
     $     ,iday,ihour,imin,imonth,isec,iyear,ihr,iftim,idate1,idate2
     $     ,ifhr,ifmn,j,imiss,izagllhsl,kk,ii,jj,jik,ijk,inum,ivcor
     $     ,iflip,imakedot,iusedotvar,ier
      real eps,ptop,yllc,xllc,yurc,xurc,xlatll,xlonll,xlatur,xlonur
     $     ,deltan,deltax,deltay,wmixrat,tvll,beta,gpll,tvavg,q,t
     $     ,tcels,p,e,esat,qsat,hlhsltot,zagllhsl
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
      real gbnds(4),ebnds(4)
      integer iebnds(4)
      real angle1,angle2,angle3 
      external dg_nrel

c     necessary after GEMPAK 5.6
      data gdattm(2)/'                    '/

      eps = 0.622
      ftimes = 0. ! array
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
      if (iargc() .gt. 0) then
         call getarg(1,filein)
         if (filein .eq. '-h' .or. filein .eq. '-help') then
            print*,'Usage: mm5togem.exe [options]'
            print *, '  where:'
            print *, '    input_file   --> name of input file'
            print *, '                     DEFAULT is input_file'
            print *, ' other options                                 ',
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
            stop 1
         endif
      endif
      open (unit=10,file=filein
     $     ,form='formatted',status='unknown',err=999)
c     read info about how different types of observations are stored
c     in the ftype variable in the program
      read(10,1005) innam
 1005 format(/,a80)
      read(10,1005) outnam
      read (10,*) dummy
      read(10,*) tempfile,newtempfile,dom_special,special_hour
      read(10,1006) ioverwrite,issig,isprs,iprecise
 1006 format(/,i1,//,i1,//,i1,//,i1)
      read(10,1007) issfp,num_pres_levs

      allocate(plevs(num_pres_levs))
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

      do iline = num_3d_outputs+1,num_3d_outputs+num_2d_outputs
         read(10,1008) i_output(iline),sig_fld(iline),prs_fld(iline)
     $        ,fld_name(iline)
      enddo
 1007 format(/,/,i1,x,i3)
 1008 format(a1,x,a12,x,a12,x,a)
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
        endif
        i=i+1
      enddo
c
c Initialize all logical variables defining output fields to true
c
      ifvv=.false.
      iftt=.false.
      ifqv=.false.
      ifqc=.false.
      ifqr=.false.
      ifww=.false.
      ifrtten=.false.
      ifcon=.false.
      ifres=.false.
c

      print *, '          MM5 TO GEMPAK CONVERTER'
      print *, ' '
      print *, 'This Program converts raw MM5 Version 1 and 2'
      print *, '         output to GEMPAK format.'
      print *, ' '
 
10    continue
c
c   Read header:
c
      open(unit=20,file=innam,status='old',form='unformatted')
      read (20,err=130,end=999) mifv1,mrfv1,mifc,mrfc
      dataform='mm5v1   '
      iprog=mifv1(1,1)
      ivert = mifv1(104,1)
      ihorz = mifv1(105,1)
      if ((mifv1(103,1).eq.0) .and. (mifv1(5,1).eq.1) .and. 
     &     (iprog.le.3)) then
         ivert = mifv1(6,1)
         ihorz = mifv1(7,1)
      endif
      kz = mifv1(101,iprog)
      if ((iprog.eq.5).or.(iprog.eq.6))
     &   kz = nint(mrfv1(101,iprog))
      goto 200
 130  rewind (20)
c
      read(20,err=140,end=180) flag
      if (flag.ne.0) then
         print*,'MM5 V3 data does not begin with a big header.'
         print*,'Stopping.'
         call exit(1)
      endif
      read(20) bhi, bhr, bhic, bhrc
c
c      do j=1,20
c         do i=1,50
c            if (bhi(i,j).ne.-999)
c     &         print 633, i,j,bhi(i,j),bhic(i,j)(1:50)
c         enddo
c         do i=1,20
c            if (abs(bhr(i,j)+999.).gt..01)
c     &         print 634, i,j,bhr(i,j),bhrc(i,j)(1:50)
c         enddo
c      enddo
c 633  format('bhi(',i2,',',i2,'): ',i6,' : ',a50)
c 634  format('bhr(',i2,',',i2,'): ',f9.2,' : ',a50)
c
      dataform='mm5v3   '
      iprog=bhi(1,1)
      ivert = bhi(16,1)
      ihorz = bhi(17,1)
      rewind (20)
      if (iprog.eq.1) then
         kz=1
      elseif (iprog.eq.2.or.iprog.eq.3) then
c
c      Number of levels needs to be obtained from the header for a
c      3D array--can't trust bhi(12,iprog) to be correct.
c
         read(20) flag
         if (flag.ne.0) then
            print*,'MM5 V3 data does not begin with a big header.'
            print*,'Stopping.'
            stop
         endif
         read(20) bhi, bhr, bhic, bhrc
 703     read(20) flag
         if (flag.eq.1) then
            read (20) ndim,start_index,end_index,time,
     &         staggering,ordering,current_date,name,
     &         units,description
            if (name.ne.'U        ') then
               read(20)
               goto 703
            else
               ntotlevels=end_index(3)-start_index(3)+1
               kz=ntotlevels-1   ! Don't include surface level in kz
            endif
         else
            print*,'Ran into a flag not =1 looking for U.'
            stop
         endif
         rewind (20)
      elseif (iprog.eq.5.or.iprog.eq.11) then
         kz = bhi(12,iprog)
      else
         print*,'For MM5V3 data, can only read output from'
         print*,'Terrain, Regrid, Rawins/Little_r, Interp, or MM5.'
         stop
      endif
      if (iprog.ge.1.and.iprog.le.3) then
c
c      Need to determine if this is "expanded domain".
c      Look for 'TERRAIN' variable.
c
         read(20) flag
         if (flag.ne.0) then
            print*,'MM5 V3 data does not begin with a big header.'
            print*,'Stopping.'
            stop
         endif
         read(20) bhi, bhr, bhic, bhrc
 733     read(20) flag
         if (flag.eq.1) then
            read (20) ndim,start_index,end_index,time,
     &         staggering,ordering,current_date,name,
     &         units,description
            if (name.ne.'TERRAIN  ') then
               read(20)
               goto 733
            else
               if (end_index(1).ne.ivert.or.end_index(2).ne.ihorz) then
c
c               Looks like this is an expanded domain.
c
                  ivert=end_index(1)-start_index(1)+1
                  ihorz=end_index(2)-start_index(2)+1
                  if (bhi(8,1).ne.1.or.bhi(13,1).ne.1.or.
     &                ivert.ne.bhi(9,1).or.ihorz.ne.bhi(10,1)) then
                     print*,'Looks like an expanded domain but ',
     &                  'some things in headers are not consistent.'
                     print*,'bhi(8,1),bhi(13,1),bhi(9,1),bhi(10,1)='
                     print*,bhi(8,1),bhi(13,1),bhi(9,1),bhi(10,1)
                     print*,'ivert,ihorz (dervied from "TERRAIN" ',
     &                      'variable header)='
                     print*,ivert,ihorz
                     stop
                  endif
                  iexpanded=1
               endif
            endif
         else
            print*,'Ran into a flag not =1 looking for TERRAIN.'
            stop
         endif
         rewind (20)
      endif
      goto 200
 140  continue
c
      print*,'The model data header is not a format'//
     &   ' that RIPDP recognizes.  Stopping.'
      stop
c
 180  print*,'Unexpected EOF reached when trying to read'
      print*,'model data header.  Stopping.'
      stop
c
 200  continue
c
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
cmmm      call ginitp(0,istat,iret)
      call ginitp(1,istat,iret)
c     from  /usr/local/gempak/gempak5.4/source/gemlib/fl/flinqr.f
c     gempak determination of file existence
      call fl_inqr(outnam,exists,dummy,iret)
      if(iret.ne.0)stop 2
      if(exists .and. ioverwrite.eq.1) then
         write(chdum,846) outnam(1:72)
 846     format('/bin/rm ',a72)
         call system(chdum)
         print *,chdum
         exists=.false.
      endif
cmmmm hmmm not sure if this would work for file that doesn't exist...
cdoubt      call dg_ofil(outnam,'',.false.,igdfln,idum,iret)
cdoubt      print *,'DG_OFIL: iret = ',iret,' idum = ',idum
cmmmmmm
c      if(exists .and. ioverwrite.eq.0) goto 999
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
      ivertl = ivert
      ihorzl = ihorz
      print *,'ivert,ihorz,kz = ',ivert,ihorz,kz
c
c     allocate memory space
c
c  Raw data 3-d fields (ivert,ihorz,kz):
      allocate(uu(ivert,ihorz,kz))
      allocate(vv(ivert,ihorz,kz))
      allocate(tt(ivert,ihorz,kz))
      allocate(qv(ivert,ihorz,kz))
      allocate(qc(ivert,ihorz,kz))
      allocate(qr(ivert,ihorz,kz))
      allocate(cice(ivert,ihorz,kz))
      allocate(sice(ivert,ihorz,kz))
      allocate(pp(ivert,ihorz,kz))
      allocate(rh(ivert,ihorz,kz))
      allocate(rtten(ivert,ihorz,kz))
c  3-d (ivert,ihorz,kz+1):
      allocate(ww(ivert,ihorz,(kz+1)))
c  3-d (ivert,ihorz,6)
      allocate(soilt(ivert,ihorz,6))
c  Raw precip 2-d fields (ivert,ihorz):
      allocate(re(ivert,ihorz))
      allocate(rc(ivert,ihorz))
      allocate(prec(ivert,ihorz))
      allocate(prec_3hr(ivert,ihorz))
      allocate(re_3hr(ivert,ihorz))
      allocate(rc_3hr(ivert,ihorz))
      allocate(last_output_hr_tot_prec(ivert,ihorz))
      allocate(last_output_hr_exp_prec(ivert,ihorz))
      allocate(last_output_hr_con_prec(ivert,ihorz))
c  Raw data 2-d fields (ivert,ihorz):
      allocate(pstx(ivert,ihorz))
      allocate(tg(ivert,ihorz))
      allocate(ter(ivert,ihorz))
      allocate(xmap(ivert,ihorz))
      allocate(dotmap(ivert,ihorz))
      allocate(dotcor(ivert,ihorz))
      allocate(sst(ivert,ihorz))
      allocate(slpcrs(ivert,ihorz))
      allocate(irst(ivert,ihorz))
      allocate(xlat(ivert,ihorz))
      allocate(xluse(ivert,ihorz))
      allocate(xsnowc(ivert,ihorz))
      allocate(xlon(ivert,ihorz))
c  Derived 3d fields (ivert,ihorz,kz):
      allocate(prs(ivert,ihorz,kz))
      allocate(dew(ivert,ihorz,kz))
      allocate(wwh(ivert,ihorz,kz))
      allocate(prsd(ivert,ihorz,kz))
      allocate(pdata(ivert,ihorz,kz))
      allocate(ttp(ivert,ihorz,kz))
      allocate(prsp(ivert,ihorz,kz))
      allocate(qvlog(ivert,ihorz,kz))
      allocate(tmp3dp(ivert,ihorz,kz))
      allocate(tmp3ds(ivert,ihorz,kz))
c     derived 3-d fields (ihorz,ivert,kz):
      allocate(unor(ihorz,ivert,kz))
      allocate(vnor(ihorz,ivert,kz))
c  ivert*ihorz*2
      allocate(hlhsl(ivert,ihorz,2))
c  Derived 2d fields (ivert,ihorz):
      allocate(pstd(ivert,ihorz))
      allocate(sfp(ivert,ihorz))
      allocate(slp(ivert,ihorz))
      allocate(dotvar(ivert,ihorz))
      allocate(dat(ivert,ihorz))
      allocate(dluse(ivert,ihorz))
      allocate(tmp2d(ivert,ihorz))
      allocate(tot_prec(ivert,ihorz))
c  Miscellaneous 2-d fields (ivert,ihorz):
      allocate(tvirt(ivert,ihorz))
c  Miscellaneous 1-d fields (kz):
      allocate(maxp(kz))
      allocate(minp(kz))
      allocate(maxpd(kz))
      allocate(minpd(kz))
      allocate(sih(kz))
c  optional 2-d fields (ivert,ihorz):
      allocate(senhtflx(ivert,ihorz))
      allocate(lathtflx(ivert,ihorz))
      allocate(lonwvrad(ivert,ihorz))
      allocate(shtwvrad(ivert,ihorz))
      allocate(pblht(ivert,ihorz))
      allocate(ustar(ivert,ihorz))
      allocate(pblregime(ivert,ihorz))
      allocate(moisture(ivert,ihorz))
      allocate(u10(ivert,ihorz))
      allocate(v10(ivert,ihorz))
      allocate(t2(ivert,ihorz))
c
c     finished allocating memory
c

c
c   Process some header information:
c
      if (dataform.eq.'mm5v1   ') then
         if(itapfrq.ne.0)then
            print*,'Output every',itapfrq,' minutes according to user.'
         else
            itapfrq= nint(mrfv1(302,6))
            if(iprog.eq.2.or.iprog.eq.3)itapfrq=mifv1(5,iprog)*60
         endif
         ptop=mrfv1(1,2)
c
c  Find the highest pressure level to use
c
         kp=num_pres_levs
         do k=num_pres_levs,1,-1
            if(ptop.gt.plevs(k))kp=k-1
         enddo
c
c     Calculate map area           ! simpler than Ken's old version using
c      refrat=mrfv1(1,1)/mrfv1(101,1)
c      yllc=mifv1(106,1)+(0)/refrat
c      yurc=mifv1(106,1)+(ivert-1.)/refrat
c      xllc=mifv1(107,1)+(0)/refrat
c      xurc=mifv1(107,1)+(ihorz-1.)/refrat
c
         yllc=mrfv1(102,1)
         xllc=mrfv1(103,1)
         yurc=mrfv1(104,1)
         xurc=mrfv1(105,1)
         call xyll(xllc,yllc,xlatll,xlonll,mrfv1(4,1),mrfv1(1,1),
     +        mrfv1(2,1),mrfv1(3,1),mifv1(2,1),mifv1(3,1))
         call xyll(xurc,yurc,xlatur,xlonur,mrfv1(4,1),mrfv1(1,1),
     +        mrfv1(2,1),mrfv1(3,1),mifv1(2,1),mifv1(3,1))
         if(mifv1(4,1).eq.1)then  !lambert conic conformal projection
            cproj = 'LCC'       !projection type (lambert conformal=LCC)
            angflg=.true.       !full map projection
            angle1=  mrfv1(6,1)   !standard latitude 1 for GEMPAK (true lat 2 MM5)
            angle2=  mrfv1(3,1)   !central longitude for GEMPAK (pole longitude in MM5)
            angle3=  mrfv1(5,1)   !standard latitude 2 for GEMPAK (true lat 1 MM5)
         elseif(mifv1(4,1).eq.2)then !polar stereographic projection
            cproj = 'STR'       !projection type (polar stereographic=STR)
            angflg=.true.       !full map projection
            angle1= mrfv1(7,1)    !latitude of point of tangency
            angle2= mrfv1(3,1)    !longitude of point of tangency
            angle3= 0.          !not used
         elseif(mifv1(4,1).eq.3)then !Mercator projection
            cproj = 'MER'       !projection type (Mercator=STR)
            angflg=.true.       !full map projection
            angle1= mrfv1(7,1)    !latitude of point of tangency
            angle2= mrfv1(3,1)    !longitude of point of tangency
            angle1= 0.          !not used
            angle3= 0.          !not used
         else
            print*,'This map projection is not supported in MM5, but'
     $           ,' could be used in GEMPAK'
            goto 999
         endif 
         deltan=2.*mrfv1(1,1)     !station spacing(deg) for barnes analysis
     &        /111.19842        !roughly twice the grid spacing
         rewind(20)
      elseif (dataform.eq.'mm5v3   ') then
         rewind(20)
         if(itapfrq.ne.0)then
            print*,'Output every',itapfrq,' minutes according to user.'
         else
            itapfrq= nint(bhr(4,12))
         endif
         ptop=bhr(2,2)/100.
c
c  Find the highest pressure level to use
c
         kp=num_pres_levs
         do k=num_pres_levs,1,-1
            if(ptop.gt.plevs(k))kp=k-1
         enddo
c
c     Calculate map area           ! simpler than Ken's old version using
c      refrat=mrfv1(1,1)/mrfv1(101,1)
c      yllc=mifv1(106,1)+(0)/refrat
c      yurc=mifv1(106,1)+(ivert-1.)/refrat
c      xllc=mifv1(107,1)+(0)/refrat
c      xurc=mifv1(107,1)+(ihorz-1.)/refrat
c
         yllc=bhr(10,1)
         xllc=bhr(11,1)
         yurc=bhr(12,1)
         xurc=bhr(13,1)
         call xyll(xllc,yllc,xlatll,xlonll,bhr(4,1),(bhr(1,1)/1000.),
     +        bhr(2,1),bhr(3,1),bhi(5,1),bhi(6,1))
         call xyll(xurc,yurc,xlatur,xlonur,bhr(4,1),(bhr(1,1)/1000.),
     +        bhr(2,1),bhr(3,1),bhi(5,1),bhi(6,1))
         if(bhi(7,1).eq.1)then  !lambert conic conformal projection
            cproj = 'LCC'       !projection type (lambert conformal=LCC)
            angflg=.true.       !full map projection
            angle1=  bhr(6,1)   !standard latitude 1 for GEMPAK (true lat 2 MM5)
            angle2=  bhr(3,1)   !central longitude for GEMPAK (pole longitude in MM5)
            angle3=  bhr(5,1)   !standard latitude 2 for GEMPAK (true lat 1 MM5)
         elseif(bhi(7,1).eq.2)then !polar stereographic projection
            cproj = 'STR'       !projection type (polar stereographic=STR)
            angflg=.true.       !full map projection
            angle1= bhr(7,1)    !latitude of point of tangency
            angle2= bhr(3,1)    !longitude of point of tangency
            angle3= 0.          !not used
         elseif(bhi(7,1).eq.3)then !Mercator projection
            cproj = 'MER'       !projection type (Mercator=STR)
            angflg=.true.       !full map projection
            angle1= bhr(7,1)    !latitude of point of tangency
            angle2= bhr(3,1)    !longitude of point of tangency
            angle1= 0.          !not used
            angle3= 0.          !not used
         else
            print*,'This map projection is not supported in MM5, but'
     $           ,' could be used in GEMPAK'
            goto 999
         endif 
         deltan=2.*(bhr(1,1)/1000.)     !station spacing(deg) for barnes analysis
     &        /111.19842        !roughly twice the grid spacing
      endif

      print *, 'I = ',ivert,' J = ',ihorz,' K = ', kz
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
        if(ioldkx.ne.ihorz.or.ioldky.ne.ivert) then
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
         CALL GR_MNAV  ( cproj, ihorz, ivert, xlatll, xlonll,
     $        xlatur, xlonur, angle1, angle2, angle3, angflg,
     +        rnvblk, iret )
         IF( iret.ne.0 ) THEN
            write(6,820) iret
 820        FORMAT('  Error building GEMPAK navigation block:',i6)
            STOP
         END IF
         print *,'built GEMPAK navigation block ok for ',cproj
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
         ighdr(1)=0
         ighdr(2)=0
         call gd_cref(outnam,navsz,rnvblk,ianlsz,anlblk,ihdrsz,maxgrd,
     &        igdfln,iret)
         if(iret.ne.0)then
            print*,'gd_cref:  problem creating output gempak file:'
     $           ,outnam
            goto 999
         endif
      endif                     ! end of creating new gempak file or reading old one

c     this was inside the loop of reading data before, why?
      call gr_snav(navsz,rnvblk,iret)      
      call dg_init(0,rnvblk,'first',iret)
      print *,'DG_INIT: return = ',iret
      if (iret .ne. 0) then
         call bail_out('error in dg_init')
      endif
      ntimes = 0
      ndates = 0
      lasttime = -1.
      lastdate = ''

 30   iend = 0
      print *,'opening old precip acuumulation tempfile for reading: ',
     $     TRIM(tempfile)
      print *,'opening new precip acuumulation tempfile writing: ',
     $     TRIM(newtempfile)
      open(unit=22,file=tempfile,form='unformatted')
      open(unit=23,file=newtempfile,form='unformatted',status='unknown')
      if (dataform.eq.'mm5v1   ') then 
         call fill_data_withv1(
     $     uu,vv,tt,qv,qc,qr,cice,sice,pp,rh,rtten,
     $     ww,soilt,
     $     re,rc,prec,
     $     pstx,tg,ter,xmap,dotmap,dotcor,sst,slpcrs,irst,
     $     xlat,xluse,xsnowc,xlon,
     $     prs,dew,wwh,prsd,pdata,ttp,prsp,qvlog,tmp3dp,tmp3ds,
     $     unor,vnor,hlhsl,
     $     pstd,sfp,slp,dotvar,dat,dluse,tmp2d,tot_prec,
     $     tvirt,
     $     maxp,minp,maxpd,minpd,sih,
     $     senhtflx,lathtflx,lonwvrad,shtwvrad,pblht,ustar,
     $     pblregime,moisture,u10,v10,t2,
     $     mrfv1,mifv1,mrfc,mifc,
     $     ivert,ihorz,kz,iprog,
     $     iend
     $     )
      else if (dataform.eq.'mm5v3   ') then 
         call fill_data_withv3(
     $     uu,vv,tt,qv,qc,qr,cice,sice,pp,rh,rtten,
     $     ww,soilt,
     $     re,rc,prec,
     $     pstx,tg,ter,xmap,dotmap,dotcor,sst,slpcrs,irst,
     $     xlat,xluse,xsnowc,xlon,
     $     prs,dew,wwh,prsd,pdata,ttp,prsp,qvlog,tmp3dp,tmp3ds,
     $     unor,vnor,hlhsl,
     $     pstd,sfp,slp,dotvar,dat,dluse,tmp2d,tot_prec,
     $     tvirt,
     $     maxp,minp,maxpd,minpd,sih,
     $     senhtflx,lathtflx,lonwvrad,shtwvrad,pblht,ustar,
     $     pblregime,moisture,u10,v10,t2,
     $     bhi,bhr,bhic,bhrc,
     $     ivert,ihorz,kz,iprog,
     $     ftimes,ntimes,lasttime,chdates,ndates,lastdate,iend
     $     )
      endif
      if(iprog.eq.2.or.iprog.eq.3)then
        ifprs=.true.
        ifsig=.false.
      endif
c
c  Return for next time
      if (iend.eq.1) goto 999
      itime = itime + 1

      if (itime.eq.1)then       ! adjust end points of land use
         do bb = 1,ihorz
            do aa = 1,ivert
c               iaabb = aa+(bb-1)*ivert
               if(aa.eq.ivert) then
                  xluse(aa,bb) = xluse(aa-1,bb)
               endif
               if(bb.eq.ihorz) then
c                  iaabbm1 = aa+(bb-2)*ivert
                  xluse(aa,bb) = xluse(aa,bb-1)
               endif
            enddo
         enddo
c     I don't think you really want to interpolate
c     land use -- categories would lose meaning and value
c         call xtodot2(xluse,dluse,ivert,ihorz,ivertl,ihorzl)
      endif
      print*,' '
c
c  Calculate date/time info and possibly prompt user
c
      if (dataform.eq.'mm5v1   ') then
         if (iprog.eq.5) then
C       INTERP OUTPUT FILE / MMINPUT FILE
            call mdate_split(mifv1(4,5)
     $           ,chdate14,iyear,imonth,iday,ihour,imin,isec)
            write(yymmdd,'(i2.2,i2.2,i2.2)') mod(iyear,100)
     $           ,imonth,iday
            ihr = ihour*100
            iftim=nint(mrfv1(1,5))
            if(mod(nint(mrfv1(1,5)),60).ne.0)then
               print *, 'Date =',mifv1(1,iprog), ' Forecast Time =',
     &              int(mrfv1(1,iprog)/60.),
     $              ' hours and',nint(mrfv1(1,iprog)
     &              - int(mrfv1(1,iprog)/60.)*60),' minutes'
            else
               print *, 'Date =',mifv1(1,iprog), ' Forecast Time =',
     &              nint(mrfv1(1,iprog)/60.0),' hours'
            endif
            print *,'mdate = ',mifv1(4,5),' -> ',chdate14,' -> ',yymmdd
         elseif (iprog.eq.6) then
C     OUTPUT FROM MM5 - - MORE TIME INFO
            idate1=mifv1(13,6)+100*mifv1(14,6)+
     $           mifv1(15,6)*10000+mifv1(16,6)*1000000
            idate2=mifv1(11,6)+100*mifv1(12,6)
            call mdate_split(mifv1(2,6)
     $           ,chdate14,iyear,imonth,iday,ihour,imin,isec)
            write(yymmdd,'(i2.2,i2.2,i2.2)') mod(iyear,100)
     $           ,imonth,iday
            ihr = ihour*100
            if(mod(itapfrq,60).ne.0)then
               write(*,'(a6,1x,i8.8,a1,i4.4,1x,a15,1x,i4,1x,a9,1x,i2
     $              ,1x,a7)')
     &              'Current Date =',idate1,'/'
     $              ,idate2,'Forecast Time =',
     &              int(mrfv1(1,iprog)/60.),
     $              'hours and',nint(mrfv1(1,iprog)
     &              - int(mrfv1(1,iprog)/60.)*60),'minutes'
            else
               write(*,'(a6,1x,i8.8,a1,i4.4,1x,a15,1x,i4,1x,a5)')
     &              'Current Date =',idate1,'/'
     $              ,idate2,'Forecast Time =',
     &              nint(mrfv1(1,iprog)/60.0),'hours'
            endif
            iftim=(itime-1)*itapfrq + int(mrfv1(1,6))
         elseif (iprog.eq.1) then
C     OUTPUT FROM PROGRAM TERRAIN
            yymmdd="010101"
            ihr=0
            ifhr=0
            ifmn=0
            write(gdattm(1),fmt="(a6,a1,i4.4,a1,i3.3,i2.2,'   ')")
     &           yymmdd,'/',ihr,'F',ifhr,ifmn
            goto 321            !can skip all the following calculations for terrain
         elseif (iprog.eq.2.or.iprog.eq.3) then
            call mdate_split(mifv1(1,iprog)
     $           ,chdate14,iyear,imonth,iday,ihour,imin,isec)
            write(yymmdd,'(i2.2,i2.2,i2.2)') mod(iyear,100)
     $           ,imonth,iday
            ihr = ihour*100
            iftim=(itime-1)*itapfrq
            print *,'mdate = ',mifv1(1,iprog),' -> '
     $           ,chdate14,' -> ',yymmdd
         elseif(iprog.eq.4.or.iprog.gt.6.or.iprog.lt.1)then
            print*,'This output is not from interp, mm5, or terrain, ',
     &           'datagrid, or rawins . . .stopping'
            goto 999
         endif
         ifhr  = nint(float(iftim)/60.)
         write(gdattm(1),fmt="(a6,a1,i4.4,a1,i3.3,'   ')")
     &        yymmdd,'/',ihr,'F',ifhr
      else
         print *,'working on dates for mm5v3, itime= ',itime
     $        ,' ntimes= ',ntimes
     $        ,' ftimes(itime)= ',ftimes(itime)
     $        ,' ftimes(ntimes)= ',ftimes(ntimes)
     $        ,' ftimes= ',ftimes
     $        ,' ndates= ',ndates
     $        ,' chdates(itime)= ',chdates(itime)
     $        ,' chdates(ntimes)= ',chdates(ntimes)
     $        ,' chdates= ',chdates
         if (iprog.ge.11) then
            write(gdattm(1),'(i2.2,i2.2,i2.2,"/",i2.2,"00F",i3.3)')
     $           (bhi(5,iprog)-bhi(5,iprog)/100*100),bhi(6,iprog)
     $           ,bhi(7,iprog),bhi(8,iprog),nint(ftimes(ntimes)/60.)
         elseif (iprog.eq.2 .or. iprog.eq.3) then
            gdattm(1) = chdates(itime)(3:4)//
     $           chdates(itime)(6:7)//
     $           chdates(itime)(9:10)//"/"//
     $           chdates(itime)(12:13)//"00F000"
         endif
c         gdattm(1) = '010618/1200F000'
      endif
      print *,'gdattm(1) = ',gdattm(1)

      tmpa='n'
      if(quiet)then             !don't prompt user
         tmpa='n'
      elseif(igemfrq.ne.0)then  !also don't prompt user
         if((iprog.eq.6.and.mod(itapfrq*(itime-1),igemfrq).eq.0).or.
     &        (iprog.eq.5.and.mod(nint(mrfv1(1,5)),   igemfrq).eq.0).or.
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

      DO 700 k = kz,1,-1        ! reverse order because sfp wants pp(i,j,kz)
         maxp(k)=0.
         minp(k)=10000.
         if (dataform.eq.'mm5v1   ') then
            sih(k) = mrfv1(101+k,iprog)
            print *,'mm5v1 got sih(',k,') = ',sih(k)
         else
            print *,'mm5v3 got sih(',k,') = ',sih(k)
         endif
         if (k.eq.kz) then
c     need all pstd's pre-calculated with standard pstx
            DO j=1,ihorz
               DO i=1,ivert
                  aa = MIN(i,ivert-1)
                  bb = MIN(j,ihorz-1)
                  cc = MAX(i-1,1)
                  dd = MAX(j-1,1)
                  pstd(i,j)= 0.25*(pstx(aa,bb)+pstx(cc,bb)
     &                 +pstx(aa,dd)+pstx(cc,dd))
               END DO
            END DO
         endif
         DO 650 j = 1,ihorz
            DO 625 i = 1,ivert
               if (k.eq.kz .and. i.le.ivert-1 .and. j.le.ihorz-1) then
C     calculate tot_prec and convert units
C     precip array.  Adjust output from CENTIMETERS to MILLIMETERS!!!!!!!
                  tot_prec(i,j) = (10.0*rc(i,j))+(10.0*re(i,j))
                  re(i,j) = 10.0*re(i,j)
                  rc(i,j) = 10.0*rc(i,j)
                  prec_3hr(i,j)=0.0
                  re_3hr(i,j)=0.0
                  rc_3hr(i,j)=0.0
               endif
C     
C     Decouple Data with original pstx
c       (LEAVE pstx units = kPa alone!)
C     
               if (ifprs.or.ifsig) then
                  if (i.le.(ivert-1) .and. j.le.(ihorz-1)) then
                     if (dataform.eq.'mm5v1   ') then
                        if (k.eq.kz) then
                           ww(i,j,(k+1)) = ww(i,j,(k+1))/pstx(i,j)
                           ww(i,j,k)   = ww(i,j,k)/pstx(i,j)
                        else
                           ww(i,j,k) = ww(i,j,k)/pstx(i,j)
                        endif
                     endif
                     wwh(i,j,k)  = (ww(i,j,k)+ww(i,j,(k+1)))/2.
                  endif
               endif
               if (k.eq.kz) then
                  if (i.le.ivert-1 .and. j.le.ihorz-1) then
                     if (dataform.eq.'mm5v1   ') then
                        tt(i,j,k) = tt(i,j,k)/pstx(i,j)
                        qv(i,j,k) = qv(i,j,k)/pstx(i,j)
                        qc(i,j,k) = qc(i,j,k)/pstx(i,j)
                        qr(i,j,k) = qr(i,j,k)/pstx(i,j)
                        cice(i,j,k) = cice(i,j,k)/pstx(i,j)
                        sice(i,j,k) = sice(i,j,k)/pstx(i,j)
                        pp(i,j,k) = 0.01*pp(i,j,k)/pstx(i,j) ! to mb
                     else
                        pp(i,j,k) = 0.01*pp(i,j,k) ! Pa to mb
                     endif
C     Find virtual temperature on lowest sigma level:
                     tvirt(i,j) = tt(i,j,k)*(1.+qv(i,j,k)/.62197)
     $                    / (1.+qv(i,j,k))
                  endif
                  if (dataform.eq.'mm5v1   ') then
                     uu(i,j,k) = uu(i,j,k)/pstd(i,j)
                     vv(i,j,k) = vv(i,j,k)/pstd(i,j)
                  endif
               else
c     convert back to original pstx
                  if (i.le.ivert-1 .and. j.le.ihorz-1) then
                     if (dataform.eq.'mm5v1   ') then
                        tt(i,j,k) = tt(i,j,k)/pstx(i,j)
                        qv(i,j,k) = qv(i,j,k)/pstx(i,j)
                        qc(i,j,k) = qc(i,j,k)/pstx(i,j)
                        qr(i,j,k) = qr(i,j,k)/pstx(i,j)
                        cice(i,j,k) = cice(i,j,k)/pstx(i,j)
                        sice(i,j,k) = sice(i,j,k)/pstx(i,j)
                        pp(i,j,k) = 0.01*pp(i,j,k)/pstx(i,j)
                     else
                        pp(i,j,k) = 0.01*pp(i,j,k)
                     endif
                  endif
c     convert back to original pstd
                  if (dataform.eq.'mm5v1   ') then
                     uu(i,j,k) = uu(i,j,k)/pstd(i,j)
                     vv(i,j,k) = vv(i,j,k)/pstd(i,j)
                  endif
               endif
               if (i.le.ivert-1 .and. j.le.ihorz-1) then
c     arrays for calculating geopotential height
                  if (k.le.kp) then
                     prsp(i,j,k) = plevs(k)
                  endif
                  qvlog(i,j,k)=log(max(qv(i,j,k),1.e-15))
               endif
C     
C     Now use pstx and pstd in mb for vapor formulas
C     Only need to do when k = kz
               if (k.eq.kz) then
                  if (i.le.ivert-1 .and. j.le.ihorz-1) then
c     converting pstx units at k = kz
ccccc                     pstx(i,j)=pstx(i,j)*10.0

                     sfp(i,j)=10.*pstx(i,j)+ptop+pp(i,j,k)
                     wmixrat=qv(i,j,k)
                     tvll=tt(i,j,k) * (1.+0.608*wmixrat )
                     if (i.eq.20 .and. j.eq.20) then
                        print *,'ptop,sih(k)',ptop,sih(k)
                        print *,'pstx,pp',(10.*pstx(i,j)),pp(i,j,k)
                     endif
                     beta=287.04/9.81*
     &                    log((10.*pstx(i,j)+ptop+pp(i,j,k))/
     &                    (sih(k)*10.*pstx(i,j)+ptop+pp(i,j,k)))
                     gpll=ter(i,j)+beta*tvll/(1.-.5*beta*0.0065)
                     tvavg = tvll + 0.0065 * (gpll - .5 * ter(i,j))

                     slp(i,j)=(10.*pstx(i,j)+ptop+pp(i,j,k))*
     &                    exp(9.81*ter(i,j)/287.04/tvavg)
                  endif
c     converting pstd units at k = kz
ccccc                  pstd(i,j)=pstd(i,j)*10.0
               endif
C     
C     The following fields require pstx/pstd in mb also
C     
               if (i.le.ivert-1 .and. j.le.ihorz-1) then

                  prs(i,j,k)=10.*pstx(i,j)*sih(k)+ptop+pp(i,j,k)
                  if(prs(i,j,k).gt.maxp(k))maxp(k)=prs(i,j,k)
                  if(prs(i,j,k).lt.minp(k))minp(k)=prs(i,j,k)

                  q = MAX(qv(i,j,k),1.e-15)
                  t = tt(i,j,k)
                  tcels = t-273.15
                  p = 10.*pstx(i,j)*sih(k)+ptop+pp(i,j,k)
                  e = q*p/(eps+(1.-eps)*q)
c     dewpoint
                  dew(i,j,k) = 5418.12/(19.84659-alog(e/6.112))
c     relative humidity
                  esat = 6.112*exp((17.67*tcels)/(tcels+243.5))
                  qsat = eps*esat/p 
                  rh(i,j,k) = q/qsat*100
               endif
 625        continue              ! i loop
 650        continue            ! j loop
         if(ifprs)then
            print *,'calculating ttp for k level ',k
            maxpd(k)=0.
            minpd(k)=10000.
c            ik = 1 + (k-1)*ivert*ihorz
            call xtodot(prs(:,:,k),dotvar,ivert,ihorz,ivertl,ihorzl)
            DO j=1,ihorz
               DO i=1,ivert
                  prsd(i,j,k)=dotvar(i,j)
                  if(prsd(i,j,k).gt.maxpd(k))maxpd(k)=prsd(i,j,k)
                  if(prsd(i,j,k).lt.minpd(k))minpd(k)=prsd(i,j,k)
               END DO           ! i loop
            END DO              ! j loop
            call pintp (tt,ttp,prs,tvirt,plevs,maxp,minp,
     $           ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,icross,3)
         endif 
 700  continue                    ! k loop

c      i5j5    = 5 + 4*ivert
c      i5j5k25 = i5j5 + 24*ivert*ihorz
c      i5j5k5  = i5j5 + 4*ivert*ihorz
      PRINT*,'uu(5,5,25) = ',uu(5,5,25) 
      PRINT*,'uu(5,5,5) = ',uu(5,5,5) 
      PRINT*,'pstx(5,5) = ',pstx(5,5) 
      PRINT*,'pstx(5,5) = ',pstx(5,5)/10. 
      PRINT*,'pstd(5,5) = ',pstd(5,5) 
      PRINT*,'vv(5,5,25) = ',vv(5,5,25) 
      PRINT*,'vv(5,5,5) = ',vv(5,5,5) 
      print*,'tt(5,5,5) = ', tt(5,5,5)
      print*,'tt(5,5,25) = ', tt(5,5,25)
      print*,'prs(5,5,5) = ', prs(5,5,5)
      print*,'prs(5,5,25) = ', prs(5,5,25)
      PRINT*,'slp(5,5) = ',slp(5,5) 
      PRINT*,'sfp(5,5) = ',sfp(5,5) 
      PRINT*,'dew(5,5,5) = ',dew(5,5,5) 
      PRINT*,'dew(5,5,25) = ',dew(5,5,25) 
      PRINT*,'rh(5,5,5) = ',rh(5,5,5) 
CCCCCC
C Calculate Rain totals for mmout files
c
c     unless this is the initialization time, compute the precip
c
      if ((dataform.eq.'mm5v1   '.and.iprog.eq.6) .or.
     $     (dataform.eq.'mm5v3   '.and.iprog.eq.11)) then
         if (dataform.eq.'mm5v1   ') then
            idomid = mifv1(101,1)
            output_hr = NINT(mrfv1(1,6)/60.0)
         else
            idomid = bhi(13,1)
            output_hr = nint(ftimes(ntimes)/60.)
         endif
         print '("idomid = ",i4," output_hr = ",i4)',idomid,output_hr
         if (output_hr.le.3 .or. (output_hr.eq.special_hour .and.
     $        idomid.eq.dom_special) ) then
c     print 1027,mifv1(101,1),output_hr,dom_special,special_hour
c 1027    format('skipping precip calc.s for domain',i2
c     $        ,' output hour = ',i3,'. dom_special = ',i2
c     $        ' and special_hour = ',i3)
            print *,'skipping precip calc.s for domain',idomid
     $           ,' output hour = ',output_hr,'. dom_special = ',
     $           dom_special,' and special_hour = ',special_hour
            continue
         else
            if (ntimes .gt. 1) then
               if (mod(output_hr,3) .eq. 0) then
                  print *,'calculating 3-hour precip by ',
     $                 output_hr,' minus ',last_output_hr
                  prec_3hr = tot_prec -
     &                 last_output_hr_tot_prec
                  re_3hr = re -
     &                 last_output_hr_exp_prec
                  rc_3hr = rc -
     &                 last_output_hr_con_prec
               else 
                  print *,'skipping output_hr ',output_hr,' since not '
     $                 ,' divisible by 3, last was ',last_output_hr
               endif
            else 
               print*,'reading the last output from tempfile'
               imiss = 0
               read(22,end=777,err=777) last_output_hr,livert,lihorz
               if (livert.ne.ivert .or. lihorz.ne.ihorz) then
                  print *,'we have a problem in the precip file'
                  print *,livert,' vs other ivert ',ivert
                  print *,lihorz,' vs other ihorz ',ihorz
                  call bail_out('error with precip tempfile dimensions')
               endif
               goto 778
 777           imiss = 1
 778           if ((output_hr-last_output_hr) .ne. 3 .or.
     $              imiss .eq. 1) then
                  print*,'WARNING: missing the last output file precip.'
     $                 //' Accumulations will be off!',last_output_hr
     $                 ,livert,lihorz,output_hr,imiss
               else
                  read(22) last_output_hr_tot_prec
                  read(22) last_output_hr_exp_prec
                  read(22) last_output_hr_con_prec
                  do j=1,ihorz-1            
                     do i=1,ivert-1
                        prec_3hr(i,j)=tot_prec(i,j)-
     &                       last_output_hr_tot_prec(i,j)
                        re_3hr(i,j)=re(i,j)-
     &                       last_output_hr_exp_prec(i,j)
                        rc_3hr(i,j)=rc(i,j)-
     &                       last_output_hr_con_prec(i,j)
                     enddo
                  enddo
               endif
            endif
c     write the current total precip grid to the newtempfile      
            close (22)
            write(23) output_hr,ivert,ihorz
            write(23) tot_prec
            write(23) re
            write(23) rc
            close(23)
         endif
      endif
      if (output_hr .ge. 3 .and. mod(output_hr,3) .eq. 0) then
         print *,'saving precip for output_hr ',output_hr
         last_output_hr_tot_prec = tot_prec
         last_output_hr_exp_prec = re
         last_output_hr_con_prec = rc
         last_output_hr = output_hr
      endif
c
c  Calculate the height above ground level of lowest half sigma level
c
c     if(.not.ifsig.and.iout.eq.1)then
      if(.not.ifsig.and.itime.eq.1)then
c         ik = 1 + (kz-2)*ivert*ihorz
         call ghtcalc(tt(:,:,kz-1),prs(:,:,kz-1),qv(:,:,kz-1),sfp,ter
     $        ,hlhsl,ivert,ihorz,2,ivert,ihorz,2,2)
         do j=1,ihorz-icross
            do i=1,ivert-icross
               hlhsl(i,j,1)=hlhsl(i,j,2)-ter(i,j)
               hlhsltot=hlhsltot+hlhsl(i,j,2)-ter(i,j)
            enddo
         enddo
         zagllhsl=hlhsltot/((ivert-1)*(ihorz-1))
         if(zagllhsl.ge.40.)then
            izagllhsl=10*nint(zagllhsl/10.)
         elseif(zagllhsl.ge.5.)then
            izagllhsl=5*nint(zagllhsl/5.)
         else
            izagllhsl=nint(zagllhsl)
         endif
      endif
c

c
c   Dump out arrays for gempak
c

 678  continue
c  
c  3-d fields first:
c
      do kk = 1,kz
         do ii = 1,ihorz
            do jj = 1,ivert
c               ijk = ii + (jj-1)*ihorz + (kk-1)*ihorz*ivert
c               jik = jj + (ii-1)*ivert + (kk-1)*ivert*ihorz
               unor(ii,jj,kk) = uu(jj,ii,kk)
               vnor(ii,jj,kk) = vv(jj,ii,kk)
            enddo
         enddo
         if (kk.eq.1) then
            ii = int((400-1)/ivert)+1
            jj = 400-(ii-1)*ivert
            jik = jj + (ii-1)*ivert
            ijk = ii + (jj-1)*ihorz
            print *,'uu(',jik,') vs unor(',ijk,') ',
     $           uu(jj,ii,kk),unor(ii,jj,kk)
         endif
c         ik = 1 + (kk-1)*ivert*ihorz
         call dg_nrel(unor(:,:,kk),vnor(:,:,kk)
     $        ,unor(:,:,kk),vnor(:,:,kk),iret)
         if (iret .ne. 0) then
            print *,'dg_nrel iret = ',iret
            call bail_out('bad grid to north-relative conversion')
         endif
      enddo
c      ii = int((400-1)/ivert)+1
c      jj = 400-(ii-1)*ivert
c      jik = jj + (ii-1)*ivert
c      ijk = ii + (jj-1)*ihorz
c      print *,'uu(',jik,') vs unor(',ijk,') ',uu(jik),unor(ijk)
c      ijk400 = ijk
      if(iprog.eq.2.or.iprog.eq.3)then !datagrid or rawins output
         do kk = 1,kz
            if (dataform.eq.'mm5v1   ') then
               ilvl(1)=mifv1(101+kk,iprog)
            else
               ilvl(1)=nint(sih(kk)/100.)
            endif
c            ik = 1 + (kk-1)*ivert*ihorz

            call write_data(igdfln,unor(:,:,kk),ihorz,ivert,ighdr,gdattm
     $           ,ilvl,1,
     $           'UWND        ',rewrite,ipktyp,iprecise,0)
            call write_data(igdfln,vnor(:,:,kk),ihorz,ivert,ighdr,gdattm
     $           ,ilvl,1,
     $           'VWND        ',rewrite,ipktyp,iprecise,0)
            call write_data(igdfln,uu(:,:,kk),ivert,ihorz,ighdr,gdattm
     $           ,ilvl,1,
     $           'UREL        ',rewrite,ipktyp,iprecise,1)
            call write_data(igdfln,vv(:,:,kk),ivert,ihorz,ighdr,gdattm
     $           ,ilvl,1,
     $           'VREL        ',rewrite,ipktyp,iprecise,1)
            call write_data(igdfln,tt(:,:,kk),ivert,ihorz,ighdr,gdattm
     $           ,ilvl,1,
     $           'TMPK        ',rewrite,ipktyp,iprecise,1)
            call write_data(igdfln,qv(:,:,kk),ivert,ihorz,ighdr,gdattm
     $           ,ilvl,1,
     $           'RELH        ',rewrite,ipktyp,iprecise,1)
            call write_data(igdfln,pp(:,:,kk),ivert,ihorz,ighdr,gdattm
     $           ,ilvl,1,
     $           'HGHT        ',rewrite,ipktyp,iprecise,1)
         enddo
         goto 321               ! plot out surface fields, same as for mmout and interp
      endif    !end of stuff for datagrid output

      do 421 inum = 1,num_3d_outputs
         if(ifprs)then
            ivcor=1
            iflip = 1
            parm=prs_fld(inum)
            print *,'calling pintp for ',parm
            if(sig_fld(inum)(1:4).eq.'UWND' .or.
     $           prs_fld(inum)(1:4).eq.'UWND')then 
               call pintp (unor,pdata,prsd,tvirt,plevs,maxpd,minpd,
     $              ihorz,ivert,kz,ihorz,ivert,kz,kp,kz,idot,0)
               imakedot = 0
               iflip = 0
            else if(sig_fld(inum)(1:4).eq.'UREL' .or.
     $              prs_fld(inum)(1:4).eq.'UREL')then 
               call pintp (uu,pdata,prsd,tvirt,plevs,maxpd,minpd,
     $              ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,idot,0)
               imakedot = 0
            else if(sig_fld(inum)(1:4).eq.'VWND' .or.
     $              sig_fld(inum)(1:4).eq.'VWND')then 
               call pintp (vnor,pdata,prsd,tvirt,plevs,maxpd,minpd,
     $              ihorz,ivert,kz,ihorz,ivert,kz,kp,kz,idot,0)
               imakedot = 0
               iflip = 0
            else if(sig_fld(inum)(1:4).eq.'VREL' .or.
     $              prs_fld(inum)(1:4).eq.'VREL')then 
               call pintp (vv,pdata,prsd,tvirt,plevs,maxpd,minpd,
     $              ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,idot,0)
               imakedot = 0
            else if(sig_fld(inum)(1:4).eq.'WWND' .or.
     $              prs_fld(inum)(1:4).eq.'WWND')then 
               call pintp (wwh,tmp3ds,prsd,tvirt,plevs,maxpd,minpd,
     $              ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,icross,0)
               imakedot = 1
            else if(sig_fld(inum)(1:4).eq.'TMPK' .or.
     $              prs_fld(inum)(1:4).eq.'TMPK')then 
c               ttp has already been calculated from tt
               imakedot = 3
            else if(sig_fld(inum)(1:4).eq.'DWPK' .or.
     $              prs_fld(inum)(1:4).eq.'DWPK')then 
               call pintp (dew,tmp3ds,prs,tvirt,plevs,maxp,minp,
     $              ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,icross,3)
               imakedot = 1
            else if(prs_fld(inum)(1:4).eq.'HGHT')then 
               imakedot = 1
               if(ighttype.eq.2)then
                  call ghtcalc(tt,prs,qv,sfp,ter,pdata
     $                 ,ivert,ihorz,kz,ivert,ihorz,kz,2)
                  call pintp(pdata,tmp3ds,prs,tvirt,plevs,maxp,minp
     $                 ,ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,icross,4)
               elseif(ighttype.eq.1)then
c     :skip following loop to do interpolation of q rather than log q
                  call pintp(qvlog,tmp3dp,prs,tvirt,plevs,maxp,minp
     $                 ,ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,icross,1)
c     :skip following loop to do interpolation of q rather than log q
                  do k=1,kp
                     do j=1,ihorz-1
                        do i=1,ivert-1
                           tmp3dp(i,j,k)=exp(tmp3dp(i,j,k))
                        enddo
                     enddo
                  enddo
                  call ghtcalc(ttp,prsp,tmp3dp,sfp,ter,tmp3ds
     $                 ,ivert,ihorz,kp,ivert,ihorz,kz,1)
               endif
            else if(sig_fld(inum)(1:4).eq.'MIXR' .or.
     $              prs_fld(inum)(1:4).eq.'MIXR')then 
               call pintp (qv,tmp3ds,prs,tvirt,plevs,maxp,minp
     $              ,ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,icross,0)
               imakedot = 2000 ! use xtodot2 and mult by 1000.
            else if(sig_fld(inum)(1:4).eq.'QCLD' .or.
     $              prs_fld(inum)(1:4).eq.'QCLD')then 
               call pintp (qc,tmp3ds,prs,tvirt,plevs,maxp,minp
     $              ,ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,icross,0)
               imakedot = 2000 ! use xtodot2 and mult by 1000.
            else if(sig_fld(inum)(1:4).eq.'QRAI' .or.
     $              prs_fld(inum)(1:4).eq.'QRAI')then 
               call pintp (qr,tmp3ds,prs,tvirt,plevs,maxp,minp
     $              ,ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,icross,0)
               imakedot = 2000 ! use xtodot2 and mult by 1000.
            else if(sig_fld(inum)(1:4).eq.'QICE' .or.
     $              prs_fld(inum)(1:4).eq.'QICE')then 
               call pintp (cice,tmp3ds,prs,tvirt,plevs,maxp,minp
     $              ,ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,icross,0)
               imakedot = 2000 ! use xtodot2 and mult by 1000.
            else if(sig_fld(inum)(1:4).eq.'QSNO' .or.
     $              prs_fld(inum)(1:4).eq.'QSNO')then 
               call pintp (sice,tmp3ds,prs,tvirt,plevs,maxp,minp
     $              ,ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,icross,0)
               imakedot = 2000 ! use xtodot2 and mult by 1000.
            else if(sig_fld(inum)(1:4).eq.'RTEN' .or.
     $              prs_fld(inum)(1:4).eq.'RTEN')then 
               call pintp (rtten,tmp3ds,prs,tvirt,plevs,maxp,minp
     $              ,ivert,ihorz,kz,ivert,ihorz,kz,kp,kz,icross,0)
               imakedot = 2    ! use xtodot2
            endif ! end of checking names of the pressure fields
            do kk = 1,kp
               ilvl(1)=int(plevs(kk))
c               ik = 1 + (kk-1)*ivert*ihorz
               if (imakedot.eq.1) then
                  call xtodot(tmp3ds(:,:,kk),pdata(:,:,kk),ivert,ihorz,
     $                 ivertl,ihorzl)
               else if (imakedot.eq.2000) then
                  do j = 1,ihorz
                     do i = 1,ivert
                        tmp3ds(i,j,kk) = 1000. * tmp3ds(i,j,kk) ! change to g/kg
                     enddo
                  enddo
                  call xtodot2(tmp3ds(:,:,kk),pdata(:,:,kk),ivert,ihorz,
     $                 ivertl,ihorzl)
               else if (imakedot.eq.2) then
                  call xtodot2(tmp3ds(:,:,kk),pdata(:,:,kk),ivert,ihorz,
     $                 ivertl,ihorzl)
               else if (imakedot.eq.3) then
                  call xtodot(ttp(:,:,kk),pdata(:,:,kk),ivert,ihorz,
     $                 ivertl,ihorzl)
               endif
               call write_data(igdfln,pdata(:,:,kk),ivert,ihorz,ighdr,
     $              gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise
     $              ,iflip)
            enddo ! end of loop through all pressure levels
         endif

         if(ifsig)then
            ivcor=4
            parm=sig_fld(inum)
            do kk = 1,kz
               ilvl(1) = nint(sih(kk)*10000.0)
c               ik = 1 + (kk-1)*ivert*ihorz
               if(sig_fld(inum)(1:4).eq.'UREL' .or.
     $              prs_fld(inum)(1:4).eq.'UREL')then 
                  if (kk.eq.1) print*,'Writing UREL on sigma levels -'
                  call write_data(igdfln,uu(:,:,kk),ivert,ihorz,ighdr,
     $                 gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,1)
                  imakedot = 0
               else if(sig_fld(inum)(1:4).eq.'VREL' .or.
     $                 prs_fld(inum)(1:4).eq.'VREL')then 
                  if (kk.eq.1) print*,'Writing VREL on sigma levels - '
                  call write_data(igdfln,vv(:,:,kk),ivert,ihorz,ighdr,
     $                 gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,1)
                  imakedot = 0
               else if(sig_fld(inum)(1:4).eq.'UWND' .or.
     $              sig_fld(inum)(1:4).eq.'UWND')then 
                  if (kk.eq.1) print*,'Writing U on sigma levels -',
     $                 uu(10,10,10),'-->',unor(10,10,10),'M/S'
                  call write_data(igdfln,unor(:,:,kk),ihorz,ivert,ighdr,
     $                 gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,0)
                  imakedot = 0
               else if(sig_fld(inum)(1:4).eq.'VWND' .or.
     $                 sig_fld(inum)(1:4).eq.'VWND')then 
                  if (kk.eq.1) print *, 'Writing V on sigma levels - '
                  call write_data(igdfln,vnor(:,:,kk),ihorz,ivert,ighdr,
     $                 gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,0)
                  imakedot = 0
               else if(sig_fld(inum)(1:4).eq.'WWND' .or.
     $                 prs_fld(inum)(1:4).eq.'WWND')then 
                  call xtodot(wwh(:,:,kk),dotvar,ivert,ihorz,ivertl
     $                 ,ihorzl)
                  if (kk.eq.1) print *, 'Writing W on sigma levels - '
                  imakedot = 1
               else if(sig_fld(inum)(1:4).eq.'TMPK' .or.
     $                 prs_fld(inum)(1:4).eq.'TMPK')then 
                  call xtodot(tt(:,:,kk),dotvar,ivert,ihorz,ivertl
     $                 ,ihorzl)
                  if (kk.eq.1) print *, 'Writing T on sigma levels - '
                  imakedot = 1
               else if(sig_fld(inum)(1:4).eq.'DWPK' .or.
     $                 prs_fld(inum)(1:4).eq.'DWPK')then 
                  call xtodot(dew(:,:,kk),dotvar,ivert,ihorz,ivertl
     $                 ,ihorzl)
                  if (kk.eq.1) print *,
     $                 'Writing dewpoint on sigma levels - '
                  imakedot = 1
               else if(sig_fld(inum)(1:4).eq.'PRES')then 
                  if (kk.eq.1) print *,
     $                 'Writing Pressure on sigma levels - '
                  call xtodot(prs(:,:,kk),dotvar,ivert,ihorz,ivertl
     $                 ,ihorzl)
                  imakedot = 1
               else if(sig_fld(inum)(1:4).eq.'QCLD' .or.
     $                 prs_fld(inum)(1:4).eq.'QCLD')then 
                  if (kk.eq.1) print *,
     $                 'Writing Cloud Water on sigma levels - '
                  call xtodot2(qc(:,:,kk),dotvar,ivert,ihorz,ivertl
     $                 ,ihorzl)
                  imakedot = 1000 ! use xtodot2 and mult by 1000.
               else if(sig_fld(inum)(1:4).eq.'MIXR' .or.
     $                 prs_fld(inum)(1:4).eq.'MIXR')then 
                  if (kk.eq.1) print *,
     $                 'Writing Mixing Ratio on sigma levels - '
                  call xtodot2(qv(:,:,kk),dotvar,ivert,ihorz,ivertl
     $                 ,ihorzl)
                  imakedot = 1000 ! use xtodot2 and mult by 1000.
               else if(sig_fld(inum)(1:4).eq.'QRAI' .or.
     $                 prs_fld(inum)(1:4).eq.'QRAI')then 
                  if (kk.eq.1) print *,
     $                 'Writing Rain Water on sigma levels - '
                  call xtodot2(qr(:,:,kk),dotvar,ivert,ihorz,ivertl
     $                 ,ihorzl)
                  imakedot = 1000 ! use xtodot2 and mult by 1000.
               else if(sig_fld(inum)(1:4).eq.'QICE' .or.
     $                 prs_fld(inum)(1:4).eq.'QICE')then 
                  if (kk.eq.1) print *,
     $                 'Writing Ice Mixing on sigma levels - '
                  call xtodot2(cice(:,:,kk),dotvar,ivert,ihorz,
     $                 ivertl,ihorzl)
                  imakedot = 1000 ! use xtodot2 and mult by 1000.
               else if(sig_fld(inum)(1:4).eq.'QSNO' .or.
     $                 prs_fld(inum)(1:4).eq.'QSNO')then 
                  if (kk.eq.1) print *,
     $                 'Writing snow water mixing ratio(g/kg) - '
                  call xtodot2(sice(:,:,kk),dotvar,ivert,ihorz,
     $                 ivertl,ihorzl)
                  imakedot = 1000 ! use xtodot2 and mult by 1000.
               else if(sig_fld(inum)(1:4).eq.'RTEN' .or.
     $                 prs_fld(inum)(1:4).eq.'RTEN')then 
                  if (kk.eq.1) print *,
     $                 'Writing Radiative Tendency(K/day) - '
                  call xtodot2(rtten(:,:,kk),dotvar,ivert,ihorz,
     $                 ivertl,ihorzl)
                  imakedot = 1
               endif
               if (imakedot.eq.1) then
                  call write_data(igdfln,dotvar,ivert,ihorz,ighdr,
     $                 gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,1)
               else if (imakedot.eq.1000) then
                  do j = 1,ihorz
                     do i = 1,ivert
                        dotvar(i,j) = 1000. * dotvar(i,j)
                     enddo
                  enddo
                  call write_data(igdfln,dotvar,ivert,ihorz,ighdr,
     $                 gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,1)
               endif
            enddo
         endif ! end of sigma if

c         if(.not.ifsig.and.iprog.ne.2)then
c            ivcor=1514227532
c            ik = 1 + (kz-1)*ivert*ihorz
c            ilvl(1) = izagllhsl
c            if(sig_fld(inum)(1:4).eq.'UWND' .or.
c     $           sig_fld(inum)(1:4).eq.'UREL' .or.
c     $           prs_fld(inum)(1:4).eq.'UWND' .or.
c     $           prs_fld(inum)(1:4).eq.'UREL')then 
c               call write_data(igdfln,uu(ik),ivert,ihorz,ighdr,
c     $              gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,1)
c.....etc.
c        endif
 421  continue                  ! end of loop for 3-d output fields
c
c  Now do surface fields
c
 321  ivcor=0
      ilvl(1) = 0
      
      do 521 inum = num_3d_outputs+1,num_3d_outputs+num_2d_outputs
         if(iprog.gt.0)then     !mm5 or interp output
            parm=prs_fld(inum)
            iusedotvar = 1
            if(sig_fld(inum)(1:4).eq.'PRES' .or.
     $           prs_fld(inum)(1:4).eq.'PRES')then 
               print *, 'Writing surface Pressure'
               call xtodot(sfp,dotvar,ivert,ihorz,ivertl,ihorzl)
            else if(sig_fld(inum)(1:4).eq.'PMSL' .or.
     $              prs_fld(inum)(1:4).eq.'PMSL')then 
c               iusedotvar = 5
               print *, 'Writing sea level Pressure'
               call xtodot(slp,dotvar,ivert,ihorz,ivertl,ihorzl)
            else if(sig_fld(inum)(1:4).eq.'SLPC' .or.
     $              prs_fld(inum)(1:4).eq.'SLPC')then 
               print *, 'Writing sea level Pressure cross points'
               call xtodot(slpcrs,dotvar,ivert,ihorz,ivertl,ihorzl)
            else if(sig_fld(inum)(1:4).eq.'T2  ') then
               print *, 'Writing T2 as though at sigma = 1- '
               parm='T2  '//sig_fld(inum)(5:12)
               call xtodot(t2,dotvar,ivert,ihorz,ivertl,ihorzl)
            else if(sig_fld(inum)(1:4).eq.'U10 ') then
               print *, 'Writing U10 as though at sigma = 1- '
               iusedotvar = 6
               parm='U10 '//sig_fld(inum)(5:12)
            else if(sig_fld(inum)(1:4).eq.'V10 ') then
               print *, 'Writing V10 as though at sigma = 1- '
               iusedotvar = 7
               parm='V10 '//sig_fld(inum)(5:12)
            else if(sig_fld(inum)(1:4).eq.'TMPK' .or.
     $              sig_fld(inum)(1:4).eq.'TGRD') then
               print *, 'Writing Ground Temp as though at sigma = 1- '
               parm='TMPK'//sig_fld(inum)(5:12)
               call xtodot(tg,dotvar,ivert,ihorz,ivertl,ihorzl)
            else if(prs_fld(inum)(1:4).eq.'TMPK' .or.
     $              prs_fld(inum)(1:4).eq.'TGRD') then
               print *, 'Writing Ground Temp as though at sigma = 1- '
               parm='TMPK'//prs_fld(inum)(5:12)
               call xtodot(tg,dotvar,ivert,ihorz,ivertl,ihorzl)
            else if(sig_fld(inum)(1:4).eq.'PTOT' .or.
     $              prs_fld(inum)(1:4).eq.'PTOT')then 
               call xtodot2(tot_prec,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing accumulated Ppt. -'
            else if(sig_fld(inum)(1:4).eq.'CTOT' .or.
     $              prs_fld(inum)(1:4).eq.'CTOT')then 
               call xtodot2(rc,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing accumulated Ppt. -'
            else if(sig_fld(inum)(1:4).eq.'ETOT' .or.
     $              prs_fld(inum)(1:4).eq.'ETOT')then 
               call xtodot2(re,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing accumulated Ppt. -'
            else if(sig_fld(inum)(1:4).eq.'P03M' .or.
     $              prs_fld(inum)(1:4).eq.'P03M')then 
               call xtodot2(prec_3hr,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing accumulated 3hr Ppt. -'
            else if(sig_fld(inum)(1:4).eq.'E03M' .or.
     $              prs_fld(inum)(1:4).eq.'E03M')then 
               call xtodot2(re_3hr,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing 3hr explicit Ppt. -'
            else if(sig_fld(inum)(1:4).eq.'C03M' .or.
     $              prs_fld(inum)(1:4).eq.'C03M')then 
               call xtodot2(rc_3hr,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing 3hr convective Ppt. -'
            else if(sig_fld(inum)(1:4).eq.'HGHT' .or.
     $              prs_fld(inum)(1:4).eq.'HGHT')then 
               print *, 'Writing Terrain - '
               call xtodot2(ter,dotvar,ivert,ihorz,ivertl,ihorzl)
            else if(sig_fld(inum)(1:4).eq.'MPDT' .or.
     $              prs_fld(inum)(1:4).eq.'MPDT')then 
               iusedotvar = 2
               print*,'Writing map scale factor on dot points -'
            else if(sig_fld(inum)(1:4).eq.'MPCR' .or.
     $              prs_fld(inum)(1:4).eq.'MPCR')then 
               call xtodot(xmap,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing map scale factor on cross points -'
            else if(sig_fld(inum)(1:3).eq.'SST' .or.
     $              prs_fld(inum)(1:3).eq.'SST')then 
               call xtodot(sst,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing sea-surface temperature-'
            else if(sig_fld(inum)(1:4).eq.'CORI' .or.
     $              prs_fld(inum)(1:4).eq.'CORI')then 
               iusedotvar = 3
               print*,'Writing map scale factor on cross points -'
            else if(sig_fld(inum)(1:4).eq.'SLAB' .or.
     $              prs_fld(inum)(1:4).eq.'SLAB')then 
               call xtodot(irst,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing infinite reservoir slab temp -'
            else if(sig_fld(inum)(1:4).eq.'DLAT' .or.
     $              prs_fld(inum)(1:4).eq.'DLAT')then 
               call xtodot(xlat,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing latitude -'
            else if(sig_fld(inum)(1:4).eq.'DLON' .or.
     $              prs_fld(inum)(1:4).eq.'DLON')then 
               call xtodot(xlon,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing longitude -'
            else if(sig_fld(inum)(1:4).eq.'LAND' .or.
     $              prs_fld(inum)(1:4).eq.'LAND')then 
               iusedotvar = 4
               print*,'Writing land use -'
            else if(sig_fld(inum)(1:4).eq.'SNOW' .or.
     $              prs_fld(inum)(1:4).eq.'SNOW')then 
               call xtodot(xsnowc,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing snow cover -'
            else if(sig_fld(inum)(1:3).eq.'PBL' .or.
     $              prs_fld(inum)(1:3).eq.'PBL')then 
               call xtodot(pblht,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing PBL height -'
            else if(sig_fld(inum)(1:4).eq.'PBLR' .or.
     $              prs_fld(inum)(1:4).eq.'PBLR')then 
               call xtodot(pblregime,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing PBL regime -'
            else if(sig_fld(inum)(1:4).eq.'SHFX' .or.
     $              prs_fld(inum)(1:4).eq.'SHFX')then 
               call xtodot(senhtflx,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing sensible heat flux -'
            else if(sig_fld(inum)(1:4).eq.'LHFX' .or.
     $              prs_fld(inum)(1:4).eq.'LHFX')then 
               call xtodot(lathtflx,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing latent heat flux -'
            else if(sig_fld(inum)(1:3).eq.'UST' .or.
     $              prs_fld(inum)(1:3).eq.'UST')then 
               call xtodot(ustar,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing friction velocity -'
            else if(sig_fld(inum)(1:4).eq.'SRAD' .or.
     $              prs_fld(inum)(1:4).eq.'SRAD')then 
               call xtodot(shtwvrad,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing shortwave radiation -'
            else if(sig_fld(inum)(1:4).eq.'LRAD' .or.
     $              prs_fld(inum)(1:4).eq.'LRAD')then 
               call xtodot(lonwvrad,dotvar,ivert,ihorz,ivertl,ihorzl)
               print*,'Writing longwave radiation -'
            else if(sig_fld(inum)(1:4).eq.'SIL1' .or.
     $              prs_fld(inum)(1:4).eq.'SIL1')then
               call xtodot(soilt(:,:,1),dotvar,ivert,ihorz,ivertl
     $              ,ihorzl)
               print*,'Writing soil temperature in layer 1 -'
            else if(sig_fld(inum)(1:4).eq.'SIL2' .or.
     $              prs_fld(inum)(1:4).eq.'SIL2')then 
c               ik = 1 + ivert*ihorz
               call xtodot(soilt(:,:,2),dotvar,ivert,ihorz,ivertl
     $              ,ihorzl)
               print*,'Writing soil temperature in layer 2 -'
            else if(sig_fld(inum)(1:4).eq.'SIL3' .or.
     $              prs_fld(inum)(1:4).eq.'SIL3')then 
c               ik = 1 + 2*ivert*ihorz
               call xtodot(soilt(:,:,3),dotvar,ivert,ihorz,ivertl
     $              ,ihorzl)
               print*,'Writing soil temperature in layer 3 -'
            else if(sig_fld(inum)(1:4).eq.'SIL4' .or.
     $              prs_fld(inum)(1:4).eq.'SIL4')then 
c               ik = 1 + 3*ivert*ihorz
               call xtodot(soilt(:,:,4),dotvar,ivert,ihorz,ivertl
     $              ,ihorzl)
               print*,'Writing soil temperature in layer 4 -'
            else if(sig_fld(inum)(1:4).eq.'SIL5' .or.
     $              prs_fld(inum)(1:4).eq.'SIL5')then 
c               ik = 1 + 4*ivert*ihorz
               call xtodot(soilt(:,:,5),dotvar,ivert,ihorz,ivertl
     $              ,ihorzl)
               print*,'Writing soil temperature in layer 5 -'
               print*,'Writing latent heat flux -'
            else if(sig_fld(inum)(1:4).eq.'SIL6' .or.
     $              prs_fld(inum)(1:4).eq.'SIL6')then 
c               ik = 1 + 5*ivert*ihorz
               call xtodot(soilt(:,:,6),dotvar,ivert,ihorz,ivertl
     $              ,ihorzl)
               print*,'Writing soil temperature in layer 6 -'
            endif
            if (iusedotvar .eq. 1) then
               call write_data(igdfln,dotvar,ivert,ihorz,ighdr,
     $              gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,1)
            else if (iusedotvar .eq. 2) then
               call write_data(igdfln,dotmap,ivert,ihorz,ighdr,
     $              gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,1)
            else if (iusedotvar .eq. 3) then
               call write_data(igdfln,dotcor,ivert,ihorz,ighdr,
     $              gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,1)
            else if (iusedotvar .eq. 4) then
               call write_data(igdfln,xluse,ivert,ihorz,ighdr,
     $              gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,1)
            else if (iusedotvar .eq. 5) then
               call write_data(igdfln,slp,ivert,ihorz,ighdr,
     $              gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,1)
            else if (iusedotvar .eq. 6) then
               call write_data(igdfln,u10,ivert,ihorz,ighdr,
     $              gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,1)
            else if (iusedotvar .eq. 7) then
               call write_data(igdfln,v10,ivert,ihorz,ighdr,
     $              gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise,1)
            endif
         endif
 521  continue

c
c  Return for next time
      if (iend.eq.1) goto 999
      goto 30
999   print *,'End'

c
c     free up the memory
c
      deallocate(uu)
      deallocate(vv)
      deallocate(tt)
      deallocate(qv)
      deallocate(qc)
      deallocate(qr)
      deallocate(cice)
      deallocate(sice)
      deallocate(pp)
      deallocate(rh)
      deallocate(rtten)
c     derived 3-d fields (ihorz,ivert,kz):
      deallocate(unor)
      deallocate(vnor)
c  3-d (ivert,ihorz,kz+1):
      deallocate(ww)
c  3-d (ivert,ihorz,6):
      deallocate(soilt)
c  Raw precip 2-d fields (ivert,ihorz):
      deallocate(re)
      deallocate(rc)
      deallocate(prec)
      deallocate(prec_3hr)
      deallocate(re_3hr)
      deallocate(rc_3hr)
      deallocate(last_output_hr_tot_prec)
      deallocate(last_output_hr_exp_prec)
      deallocate(last_output_hr_con_prec)
c  Raw data 2-d fields (ivert,ihorz):
      deallocate(pstx)
      deallocate(tg)
      deallocate(ter)
      deallocate(xmap)
      deallocate(dotmap)
      deallocate(dotcor)
      deallocate(irst)
      deallocate(xlat)
      deallocate(xluse)
      deallocate(xsnowc)
      deallocate(xlon)
c  Derived 3d fields (ivert,ihorz,kz):
      deallocate(prs)
      deallocate(dew)
      deallocate(wwh)
      deallocate(prsd)
      deallocate(pdata)
      deallocate(ttp)
      deallocate(prsp)
      deallocate(qvlog)
      deallocate(tmp3dp)
      deallocate(tmp3ds)
c     ivert*ihorz*2
      deallocate(hlhsl)
c  Derived 2d fields (ivert,ihorz):
      deallocate(pstd)
      deallocate(sfp)
      deallocate(slp)
      deallocate(dotvar)
      deallocate(dat)
      deallocate(dluse)
      deallocate(tmp2d)
      deallocate(tot_prec)
c  Miscellaneous 2-d fields (ivert,ihorz):
      deallocate(tvirt)
c     Miscellaneous 1-d fields (kz):
      deallocate(plevs)
      deallocate(maxp)
      deallocate(minp)
      deallocate(maxpd)
      deallocate(minpd)
      deallocate(sih)
c     optional 2-d fields (ivert,ihorz):
      deallocate(senhtflx)
      deallocate(lathtflx)
      deallocate(lonwvrad)
      deallocate(shtwvrad)
      deallocate(pblht)
      deallocate(ustar)
      deallocate(pblregime)
      deallocate(moisture)
      deallocate(u10)
      deallocate(v10)
      deallocate(t2)

      call gd_clos(igdfln,iret)
      call gendp(1,iret)
      if  ( iret .ne. 0 )  call er_wmsg ( 'GEMPLT', iret, ' ', ier )
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

      subroutine xtodot2(xslab,dotslab,ivert,ihorz,maxivert,maxihorz)
c  This subroutine interpolates to dot points, but along edges of
c  dotpoint domain, only uses the same values as the crosspoint edge.
c  The result is that fields like precip can not obtain negative values
c  as a result of extrapolation
      dimension xslab(ivert,ihorz),dotslab(ivert,ihorz)
      integer aa,bb,cc,dd
      do i=1,maxivert
      do j=1,maxihorz
        aa = min(i,maxivert-1)
        bb = min(j,maxihorz-1)
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
      subroutine xtodot(xslab,dotslab,ivert,ihorz,maxivert,maxihorz)
c
      dimension xslab(ivert,ihorz), dotslab(ivert,ihorz)
c
c   Interpolate in the interior.
c
      do 100 j=2,maxihorz-1
      do 100 i=2,maxivert-1
         dotslab(i,j)=.25*(xslab(i-1,j-1)+xslab(i,j-1)+xslab(i-1,j)+
     &      xslab(i,j))
  100    continue
c
c   Extrapolate out to top and bottom edges.
c
      do 200 j=2,maxihorz-1
         dotslab(1,j)=(3.*(xslab(1,j-1)+xslab(1,j))-
     &      (xslab(2,j-1)+xslab(2,j)))/4.
         dotslab(maxivert,j)=(3.*(xslab(maxivert-1,j-1)+
     $        xslab(maxivert-1,j))-
     &      (xslab(maxivert-2,j-1)+xslab(maxivert-2,j)))/4.
  200 continue
c
c   Extrapolate out to left and right edges.
c
      do 300 i=2,maxivert-1
         dotslab(i,1)=(3.*(xslab(i-1,1)+xslab(i,1))-
     &      (xslab(i-1,2)+xslab(i,2)))/4.
         dotslab(i,maxihorz)=(3.*(xslab(i-1,maxihorz-1)+
     $        xslab(i,maxihorz-1))-
     &      (xslab(i-1,maxihorz-2)+xslab(i,maxihorz-2)))/4.
  300 continue
c
c   Extrapolate out to corners.
c
      dotslab(1,1)=(3.*xslab(1,1)-xslab(2,2))/2.
      dotslab(maxivert,1)=(3.*xslab(maxivert-1,1)-
     $     xslab(maxivert-2,2))/2.
      dotslab(1,maxihorz)=(3.*xslab(1,maxihorz-1)-
     $     xslab(2,maxihorz-2))/2.
      dotslab(maxivert,maxihorz)=(3.*xslab(maxivert-1,maxihorz-1)-
     &   xslab(maxivert-2,maxihorz-2))/2.
      return
      end



      subroutine read_3d(iunit,scr3d,chdum,ni,nj,nk,itime,ifchoose)
      implicit none
      integer itime,ni,nj,nk,i,j,k,iunit
      real scr3d(ni,nj,nk)
      character chdum*80, tmpa*3
      logical ifchoose

      print*,'Reading ', chdum(1:67)
      read (iunit) (((scr3d(i,j,k),i=1,ni),j=1,nj),k=1,nk)
      if(itime.eq.1 .and. ifchoose)then
         print*,'Would you like this field written out(y/n)?'
         read(5, '(A)') tmpa
         if(tmpa(1:1).eq.'y') print *,'not yet implemented, sorry'
      endif
      
      return
      end
      
      subroutine read_2d(iunit,scr2d,chdum,ni,nj,itime,ifchoose)
      implicit none
      integer itime,ni,nj,iunit,i,j
      real scr2d(ni,nj)
      character chdum*80, tmpa*3
      logical ifchoose

      print*,'Reading ', chdum(1:67)
      read (iunit) ((scr2d(i,j),i=1,ni),j=1,nj)
      if(itime.eq.1 .and. ifchoose)then
         print*,'Would you like this field written out(y/n)?'
         read(5, '(A)') tmpa
         if(tmpa(1:1).eq.'y') print *,'not yet implemented, sorry'
      endif
      
      return
      end
      
      subroutine write_data (igdfln,grid,igx,igy,ighdr,
     $     gdattm,ilvl,ivcor,parm,rewrite,ipktyp,iprecise
     $     ,iflip)
      parameter (nhdrsz=2)
      integer ighdr(nhdrsz),ilvl(2)
      real grid(igx,igy)
      real, allocatable, dimension(:,:) :: flip
      character gdattm(2)*20,parm*12
      logical rewrite

      if (ivcor .eq. 1) then
         print *,'Writing ',parm,' on pressure levels',
     $        grid(igy/2,igx/2),'M/S'
      endif
      if (iflip.eq.1) then
         allocate(flip(igy,igx))
      else
         allocate(flip(igx,igy))
      endif
      do j = 1,igy
         do i = 1,igx
c            ij = i + (j-1)*igx
c            ji = j + (i-1)*igy
            if (iflip.eq.1) then
               flip(j,i) = grid(i,j)
            else
               flip(i,j) = grid(i,j)
            endif
         enddo
      enddo
      if (iflip.eq.1) then
         call gd_wpgd(igdfln,flip,igy,igx,ighdr,gdattm,ilvl,ivcor,
     $        parm,rewrite,ipktyp,iprecise,iret)
      else
         call gd_wpgd(igdfln,flip,igx,igy,ighdr,gdattm,ilvl,ivcor,
     $        parm,rewrite,ipktyp,iprecise,iret)
      endif
      if(iret.ne.0) then
         print *,'gd_wpgd: iret=',iret,' Stopping . . .'
         call bail_out('error in write_data')
      endif
      deallocate(flip)

c          call gd_wpgd(igdfln,dat,igx,igy,ighdr,gdattm,ilvl,ivcor,
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
      return
c      ivert -> igx
c      ihorz -> igy
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

      subroutine mdate_split(mdate
     $     ,yyyymmddhhnnss,iyear,imonth,iday,ihour,imin,isec)

      character *(*) yyyymmddhhnnss

      if(mdate.lt.1e9)then
c     
c     ... decompose mdate into year, month, day, hour and minute
c     
         iyear  = mdate/1000000
         imonth = (mdate-iyear*1000000)/10000
         iday   = (mdate-iyear*1000000-imonth*10000)/100
         ihour  = (mdate-iyear*1000000-imonth*10000-iday*100)
         imin   = 0
      elseif(mdate.gt.1e9)then
c     
c     ... decompose mdate into year, month, day, hour and minute
c     
         iyear  = mdate/100000000
         imonth = (mdate-iyear*100000000)/1000000
         iday   = (mdate-iyear*100000000-imonth*1000000)/10000
         ihour  = (mdate-iyear*100000000-imonth*1000000-iday*
     +        10000)/100
         imin   = (mdate-iyear*100000000-imonth*1000000-iday*
     +        10000-ihour*100)
c     
      endif
      isec  = 0
      iyear = iyear + 1900
      write(yyyymmddhhnnss,'(i4.4,i2.2,i2.2,i2.2,i2.2,i2.2)')
     $     iyear,imonth,iday,ihour,imin,isec
      return
      end

      subroutine fill_data_withv1(
     $     uu,vv,tt,qv,qc,qr,cice,sice,pp,rh,rtten,
     $     ww,soilt,
     $     re,rc,prec,
     $     pstx,tg,ter,xmap,dotmap,dotcor,sst,slpcrs,irst,
     $     xlat,xluse,xsnowc,xlon,
     $     prs,dew,wwh,prsd,pdata,ttp,prsp,qvlog,tmp3dp,tmp3ds,
     $     unor,vnor,hlhsl,
     $     pstd,sfp,slp,dotvar,dat,dluse,tmp2d,tot_prec,
     $     tvirt,
     $     maxp,minp,maxpd,minpd,sih,
     $     senhtflx,lathtflx,lonwvrad,shtwvrad,pblht,ustar,
     $     pblregime,moisture,u10,v10,t2,
     $     mrfv1,mifv1,mrfc,mifc,
     $     ivert,ihorz,kz,iprog,
     $     iend
     $     )
      implicit none
      integer ivert,ihorz,kz,iprog,iend
      real uu(ivert,ihorz,kz)
      real vv(ivert,ihorz,kz)
      real tt(ivert,ihorz,kz)
      real qv(ivert,ihorz,kz)
      real qc(ivert,ihorz,kz)
      real qr(ivert,ihorz,kz)
      real cice(ivert,ihorz,kz)
      real sice(ivert,ihorz,kz)
      real pp(ivert,ihorz,kz)
      real rh(ivert,ihorz,kz)
      real rtten(ivert,ihorz,kz)
c  3-d (ivert,ihorz,kz+1):
      real ww(ivert,ihorz,(kz+1))
c  3-d (ivert,ihorz,6)
      real soilt(ivert,ihorz,6)
c  Raw precip 2-d fields (ivert,ihorz):
      real re(ivert,ihorz)
      real rc(ivert,ihorz)
      real prec(ivert,ihorz)
      real prec_3hr(ivert,ihorz)
      real re_3hr(ivert,ihorz)
      real rc_3hr(ivert,ihorz)
c  Raw data 2-d fields (ivert,ihorz):
      real pstx(ivert,ihorz)
      real tg(ivert,ihorz)
      real ter(ivert,ihorz)
      real xmap(ivert,ihorz)
      real dotmap(ivert,ihorz)
      real dotcor(ivert,ihorz)
      real sst(ivert,ihorz)
      real slpcrs(ivert,ihorz)
      real irst(ivert,ihorz)
      real xlat(ivert,ihorz)
      real xluse(ivert,ihorz)
      real xsnowc(ivert,ihorz)
      real xlon(ivert,ihorz)
c  Derived 3d fields (ivert,ihorz,kz):
      real prs(ivert,ihorz,kz)
      real dew(ivert,ihorz,kz)
      real wwh(ivert,ihorz,kz)
      real prsd(ivert,ihorz,kz)
      real pdata(ivert,ihorz,kz)
      real ttp(ivert,ihorz,kz)
      real prsp(ivert,ihorz,kz)
      real qvlog(ivert,ihorz,kz)
      real tmp3dp(ivert,ihorz,kz)
      real tmp3ds(ivert,ihorz,kz)
c     derived 3-d fields (ihorz,ivert,kz):
      real unor(ihorz,ivert,kz)
      real vnor(ihorz,ivert,kz)
c  ivert*ihorz*2
      real hlhsl(ivert,ihorz,2)
c  Derived 2d fields (ivert,ihorz):
      real pstd(ivert,ihorz)
      real sfp(ivert,ihorz)
      real slp(ivert,ihorz)
      real dotvar(ivert,ihorz)
      real dat(ivert,ihorz)
      real dluse(ivert,ihorz)
      real tmp2d(ivert,ihorz)
      real tot_prec(ivert,ihorz)
c  Miscellaneous 2-d fields (ivert,ihorz):
      real tvirt(ivert,ihorz)
c  Miscellaneous 1-d fields (kz):
      real maxp(kz)
      real minp(kz)
      real maxpd(kz)
      real minpd(kz)
      real sih(kz)
c  optional 2-d fields (ivert,ihorz):
      real senhtflx(ivert,ihorz)
      real lathtflx(ivert,ihorz)
      real lonwvrad(ivert,ihorz)
      real shtwvrad(ivert,ihorz)
      real pblht(ivert,ihorz)
      real ustar(ivert,ihorz)
      real pblregime(ivert,ihorz)
      real moisture(ivert,ihorz)
      real u10(ivert,ihorz)
      real v10(ivert,ihorz)
      real t2(ivert,ihorz)
      real mrfv1(1000,20)
      integer mifv1(1000,20)
      character*80 mrfc(1000,20),mifc(1000,20)
      integer idry,imoist,inhyd,iice,inav,iiceg,num3d,num2d,isoil,
     $     irddim,ik,ifld,itime
      logical ifchoose
      real ref_temp

      read (20,END=999) mifv1,mrfv1,mifc,mrfc
      if (iprog.eq.6) then
         idry = mifv1(3,6)
         imoist = mifv1(4,6)
         inhyd = mifv1(5,6)
         ref_temp = mrfv1(3,6)
         iice = mifv1(7,6)
         inav = mifv1(8,6)
         iiceg = mifv1(9,6)
         num3d = mifv1(201,6)
         num2d = mifv1(202,6)
         isoil = mifv1(355,6)
         irddim = mifv1(201,6) -
     $        (3+(1-idry)+(1-idry)*2*(1-mod(imoist,2))
     $        +2*iice*(1-mod(imoist,2))+2*inhyd+inav+2*iiceg*(1-
     $        mod(imoist,2)))
      endif

c
c   Read data
c
      num3d=mifv1(201,iprog)
      num2d=mifv1(202,iprog)
c  3-D fields:
      do 205 ifld=205,204+num3d
         if(mifc(ifld,iprog)(1:8).eq.'U       ')then
C     1) u component wind
            call read_3d(20,uu,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'V       ')then
C     2) v component wind
            call read_3d(20,vv,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'T       ')then
C     3) temperature
            call read_3d(20,tt,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'Q       '.or.
     &           mifc(ifld,iprog)(1:8).eq.'RH      ')then
C     4) qv or rh
            call read_3d(20,qv,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'CLW     ')then
C     5) qv cloud
            call read_3d(20,qc,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'RNW     ')then
C     6) qv rain
            call read_3d(20,qr,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'ICE     ')then
C     3.25) cloud ice
            call read_3d(20,cice,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'SNOW    ')then
C     3) snow ice
            call read_3d(20,sice,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'GRAW    ' .or.
     $           mifc(ifld,iprog)(1:8).eq.'GRAUPEL ')then
c     graupel
            call read_3d(20,tmp3ds,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'NICEPART' .or.
     $           mifc(ifld,iprog)(1:8).eq.'NCI     ')then
c     number concentration
            call read_3d(20,tmp3ds,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'TKE     ')then
C     7) TURBULENT KINETIC ENERGY
            call read_3d(20,tmp3ds,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'RAD TEND' .or.
     $           mifc(ifld,iprog)(1:8).eq.'RTTEN   ')then
c     8) ATMOSPHERIC RADIATION TENDENCY
            call read_3d(20,rtten,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'W       ')then
C     8) VERTICAL WIND COMPONENT
            call read_3d(20,ww,mifc(ifld,iprog),ivert,ihorz,(kz+1),
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'PP      '.or.
     &           mifc(ifld,iprog)(1:8).eq.'H       ')then
C     9) PRESSURE PERTURBATION
            call read_3d(20,pp,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         else
C     UNKNOWN FIELD 
            call read_3d(20,tmp3ds,mifc(ifld,iprog),ivert,ihorz,kz,
     $           itime,ifchoose)
         endif
 205  continue
C  2-D fields:
      do 206 ifld=205+num3d,204+num3d+num2d
         if(mifc(ifld,iprog)(1:8).eq.'PSTARCRS')then
C     1) PSTAR in kPa
            call read_2d(20,pstx,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'GROUND T')then
C     1) GROUND TEMPERATURE
            call read_2d(20,tg,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'RAIN CON')then
C     1) ACCUMULATED CONVECTIVE PRECIPITATION
            call read_2d(20,rc,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'RAIN NON')then
C     2) ACCUMULATED RESOLVED PRECIPITATION
            call read_2d(20,re,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'TERRAIN ')then
C     1) TERRAIN
            call read_2d(20,ter,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'MAPFACCR')then
C     2) MAP FACTOR CROSS POINTS
            call read_2d(20,xmap,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'MAPFACDT')then
C     3) MAP FACTOR DOT POINTS
            call read_2d(20,dotmap,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:7).eq.'TSEASFC')then
C     SEA-SURFACE TEMPERATURE
            call read_2d(20,sst,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'CORIOLIS')then
C     4) CORIOLIS
            call read_2d(20,dotcor,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'RES TEMP')then
C     5) (OPTIONAL) RESERVOIR TEMPERATURE
            call read_2d(20,irst,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'LATITCRS')then
C     6) LATITUDE
            call read_2d(20,xlat,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'LONGICRS')then
C     7) LONGITUDE
            call read_2d(20,xlon,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'LAND USE')then
C     8) LAND USE CATEGORY
            call read_2d(20,xluse,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'SNOWCOVR')then
C     9) SNOW COVER
            call read_2d(20,xsnowc,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:18).eq.'SENSIBLE HEAT FLUX' .or.
     $           mifc(ifld,iprog)(1:8).eq.'SHFLUX  ')then
C     12) SENSIBLE HEAT FLUX
            call read_2d(20,senhtflx,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:16).eq.'LATENT HEAT FLUX' .or.
     $           mifc(ifld,iprog)(1:8).eq.'LHFLUX  ')then
C     13) LATENT HEAT FLUX
            call read_2d(20,lathtflx,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'LONGWAVE' .or.
     $           mifc(ifld,iprog)(1:8).eq.'LWDOWN  ')then
C     16) LONGWAVE RADIATION
            call read_2d(20,lonwvrad,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:9).eq.'SHORTWAVE' .or.
     $           mifc(ifld,iprog)(1:8).eq.'SWDOWN  ')then
C     15) SHORTWAVE RADIATION
            call read_2d(20,shtwvrad,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:9).eq.'PLANETARY' .or.
     $           mifc(ifld,iprog)(1:8).eq.'PBL HGT ')then
C     10) PLANETARY BOUNDARY LAYER HEIGHT
            call read_2d(20,pblht,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'FRICTION' .or.
     $           mifc(ifld,iprog)(1:8).eq.'UST     ')then
C     14) FRICTION VELOCITY
            call read_2d(20,ustar,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:10).eq.'PBL REGIME' .or.
     $           mifc(ifld,iprog)(1:8).eq.'REGIME  ')then
C     11) PBL REGIME
            call read_2d(20,pblregime,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:21).eq.'MOISTURE AVAILABILITY')then
C     18) MOISTURE AVAILABILITY
            call read_2d(20,moisture,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:7).eq.'SOIL T ')then
C     20-25) SOIL TEMPERATURES
            read(mifc(ifld,iprog)(8:8),'(i1)') ik
c            ik = 1 + (ik-1)*ivert*ihorz
            print *,'reading soilt at ',mifc(ifld,iprog)(8:8),' = ',ik
            call read_2d(20,soilt(:,:,ik),mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'PSEALVLD')then
C     SEA-LEVEL PRESSURE FROM DATAGRID
            call read_2d(20,slp,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         elseif(mifc(ifld,iprog)(1:8).eq.'PSEALVLC')then
C     SEA-LEVEL PRESSURE FROM DATAGRID CROSS POINT
            call read_2d(20,slpcrs,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         else
C     UNKNOWN FIELD
            call read_2d(20,tmp2d,mifc(ifld,iprog),ivert,ihorz,
     $           itime,ifchoose)
         endif                  ! enf of scanning the field names
 206  continue                  ! end of looping through the 2-d fields
      return
 999  print *,'end of file'
      iend = 1
      return
      end


      subroutine fill_data_withv3(
     $     uu,vv,tt,qv,qc,qr,cice,sice,pp,rh,rtten,
     $     ww,soilt,
     $     re,rc,prec,
     $     pstx,tg,ter,xmap,dotmap,dotcor,sst,slpcrs,irst,
     $     xlat,xluse,xsnowc,xlon,
     $     prs,dew,wwh,prsd,pdata,ttp,prsp,qvlog,tmp3dp,tmp3ds,
     $     unor,vnor,hlhsl,
     $     pstd,sfp,slp,dotvar,dat,dluse,tmp2d,tot_prec,
     $     tvirt,
     $     maxp,minp,maxpd,minpd,sih,
     $     senhtflx,lathtflx,lonwvrad,shtwvrad,pblht,ustar,
     $     pblregime,moisture,u10,v10,t2,
     $     bhi,bhr,bhic,bhrc,
     $     ivert,ihorz,kz,iprog,
     $     ftimes,ntimes,lasttime,chdates,ndates,lastdate,iend
     $     )
      implicit none
      integer ivert,ihorz,kz,iprog,iend
      real uu(ivert,ihorz,kz)
      real vv(ivert,ihorz,kz)
      real tt(ivert,ihorz,kz)
      real qv(ivert,ihorz,kz)
      real qc(ivert,ihorz,kz)
      real qr(ivert,ihorz,kz)
      real cice(ivert,ihorz,kz)
      real sice(ivert,ihorz,kz)
      real pp(ivert,ihorz,kz)
      real rh(ivert,ihorz,kz)
      real rtten(ivert,ihorz,kz)
c  3-d (ivert,ihorz,kz+1):
      real ww(ivert,ihorz,(kz+1))
c  3-d (ivert,ihorz,6)
      real soilt(ivert,ihorz,6)
c  Raw precip 2-d fields (ivert,ihorz):
      real re(ivert,ihorz)
      real rc(ivert,ihorz)
      real prec(ivert,ihorz)
c  Raw data 2-d fields (ivert,ihorz):
      real pstx(ivert,ihorz)
      real tg(ivert,ihorz)
      real ter(ivert,ihorz)
      real xmap(ivert,ihorz)
      real dotmap(ivert,ihorz)
      real dotcor(ivert,ihorz)
      real sst(ivert,ihorz)
      real slpcrs(ivert,ihorz)
      real irst(ivert,ihorz)
      real xlat(ivert,ihorz)
      real xluse(ivert,ihorz)
      real xsnowc(ivert,ihorz)
      real xlon(ivert,ihorz)
c  Derived 3d fields (ivert,ihorz,kz):
      real prs(ivert,ihorz,kz)
      real dew(ivert,ihorz,kz)
      real wwh(ivert,ihorz,kz)
      real prsd(ivert,ihorz,kz)
      real pdata(ivert,ihorz,kz)
      real ttp(ivert,ihorz,kz)
      real prsp(ivert,ihorz,kz)
      real qvlog(ivert,ihorz,kz)
      real tmp3dp(ivert,ihorz,kz)
      real tmp3ds(ivert,ihorz,kz)
c     derived 3-d fields (ihorz,ivert,kz):
      real unor(ihorz,ivert,kz)
      real vnor(ihorz,ivert,kz)
c  ivert*ihorz*2
      real hlhsl(ivert,ihorz,2)
c  Derived 2d fields (ivert,ihorz):
      real pstd(ivert,ihorz)
      real sfp(ivert,ihorz)
      real slp(ivert,ihorz)
      real dotvar(ivert,ihorz)
      real dat(ivert,ihorz)
      real dluse(ivert,ihorz)
      real tmp2d(ivert,ihorz)
      real tot_prec(ivert,ihorz)
c  Miscellaneous 2-d fields (ivert,ihorz):
      real tvirt(ivert,ihorz)
c  Miscellaneous 1-d fields (kz):
      real maxp(kz)
      real minp(kz)
      real maxpd(kz)
      real minpd(kz)
      real sih(kz)
c  optional 2-d fields (ivert,ihorz):
      real senhtflx(ivert,ihorz)
      real lathtflx(ivert,ihorz)
      real lonwvrad(ivert,ihorz)
      real shtwvrad(ivert,ihorz)
      real pblht(ivert,ihorz)
      real ustar(ivert,ihorz)
      real pblregime(ivert,ihorz)
      real moisture(ivert,ihorz)
      real u10(ivert,ihorz)
      real v10(ivert,ihorz)
      real t2(ivert,ihorz)
      integer idry,imoist,inhyd,iice,inav,iiceg,num3d,num2d,isoil
     $     irddim,i,j,k,icd
      real ref_temp
      integer bhi(50,20), flag, ik
      real bhr(20,20)
      character*80 bhic(50,20),bhrc(20,20)
      integer :: ndim, ntimes, ndates
      real :: time, sample, ftimes(60), lasttime
      character chdates(60)*24, lastdate*24
      integer, dimension(4) :: start_index, end_index
      character (len= 4) :: staggering
      character (len= 4) :: ordering
      character (len=24) :: start_date
      character (len=24) :: current_date
      character (len= 9) :: name
      character (len=25) :: units
      character (len=46) :: description

c      num_1d_outputs = 0
c      num_2d_outputs = 0
c      num_3d_outputs = 0
      print *,'filling v3 for iprog ',iprog
 421  read(20,end=422,err=422) flag
      if (flag.eq.0) then
        read(20,err=422) bhi, bhr, bhic, bhrc
      else if (flag.eq.1) then
         read (20,err=422) ndim,start_index,end_index,time,
     &      staggering,ordering,current_date,name,
     &      units,description
         if (time.ne.lasttime) then
            ntimes = ntimes + 1
            ftimes(ntimes) = time
         endif
         if (current_date.ne.lastdate) then
            ndates = ndates + 1
            chdates(ndates) = current_date
         endif
         lastdate = current_date
         lasttime = time
         print *,'got name ',name,description,current_date,ndim,time
c         print*
c         print*,'ndim=',ndim
c         print*,'start_index=',start_index
c         print*,'end_index=',end_index
c         print*,'time=',time
c         print*,'staggering=',staggering
c         print*,'ordering=',ordering
c         print*,'current_date=',current_date
c         print*,'name=',name
c         print*,'units=',units
c         print*,'description=',description
c         if (ndim == 1) then
c            allocate(data(end_index(1), 1, 1, 1))
c            num_1d_outputs = num_1d_outputs + 1
c         elseif (ndim == 2) then
c            allocate(data(end_index(1), end_index(2), 1, 1))
c            num_2d_outputs = num_2d_outputs + 1
c         elseif (ndim == 3) then
c            allocate(data(end_index(1), end_index(2), end_index(3), 1))
c            num_3d_outputs = num_3d_outputs + 1
c        endif
c
c        read(iunit) data

         if (name.eq.'PSEALVLD ') then
            read(20) slp
         elseif (name.eq.'LATITDOT '.or.
     &           name.eq.'HSFC     '.or.
     &           name.eq.'LONGIDOT ') then
            read(20) tmp2d
         elseif (name.eq.'LATITCRS ') then
            read(20) xlat
         elseif (name.eq.'LONGICRS ') then
            read(20) xlon
         elseif (name.eq.'RES TEMP '.or.
     &           name.eq.'RES TMP  ') then
            read(20) irst
         elseif (name.eq.'U        ') then
            read(20) uu
         elseif (name.eq.'V        ') then
            read(20) vv
         elseif (name.eq.'T        ') then
            read(20) tt
         elseif (name.eq.'Q        ') then
            read(20) qv
         elseif (name.eq.'RH       ') then
            print *,'reading rh into qv'
            read(20) qv
         elseif (name.eq.'CLW      ') then
            read(20) qc
         elseif (name.eq.'RNW      ') then
            read(20) qr
         elseif (name.eq.'ICE      ') then
            read(20) cice
         elseif (name.eq.'SNOW     ') then
            read(20) sice
         elseif (name.eq.'GRAUPEL  ') then
            read(20) tmp3ds
         elseif (name.eq.'NCI      ') then
            read(20) tmp3ds
         elseif (name.eq.'W        ') then
            read(20) ww
         elseif (name.eq.'PP       ') then
            read(20) pp
         elseif (name.eq.'H        ') then
            print *,'reading H into pp'
            if (iprog.ne.2.and.iprog.ne.3) then
               print*,'RIP only expects geop. hgt. in prs-level data.'
               print*,'Can''t handle it here.  Stopping.'
               stop
            endif
            read(20) pp
         elseif (name.eq.'PSTARCRS ') then
            read(20) pstx
            pstx = pstx * 0.001 ! Pascals to kPa same units as V2
         elseif (name.eq.'PSEALVLC ') then
            read(20) slpcrs
         elseif (name.eq.'GROUND T ') then
            read(20) tg
         elseif (name.eq.'TSEASFC  ') then
            read(20) sst
         elseif (name.eq.'RAIN CON ') then
            read(20) rc
            print *,'filled v3 rc'
         elseif (name.eq.'RAIN NON ') then
            read(20) re
            print *,'filled v3 re'
         elseif (name.eq.'TERRAIN  ') then
            read(20) ter
         elseif (name.eq.'MAPFACCR ') then
            read(20) xmap
         elseif (name.eq.'MAPFACDT '.or.name.eq.'MAPFADOT ') then
            read(20) dotmap
         elseif (name.eq.'CORIOLIS ') then
            read(20) dotcor
         elseif (name.eq.'LAND USE ') then
            read(20) xluse
         elseif (name.eq.'SNOWCOVR ') then
            read(20) xsnowc
         elseif (name.eq.'PBL HGT  ') then
            read(20) pblht
         elseif (name.eq.'REGIME   ') then
            read(20) pblregime
         elseif (name.eq.'SHFLUX   ') then
            read(20) senhtflx
         elseif (name.eq.'LHFLUX   ') then
            read(20) lathtflx
         elseif (name.eq.'MAVAIL   ') then
            read(20) moisture
         elseif (name.eq.'UST      ') then
            read(20) ustar
         elseif (name.eq.'SWDOWN   ') then
            read(20) shtwvrad
         elseif (name.eq.'LWDOWN   ') then
            read(20) lonwvrad
         elseif (name.eq.'U10      '.or.
     $           name(1:5).eq.'USFC ') then
            read(20) u10
            print *,'filling u10 array with ',name
         elseif (name.eq.'V10      '.or.
     $           name(1:5).eq.'VSFC ') then
            read(20) v10
            print *,'filling v10 array with ',name
         elseif (name.eq.'T2       '.or.
     $           name(1:5).eq.'TSFC ') then
            print *,'filling t2 array with ',name
            read(20) t2
         elseif (name.eq.'SIGMAH   ') then
            read(20) sih
            print *,'got sih = ',sih
         elseif (name.eq.'PRESSURE ') then
            if (iprog.eq.2 .or. iprog.eq.3) then
               read(20) sih
               print *,' *** WARNING putting pressure into sih'
               print *,'iprog = ',iprog,' should be 2 or 3'
c               sih = sih*.01
               print *,' *** CHECK OUTPUT, I have not tested ***',sih
            endif
         elseif (name(1:6).eq.'SOIL T') then
            read(name(8:8),'(i1)') ik
            print *,'reading soilt at ',name(8:8),' = ',ik
            read(20) soilt(:,:,ik)
         elseif (ordering.eq.'YXS '.or.ordering.eq.'YXW ') then
c
c         unknown 3d field
c
c         Note: V3 header: name*9, units*25, description*46
c
            if (iprog.eq.2.or.iprog.eq.3) then
               read(20) tmp2d,
     &              (((tmp3ds(i,j,k),i=1,ivert),j=1,ihorz),k=kz,1,-1)
            else
               read(20) tmp3ds
            endif
            icd=1
            if (staggering.eq.'D   ') icd=0
            if (end_index(3).eq.kz+1) then ! interpolate to half sigma levels
               do k=1,kz
                  do j=1,ihorz-icd
                     do i=1,ivert-icd
                        tmp3ds(i,j,k)=.5*(tmp3ds(i,j,k)+tmp3ds(i,j,k+1))
                     enddo
                  enddo
               enddo
               end_index(3)=kz
            endif
c
c         Write out surface part from pressure level data.
c
         else
            read(20)
            print*,'   Discarding ',name,
     &             '   because it is not 2- or 3-d data.'
         endif
 353     continue
         print*,'Processing MM5 variable ',name
      elseif (flag.eq.2) then
         return
      else
         print*,'Ran into a flag not =1 or 2 while reading data. = '
     $        ,flag
         iend = 1
         return
      endif
      goto 421
c
 422  iend = 1
      return
      end
