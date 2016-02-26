#	Top-level Makefile for rewrite_domain3_header

#	Macros, these should be generic for all machines
.IGNORE:

MAKE	=	make

default:
	uname -s > .tmpfile
	@grep CRAY .tmpfile ; \
	if [ $$? = 0 ] ; then echo "Compiling for Cray"							; \
		( $(MAKE) -f Makefile.Cray )								; \
	else grep OSF .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for Compaq"						; \
		( $(MAKE) -f Makefile.OSF1 )								; \
	else grep SunOS .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for SunOS"						; \
		( $(MAKE) -f Makefile.SunOS )								; \
	else grep Linux .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for Linux"						; \
		( $(MAKE) -f Makefile.Linux linuxfile all)								; \
	else echo "Do not know how to compile for the `cat .tmpfile` machine." 				; \
	fi ; \
	fi ; \
	fi ; \
	fi

include Makefile.common


