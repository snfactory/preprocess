######################################################################
## Filename:      Makefile.am
## Version:       1.0
## Description:   Main AutoMake makefile
## Author:        Yannick Copin <y.copin@ipnl.in2p3.fr>
## Created at:    Tue Mar 25 16:21:44 2003
######################################################################

SUBDIRS=pkg

EXTRA_DIST=add_defs

AM_MAKEFLAGS = CC=$(CC) F77=$(F77) F_OPT="$(F_OPT) -ftypeless-boz" LD=$(CC) \
	ADD_LIBS="-lc @FLIBS@" RANLIB=$(RANLIB) C_OPT=$(C_OPT)

# needs to be sure makefiles are sync with add_defs/makedefs
Makefile : $(srcdir)/add_defs/makedefs 

$(srcdir)/add_defs/makedefs : $(srcdir)/add_defs/makedefs.in $(top_builddir)/config.status $(BUILT_SOURCES)
	cd $(top_builddir) \
	  && CONFIG_FILES=$@ CONFIG_HEADERS= $(SHELL) ./config.status
