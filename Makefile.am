######################################################################
## Filename:      Makefile.am
## Version:       1.0
## Description:   Main AutoMake makefile
## Author:        Yannick Copin <y.copin@ipnl.in2p3.fr>
## Created at:    Tue Mar 25 16:21:44 2003
######################################################################

SUBDIRS=pkg doc

EXTRA_DIST=add_defs

AM_MAKEFLAGS = CC=$(CC) F77=$(F77) F_OPT="$(F_OPT) -ftypeless-boz" LD=$(CC) \
	ADD_LIBS="-lc @FLIBS@" RANLIB=$(RANLIB) C_OPT=$(C_OPT)

# needs to be sure makefiles are sync with add_defs/makedefs and configure

# apparently missing from automake output
# does not work with old automakes 
# $(SHELL) ./config.status $@ $(am__depfiles_maybe)

Makefile : $(top_builddir)/config.status $(srcdir)/Makefile.in  
	cd $(top_builddir) && CONFIG_FILES=$@ CONFIG_HEADERS= CONFIG_LINKS= $(SHELL) ./config.status


# trick to get config-clean-files checked somwhere in the build process
#Makefile : $(CONFIG_CLEAN_FILES)
# if set, the makefile target from pkg/preprocess/source will not be 
# executed properly
# it implies writing a correct makefile rebuild scheme for all subdirs

all : all-am
