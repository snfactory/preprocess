######################################################################
## Filename:      makefile.in
## Version:       
## Description:   
## Author:        Yannick Copin <y.copin@ipnl.in2p3.fr>
## Created at:    Tue Mar 25 16:21:56 2003
## $Id: makefile.in,v 1.30 2012/03/13 09:45:58 gangler Exp $
######################################################################

topdir = ..
include $(topdir)/makedefs

# ============================================================

# Local libraries

LOCLIB = $(O)/libpreprocess.a

# Targets

WSRC = bichip_build.cxx\
       preprocess_bias.cxx\

SRC =	 preprocess.cxx\
	 clean_cosmic.cxx \
	 fit_image.cxx\
	 remove_dark.cxx\
	 remove_bias_dark.cxx\
	 fit_dark.cxx\
	 guess_dark.cxx\
	 build_dark.cxx\
         quick_median.cxx\
	 ccd_rem_cosmic.cxx\
	 split_image.cxx\
	 fft_image.cxx\
	 copy_image.cxx\
	 copy_image_slice.cxx\
	 stack_image.cxx\
	 stack_bichip.cxx\
	 preprocess_dark.cxx\
         filter_image.cxx\
	 add_image.cxx\
	 image_flip.cxx\
	 dark_fit.cxx\
	 stack_dark.cxx


OPTSRC = substract_overscan.cxx\
	 substract_dark.cxx\
	 substract_flat.cxx\
	 build_flat.cxx\
	 add_poisson.cxx\
	 bichip_assemble.cxx\
	 bichip_sub_bias.cxx\
	 bichip_overscan.cxx\
	 bichip_variance.cxx\
	 stack_bichip_gain.cxx\
	 wr_desc_cat.cxx\
	 bichip_build.cxx\
	 hack_B_headers.cxx\
	 test_extension.cxx\
	 rickify.cxx\
	 lrerror.cxx\
	 preprocess_bias.cxx\
	 preprocess_flat.cxx\
	 preprocess_overscan.cxx\
	 quick_look.cxx\
#	 odd_even_correct.cxx\  # non-maintained


EXE = $(patsubst %.cxx, $(E)/%, $(SRC))
#OPTEXE = $(patsubst %.cxx, $(E)/%, $(OPTSRC))
#WEXE = $(patsubst %.cxx, $(E)/%, $(WSRC))

#OBJ = $(patsubst %.cxx, $(O)/%.o, $(SRC) )
#OBJ = $(patsubst %.cxx, $(O)/%.o, $(SRC) $(OPTSRC) )

# ============================================================

# Objects C++
$(O)/%.o : %.cxx $(HEADERS)
	$(CXXCOMPILE)

# Executables (beware: the % syntax might not be standard)
$(E)/% : $(O)/%.o $(LOCLIB) $(IOLIB)
	$(CXXLINK)
	$(STRIP) $@

all : dirbin dirobj $(EXE)

.PHONY : all
