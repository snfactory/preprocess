topdir = ..
include $(topdir)/makedefs
include $(topdir)/makedefs.local

# ============================================================

# Targets

SRC = utils.c

SRCC = bichip.cxx\
	catorfile.cxx\
	image.cxx \
	imagesnifs.cxx \
	section.cxx \
	sectionlist.cxx \
	analyser.cxx \
	ioplain.cxx \
	ioslice.cxx \
	imagestack.cxx \
	kombinator.cxx \
	kombinatorfit.cxx \
	imagestacksnifs.cxx \
	valuegetter.cxx \
	preprocessor.cxx \
	algocams.cxx\
	overscan.cxx\
	filter.cxx\
	darkmodel.cxx\

COBJ = $(patsubst %.c, $(O)/%.o, $(SRC))
CCOBJ = $(patsubst %.cxx, $(O)/%.o, $(SRCC))
OBJ = $(COBJ) $(CCOBJ)

LIB = $(O)/libpreprocess.a 

# ============================================================

OUT : dirobj $(LIB)

# Objects
$(COBJ) : $(O)/%.o : %.c makefile
	$(C_COMPILE)

# Objects C++
$(CCOBJ) : $(O)/%.o : %.cxx makefile
	$(CXXCOMPILE)

# Library
$(LIB) : $(OBJ)
	$(AR) rv $(LIB) $(OBJ)
	$(RANLIB) $@

.PHONY : OUT all
