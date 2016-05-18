topdir = ..
include $(topdir)/makedefs

CSRC = utils.c
CXXSRC = bichip.cxx \
	catorfile.cxx \
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
	algocams.cxx \
	overscan.cxx \
	filter.cxx \
	darkmodel.cxx

COBJ = $(patsubst %.c, $(O)/%.o, $(CSRC))
CXXOBJ = $(patsubst %.cxx, $(O)/%.o, $(CXXSRC))
ALLOBJ = $(COBJ) $(CXXOBJ)

all : $(LIB)

$(LIB) : $(ALLOBJ)
	$(AR) rv $(LIB) $(ALLOBJ)
	$(RANLIB) $@

$(COBJ) : $(O)/%.o : %.c $(HEADERS) $(O)
	$(CCOMPILE)

$(CXXOBJ) : $(O)/%.o : %.cxx $(HEADERS) $(O)
	$(CXXCOMPILE)

$(O) :
	mkdir $(O)

clean :
	$(RM) $(ALLOBJ) $(LIB)

.PHONY : all clean
