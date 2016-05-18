topdir = .
include $(topdir)/makedefs

all : libsrc src

libsrc :
	$(MAKE) -C libsrc

src : libsrc
	$(MAKE) -C src

# rmdir $(O) is here because it is written to from both src and libsrc
clean :
	cd libsrc && $(MAKE) clean
	cd src && $(MAKE) clean
	-rmdir $(O)

.PHONY : libsrc src all
