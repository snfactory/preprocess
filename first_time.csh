#/bin/csh -f
aclocal
autoconf
automake -a --foreign
./configure 
make
