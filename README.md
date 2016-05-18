preprocess
==========

*CCD-level processing in the Nearby Supernova Factory pipeline*

Dependencies:

- wcslib
- cfitsio
- gsl
- [ifuio](http://github.com/snfactory/ifuio)

## Build

It is assumed that the `ifuio` libaries are not installed, so we need
to tell make where to find the static libaries and header files. This
is done by creating a `make.user` file with contents:

```
IFUIO_LIBS = /path/to/ifuio/obj/libio.a /path/to/ifuio/obj/libgen.a
IFUIO_INC = /path/to/ifuio/include
```

Be sure to use *absolute paths*, as these variables are used when
running make from subdirectories.

Once the `make.user` file is created, type `make`.

**TODO:** automated install

## License & author

The license for the code here is MIT, but as two dependencies (gsl and ifuio)
are GPL-licensed, the use of this code is governed by the GPL.

The primary author of the code is Manu Ganglier.