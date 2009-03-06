SHELL = /bin/sh

subdirs = src

ifndef ASP_BIN_DIR
ASP_BIN_DIR = ./bin
endif

ASP_INCLUDE_DIR = ./include

#### Set directories here for fitsio*.h include file and  ####
#### libcfitsio.a library files, respectively.            ####
#FITS_INCLUDE_DIR = /usr/include
#FITS_LIB_DIR = /usr/lib
#FITS_INCLUDE_DIR = /sw/include
#FITS_LIB_DIR = /sw/lib
FITS_INCLUDE_DIR = /astro/mahler/gonzalez/include
FITS_LIB_DIR = /astro/mahler/gonzalez/lib

CC = gcc 
CFLAGS = -Wall -O -c -g -DASP -I$(ASP_INCLUDE_DIR) -I$(FITS_INCLUDE_DIR) 
CFLAGS_NO_OPT = -Wall -g -DASP -I$(ASP_INCLUDE_DIR) -I$(FITS_INCLUDE_DIR) 
CFLAGS_OPT2 = -Wall -O3 -g -DASP -I$(ASP_INCLUDE_DIR) -I$(FITS_INCLUDE_DIR) 
F77 = g77
F77FLAGS = -c -I$(ASP_INCLUDE_DIR)



#
# Object dependencies:
#-----------------------


all:
	@for dir in ${subdirs}; do \
	  (cd $$dir && $(MAKE) all) \
	  || case "$(MFLAGS)" in *k*) fail=yes;; *) exit 1;; esac; \
	done && test -z "$$fail"

$(ASP_BIN_DIR):
	mkdir $(ASP_BIN_DIR)


clean:
	/bin/rm -f *~
	@for dir in ${subdirs}; do \
	  (cd $$dir && $(MAKE) clean) \
	  || case "$(MFLAGS)" in *k*) fail=yes;; *) exit 1;; esac; \
	done && test -z "$$fail"

#distclean: clean
#	/bin/rm -f Makefile config.h config.status config.cache config.log
#	@for dir in ${subdirs}; do \
#	  (cd $$dir && $(MAKE) distclean) \
#	  || case "$(MFLAGS)" in *k*) fail=yes;; *) exit 1;; esac; \
#	done && test -z "$$fail"

