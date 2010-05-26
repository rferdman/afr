SHELL = /bin/sh

subdirs = src

ifndef ASP_BIN_DIR
ASP_BIN_DIR = ./bin
endif

ASP_INCLUDE_DIR = ./include

### SET TO YOUR PREFERRED BIN DIRECTORY FOR INSTALLING AFR EXECUTABLES ###
INSTALL_BIN_DIR = /usr/local/bin


### SPECIFY CFITSIO INCLUDE FILE AND LIBRARY DIRECTORIES HERE ###
FITS_INCLUDE_DIR = /usr/local/include
FITS_LIB_DIR = /usr/local/lib

### IF USING GRFORTRAN, UNCOMMENT THE FOLLWOING ###
LCFITSIO = -L${FITS_LIB_DIR} -lcfitsio -lm 
 
### IF USING G77, UNCOMMENT THE FOLLOWING:
#LCFITSIO = -L${FITS_LIB_DIR} -lcfitsio -lm -lg2c


### SPECIFY PGPLOT INCLUDE AND LIBRARY DIRECTORIES HERE ###
PGPLOT_INCLUDE_DIR = /usr/local/pgplot
PGPLOT_LIB_DIR = /usr/local/pgplot

### SPECIFY X11 LIBRARY DIRECTORY HERE ###
X11_LIB_DIR = /usr/X11R6/lib

LCPGPLOT = -L$(PGPLOT_LIB_DIR) -lcpgplot -lpgplot  -L$(X11_LIB_DIR) -lX11 -lpng # -L/sw/lib -lg2c -pg


### C COMPILER AND FLAGS ###
CC = gcc
CFLAGS = -Wall -O -c -g -DASP -I$(ASP_INCLUDE_DIR) -I$(FITS_INCLUDE_DIR) -I$(PGPLOT_INCLUDE_DIR)
# CFLAGS_NO_OPT = -Wall -g -DASP -I$(ASP_INCLUDE_DIR) -I$(FITS_INCLUDE_DIR) -I$(PGPLOT_INCLUDE_DIR)
#CFLAGS_OPT2 = -Wall -O3 -g -DASP -I$(ASP_INCLUDE_DIR) -I$(FITS_INCLUDE_DIR)  -I$(PGPLOT_INCLUDE_DIR)


### GFORTRAN FLAGS ###
F77 = gfortran
F77FLAGS = -c -Wall -ffixed-line-length-none -fPIC -O -I$(ASP_INCLUDE_DIR)

### G77 FLAGS ###
# F77 = g77
# F77FLAGS = -c -Wall -I$(ASP_INCLUDE_DIR)



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

install: $(INSTALL_BIN_DIR)
	cd $(ASP_BIN_DIR); ln -sf `pwd`/ASP* $(INSTALL_BIN_DIR)