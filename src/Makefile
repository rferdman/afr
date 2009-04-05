SHELL = /bin/sh

ifndef ASP_BIN_DIR
ASP_BIN_DIR = ../bin
endif

ASP_INCLUDE_DIR = ../include

#### Set directories here for fitstio*.h include file and ####
#### libcfitsio.a library files, respectively.            ####
#FITS_INCLUDE_DIR = /usr/include
#FITS_LIB_DIR = /usr/lib
FITS_INCLUDE_DIR = /sw/include
FITS_LIB_DIR = /sw/lib
#FITS_INCLUDE_DIR = /astro/mahler/gonzalez/include
#FITS_LIB_DIR = /astro/mahler/gonzalez/lib


CC = gcc
CFLAGS = -Wall -O -c -g -DASP -I$(ASP_INCLUDE_DIR) -I$(FITS_INCLUDE_DIR) 
CFLAGS_NO_OPT = -Wall -g -DASP -I$(ASP_INCLUDE_DIR) -I$(FITS_INCLUDE_DIR) 
CFLAGS_OPT2 = -Wall -O3 -g -DASP -I$(ASP_INCLUDE_DIR) -I$(FITS_INCLUDE_DIR) 
F77 = g77
F77FLAGS = -c -I$(ASP_INCLUDE_DIR)

ALLEXECS = ASPFitsReader ASPStokes ASPCal ASPThetaBB ASPBinDown ASPRotate \
	ASPToa ASPMatch ASPHead ASPAdd ASPRotateAsc


# Object file names:
#-------------------

OBJS1 = ASPFitsReader.o CmdLine.o GetOptions.o GetChans.o GetOmit.o Sum.o \
	Freq2Chan.o InitPars.o ReadCal.o GetCalData.o ReadASPHdr.o \
	GetPoly.o cprofc.o ReadASPData.o MakeStokes.o PhaseShift.o PhaseCalc.o \
	BinDown.o MakePol.o WriteStokes.o DutyLookup.o BMask.o Sort.o Zero.o \
	MaxMin.o FindPeak.o RotateProf.o FitAngle.o GetChi.o \
	FitThetaBB.o FitMueller.o WrtASPHdr.o WrtASPStokes.o \
	Median.o

OBJS2 = ffft.o fftfit.o cprof.o fccf.o brent.o 

STOKESOBJS = ASPStokes.o StokesCmdLine.o ReadASPHdr.o ReadASPStokes.o BMask.o \
	Sort.o FindPeak.o DutyLookup.o GetChi.o

CALOBJS = ASPCal.o CalCmdLine.o InitPars.o ReadASPHdr.o GetCalOpt.o \
	GetCalData.o ReadASPData.o Zero.o MaxMin.o Sort.o Median.o

THETAOBJS = ASPThetaBB.o ThetaCmdLine.o InitPars.o ReadASPHdr.o GetCalData.o \
	ReadASPData.o MaxMin.o Sort.o Median.o Zero.o  

BINDOWNOBJS = ASPBinDown.o BinDownCmdLine.o BinDown.o MakePol.o Zero.o \
	WriteStokes.o  DutyLookup.o Sort.o FindPeak.o BMask.o 

ROTATEOBJS2 = ASPRotateAsc.o RotateCmdLine.o RotateProf.o cprofc.o MakePol.o \
	Zero.o WriteStokes.o DutyLookup.o Sort.o FindPeak.o BMask.o ffft.o

ROTATEOBJS = ASPRotate.o RotateCmdLine.o RotateProf.o ReadASPHdr.o \
	ReadASPStokes.o BMask.o Sort.o FindPeak.o DutyLookup.o cprofc.o \
	MakePol.o Zero.o WrtASPHdr.o WrtASPStokes.o ffft.o

MATCHOBJS = ASPMatch.o MatchCmdLine.o ReadASPAsc.o RotateProf.o cprofc.o \
	MakePol.o Zero.o WriteStokes.o DutyLookup.o Sort.o FindPeak.o BMask.o 

HEADOBJS = ASPHead.o HeadCmdLine.o

ADDOBJS = ASPAdd.o AddCmdLine.o ReadASPHdr.o ReadASPStokes.o Sum.o MakePol.o \
	Zero.o WriteStokes.o DutyLookup.o Sort.o FindPeak.o BMask.o

TOAOBJS = ASPToa.o ToaCmdLine.o ReadASPAsc.o ReadASPHdr.o ReadASPStokes.o \
	cprofc.o

%.o: %.f
	$(F77) $(F77FLAGS) $<

%.o: %.c
	$(CC) $(CFLAGS) $<


# Object dependencies:
#-----------------------

ASPFitsReader: ${OBJS1} ${OBJS2} $(ASP_BIN_DIR)
	${CC} -o $(ASP_BIN_DIR)/ASPFitsReader ${OBJS1} ${OBJS2} \
	-L${FITS_LIB_DIR} -lcfitsio -lm -lg2c

ASPStokes: ${STOKESOBJS} $(ASP_BIN_DIR)
	${CC} -o $(ASP_BIN_DIR)/ASPStokes ${STOKESOBJS} \
	-L${FITS_LIB_DIR} -lcfitsio -lm -lg2c

ASPCal: ${CALOBJS} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPCal ${CALOBJS} -L${FITS_LIB_DIR} \
	-lcfitsio -lm -lg2c

ASPThetaBB: ${THETAOBJS} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPThetaBB ${THETAOBJS} \
	-L${FITS_LIB_DIR} -lcfitsio -lm -lg2c 

ASPBinDown: ${BINDOWNOBJS} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPBinDown ${BINDOWNOBJS} \
	-L${FITS_LIB_DIR} -lcfitsio -lm -lg2c 

ASPRotateAsc: ${ROTATEOBJS2} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPRotateAsc ${ROTATEOBJS2} \
	-L${FITS_LIB_DIR} -lcfitsio -lm -lg2c 

ASPRotate: ${ROTATEOBJS} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPRotate ${ROTATEOBJS} \
	-L${FITS_LIB_DIR} -lcfitsio -lm -lg2c 

ASPMatch: ${MATCHOBJS} ${OBJS2} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPMatch ${MATCHOBJS} ${OBJS2} \
	-L${FITS_LIB_DIR} -lcfitsio -lm -lg2c 

ASPHead: ${HEADOBJS} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPHead ${HEADOBJS} \
	-L${FITS_LIB_DIR} -lcfitsio -lm -lg2c 

ASPAdd: ${ADDOBJS} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPAdd ${ADDOBJS} \
	-L${FITS_LIB_DIR} -lcfitsio -lm -lg2c 

ASPToa: ${TOAOBJS} ${OBJS2} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPToa ${TOAOBJS} ${OBJS2} \
	-L${FITS_LIB_DIR} -lcfitsio -lm -lg2c 

CmdLine.c CmdLine.h: CmdLine.cli
	clig -d CmdLine.cli
	mv CmdLine.h $(ASP_INCLUDE_DIR)

CalCmdLine.c CalCmdLine.h: CalCmdLine.cli
	clig -d CalCmdLine.cli
	mv CalCmdLine.h $(ASP_INCLUDE_DIR)

StokesCmdLine.c StokesCmdLine.h: StokesCmdLine.cli
	clig -d StokesCmdLine.cli
	mv StokesCmdLine.h $(ASP_INCLUDE_DIR)

ThetaCmdLine.c ThetaCmdLine.h: ThetaCmdLine.cli
	clig -d ThetaCmdLine.cli
	mv ThetaCmdLine.h $(ASP_INCLUDE_DIR)

BinDownCmdLine.c BinDownCmdLine.h: BinDownCmdLine.cli
	clig -d BinDownCmdLine.cli
	mv BinDownCmdLine.h $(ASP_INCLUDE_DIR)

RotateCmdLine.c RotateCmdLine.h: RotateCmdLine.cli
	clig -d RotateCmdLine.cli
	mv RotateCmdLine.h $(ASP_INCLUDE_DIR)

MatchCmdLine.c MatchCmdLine.h: MatchCmdLine.cli
	clig -d MatchCmdLine.cli
	mv MatchCmdLine.h $(ASP_INCLUDE_DIR)

HeadCmdLine.c HeadCmdLine.h: HeadCmdLine.cli
	clig -d HeadCmdLine.cli
	mv HeadCmdLine.h $(ASP_INCLUDE_DIR)

AddCmdLine.c AddCmdLine.h: AddCmdLine.cli
	clig -d AddCmdLine.cli
	mv AddCmdLine.h $(ASP_INCLUDE_DIR)

ToaCmdLine.c ToaCmdLine.h: ToaCmdLine.cli
	clig -d ToaCmdLine.cli
	mv ToaCmdLine.h $(ASP_INCLUDE_DIR)

ASPFitsReader.o: ASPFitsReader.c CmdLine.h

ASPStokes.o: ASPStokes.c StokesCmdLine.h

ASPCal.o: ASPCal.c CalCmdLine.h

ASPThetaBB.o: ASPThetaBB.c ThetaCmdLine.h

ASPBinDown.o: ASPBinDown.c BinDownCmdLine.h

ASPRotateAsc.o: ASPRotateAsc.c RotateCmdLine.h

ASPRotate.o: ASPRotate.c RotateCmdLine.h

ASPMatch.o: ASPMatch.c MatchCmdLine.h

ASPHead.o: ASPHead.c HeadCmdLine.h

ASPAdd.o: ASPAdd.c AddCmdLine.h

ASPToa.o: ASPToa.c ToaCmdLine.h

$(ASP_BIN_DIR):
	mkdir $(ASP_BIN_DIR)

all: $(ALLEXECS)


afrclean:
	/bin/rm -f ${OBJS1} ${OBJS2} *~ PI* $(ASP_BIN_DIR)/ASPFitsReader rdb.* startup.* core* CmdLine.[ch] 

calclean:
	/bin/rm -f ${CALOBJS} $(ASP_BIN_DIR)/ASPCal core* *~ CalCmdLine.[ch]

thetaclean:
	/bin/rm -f ${THETAOBJS} $(ASP_BIN_DIR)/ASPThetaBB core* *~ ThetaCmdLine.[ch]

stokesclean:
	/bin/rm -f ${STOKESOBJS} $(ASP_BIN_DIR)/ASPStokes core* *~ StokesCmdLine.[ch]

bindownclean:
	/bin/rm -f ${BINDOWNOBJS} $(ASP_BIN_DIR)/ASPBinDown core* *~ BinDownCmdLine.[ch]

rotclean2:
	/bin/rm -f ${ROTATEOBJS2} $(ASP_BIN_DIR)/ASPRotateAsc core* *~ RotateCmdLine.[ch]

rotclean:
	/bin/rm -f ${ROTATEOBJS} $(ASP_BIN_DIR)/ASPRotate core* *~ RotateCmdLine.[ch]

matchclean:
	/bin/rm -f ${MATCHOBJS} ${OBJS2} $(ASP_BIN_DIR)/ASPMatch core* *~ MatchCmdLine.[ch]

headclean:
	/bin/rm -f ${HEADOBJS} $(ASP_BIN_DIR)/ASPHead core* *~ HeadCmdLine.[ch]

addclean:
	/bin/rm -f ${ADDOBJS} $(ASP_BIN_DIR)/ASPAdd core* *~ AddCmdLine.[ch]

toaclean:
	/bin/rm -f ${TOAOBJS} $(ASP_BIN_DIR)/ASPToa core* *~ ToaCmdLine.[ch]

includeclean:
	/bin/rm -f $(ASP_INCLUDE_DIR)/*CmdLine.h

oneclean:
	/bin/rm -f *.1

clean:
	make afrclean calclean thetaclean stokesclean bindownclean rotclean rotclean2 matchclean headclean addclean toaclean includeclean oneclean

#distclean:
#	/bin/rm -f Makefile config.h config.status config.cache config.log

