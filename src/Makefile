SHELL = /bin/sh

ifndef ASP_BIN_DIR
ASP_BIN_DIR = ../bin
endif


ASP_INCLUDE_DIR = ../include


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

### SPECIFY X11 LIBRARY DIRECTASPORY HERE ###
X11_LIB_DIR = /usr/X11R6/lib

LCPGPLOT = -L$(PGPLOT_LIB_DIR) -lcpgplot -lpgplot  -L$(X11_LIB_DIR) -lX11 -lpng -pg -lgfortran # -L/sw/lib -lg2c -pg


### C COMPILER AND FLAGS ###
CC = gcc
CFLAGS = -Wall -c -g -DASP -I$(ASP_INCLUDE_DIR) -I$(FITS_INCLUDE_DIR) -I$(PGPLOT_INCLUDE_DIR)
# CFLAGS_NO_OPT = -Wall -g -DASP -I$(ASP_INCLUDE_DIR) -I$(FITS_INCLUDE_DIR) -I$(PGPLOT_INCLUDE_DIR)
#CFLAGS_OPT2 = -Wall -O3 -g -DASP -I$(ASP_INCLUDE_DIR) -I$(FITS_INCLUDE_DIR)  -I$(PGPLOT_INCLUDE_DIR)


### GFORTRAN FLAGS ###
F77 = gfortran
F77FLAGS = -c -Wall -ffixed-form -ffixed-line-length-none -fPIC -std=legacy -O -I$(ASP_INCLUDE_DIR)

### G77 FLAGS ###
# F77 = g77
# F77FLAGS = -c -Wall -I$(ASP_INCLUDE_DIR)



ALLEXECS = ASPFitsReader ASPStokes ASPCal ASPThetaBB ASPBinDown ASPRotate \
	ASPZap ASPToa ASPMatch ASPHead ASPAdd ASPDiff ASPTemplate ASPRotateAsc \
	GetPhases GetWidths


# Object file names:
#-------------------

OBJS1 = ASPFitsReader.o CmdLine.o AllocRunMode.o GetOptions.o GetChans.o \
	GetOmit.o Sum.o Freq2Chan.o InitPars.o ReadCal.o \
	ReadHdr.o ReadASPHdr.o ReadPSRFITSHdr.o GetTelescope.o GetPoly.o \
	cprofc.o ReadData.o ReadASPData.o ReadPSRFITSData.o \
	Dedisperse.o MakeStokes.o PhaseShift.o PhaseCalc.o BinDown.o \
	MakePol.o WriteStokes.o DutyLookup.o BMask.o Sort.o Zero.o MaxMin.o \
	FindPeak.o RotateProf.o FitAngle.o GetChi.o GetLST.o FitThetaBB.o \
	FitMueller.o WrtASPHdr.o WrtASPStokes.o Median.o

OBJS2 = ffft.o fftfit.o cprof.o fccf.o brent.o 

STOKESOBJS = ASPStokes.o StokesCmdLine.o ReadASPHdr.o GetTelescope.o \
	ReadASPStokes.o BMask.o Sort.o FindPeak.o DutyLookup.o GetChi.o \
	GetLST.o

CALOBJS = ASPCal.o CalCmdLine.o InitPars.o ReadHdr.o ReadASPHdr.o \
	ReadPSRFITSHdr.o GetTelescope.o AllocRunMode.o GetCalOpt.o BMask.o \
	FindPeak.o DutyLookup.o MakePol.o Zero.o GetCalData.o ReadData.o \
	ReadASPData.o ReadPSRFITSData.o PhaseCalc.o MaxMin.o Sort.o Median.o

THETAOBJS = ASPThetaBB.o ThetaCmdLine.o InitPars.o ReadHdr.o ReadASPHdr.o \
	ReadPSRFITSHdr.o GetTelescope.o AllocRunMode.o GetCalData.o \
	ReadData.o ReadASPData.o ReadPSRFITSData.o PhaseCalc.o MaxMin.o \
	Sort.o Median.o Zero.o Freq2Chan.o

BINDOWNOBJS = ASPBinDown.o BinDownCmdLine.o AllocRunMode.o BinDown.o \
	MakePol.o Zero.o WriteStokes.o  DutyLookup.o Sort.o FindPeak.o BMask.o 

ROTATEOBJS2 = ASPRotateAsc.o RotateCmdLine.o AllocRunMode.o RotateProf.o \
	cprofc.o MakePol.o Zero.o WriteStokes.o DutyLookup.o Sort.o \
	FindPeak.o BMask.o ffft.o 

ROTATEOBJS = ASPRotate.o RotateCmdLine.o AllocRunMode.o RotateProf.o \
	ReadASPHdr.o ReadASPStokes.o BMask.o Sort.o FindPeak.o DutyLookup.o \
	cprofc.o MakePol.o Zero.o WrtASPHdr.o WrtASPStokes.o ffft.o

ZAPOBJS = ASPZap.o ZapCmdLine.o AllocRunMode.o ReadASPAsc.o ReadHdr.o \
	ReadASPHdr.o ReadPSRFITSHdr.o ReadData.o ReadASPData.o \
	ReadPSRFITSData.o Dedisperse.o GetPoly.o PhaseShift.o MakeStokes.o \
	BMask.o Sort.o FindPeak.o \
	DutyLookup.o MakePol.o Zero.o MaxMin.o ArrayZero.o BinDown.o cprofc.o \
	RotateProf.o PhaseCalc.o

MATCHOBJS = ASPMatch.o MatchCmdLine.o AllocRunMode.o ReadASPAsc.o \
	RotateProf.o cprofc.o MakePol.o Zero.o WriteStokes.o DutyLookup.o \
	Sort.o FindPeak.o BMask.o 

HEADOBJS = ASPHead.o HeadCmdLine.o

ADDOBJS = ASPAdd.o AddCmdLine.o ReadASPHdr.o AllocRunMode.o ReadASPStokes.o \
	ArrayZero.o MakePol.o Zero.o WriteStokes.o DutyLookup.o Sort.o \
	FindPeak.o BMask.o

DIFFOBJS = ASPDiff.o DiffCmdLine.o AllocRunMode.o RotateProf.o ReadASPAsc.o \
	ReadASPHdr.o ReadASPStokes.o BMask.o Sort.o FindPeak.o DutyLookup.o \
	cprofc.o MakePol.o Zero.o ArrayZero.o WrtASPHdr.o WrtASPStokes.o

TEMPLATEOBJS = ASPTemplate.o TemplateCmdLine.o ReadASPHdr.o AllocRunMode.o \
	ReadASPAsc.o cprofc.o ReadASPStokes.o ArrayZero.o RotateProf.o \
	MakePol.o Zero.o Median.o WriteStokes.o DutyLookup.o Sort.o \
	FindPeak.o BMask.o

TOAOBJS = ASPToa.o ToaCmdLine.o ReadASPAsc.o ReadASPHdr.o ReadASPStokes.o \
	Sum.o cprofc.o

PHASEOBJS = GetPhases.o PhaseCmdLine.o GetPoly.o PhaseCalc.o # FileExists.o

WIDTHOBJS = GetWidths.o WidthCmdLine.o ReadASPAsc.o FindPeak.o DutyLookup.o \
	BMask.o Sort.o MaxMin.o 

TESTOBJS = TestPSRFITS.o TestCmdLine.o InitPars.o ReadHdr.o ReadASPHdr.o \
	ReadPSRFITSHdr.o GetTelescope.o AllocRunMode.o ReadData.o \
	ReadASPData.o ReadPSRFITSData.o Dedisperse.o Zero.o PhaseCalc.o \
	RotateProf.o cprofc.o plot_profile.o 

%.o: %.f
	$(F77) $(F77FLAGS) $<

%.o: %.c
	$(CC) $(CFLAGS) $<


# Object dependencies:
#-----------------------

ASPFitsReader: ${OBJS1} ${OBJS2} $(ASP_BIN_DIR)
	${CC} -o $(ASP_BIN_DIR)/ASPFitsReader ${OBJS1} ${OBJS2} \
	${LCFITSIO}  

ASPStokes: ${STOKESOBJS} $(ASP_BIN_DIR)
	${CC} -o $(ASP_BIN_DIR)/ASPStokes ${STOKESOBJS} \
	${LCFITSIO}  

ASPCal: ${CALOBJS} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPCal ${CALOBJS} -L${FITS_LIB_DIR} \
	-lcfitsio -lm  

ASPThetaBB: ${THETAOBJS} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPThetaBB ${THETAOBJS} \
	${LCFITSIO}   

ASPBinDown: ${BINDOWNOBJS} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPBinDown ${BINDOWNOBJS} \
	${LCFITSIO}   

ASPRotateAsc: ${ROTATEOBJS2} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPRotateAsc ${ROTATEOBJS2} \
	${LCFITSIO}   

ASPRotate: ${ROTATEOBJS} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPRotate ${ROTATEOBJS} \
	${LCFITSIO}   

ASPZap: ${ZAPOBJS} ${OBJS2} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPZap ${ZAPOBJS} ${OBJS2} \
	${LCFITSIO} $(LCPGPLOT)

ASPMatch: ${MATCHOBJS} ${OBJS2} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPMatch ${MATCHOBJS} ${OBJS2} \
	${LCFITSIO}   

ASPHead: ${HEADOBJS} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPHead ${HEADOBJS} \
	${LCFITSIO}   

ASPAdd: ${ADDOBJS} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPAdd ${ADDOBJS} \
	${LCFITSIO}   

ASPDiff: ${DIFFOBJS} ${OBJS2} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPDiff ${DIFFOBJS} ${OBJS2}  \
	${LCFITSIO}   

ASPTemplate: ${TEMPLATEOBJS} ${OBJS2} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPTemplate ${TEMPLATEOBJS} ${OBJS2} \
	${LCFITSIO}   

ASPToa: ${TOAOBJS} ${OBJS2} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/ASPToa ${TOAOBJS} ${OBJS2} \
	${LCFITSIO}  

GetPhases: ${PHASEOBJS} ${SAP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/GetPhases ${PHASEOBJS} \
	$(LCPGPLOT)

GetWidths: ${WIDTHOBJS} ${SAP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/GetWidths ${WIDTHOBJS} \
	$(LCPGPLOT)

TestPSRFITS: ${TESTOBJS} ${OBJS2} ${ASP_BIN_DIR}
	${CC} -o $(ASP_BIN_DIR)/TestPSRFITS ${TESTOBJS} ${OBJS2} \
	${LCFITSIO}   $(LCPGPLOT)

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

ZapCmdLine.c ZapCmdLine.h: ZapCmdLine.cli
	clig -d ZapCmdLine.cli
	mv ZapCmdLine.h $(ASP_INCLUDE_DIR)

MatchCmdLine.c MatchCmdLine.h: MatchCmdLine.cli
	clig -d MatchCmdLine.cli
	mv MatchCmdLine.h $(ASP_INCLUDE_DIR)

HeadCmdLine.c HeadCmdLine.h: HeadCmdLine.cli
	clig -d HeadCmdLine.cli
	mv HeadCmdLine.h $(ASP_INCLUDE_DIR)

AddCmdLine.c AddCmdLine.h: AddCmdLine.cli
	clig -d AddCmdLine.cli
	mv AddCmdLine.h $(ASP_INCLUDE_DIR)

DiffCmdLine.c DiffCmdLine.h: DiffCmdLine.cli
	clig -d DiffCmdLine.cli
	mv DiffCmdLine.h $(ASP_INCLUDE_DIR)

TemplateCmdLine.c TemplateCmdLine.h: TemplateCmdLine.cli
	clig -d TemplateCmdLine.cli
	mv TemplateCmdLine.h $(ASP_INCLUDE_DIR)

ToaCmdLine.c ToaCmdLine.h: ToaCmdLine.cli
	clig -d ToaCmdLine.cli
	mv ToaCmdLine.h $(ASP_INCLUDE_DIR)

PhaseCmdLine.c PhaseCmdLine.h: PhaseCmdLine.cli
	clig -d PhaseCmdLine.cli
	mv PhaseCmdLine.h $(ASP_INCLUDE_DIR)

WidthCmdLine.c WidthCmdLine.h: WidthCmdLine.cli
	clig -d WidthCmdLine.cli
	mv WidthCmdLine.h $(ASP_INCLUDE_DIR)

TestCmdLine.c TestCmdLine.h: TestCmdLine.cli
	clig -d TestCmdLine.cli
	mv TestCmdLine.h $(ASP_INCLUDE_DIR)

ASPFitsReader.o: ASPFitsReader.c CmdLine.h

ASPStokes.o: ASPStokes.c StokesCmdLine.h

ASPCal.o: ASPCal.c CalCmdLine.h

ASPThetaBB.o: ASPThetaBB.c ThetaCmdLine.h

ASPBinDown.o: ASPBinDown.c BinDownCmdLine.h

ASPRotateAsc.o: ASPRotateAsc.c RotateCmdLine.h

ASPRotate.o: ASPRotate.c RotateCmdLine.h

ASPZap.o: ASPZap.c ZapCmdLine.h

ASPMatch.o: ASPMatch.c MatchCmdLine.h

ASPHead.o: ASPHead.c HeadCmdLine.h

ASPAdd.o: ASPAdd.c AddCmdLine.h

ASPDiff.o: ASPDiff.c DiffCmdLine.h

ASPTemplate.o: ASPTemplate.c TemplateCmdLine.h

ASPToa.o: ASPToa.c ToaCmdLine.h

GetPhases.o: GetPhases.c PhaseCmdLine.h

GetWidths.o: GetWidths.c WidthCmdLine.h

TestPSRFITS.o: TestPSRFITS.c TestCmdLine.h

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

zapclean:
	/bin/rm -f ${ZAPOBJS} $(ASP_BIN_DIR)/ASPZap core* *~ ZapCmdLine.[ch]

matchclean:
	/bin/rm -f ${MATCHOBJS} ${OBJS2} $(ASP_BIN_DIR)/ASPMatch core* *~ MatchCmdLine.[ch]

headclean:
	/bin/rm -f ${HEADOBJS} $(ASP_BIN_DIR)/ASPHead core* *~ HeadCmdLine.[ch]

addclean:
	/bin/rm -f ${ADDOBJS} $(ASP_BIN_DIR)/ASPAdd core* *~ AddCmdLine.[ch]

diffclean:
	/bin/rm -f ${DIFFOBJS} $(ASP_BIN_DIR)/ASPDiff core* *~ DiffCmdLine.[ch]

templateclean:
	/bin/rm -f ${TEMPLATEOBJS} ${OBJS2} $(ASP_BIN_DIR)/ASPTemplate core* *~ TemplateCmdLine.[ch]

toaclean:
	/bin/rm -f ${TOAOBJS} $(ASP_BIN_DIR)/ASPToa core* *~ ToaCmdLine.[ch]

phaseclean:
	/bin/rm -f ${PHASEOBJS} $(ASP_BIN_DIR)/GetPhases core* *~ PhaseCmdLine.[ch]

widthclean:
	/bin/rm -f ${WIDTHOBJS} $(ASP_BIN_DIR)/GetWidths core* *~ WidthCmdLine.[ch]

testclean:
	/bin/rm -f ${TESTOBJS} $(ASP_BIN_DIR)/TestPSRFITS core* *~ TestCmdLine.[ch]

includeclean:
	/bin/rm -f $(ASP_INCLUDE_DIR)/*CmdLine.h

oneclean:
	/bin/rm -f *.1

clean:
	make afrclean calclean thetaclean stokesclean bindownclean rotclean rotclean2 zapclean matchclean headclean addclean diffclean templateclean toaclean includeclean oneclean widthclean phaseclean

#distclean:
#	/bin/rm -f Makefile config.h config.status config.cache config.log

