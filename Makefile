
CC=     g++

CFLAGS= -O3 -D_USE_64 -msse4.2 -funroll-loops -fomit-frame-pointer

LFLAGS= -std=c++11 -DNDEBUG -lz -lm -lpthread -I . \
        -I ./vcflib/tabixpp/ -I ./vcflib/tabixpp/htslib/ -I ./vcflib/smithwaterman/ -I ./vcflib/multichoose/ -I ./vcflib/filevercmp/ -I ./vcflib/src/ \
        -L ./vcflib/ -L ./vcflib/tabixpp/htslib/ -lvcflib -lhts -Wl,-rpath=$(PWD)/vcflib/ -Wl,-rpath=$(PWD)/vcflib/tabixpp/htslib/

EXE=    edso

SRC=    main.cpp

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .cpp .o

OBJ=    $(SRC:.cpp=.o)

.cpp.o:
	$(CC) $(CFLAGS) -c $(LFLAGS) $<

all:    $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ): $(MF) $(HD)

clean:
	rm -f $(OBJ) $(EXE) *~

clean-all:
	rm -f $(OBJ) $(EXE) *~
	rm -rf vcflib
