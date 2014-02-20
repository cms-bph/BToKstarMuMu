GCCFLAGS=-Wall -g
ROOTFLAGS=`root-config --cflags --libs`

TOOLSFLAGS=$(GCCFLAGS) -I$(ROOTSYS)/include -c 
SELECTORFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lProof
FIGURESFLAGS=$(GCCFLAGS) $(ROOTFLAGS) 
FITFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lRooFit -lRooFitCore -lMinuit

all: sel fig fit  
	@echo "Full build successful."


tools.o: tools.cc tools.h 
	g++ $(TOOLSFLAGS) $< 

sel: SingleBuToKstarMuMuSelector.cc 
	g++ $(SELECTORFLAGS) $< -o $@

fig: SingleBuToKstarMuMuFigures.cc tools.o 
	g++ $(FIGURESFLAGS) $^ -o $@

fit: SingleBuToKstarMuMuFitters.cc tools.o 
	g++ $(FITFLAGS) $^ -o $@	

sel2011: SingleBdToKstarMuMuSelector_xCheck2011.cc
	g++ ${SELECTORFLAGS} $< -o $@

test: testMinuitFitter.cc
	g++ ${FITFLAGS} $< -o $@


clean:
	rm -rf tools.o sel fig fit 
