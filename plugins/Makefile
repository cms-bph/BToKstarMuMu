GCCFLAGS=-Wall -g
ROOTFLAGS=`root-config --cflags --libs`

SELECTORFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lProof
FIGURESFLAGS=$(GCCFLAGS) $(ROOTFLAGS) 

all: SingleBuToKstarMuMuSelector SingleBuToKstarMuMuFigures
	@echo "Full build successful."

SingleBuToKstarMuMuSelector: SingleBuToKstarMuMuSelector.cc 
	g++ $(SELECTORFLAGS) $< -o $@

SingleBuToKstarMuMuFigures: SingleBuToKstarMuMuFigures.cc 
	g++ $(FIGURESFLAGS) $< -o $@

clean:
	rm -rf SingleBuToKstarMuMuSelector SingleBuToKstarMuMuFigures 
