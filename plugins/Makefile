GCCFLAGS=-Wall -g
ROOTFLAGS=`root-config --cflags --libs`

SELECTORFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lProof
FIGURESFLAGS=$(GCCFLAGS) $(ROOTFLAGS) 
FITFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lRooFit -lRooFitCore 

all: sel fig fit  
	@echo "Full build successful."

sel: SingleBuToKstarMuMuSelector.cc 
	g++ $(SELECTORFLAGS) $< -o $@

fig: SingleBuToKstarMuMuFigures.cc 
	g++ $(FIGURESFLAGS) $< -o $@

fit: SingleBuToKstarMuMuFitters.cc 
	g++ $(FITFLAGS) $< -o $@	

clean:
	rm -rf sel fig fit 
