GCCFLAGS=-Wall -g
ROOTFLAGS=`root-config --cflags --libs`

SELECTORFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lProof
FIGURESFLAGS=$(GCCFLAGS) $(ROOTFLAGS) 

all: sel fig presel
	@echo "Full build successful."

sel: SingleBuToKstarMuMuSelector.cc 
	g++ $(SELECTORFLAGS) $< -o $@

fig: SingleBuToKstarMuMuFigures.cc 
	g++ $(FIGURESFLAGS) $< -o $@

presel: draw_preselvars.cc 
	g++ $(FIGURESFLAGS) $< -o $@

clean:
	rm -rf SingleBuToKstarMuMuSelector SingleBuToKstarMuMuFigures 
