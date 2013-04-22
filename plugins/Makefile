#BINDIR=$(HOME)/local/bin
FLAGS=-Wall
FLAGS+=`root-config --cflags --libs`
FLAGS+=-lProof

all: SingleBuToKstarMuMuSelector
	@echo "Full build successful."

SingleBuToKstarMuMuSelector: SingleBuToKstarMuMuSelector.cc 
	g++ $(FLAGS) $< -o $@

clean:
	rm -rf SingleBuToKstarMuMuSelector 
