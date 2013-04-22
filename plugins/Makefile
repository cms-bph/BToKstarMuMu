BINDIR=$(HOME)/local/bin
FLAGS=-Wall
FLAGS+=`root-config --cflags --libs`
FLAGS+=-lProof


SingleBuToKstarMuMuSelector: SingleBuToKstarMuMuSelector.cc 
	g++ $(FLAGS) SingleBuToKstarMuMuSelector.cc -o $(BINDIR)/SingleBuToKstarMuMuSelector

clean:
	rm -rf $(BINDIR)/SingleBuToKstarMuMuSelector 