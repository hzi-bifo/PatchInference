ANTI = AntiPatch
SRC = $(ANTI)/src
MAX = $(SRC)/maxflow
BIN = $(ANTI)/bin
GPP = g++

antipatch: $(BIN)/graph.o $(BIN)/maxflow.o $(BIN)/Table.o $(BIN)/SubPatcher.o $(BIN)/Residue.o $(BIN)/AntiPatch.o $(BIN)/Printer.o $(SRC)/main.cpp clean
	$(GPP) -o antipatch $(SRC)/main.cpp $(BIN)/graph.o $(BIN)/maxflow.o $(BIN)/Table.o $(BIN)/SubPatcher.o $(BIN)/Residue.o $(BIN)/AntiPatch.o $(BIN)/Printer.o

$(BIN)/graph.o: $(MAX)/graph.h $(MAX)/graph.cpp
	mkdir -p $(BIN)
	$(GPP) -c $(MAX)/graph.h $(MAX)/graph.cpp
	mv graph.o $(BIN)/graph.o

$(BIN)/maxflow.o: $(MAX)/maxflow.cpp
	g++ -c $(MAX)/maxflow.cpp
	mv maxflow.o $(BIN)/maxflow.o

$(BIN)/%.o: $(SRC)/%.h $(SRC)/%.cpp 
	$(GPP) -c $(SRC)/$*.h $(SRC)/$*.cpp
	mv $*.o $(BIN)/$*.o

clean:
	rm -f $(MAX)/*.gch
	rm -f $(SRC)/*.gch

dist-clean: clean
	rm -rf $(BIN)
