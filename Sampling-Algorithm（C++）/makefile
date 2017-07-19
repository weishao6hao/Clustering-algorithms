NAME = sampling
COMPILER= g++

BIN=bin
INCLUDE = -I.
CCFLAGS=$(INCLUDE) -Wall -O2

OBJECTS= $(BIN)/ProjetPro.o $(BIN)/Dides.o $(BIN)/Dendis.o $(BIN)/ADN_Esclave1.o  $(BIN)/ADN_Esclave2.o $(BIN)/ADN_Esclave3.o $(BIN)/ADN_Esclave4.o $(BIN)/ADN_util.o

$(BIN)/$(NAME) : $(OBJECTS)  
	$(COMPILER) -o $(BIN)/$(NAME)  $(OPT) $(OBJECTS) -lm 

$(BIN)/ProjetPro.o : ProjetPro.cpp ADN_Esclave1.h ADN_Esclave2.h ADN_Esclave3.h ADN_Esclave4.h ADN_util.h
	mkdir -p $(BIN)
	$(COMPILER) -c $(CCFLAGS) -o $(BIN)/ProjetPro.o ProjetPro.cpp

$(BIN)/Dides.o : Dides.cpp ADN_Esclave1.h ADN_Esclave2.h ADN_util.h DIDES.h
	$(COMPILER) -c $(CCFLAGS) -o $(BIN)/Dides.o Dides.cpp

$(BIN)/Dendis.o : Dendis.cpp ADN_Esclave1.h ADN_Esclave2.h ADN_util.h DENDIS.h
	$(COMPILER) -c $(CCFLAGS) -o $(BIN)/Dendis.o Dendis.cpp

$(BIN)/ADN_Esclave1.o : ADN_Esclave1.cpp ADN_Esclave1.h ADN_util.h 
	$(COMPILER) -c $(CCFLAGS) -o $(BIN)/ADN_Esclave1.o ADN_Esclave1.cpp

$(BIN)/ADN_Esclave2.o : ADN_Esclave2.cpp ADN_Esclave1.h ADN_Esclave2.h ADN_util.h
	$(COMPILER) -c $(CCFLAGS) -o $(BIN)/ADN_Esclave2.o ADN_Esclave2.cpp

$(BIN)/ADN_Esclave3.o : ADN_Esclave3.cpp ADN_Esclave1.h ADN_Esclave2.h ADN_Esclave3.h ADN_util.h
	$(COMPILER) -c $(CCFLAGS) -o $(BIN)/ADN_Esclave3.o ADN_Esclave3.cpp

$(BIN)/ADN_Esclave4.o : ADN_Esclave4.cpp ADN_Esclave1.h ADN_Esclave2.h ADN_Esclave3.h ADN_Esclave4.h ADN_util.h
	$(COMPILER) -c $(CCFLAGS) -o $(BIN)/ADN_Esclave4.o ADN_Esclave4.cpp

$(BIN)/ADN_util.o : ADN_util.cpp ADN_Esclave1.h ADN_util.h 
	$(COMPILER) -c $(CCFLAGS) -o $(BIN)/ADN_util.o ADN_util.cpp

clean : 
	rm -f $(OBJECTS) $(BIN)/$(NAME)*
