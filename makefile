# 02-03-2018
# billy

EXE = hid # Executable
OBJ = hidrogen.o prova_hidrogen.o # Objectes
COM = gcc # Compilador
LIB = -lm # Llibreries (-l*, -L*, -I*)
MAC = -D_GNU_SOURCE # Macros (-D*)
AVS = -W -Wall -Wextra -ansi -pedantic # Avisos
OPT = -march=native -O2 -pipe # Optimitzacio
DEP = -g -DDEBUG # Depuracio, no recomanable junt amb $(OPT)
OPC = $(DEP) $(AVS) $(MAC) -std=c11 # Opcions
DIR = /usr/local/bin # Directori per a instalar

all: $(EXE)

$(EXE): $(OBJ)
	@echo Enllaçant $(OBJ) per a crear $(EXE)
	$(COM) $(LIB) $(OBJ) -o $@

hidrogen.o: hidrogen.c hidrogen.h
	@echo Compilant $<
	$(COM) $(OPC) -c $<

prova_hidrogen.o: prova_hidrogen.c hidrogen.h
	@echo Compilant $<
	$(COM) $(OPC) -c $<

exe: all
	./$(EXE)

install: all
	mkdir -p $(DIR)
	cp $(EXE) $(DIR)
	chown root:staff $(EXE)

clean:
	@echo Netejant...
	rm -rf $(EXE) $(OBJ) *~
