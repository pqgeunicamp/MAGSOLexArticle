
# compilador
COMPILADOR = g++ -std=c++2a

# arquivos .cpp da pasta src
CPPS = $(wildcard src/*.cpp)

# arquivos .cpp úteis do projeto RECMOL
CPPSRECMOL = \
             src/repomme/Mistura.cpp \
			 src/repomme/MG.cpp \
			 src/repomme/Randomic.cpp \

# nome e caminho do binário gerado
BIN = bin/MAGSOLexExamples 
# -Iinclude 

# flags
FLAGS = -Wall\
        -g\
		-Wno-unknown-pragmas\
		-O0\
		-Wno-unused-variable\
		-Wno-unused-but-set-variable


# includes: variável que pode indicar o caminho de pastas que contém arquivos .h ou .hpp do sistema ou do usuário
INC = -I/usr/local/include/eigen3

# executando ações definidas
all: info_in command info_out

info_in:
	@echo "Compilando codigo ... Aguarde ..."

command:
	@$(COMPILADOR) $(INC) $(FLAGS) $(CPPS) $(CPPSRECMOL) -o $(BIN) 

info_out:
	@echo "Binario construido com sucesso"
