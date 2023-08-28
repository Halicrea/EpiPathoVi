#######################################################
#           Editable options   				          #
#######################################################
PROGR_NAME	= EpiVib
CREATE_SEQ	= Finding_context_of_motif
SUBMIT_NAME	= Bioinfo_EpiVib_GitLab_M2-ElynaB.tar.gz
#######################################################
#           Folders 					       	      #
#######################################################
SRC_DIR		= ./2_src
BIN_DIR		= ./0_bin
OBJECT_DIR	= ./1_obj

SRC_LIST	= $(wildcard $(SRC_DIR)/*.cpp)
OBJ_LIST	= $(BIN_DIR)/$(notdir $(SRC_LIST:.cpp=.o))

BIN			= 0_bin/$(CREATE_SEQ)
OBJS_GEN	= 
OBJS_ALIGN	= 

#######################################################
#           Compilation for the different files  	  #
#######################################################


all: $(BIN)

## Compilation du main
$(BIN): 2_src/Finding_context_of_motif.cpp $(OBJS_GEN)
	g++ 2_src/Finding_context_of_motif.cpp $(OBJS_GEN) -o $@


## Compilation des fichiers de fonctions
obj/%.o: 2_src/%.cpp
	g++ -c $< -o $@

## Cleaning the folders with compiled files
clean:
	rm -r bin/* obj/*

submit:
	rm $(SUBMIT_NAME)
	tar -czvf $(SUBMIT_NAME) *