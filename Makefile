.PHONY = all clean

RM := rm -vf
MAKE := make
TARGET := fit_potential.x
PETSCDIR := $(HOME)/petsc

OUTPUT_DIR := output
BIN_DIR := bin
SRC_DIR := src

BUILD_DIRS := $(shell find $(SRC_DIR) -mindepth 1 -maxdepth 2 -type d 2>/dev/null | sed 's/$(SRC_DIR)/build/g')

# Detect OS
UNAME := $(shell uname)
SED := sed

ifeq ($(UNAME), Darwin)
    SED := gsed  # Use GNU sed on macOS if installed via Homebrew
endif

print_os:
	@echo "Operating System: $(UNAME). Using $(SED) for sed commands."

all: create_folders
	@$(MAKE) -j -C build all

create_folders: 
	@mkdir -p $(BUILD_DIRS) $(BIN_DIR) $(OUTPUT_DIR)

config:
	@echo "Configuring for PETSc..."
	@echo "Using PETSc directory: $(PETSCDIR), if this is not correct, please set it in the Makefile (for info type 'make help')"
	$(SED) -i 's|-include .*|-include $(PETSCDIR)/petscdir.mk|' build/Makefile

run: all
	./bin/$(TARGET) -tao_monitor_short -tao_max_it 10000 -tao_type pounders -tao_gatol 1.e-8

compare:
	@echo  "@with g0 \n\
	@		  title \"Compare fitted solution to AV18 solution for 1P1\" \n\
	@		  xaxis  label \"k\\S2\\N (fm\\S-2\\N)\" \n\
	@		  yaxis  label \"k\\S3\\N cot delta (fm\\S-3\\N)\" \n\
	@     s0 legend \"AV18\" \n\
	@     s1 legend \"fitted EFT-pless\" \n\
	@     legend on \n\
	@		  legend 0.8, 0.8 \n\
	@     s0 line linewidth 2 \n\
	@     s1 line linewidth 2 \n\
	@     s1 line linestyle 3 \n\
	"> config.agr
	@sed -i 's/delta/\\xd\\f{}/g' config.agr
	xmgrace $(OUTPUT_DIR)/*.dat config.agr -saveall $(OUTPUT_DIR)/compare.agr
	@$(RM) -v config.agr


clean:
	@$(MAKE) -C build clean

del_out:
	@$(RM) -v output/*.dat

help: 
	@echo "Makefile for fit_potential"
	@echo "Usage:"
	@echo "  make all        - Build the project"
	@echo "  make config     - Configure the project for PETSc, pass the path to PETSc using \"make config PETSCDIR=/path/to/petsc\""
	@echo "  make run        - Run the project"
	@echo "  make compare    - Compare fitted solution to AV18 solution"
	@echo "  make clean      - Clean the build directory"
	@echo "  make del_out    - Delete output files"
	@echo "  make help       - Show this help message"