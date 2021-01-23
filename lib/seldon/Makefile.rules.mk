# Makefile rules
#
ifeq ($(USE_BLAS),YES)
  BLAS_MESSAGE="Utilisation de Blas"
else
  BLAS_MESSAGE="Blas non utilise"
endif

ifeq ($(USE_ARPACK),YES)
  ARPACK_MESSAGE="Utilisation de Arpack"
else
  ARPACK_MESSAGE="Arpack non utilise"
endif

ifeq ($(USE_MUMPS),YES)
  MUMPS_MESSAGE="Utilisation de Mumps"
else
  MUMPS_MESSAGE="Mumps non utilise"
endif

ifeq ($(USE_PARDISO),YES)
  PARDISO_MESSAGE="Utilisation de Pardiso"
else
  PARDISO_MESSAGE="Pardiso non utilise"
endif

ifeq ($(USE_PASTIX),YES)
  PASTIX_MESSAGE="Utilisation de Pastix"
else
  PASTIX_MESSAGE="Pastix non utilise"
endif

ifeq ($(USE_CHOLMOD),YES)
  CHOLMOD_MESSAGE="Utilisation de Cholmod"
else
  CHOLMOD_MESSAGE="Cholmod non utilise"
endif

ifeq ($(OPTIMIZATION),YES)
  OPTIM_MESSAGE="Compilation en mode optimise"
else
  OPTIM_MESSAGE="Compilation en mode debogage"
endif

ifeq ($(USE_MPI),YES)
  MPI_MESSAGE="Compilation en mode parallele"
else
  MPI_MESSAGE="Compilation en mode sequentiel"
endif

ifeq ($(USE_VERBOSE),YES)
  VERBOSE=
else
  VERBOSE=@
endif

# Default executable
regular: lib $(BIN)

# Semi-external libraries
lib: messages lib/etime.o

messages :
	@echo $(BLAS_MESSAGE)
	@echo $(ARPACK_MESSAGE)
	@echo $(MUMPS_MESSAGE)
	@echo $(PASTIX_MESSAGE)
	@echo $(PARDISO_MESSAGE)
	@echo $(CHOLMOD_MESSAGE)
	@echo
	@echo $(MPI_MESSAGE)
	@echo $(OPTIM_MESSAGE)
	@echo

lib/etime.o : lib/etime.f
	$(VERBOSE)$(FCC) -c $< -o $@ -O3 -fPIC

include Makefile.lib.mk

# Default binary
$(BIN) : lib $(LSTBIN) $(LIB_SELDON)
	@echo "Compilation de la regle $@"
	@echo "a partir de $(LSTBIN)"
	$(VERBOSE)$(CC) $(FLAGS) $(LSTBIN) -c -o $(BIN).o $(INCLUDE) \
	  $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) $(BIN).o -o $(BIN) $(AGRESSIVE_OPTIM) $(LIB)
	$(VERBOSE)mv $(BIN) .
	@echo "Regle $@ compilee avec succes."

# deleting object files
clean :
	rm lib/*.o lib/*.a lib/Compil/Seldon/*.o

.PHONY: clean
