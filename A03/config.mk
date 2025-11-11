# Supported: GCC, CLANG, ICC
TAG ?= ICX
ENABLE_OPENMP ?= false

# used to compile the NON-blocking version of code 
COMM_NB ?= true

ifeq ($(COMM_NB),true)
	OPTIONS += -DNB_COMMUICATION
endif

# uncomment while benchmarking to prevent MPI_CALL check
OPTIONS +=  -DNDEBUG

#Feature options
OPTIONS +=  -DARRAY_ALIGNMENT=64
#OPTIONS +=  -DVERBOSE_AFFINITY
#OPTIONS +=  -DVERBOSE_DATASIZE
#OPTIONS +=  -DVERBOSE_TIMER
