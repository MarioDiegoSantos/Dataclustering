#---------------------------------------------------------------------
#
#   Compation Settings and Paths
#
#---------------------------------------------------

CXX =  /usr/local/cuda-10.1/bin/nvcc
CXXFLAGS = -O3 -DNDEBUG -arch=sm_35 -rdc=true -ccbin g++-7 --compiler-options "-Wall -fopenmp" -c
CXXDFLAGS = -O0 -g -arch=sm_35 -rdc=true -ccbin g++-7 --compiler-options "-Wall" -c 
INCLUDE = -I/usr/local/cuda-10.1/lib64   -Icub-1.8.0
LFFLAGS = -O3 -ccbin g++-7 -lgomp -arch=sm_35
LFFDFLAGS = -O0 -ccbin g++-7 -lgomp -arch=sm_35

# CXX = gcc
# CXXFLAGS = -O3 -DNDEBUG -Wall -fopenmp -c
# CXXDFLAGS = -O0 -g -Wall -c -fopenmp
# INCLUDE = -Icub-1.8.0
# LFFLAGS = -O3 -lgomp -lpython2.7 -lm
# LFFDFLAGS = -O0 -lgomp -lpython2.7 -lm


#---------------------------------------------------------------------
#
#   SRC REPS
#
#----------------------------------------------------------------------

SRC_REP = src

#---------------------------------------------------------------------
#
#   SRC FILES
#
#---------------------------------------------------------------------

SRC = $(wildcard $(SRC_REP)/*.cpp) \
		$(wildcard $(SRC_REP)/*.cu)	

#---------------------------------------------------------------------
#
#   OBJ REPS
#
#---------------------------------------------------------------------

OBJ_REP = obj

#---------------------------------------------------------------------
#
#   OBJ files
#
#---------------------------------------------------------------------

OBJ_TEMP = $(addprefix $(OBJ_REP)/,$(notdir $(patsubst %.cu, %.o, $(SRC))))
OBJ = $(addprefix $(OBJ_REP)/,$(notdir $(patsubst %.cpp, %.o, $(OBJ_TEMP))))

print-%  : ; @echo $* = $($*)

#---------------------------------------------------------------------
#
#   Dependencies	 
#
#---------------------------------------------------------------------

DEP := $(OBJ:.o=.d)

#---------------------------------------------------------------------
#
#   Targets
#
#---------------------------------------------------------------------

cgrasp: $(OBJ) 
	$(CXX) $^ $(LFFLAGS) -o cgrasp

debug: CXXFLAGS=$(CXXDFLAGS)
debug: LFFLAGS=$(LFFDFLAGS)
debug: cgrasp

clean:
	rm -rf obj
	rm -f core	
	rm -f cgrasp

#---------------------------------------------------------------------
#
#   General Rules
#
#---------------------------------------------------------------------

$(OBJ_REP)/%.o: $(SRC_REP)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(INCLUDE) $(CXXFLAGS) $< -o $@

$(OBJ_REP)/%.o: $(SRC_REP)/%.cu
	@mkdir -p $(@D)
	$(CXX) $(INCLUDE) $(CXXFLAGS) $< -o $@
	

-include $(DEP)
