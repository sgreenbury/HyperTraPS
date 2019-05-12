SOURCEDIR=src/cpp/
BINDIR=bin/
OBJECTDIR=obj/

CXX = g++
CXXFLAGS = -O3 -g -Wall -std=c++11 -pthread

LIBS =

MAIN_MCMC = run_mcmc_sampler.cpp
EXECUTABLE_MCMC = RUN_MCMC_SAMPLER
SOURCES_MCMC = $(MAIN_MCMC) utilities.cpp transition.cpp state.cpp trajectory.cpp updater_fast.cpp likelihood_sampler_thread.cpp mcmc_sampler.cpp
OBJECTS_MCMC = $(SOURCES_MCMC:.cpp=.o)
CSOURCES_MCMC=$(addprefix $(SOURCEDIR),$(SOURCES_MCMC))
COBJECTS_MCMC=$(addprefix $(OBJECTDIR),$(OBJECTS_MCMC))
CEXECUTABLE_MCMC=$(addprefix $(BINDIR),$(EXECUTABLE_MCMC))


MAIN_IC = run_information_criterion.cpp
EXECUTABLE_IC = RUN_INFORMATION_CRITERION
SOURCES_IC = $(MAIN_IC) utilities.cpp transition.cpp state.cpp trajectory.cpp updater_fast.cpp likelihood_sampler_thread.cpp
OBJECTS_IC = $(SOURCES_IC:.cpp=.o)
CSOURCES_IC=$(addprefix $(SOURCEDIR),$(SOURCES_IC))
COBJECTS_IC=$(addprefix $(OBJECTDIR),$(OBJECTS_IC))
CEXECUTABLE_IC=$(addprefix $(BINDIR),$(EXECUTABLE_IC))

MAIN_PW = run_pw.cpp
EXECUTABLE_PW = RUN_PW
SOURCES_PW = $(MAIN_PW) utilities.cpp transition.cpp state.cpp trajectory.cpp updater_fast.cpp graph.cpp
OBJECTS_PW = $(SOURCES_PW:.cpp=.o)
CSOURCES_PW=$(addprefix $(SOURCEDIR),$(SOURCES_PW))
COBJECTS_PW=$(addprefix $(OBJECTDIR),$(OBJECTS_PW))
CEXECUTABLE_PW=$(addprefix $(BINDIR),$(EXECUTABLE_PW))

all: mcmc pw ic

mcmc: $(CSOURCES_MCMC) $(CEXECUTABLE_MCMC)

$(CEXECUTABLE_MCMC): $(COBJECTS_MCMC)
	$(CXX) $(CXXFLAGS) $(LIBS) $(COBJECTS_MCMC) -o $@

pw: $(CSOURCES_PW) $(CEXECUTABLE_PW)

$(CEXECUTABLE_PW): $(COBJECTS_PW)
	$(CXX) $(CXXFLAGS) $(LIBS) $(COBJECTS_PW) -o $@

ic: $(CSOURCES_IC) $(CEXECUTABLE_IC)

$(CEXECUTABLE_IC): $(COBJECTS_IC)
	$(CXX) $(CXXFLAGS) $(LIBS) $(COBJECTS_IC) -o $@

$(OBJECTDIR)%.o: $(SOURCEDIR)%.cpp
	$(CXX) $(CXXFLAGS) $(LIBS) $< -c -o $@

clean:
	rm -f $(CEXECUTABLE_MCMC)
	rm -f $(CEXECUTABLE_IC)
	rm -f $(CEXECUTABLE_PW)
	rm -f $(addprefix $(OBJECTDIR),*.o)
