EXECSOURCES := src/tw-heuristic.cpp
SOURCES := src/minimum_degree_heuristic.cpp

###################################################

SRCDIR := src
OBJDIR := obj
BINDIR := bin

OBJECTS := $(addprefix $(OBJDIR)/,$(notdir $(SOURCES:.cpp=.o)))
EXECOBJECTS := $(addprefix $(OBJDIR)/,$(notdir $(EXECSOURCES:.cpp=.o)))
BINARIES := $(addprefix $(BINDIR)/,$(notdir $(EXECSOURCES:.cpp=)))

CXX = $(shell which clang++ || which g++)
LDFLAGS ?=

# CXX options
CXXFLAGS  = -std=c++11 -Isrc -fno-exceptions -MMD -MP

# Warning/error flags
CXXFLAGS += -Werror -Wall -Wpedantic -Wcast-align -Wcast-qual -Wold-style-cast
CXXFLAGS += -Wno-c++98-compat -Wno-global-constructors -Wno-padded 
CXXFLAGS += -Wredundant-decls -Wshadow -Wmissing-include-dirs
CXXFLAGS += -Wno-disabled-macro-expansion

DEBUG ?= 0
ifeq ($(DEBUG), 1)
  CXXFLAGS += -DDEBUG -O0 -fno-omit-frame-pointer -g
else
  CXXFLAGS += -DNDEBUG -O3 -fomit-frame-pointer -march=native
  CXXFLAGS += -Wdisabled-optimization
endif

VALIDATE_TD ?= 0
ifeq ($(VALIDATE_TD), 1)
  CXXFLAGS += -DVALIDATE_TD
endif

CMD = $(CXX) $(CXXFLAGS)
 
all: $(OBJECTS) $(EXECOBJECTS) $(BINARIES) symlink
.PHONY: all force

CLEANFLAGS = $(CMD) $(LDFLAGS)
cleanfile: force
	@echo '$(CLEANFLAGS)' | cmp -s - $@ || echo '$(CLEANFLAGS)' > cleanfile

$(BINARIES) : $(BINDIR)/% : $(OBJDIR)/%.o $(OBJECTS)
	$(CMD) -o $@ $^ $(LDFLAGS) 

$(OBJECTS) : | $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR) $(BINDIR)

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp cleanfile
	$(CMD) -o $@ -c $<

symlink: bin/tw-heuristic
	@rm -f tw-heuristic
	@ln -s bin/tw-heuristic tw-heuristic

clean:
	rm -f $(BINARIES) *.d $(OBJDIR)/*.d *.o $(OBJDIR)/*.o

test:
	bin/test_all

kill:
	@kill -9 `pidof tw-heuristic`

-include $(wildcard $(OBJDIR)/*.d)

