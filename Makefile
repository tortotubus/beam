# —————————————————————————————————————————————
# Variables
# —————————————————————————————————————————————
CXX      := g++
AR       := ar
RM       := rm -rf

# where we put headers and libs
INCDIR   := include
LIBDIR   := lib
OBJDIR   := obj

# compilation flags
CXXFLAGS := -std=c++11 -fPIC -I$(INCDIR) -Wall -Wextra

# linker flags for external numeric libs
LIBS     := -llapacke -llapack -lcblas

# library name
LIB      := beam
STATIC   := $(LIBDIR)/lib$(LIB).a
SHARED   := $(LIBDIR)/lib$(LIB).so

# source / object lists
SRCS     := $(wildcard src/*.cpp)
OBJS     := $(patsubst src/%.cpp,$(OBJDIR)/%.o,$(SRCS))

# —————————————————————————————————————————————
# Top-level targets
# —————————————————————————————————————————————
.PHONY: all clean install

all: $(STATIC) $(SHARED)

clean:
	rm -rf lib obj

# —————————————————————————————————————————————
# Build object files
# —————————————————————————————————————————————
$(OBJDIR)/%.o: src/%.cpp | $(OBJDIR)
	@echo "⏳  Compiling $<"
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

# —————————————————————————————————————————————
# Build static library
# —————————————————————————————————————————————
$(STATIC): $(OBJS) | $(LIBDIR)
	@echo "📦  Archiving $@"
	$(AR) rcs $@ $(OBJS)

# —————————————————————————————————————————————
# Build shared library
# —————————————————————————————————————————————
$(SHARED): $(OBJS) | $(LIBDIR)
	@echo "🔗  Linking $@"
	$(CXX) -shared -o $@ $(OBJS) $(LIBS)

$(LIBDIR):
	mkdir -p $(LIBDIR)

# —