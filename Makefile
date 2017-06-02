
# PROJECT crystPEstWithAdolcIpopt
# Created on:    Thu, Jun  1, 2017  2:30:44 PM
#TODO: ....
#
##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change LIB/INC/CFLAGS entries below.   #
##########################################################################

CC := g++ # This is the main compiler
SRCDIR := src
BUILDDIR := build
TARGET := bin/ipoptFiniteDif
 
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -pipe -DNDEBUG -Wparentheses -Wreturn-type \
			-Wconversion -Wno-unknown-pragmas -Wno-long-long -DIPOPT_BUILD \
			-Wcast-qual -Wall -Wpointer-arith -Wwrite-strings \
			-g
			#-O3
LIB := `PKG_CONFIG_PATH=/mingw64/lib64/pkgconfig:/mingw64/lib/pkgconfig:/mingw64/share/pkgconfig:/mingw64/lib/pkgconfig:/mingw64/share/pkgconfig pkg-config --libs ipopt` \
		-L /c/Users/Marcellos/lib/ADOL-C-2.6.3/ADOL-C/.libs -ladolc \
		-L /c/Users/Marcellos/OneDrive/Documents/Projects/DSC_Docs/cpp_files_shared/ \
		-lmyLib \
		-lgsl -lgslcblas \
		-lhdf5 -lhdf5_cpp \
		-larmadillo
INC := `PKG_CONFIG_PATH=/mingw64/lib64/pkgconfig:/mingw64/lib/pkgconfig:/mingw64/share/pkgconfig:/mingw64/lib/pkgconfig:/mingw64/share/pkgconfig pkg-config --cflags ipopt` \
		-I include \
		-I /c/Users/Marcellos/OneDrive/Documents/Projects/DSC_Docs/cpp_files_shared/


##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile Ipopt.   #
##########################################################################

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

# Tests
tester:
	$(CC) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester

# Spikes
ticket:
	$(CC) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean

