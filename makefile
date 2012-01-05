# Makefile for cswt library


# ======== OPTIONS ========

USEPGPLOT = no


# ======== COMPILER ========

FC      = gfortran
#FC      = f95
#FC      = g95

ifneq ($(USEPGPLOT),yes)
  OPTPGPLOT     = -DNO_PGPLOT
endif
OPT = $(OPTPGPLOT) -m64 -O3 \
      -DCSWT_VERSION=\"1.0b2\" -DCSWT_BUILD=\"`svnversion -n .`\" 


# ======== LINKS ========

PROGDIR    = ..

HPIXDIR = $(PROGDIR)/Healpix
HPIXLIB = $(HPIXDIR)/lib
HPIXLIBNM= healpix
HPIXINC = $(HPIXDIR)/include

S2DIR  = $(PROGDIR)/s2
S2LIB  = $(S2DIR)/lib
S2LIBNM= s2
S2INC  = $(S2DIR)/include
S2SRC  = $(S2DIR)/src/mod
S2PROG = $(S2DIR)/src/prog
S2BIN  = $(S2DIR)/bin
S2DOC  = $(S2DIR)/doc

CSWTDIR = $(PROGDIR)/fastcswt
CSWTSRC = $(CSWTDIR)/src/mod
CSWTPROG= $(CSWTDIR)/src/prog
CSWTINC = $(CSWTDIR)/include
CSWTBIN = $(CSWTDIR)/bin
CSWTLIB = $(CSWTDIR)/lib
CSWTDOC = $(CSWTDIR)/doc
CSWTLIBNM= fastcswt

CFITSIOLIB   = $(PROGDIR)/cfitsio/lib
CFITSIOLIBNM = cfitsio

PGPLOTLIB    = $(PROGDIR)/pgplot
PGPLOTLIBNM  = pgplot
X11LIB       = /usr/X11R6/lib
X11LIBNM     = X11

FFTWLIB      = $(PROGDIR)/fftw-2.1.5/lib
FFTWLIBNM    = fftw
RFFTWLIBNM   = rfftw


# ======== FFFLAGS ========

FFLAGS  = -I$(HPIXINC) -I$(S2INC) -I$(CSWTINC) -I.


# ======== LDFLAGS ========

ifeq ($(USEPGPLOT),yes)
  LDFLAGSPGPLOT = -L$(PGPLOTLIB) -L$(X11LIB) \
                  -l$(PGPLOTLIBNM) -l$(X11LIBNM)
endif

LDFLAGS = -L$(CSWTLIB) -l$(CSWTLIBNM) \
          -L$(S2LIB) -l$(S2LIBNM) \
          -L$(HPIXLIB) -l$(HPIXLIBNM) \
          -L$(CFITSIOLIB) -l$(CFITSIOLIBNM) \
          -L$(FFTWLIB) -l$(RFFTWLIBNM) -l$(FFTWLIBNM) \
          $(LDFLAGSPGPLOT)


# ======== PPFLAGS ========

ifeq ($(FC),f95)
  PPFLAGS = -fpp $(OPT)
else ifeq ($(FC),g95)
  PPFLAGS = -cpp $(OPT)
else ifeq ($(FC),gfortran)
  PPFLAGS = -cpp $(OPT)
endif


# ======== OBJECT FILES TO MAKE ========

CSWTOBJ = $(CSWTINC)/cswt_error_mod.o \
          $(CSWTINC)/cswt_tmpl_mod.o  \
          $(CSWTINC)/cswt_swav_mod.o  \
          $(CSWTINC)/cswt_tr_mod.o


# ======== MAKE RULES ========

default: all

all:     lib prog

lib:     $(CSWTLIB)/lib$(CSWTLIBNM).a

prog:    $(CSWTBIN)/cswt_analysis            \
         $(CSWTBIN)/cswt_tr2sky              \
         $(CSWTBIN)/cswt_tr_nsigma           \
         $(CSWTBIN)/cswt_mask_gen            \
         $(CSWTBIN)/cswt_mask_copy           \
         $(CSWTBIN)/cswt_mask_apply          \
         $(CSWTBIN)/cswt_mask_invert         \
         $(CSWTBIN)/cswt_mask_nonzero        \
         $(CSWTBIN)/cswt_mask_nonzero_weight \
         $(CSWTBIN)/cswt_plot_swav           \
         $(CSWTBIN)/cswt_swav_azbandlim      \
         $(CSWTBIN)/cswt_about


$(CSWTINC)/%.o: $(CSWTSRC)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 
	mv *.mod $(CSWTINC)

$(CSWTINC)/%.o: $(CSWTPROG)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 


# Library

$(CSWTLIB)/lib$(CSWTLIBNM).a: $(CSWTOBJ)
	ar -r $(CSWTLIB)/lib$(CSWTLIBNM).a $(CSWTOBJ)


# Documentation

docs:
	f90doc_fpp $(CSWTSRC)/*.f90
	f90doc_fpp $(CSWTPROG)/*.f90
	ln_multi $(S2DOC)/s2_*
	ln_multi $(S2DOC)/index_s2.html
	mv *.html $(CSWTDOC)/.
	addstyle $(CSWTDOC)/cswt_*

cleandocs:
	rm -f $(CSWTDOC)/cswt_*.html
	rm -f $(CSWTDOC)/s2_*.html
	rm -f $(CSWTDOC)/index_s2.html


# Cleaning up

clean:	tidy
	rm -f $(CSWTINC)/*.mod
	rm -f $(CSWTINC)/*.o
	rm -f $(CSWTLIB)/lib$(CSWTLIBNM).a
	rm -f $(CSWTBIN)/*

tidy:	
	rm -f *.mod
	rm -f $(CSWTSRC)/*~
	rm -f $(CSWTPROG)/*~


# Module dependencies

$(CSWTINC)/cswt_error_mod.o: $(CSWTSRC)/cswt_error_mod.f90
$(CSWTINC)/cswt_tmpl_mod.o:  $(CSWTSRC)/cswt_tmpl_mod.f90  \
                               $(CSWTINC)/cswt_error_mod.o
$(CSWTINC)/cswt_swav_mod.o:  $(CSWTSRC)/cswt_swav_mod.f90  \
                               $(CSWTINC)/cswt_tmpl_mod.o  \
                               $(CSWTINC)/cswt_error_mod.o
$(CSWTINC)/cswt_tr_mod.o:    $(CSWTSRC)/cswt_tr_mod.f90    \
                               $(CSWTINC)/cswt_tmpl_mod.o  \
                               $(CSWTINC)/cswt_error_mod.o \
                               $(CSWTINC)/cswt_swav_mod.o


# Program dependencies and compilation

$(CSWTINC)/cswt_analysis.o:            $(CSWTPROG)/cswt_analysis.f90 lib
$(CSWTBIN)/cswt_analysis:              $(CSWTINC)/cswt_analysis.o
	$(FC) -o $(CSWTBIN)/cswt_analysis $(CSWTINC)/cswt_analysis.o \
	$(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_tr2sky.o:              $(CSWTPROG)/cswt_tr2sky.f90 lib
$(CSWTBIN)/cswt_tr2sky:                $(CSWTINC)/cswt_tr2sky.o
	$(FC) -o $(CSWTBIN)/cswt_tr2sky $(CSWTINC)/cswt_tr2sky.o \
	$(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_tr_nsigma.o:           $(CSWTPROG)/cswt_tr_nsigma.f90 lib
$(CSWTBIN)/cswt_tr_nsigma:             $(CSWTINC)/cswt_tr_nsigma.o
	$(FC) -o $(CSWTBIN)/cswt_tr_nsigma $(CSWTINC)/cswt_tr_nsigma.o \
	$(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_mask_gen.o:            $(CSWTPROG)/cswt_mask_gen.f90 lib
$(CSWTBIN)/cswt_mask_gen:              $(CSWTINC)/cswt_mask_gen.o
	$(FC) -o $(CSWTBIN)/cswt_mask_gen $(CSWTINC)/cswt_mask_gen.o \
	$(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_mask_copy.o:           $(CSWTPROG)/cswt_mask_copy.f90 lib
$(CSWTBIN)/cswt_mask_copy:             $(CSWTINC)/cswt_mask_copy.o
	$(FC) -o $(CSWTBIN)/cswt_mask_copy $(CSWTINC)/cswt_mask_copy.o \
	$(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_mask_apply.o:          $(CSWTPROG)/cswt_mask_apply.f90 lib
$(CSWTBIN)/cswt_mask_apply:            $(CSWTINC)/cswt_mask_apply.o
	$(FC) -o $(CSWTBIN)/cswt_mask_apply $(CSWTINC)/cswt_mask_apply.o \
	$(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_mask_invert.o:         $(CSWTPROG)/cswt_mask_invert.f90 lib
$(CSWTBIN)/cswt_mask_invert:           $(CSWTINC)/cswt_mask_invert.o
	$(FC) -o $(CSWTBIN)/cswt_mask_invert $(CSWTINC)/cswt_mask_invert.o \
	$(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_mask_nonzero.o:        $(CSWTPROG)/cswt_mask_nonzero.f90 lib
$(CSWTBIN)/cswt_mask_nonzero:          $(CSWTINC)/cswt_mask_nonzero.o
	$(FC) -o $(CSWTBIN)/cswt_mask_nonzero $(CSWTINC)/cswt_mask_nonzero.o \
	$(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_mask_nonzero_weight.o: $(CSWTPROG)/cswt_mask_nonzero_weight.f90 lib
$(CSWTBIN)/cswt_mask_nonzero_weight:   $(CSWTINC)/cswt_mask_nonzero_weight.o
	$(FC) -o $(CSWTBIN)/cswt_mask_nonzero_weight \
	$(CSWTINC)/cswt_mask_nonzero_weight.o $(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_plot_swav.o:           $(CSWTPROG)/cswt_plot_swav.f90 lib
$(CSWTBIN)/cswt_plot_swav:	       $(CSWTINC)/cswt_plot_swav.o
	$(FC) -o $(CSWTBIN)/cswt_plot_swav \
	$(CSWTINC)/cswt_plot_swav.o $(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_swav_azbandlim.o:      $(CSWTPROG)/cswt_swav_azbandlim.f90 lib
$(CSWTBIN)/cswt_swav_azbandlim:	       $(CSWTINC)/cswt_swav_azbandlim.o
	$(FC) -o $(CSWTBIN)/cswt_swav_azbandlim \
	$(CSWTINC)/cswt_swav_azbandlim.o $(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_test_swav.o:           $(CSWTPROG)/cswt_test_swav.f90 lib
$(CSWTBIN)/cswt_test_swav:	       $(CSWTINC)/cswt_test_swav.o
	$(FC) -o $(CSWTBIN)/cswt_test_swav \
	$(CSWTINC)/cswt_test_swav.o $(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_test_tr.o:             $(CSWTPROG)/cswt_test_tr.f90 lib
$(CSWTBIN)/cswt_test_tr:               $(CSWTINC)/cswt_test_tr.o
	$(FC) -o $(CSWTBIN)/cswt_test_tr \
	$(CSWTINC)/cswt_test_tr.o $(LDFLAGS) $(PPFLAGS)

$(CSWTINC)/cswt_about.o:             $(CSWTPROG)/cswt_about.f90 lib
$(CSWTBIN)/cswt_about:               $(CSWTINC)/cswt_about.o
	$(FC) -o $(CSWTBIN)/cswt_about \
	$(CSWTINC)/cswt_about.o $(LDFLAGS) $(PPFLAGS)



