#!/usr/bin/make -f

PREFIX ?= /usr
LIBDIR ?= lib
LV2DIR ?= $(PREFIX)/$(LIBDIR)/lv2

OPTIMIZATIONS ?= -msse -msse2 -mfpmath=sse -ffast-math -fomit-frame-pointer -O3 -fno-finite-math-only

#OPTIMIZATIONS ?=

LDFLAGS ?= -Wl,--as-needed
CXXFLAGS ?= $(OPTIMIZATIONS) -Wall
CFLAGS ?= $(OPTIMIZATIONS) -Wall

###############################################################################
BUNDLE = zameq2.lv2

CXXFLAGS += -fPIC -DPIC
CFLAGS += -fPIC -DPIC

UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
  LIB_EXT=.dylib
  LDFLAGS += -dynamiclib
else
  LDFLAGS += -shared -Wl,-Bstatic -Wl,-Bdynamic
  LIB_EXT=.so
endif


ifeq ($(shell pkg-config --exists lv2 lv2core lv2-plugin || echo no), no)
  $(error "LV2 SDK was not found")
else
  LV2FLAGS=`pkg-config --cflags --libs lv2 lv2core lv2-plugin`
endif

ifeq ($(shell pkg-config --exists lv2-gui || echo no), no)
  $(error "LV2-GUI is required ")
else
  LV2GUIFLAGS=`pkg-config --cflags --libs lv2-gui lv2 lv2core lv2-plugin`
endif


$(BUNDLE): manifest.ttl zameq2.ttl zameq2$(LIB_EXT)
#zameq2_gui$(LIB_EXT)
	rm -rf $(BUNDLE)
	mkdir $(BUNDLE)
	cp manifest.ttl zameq2.ttl zameq2$(LIB_EXT) $(BUNDLE)

zameq2$(LIB_EXT): zameq2.c
	$(CXX) -o zameq2$(LIB_EXT) \
		$(CXXFLAGS) \
		zameq2.c \
		$(LV2FLAGS) $(LDFLAGS)

#zamcomp_gui$(LIB_EXT): zamcomp_gui.cpp zamcomp.peg
#	$(CXX) -o zamcomp_gui$(LIB_EXT) \
#		$(CXXFLAGS) \
#		zamcomp_gui.cpp \
#		$(LV2GUIFLAGS) $(LDFLAGS)

zameq2.peg: zameq2.ttl
	lv2peg zameq2.ttl zameq2.peg

install: $(BUNDLE)
	install -d $(DESTDIR)$(LV2DIR)/$(BUNDLE)
	install -t $(DESTDIR)$(LV2DIR)/$(BUNDLE) $(BUNDLE)/*
#	install zameq2_gui$(LIB_EXT) $(DESTDIR)$(LV2DIR)/$(BUNDLE)

uninstall:
	rm -rf $(DESTDIR)$(LV2DIR)/$(BUNDLE)

clean:
	rm -rf $(BUNDLE) zameq2$(LIB_EXT) zameq2_gui$(LIB_EXT) zameq2.peg

.PHONY: clean install uninstall
