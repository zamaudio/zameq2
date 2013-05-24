//-----------------------------------------------------
// name: "ZamEQ2"
// author: "Damien Zammit"
// copyright: "2013"
// version: "1.00"
// license: "GPLv2"
//
// Code generated with Faust 0.9.61 (http://faust.grame.fr)
//-----------------------------------------------------
/* link with  */
#include <math.h>
#include <cmath>
template <int N> inline float faustpower(float x) 		{ return powf(x,N); } 
template <int N> inline double faustpower(double x) 	{ return pow(x,N); }
template <int N> inline int faustpower(int x) 			{ return faustpower<N/2>(x) * faustpower<N-N/2>(x); } 
template <> 	 inline int faustpower<0>(int x) 		{ return 1; }
template <> 	 inline int faustpower<1>(int x) 		{ return x; }
/************************************************************************
 ************************************************************************
    FAUST Architecture File
    Copyright (C) 2009-2011 Albert Graef <Dr.Graef@t-online.de>
    ---------------------------------------------------------------------
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2.1 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with the GNU C Library; if not, write to the Free
    Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
    02111-1307 USA.
 ************************************************************************
 ************************************************************************/

/* LV2 architecture for Faust. */

#include <stdlib.h>
#include <math.h>
#include <list>
#include <map>

using namespace std;

// On Intel set FZ (Flush to Zero) and DAZ (Denormals Are Zero)
// flags to avoid costly denormals
#ifdef __SSE__
    #include <xmmintrin.h>
    #ifdef __SSE2__
        #define AVOIDDENORMALS _mm_setcsr(_mm_getcsr() | 0x8040)
    #else
        #define AVOIDDENORMALS _mm_setcsr(_mm_getcsr() | 0x8000)
    #endif
#else
  #define AVOIDDENORMALS
#endif

typedef pair<const char*,const char*> strpair;

struct Meta
{
  list< strpair > data;
  void declare (const char* key, const char* value)
  { data.push_back(strpair(key, value)); }
};

//-------------------------------------------------------------------
// Generic min and max using c++ inline
//-------------------------------------------------------------------

inline int 	max (unsigned int a, unsigned int b) { return (a>b) ? a : b; }
inline int 	max (int a, int b)		{ return (a>b) ? a : b; }

inline long 	max (long a, long b) 		{ return (a>b) ? a : b; }
inline long 	max (int a, long b) 		{ return (a>b) ? a : b; }
inline long 	max (long a, int b) 		{ return (a>b) ? a : b; }

inline float 	max (float a, float b) 		{ return (a>b) ? a : b; }
inline float 	max (int a, float b) 		{ return (a>b) ? a : b; }
inline float 	max (float a, int b) 		{ return (a>b) ? a : b; }
inline float 	max (long a, float b) 		{ return (a>b) ? a : b; }
inline float 	max (float a, long b) 		{ return (a>b) ? a : b; }

inline double 	max (double a, double b) 	{ return (a>b) ? a : b; }
inline double 	max (int a, double b) 		{ return (a>b) ? a : b; }
inline double 	max (double a, int b) 		{ return (a>b) ? a : b; }
inline double 	max (long a, double b) 		{ return (a>b) ? a : b; }
inline double 	max (double a, long b) 		{ return (a>b) ? a : b; }
inline double 	max (float a, double b) 	{ return (a>b) ? a : b; }
inline double 	max (double a, float b) 	{ return (a>b) ? a : b; }


inline int	min (int a, int b)		{ return (a<b) ? a : b; }

inline long 	min (long a, long b) 		{ return (a<b) ? a : b; }
inline long 	min (int a, long b) 		{ return (a<b) ? a : b; }
inline long 	min (long a, int b) 		{ return (a<b) ? a : b; }

inline float 	min (float a, float b) 		{ return (a<b) ? a : b; }
inline float 	min (int a, float b) 		{ return (a<b) ? a : b; }
inline float 	min (float a, int b) 		{ return (a<b) ? a : b; }
inline float 	min (long a, float b) 		{ return (a<b) ? a : b; }
inline float 	min (float a, long b) 		{ return (a<b) ? a : b; }

inline double 	min (double a, double b) 	{ return (a<b) ? a : b; }
inline double 	min (int a, double b) 		{ return (a<b) ? a : b; }
inline double 	min (double a, int b) 		{ return (a<b) ? a : b; }
inline double 	min (long a, double b) 		{ return (a<b) ? a : b; }
inline double 	min (double a, long b) 		{ return (a<b) ? a : b; }
inline double 	min (float a, double b) 	{ return (a<b) ? a : b; }
inline double 	min (double a, float b) 	{ return (a<b) ? a : b; }

// abs is now predefined
//template<typename T> T abs (T a)		{ return (a<T(0)) ? -a : a; }

inline int	lsr (int x, int n)		{ return int(((unsigned int)x) >> n); }

/******************************************************************************
*******************************************************************************

							       VECTOR INTRINSICS

*******************************************************************************
*******************************************************************************/

//inline void *aligned_calloc(size_t nmemb, size_t size) { return (void*)((unsigned)(calloc((nmemb*size)+15,sizeof(char)))+15 & 0xfffffff0); }
//inline void *aligned_calloc(size_t nmemb, size_t size) { return (void*)((size_t)(calloc((nmemb*size)+15,sizeof(char)))+15 & ~15); }


/******************************************************************************
*******************************************************************************

			ABSTRACT USER INTERFACE

*******************************************************************************
*******************************************************************************/

class UI
{
  bool	fStopped;
public:

  UI() : fStopped(false) {}
  virtual ~UI() {}

  virtual void addButton(const char* label, float* zone) = 0;
  virtual void addCheckButton(const char* label, float* zone) = 0;
  virtual void addVerticalSlider(const char* label, float* zone, float init, float min, float max, float step) = 0;
  virtual void addHorizontalSlider(const char* label, float* zone, float init, float min, float max, float step) = 0;
  virtual void addNumEntry(const char* label, float* zone, float init, float min, float max, float step) = 0;

  virtual void addHorizontalBargraph(const char* label, float* zone, float min, float max) = 0;
  virtual void addVerticalBargraph(const char* label, float* zone, float min, float max) = 0;

  virtual void openTabBox(const char* label) = 0;
  virtual void openHorizontalBox(const char* label) = 0;
  virtual void openVerticalBox(const char* label) = 0;
  virtual void closeBox() = 0;

  virtual void run() = 0;

  void stop() { fStopped = true; }
  bool stopped() { return fStopped; }

  virtual void declare(float* zone, const char* key, const char* value) {}
};

/***************************************************************************
   LV2 UI interface
 ***************************************************************************/

enum ui_elem_type_t {
  UI_BUTTON, UI_CHECK_BUTTON,
  UI_V_SLIDER, UI_H_SLIDER, UI_NUM_ENTRY,
  UI_V_BARGRAPH, UI_H_BARGRAPH,
  UI_END_GROUP, UI_V_GROUP, UI_H_GROUP, UI_T_GROUP
};

struct ui_elem_t {
  ui_elem_type_t type;
  const char *label;
  int port;
  float *zone;
  void *ref;
  float init, min, max, step;
};

class LV2UI : public UI
{
public:
  int nelems, nports;
  ui_elem_t *elems;
  map< int, list<strpair> > metadata;

  LV2UI();
  virtual ~LV2UI();

protected:
  void add_elem(ui_elem_type_t type, const char *label = NULL);
  void add_elem(ui_elem_type_t type, const char *label, float *zone);
  void add_elem(ui_elem_type_t type, const char *label, float *zone,
		float init, float min, float max, float step);
  void add_elem(ui_elem_type_t type, const char *label, float *zone,
		float min, float max);

public:
  virtual void addButton(const char* label, float* zone);
  virtual void addCheckButton(const char* label, float* zone);
  virtual void addVerticalSlider(const char* label, float* zone, float init, float min, float max, float step);
  virtual void addHorizontalSlider(const char* label, float* zone, float init, float min, float max, float step);
  virtual void addNumEntry(const char* label, float* zone, float init, float min, float max, float step);

  virtual void addHorizontalBargraph(const char* label, float* zone, float min, float max);
  virtual void addVerticalBargraph(const char* label, float* zone, float min, float max);

  virtual void openTabBox(const char* label);
  virtual void openHorizontalBox(const char* label);
  virtual void openVerticalBox(const char* label);
  virtual void closeBox();

  virtual void run();

  virtual void declare(float* zone, const char* key, const char* value);
};

LV2UI::LV2UI()
{
  nelems = nports = 0;
  elems = NULL;
}

LV2UI::~LV2UI()
{
  if (elems) free(elems);
}

void LV2UI::declare(float* zone, const char* key, const char* value)
{
  map< int, list<strpair> >::iterator it = metadata.find(nelems);
  if (it != metadata.end())
    it->second.push_back(strpair(key, value));
  else
    metadata[nelems] = list<strpair>(1, strpair(key, value));
}

inline void LV2UI::add_elem(ui_elem_type_t type, const char *label)
{
  ui_elem_t *elems1 = (ui_elem_t*)realloc(elems, (nelems+1)*sizeof(ui_elem_t));
  if (elems1)
    elems = elems1;
  else
    return;
  elems[nelems].type = type;
  elems[nelems].label = label;
  elems[nelems].port = -1;
  elems[nelems].zone = NULL;
  elems[nelems].ref = NULL;
  elems[nelems].init = 0.0;
  elems[nelems].min = 0.0;
  elems[nelems].max = 0.0;
  elems[nelems].step = 0.0;
  nelems++;
}

inline void LV2UI::add_elem(ui_elem_type_t type, const char *label, float *zone)
{
  ui_elem_t *elems1 = (ui_elem_t*)realloc(elems, (nelems+1)*sizeof(ui_elem_t));
  if (elems1)
    elems = elems1;
  else
    return;
  elems[nelems].type = type;
  elems[nelems].label = label;
  elems[nelems].port = nports++;
  elems[nelems].zone = zone;
  elems[nelems].ref = NULL;
  elems[nelems].init = 0.0;
  elems[nelems].min = 0.0;
  elems[nelems].max = 0.0;
  elems[nelems].step = 0.0;
  nelems++;
}

inline void LV2UI::add_elem(ui_elem_type_t type, const char *label, float *zone,
			     float init, float min, float max, float step)
{
  ui_elem_t *elems1 = (ui_elem_t*)realloc(elems, (nelems+1)*sizeof(ui_elem_t));
  if (elems1)
    elems = elems1;
  else
    return;
  elems[nelems].type = type;
  elems[nelems].label = label;
  elems[nelems].port = nports++;
  elems[nelems].zone = zone;
  elems[nelems].ref = NULL;
  elems[nelems].init = init;
  elems[nelems].min = min;
  elems[nelems].max = max;
  elems[nelems].step = step;
  nelems++;
}

inline void LV2UI::add_elem(ui_elem_type_t type, const char *label, float *zone,
			     float min, float max)
{
  ui_elem_t *elems1 = (ui_elem_t*)realloc(elems, (nelems+1)*sizeof(ui_elem_t));
  if (elems1)
    elems = elems1;
  else
    return;
  elems[nelems].type = type;
  elems[nelems].label = label;
  elems[nelems].port = nports++;
  elems[nelems].zone = zone;
  elems[nelems].ref = NULL;
  elems[nelems].init = 0.0;
  elems[nelems].min = min;
  elems[nelems].max = max;
  elems[nelems].step = 0.0;
  nelems++;
}

void LV2UI::addButton(const char* label, float* zone)
{ add_elem(UI_BUTTON, label, zone); }
void LV2UI::addCheckButton(const char* label, float* zone)
{ add_elem(UI_CHECK_BUTTON, label, zone); }
void LV2UI::addVerticalSlider(const char* label, float* zone, float init, float min, float max, float step)
{ add_elem(UI_V_SLIDER, label, zone, init, min, max, step); }
void LV2UI::addHorizontalSlider(const char* label, float* zone, float init, float min, float max, float step)
{ add_elem(UI_H_SLIDER, label, zone, init, min, max, step); }
void LV2UI::addNumEntry(const char* label, float* zone, float init, float min, float max, float step)
{ add_elem(UI_NUM_ENTRY, label, zone, init, min, max, step); }

void LV2UI::addHorizontalBargraph(const char* label, float* zone, float min, float max)
{ add_elem(UI_H_BARGRAPH, label, zone, min, max); }
void LV2UI::addVerticalBargraph(const char* label, float* zone, float min, float max)
{ add_elem(UI_V_BARGRAPH, label, zone, min, max); }

void LV2UI::openTabBox(const char* label)
{ add_elem(UI_T_GROUP, label); }
void LV2UI::openHorizontalBox(const char* label)
{ add_elem(UI_H_GROUP, label); }
void LV2UI::openVerticalBox(const char* label)
{ add_elem(UI_V_GROUP, label); }
void LV2UI::closeBox()
{ add_elem(UI_END_GROUP); }

void LV2UI::run() {}

/******************************************************************************
*******************************************************************************

			    FAUST DSP

*******************************************************************************
*******************************************************************************/

//----------------------------------------------------------------
//  abstract definition of a signal processor
//----------------------------------------------------------------

class dsp {
 protected:
  int fSamplingFreq;
 public:
  // internal freelist for custom voice allocation
  dsp *prev, *next;
  dsp() {}
  virtual ~dsp() {}
  virtual int getNumInputs() = 0;
  virtual int getNumOutputs() = 0;
  virtual void buildUserInterface(UI* interface) = 0;
  virtual void init(int samplingRate) = 0;
  virtual void compute(int len, float** inputs, float** outputs) = 0;
};

//----------------------------------------------------------------------------
//  FAUST generated signal processor
//----------------------------------------------------------------------------

#ifndef FAUSTFLOAT
#define FAUSTFLOAT float
#endif  

typedef long double quad;

#ifndef FAUSTCLASS 
#define FAUSTCLASS zameq2
#endif

class zameq2 : public dsp {
  private:
	FAUSTFLOAT 	fslider0;
	double 	fRec0[2];
	FAUSTFLOAT 	fslider1;
	double 	fRec1[2];
	int 	iConst0;
	double 	fConst1;
	FAUSTFLOAT 	fslider2;
	double 	fRec2[2];
	FAUSTFLOAT 	fslider3;
	double 	fRec4[2];
	FAUSTFLOAT 	fslider4;
	double 	fRec5[2];
	FAUSTFLOAT 	fslider5;
	double 	fRec6[2];
	FAUSTFLOAT 	fslider6;
	double 	fRec8[2];
	FAUSTFLOAT 	fslider7;
	double 	fRec9[2];
	int 	iConst2;
	double 	fConst3;
	FAUSTFLOAT 	fslider8;
	double 	fRec10[2];
	double 	fConst4;
	double 	fConst5;
	double 	fConst6;
	FAUSTFLOAT 	fslider9;
	double 	fRec12[2];
	FAUSTFLOAT 	fslider10;
	double 	fRec13[2];
	FAUSTFLOAT 	fslider11;
	double 	fRec14[2];
	double 	fRec15[3];
	double 	fRec11[3];
	double 	fRec7[3];
	double 	fRec3[3];
	FAUSTFLOAT 	fslider12;
	double 	fRec16[2];
	double 	fRec18[3];
	double 	fRec17[3];
	FAUSTFLOAT 	fcheckbox0;
	FAUSTFLOAT 	fslider13;
	double 	fRec19[2];
	double 	fRec21[3];
	double 	fRec20[3];
	FAUSTFLOAT 	fcheckbox1;
	FAUSTFLOAT 	fslider14;
	double 	fRec22[2];
  public:
	static void metadata(Meta* m) 	{ 
		m->declare("name", "ZamEQ2");
		m->declare("author", "Damien Zammit");
		m->declare("copyright", "2013");
		m->declare("version", "1.00");
		m->declare("license", "GPLv2");
		m->declare("math.lib/name", "Math Library");
		m->declare("math.lib/author", "GRAME");
		m->declare("math.lib/copyright", "GRAME");
		m->declare("math.lib/version", "1.0");
		m->declare("math.lib/license", "LGPL with exception");
		m->declare("filter.lib/name", "Faust Filter Library");
		m->declare("filter.lib/author", "Julius O. Smith (jos at ccrma.stanford.edu)");
		m->declare("filter.lib/copyright", "Julius O. Smith III");
		m->declare("filter.lib/version", "1.29");
		m->declare("filter.lib/license", "STK-4.3");
		m->declare("filter.lib/reference", "https://ccrma.stanford.edu/~jos/filters/");
		m->declare("music.lib/name", "Music Library");
		m->declare("music.lib/author", "GRAME");
		m->declare("music.lib/copyright", "GRAME");
		m->declare("music.lib/version", "1.0");
		m->declare("music.lib/license", "LGPL with exception");
	}

	virtual int getNumInputs() 	{ return 1; }
	virtual int getNumOutputs() 	{ return 1; }
	static void classInit(int samplingFreq) {
	}
	virtual void instanceInit(int samplingFreq) {
		fSamplingFreq = samplingFreq;
		fslider0 = 0.0;
		for (int i=0; i<2; i++) fRec0[i] = 0;
		fslider1 = 1e+02;
		for (int i=0; i<2; i++) fRec1[i] = 0;
		iConst0 = min(192000, max(1, fSamplingFreq));
		fConst1 = (2764.601535159018 / double(iConst0));
		fslider2 = 1.0;
		for (int i=0; i<2; i++) fRec2[i] = 0;
		fslider3 = 0.0;
		for (int i=0; i<2; i++) fRec4[i] = 0;
		fslider4 = 4e+01;
		for (int i=0; i<2; i++) fRec5[i] = 0;
		fslider5 = 1.0;
		for (int i=0; i<2; i++) fRec6[i] = 0;
		fslider6 = 0.0;
		for (int i=0; i<2; i++) fRec8[i] = 0;
		fslider7 = 88.0;
		for (int i=0; i<2; i++) fRec9[i] = 0;
		iConst2 = faustpower<2>(iConst0);
		fConst3 = (7643021.648203599 / double(iConst2));
		fslider8 = 0.0;
		for (int i=0; i<2; i++) fRec10[i] = 0;
		fConst4 = (1910755.4120508998 / double(iConst2));
		fConst5 = (1382.300767579509 / double(iConst0));
		fConst6 = (2.2e+02 / double(iConst0));
		fslider9 = 0.0;
		for (int i=0; i<2; i++) fRec12[i] = 0;
		fslider10 = 64.0;
		for (int i=0; i<2; i++) fRec13[i] = 0;
		fslider11 = 0.0;
		for (int i=0; i<2; i++) fRec14[i] = 0;
		for (int i=0; i<3; i++) fRec15[i] = 0;
		for (int i=0; i<3; i++) fRec11[i] = 0;
		for (int i=0; i<3; i++) fRec7[i] = 0;
		for (int i=0; i<3; i++) fRec3[i] = 0;
		fslider12 = 7e+01;
		for (int i=0; i<2; i++) fRec16[i] = 0;
		for (int i=0; i<3; i++) fRec18[i] = 0;
		for (int i=0; i<3; i++) fRec17[i] = 0;
		fcheckbox0 = 0.0;
		fslider13 = 2e+01;
		for (int i=0; i<2; i++) fRec19[i] = 0;
		for (int i=0; i<3; i++) fRec21[i] = 0;
		for (int i=0; i<3; i++) fRec20[i] = 0;
		fcheckbox1 = 0.0;
		fslider14 = 0.0;
		for (int i=0; i<2; i++) fRec22[i] = 0;
	}
	virtual void init(int samplingFreq) {
		classInit(samplingFreq);
		instanceInit(samplingFreq);
	}
	virtual void buildUserInterface(UI* interface) {
		interface->openVerticalBox("zameq2");
		interface->addHorizontalSlider("0 Boost Lowshelf (dB)", &fslider3, 0.0, -2e+01, 2e+01, 0.1);
		interface->declare(&fslider4, "unit", "PK");
		interface->addHorizontalSlider("0 Freq Lowshelf", &fslider4, 4e+01, 1.0, 1e+02, 0.01);
		interface->addHorizontalSlider("0 Slope Lowshelf", &fslider5, 1.0, 1.0, 1.5, 0.1);
		interface->addHorizontalSlider("1 Boost (dB)", &fslider9, 0.0, -2e+01, 2e+01, 0.1);
		interface->declare(&fslider10, "unit", "PK");
		interface->addHorizontalSlider("1 Freq", &fslider10, 64.0, 1.0, 1e+02, 0.01);
		interface->addHorizontalSlider("1 Q (dB)", &fslider11, 0.0, -1e+01, 4e+01, 0.1);
		interface->addHorizontalSlider("2 Boost (dB)", &fslider6, 0.0, -2e+01, 2e+01, 0.1);
		interface->declare(&fslider7, "unit", "PK");
		interface->addHorizontalSlider("2 Freq", &fslider7, 88.0, 1.0, 1e+02, 0.01);
		interface->addHorizontalSlider("2 Q (dB)", &fslider8, 0.0, -1e+01, 4e+01, 0.1);
		interface->addHorizontalSlider("3 Boost Highshelf (dB)", &fslider0, 0.0, -2e+01, 2e+01, 0.1);
		interface->declare(&fslider1, "unit", "PK");
		interface->addHorizontalSlider("3 Freq Highshelf", &fslider1, 1e+02, 1.0, 1e+02, 0.01);
		interface->addHorizontalSlider("3 Slope Highshelf", &fslider2, 1.0, 1.0, 1.5, 0.1);
		interface->declare(&fslider12, "unit", "PK");
		interface->addHorizontalSlider("4 Freq Lowpass", &fslider12, 7e+01, 1.0, 1e+02, 0.01);
		interface->declare(&fslider13, "unit", "PK");
		interface->addHorizontalSlider("5 Freq Highpass", &fslider13, 2e+01, 1.0, 1e+02, 0.01);
		interface->addVerticalSlider("6 Master Trim (dB)", &fslider14, 0.0, -12.0, 12.0, 0.1);
		interface->addCheckButton("HP", &fcheckbox1);
		interface->addCheckButton("LP", &fcheckbox0);
		interface->closeBox();
	}
	virtual void compute (int count, FAUSTFLOAT** input, FAUSTFLOAT** output) {
		double 	fSlow0 = (0.010000000000000009 * fslider0);
		double 	fSlow1 = (0.010000000000000009 * fslider1);
		double 	fSlow2 = (0.010000000000000009 * fslider2);
		double 	fSlow3 = (0.010000000000000009 * fslider3);
		double 	fSlow4 = (0.010000000000000009 * fslider4);
		double 	fSlow5 = (0.010000000000000009 * fslider5);
		double 	fSlow6 = (0.010000000000000009 * fslider6);
		double 	fSlow7 = (0.010000000000000009 * fslider7);
		double 	fSlow8 = (0.010000000000000009 * fslider8);
		double 	fSlow9 = (0.010000000000000009 * fslider9);
		double 	fSlow10 = (0.010000000000000009 * fslider10);
		double 	fSlow11 = (0.010000000000000009 * fslider11);
		double 	fSlow12 = (0.010000000000000009 * fslider12);
		int 	iSlow13 = int((fcheckbox0 == 1.0));
		double 	fSlow14 = (0.010000000000000009 * fslider13);
		int 	iSlow15 = int((fcheckbox1 == 1.0));
		double 	fSlow16 = (0.010000000000000009 * fslider14);
		FAUSTFLOAT* input0 = input[0];
		FAUSTFLOAT* output0 = output[0];
		for (int i=0; i<count; i++) {
			fRec0[0] = (fSlow0 + (0.99 * fRec0[1]));
			double fTemp0 = sqrt(pow(1e+01,(0.05 * fRec0[0])));
			fRec1[0] = (fSlow1 + (0.99 * fRec1[1]));
			double fTemp1 = (fConst1 * pow(2.0,(0.08333333333333333 * (fRec1[0] - 49.0))));
			fRec2[0] = (fSlow2 + (0.99 * fRec2[1]));
			double fTemp2 = ((sqrt((2.0 + (((1.0 / fRec2[0]) - 1.0) * (fTemp0 + (1.0 / fTemp0))))) * sin(fTemp1)) * sqrt(fTemp0));
			double fTemp3 = (fTemp0 + fTemp2);
			double fTemp4 = cos(fTemp1);
			double fTemp5 = (fTemp4 * (1.0 - fTemp0));
			double fTemp6 = (1.0 + (fTemp5 + fTemp3));
			double fTemp7 = (1.0 + fTemp0);
			int iTemp8 = int((fTemp6 == 0.0));
			fRec4[0] = (fSlow3 + (0.99 * fRec4[1]));
			double fTemp9 = sqrt(pow(1e+01,(0.05 * fRec4[0])));
			fRec5[0] = (fSlow4 + (0.99 * fRec5[1]));
			double fTemp10 = (fConst1 * pow(2.0,(0.08333333333333333 * (fRec5[0] - 49.0))));
			fRec6[0] = (fSlow5 + (0.99 * fRec6[1]));
			double fTemp11 = ((sqrt((2.0 + (((1.0 / fRec6[0]) - 1.0) * (fTemp9 + (1.0 / fTemp9))))) * sin(fTemp10)) * sqrt(fTemp9));
			double fTemp12 = (fTemp9 + fTemp11);
			double fTemp13 = cos(fTemp10);
			double fTemp14 = (fTemp13 * (fTemp9 - 1.0));
			double fTemp15 = (1.0 + (fTemp14 + fTemp12));
			double fTemp16 = (1.0 + fTemp9);
			int iTemp17 = int((fTemp15 == 0.0));
			fRec8[0] = (fSlow6 + (0.99 * fRec8[1]));
			double fTemp18 = exp((0.1151292546497023 * fRec8[0]));
			int iTemp19 = int((fRec8[0] == 0.0));
			double fTemp20 = faustpower<2>(((iTemp19)?1.0:((int((fRec8[0] < 0.0)))?(1.4125375446227544 * fTemp18):(0.7079457843841379 * fTemp18))));
			double fTemp21 = faustpower<2>(fTemp18);
			double fTemp22 = fabs((fTemp21 - fTemp20));
			fRec9[0] = (fSlow7 + (0.99 * fRec9[1]));
			double fTemp23 = pow(2.0,(0.08333333333333333 * (fRec9[0] - 49.0)));
			double fTemp24 = faustpower<2>(fTemp23);
			double fTemp25 = faustpower<2>(((fConst3 * fTemp24) - 9.869604401089358));
			fRec10[0] = (fSlow8 + (0.99 * fRec10[1]));
			double fTemp26 = exp((0.1151292546497023 * fRec10[0]));
			double fTemp27 = (fTemp22 * faustpower<2>(fTemp26));
			double fTemp28 = fabs((fTemp20 - 1.0));
			double fTemp29 = sqrt(((fTemp25 + (fConst4 * (((fTemp21 * fTemp24) * fTemp28) / fTemp27))) / ((fConst4 * ((fTemp24 * fTemp28) / fTemp27)) + fTemp25)));
			double fTemp30 = faustpower<2>(fTemp29);
			double fTemp31 = fabs((fTemp20 - fTemp30));
			double fTemp32 = (fabs((fTemp20 - fTemp29)) - sqrt((fTemp28 * fTemp31)));
			double fTemp33 = fabs((fTemp21 - 1.0));
			double fTemp34 = fabs((fTemp21 - fTemp30));
			double fTemp35 = (fabs((fTemp21 - fTemp29)) - sqrt((fTemp34 * fTemp33)));
			double fTemp36 = sqrt((fTemp34 / fTemp33));
			double fTemp37 = faustpower<2>(tan((fConst5 * fTemp23)));
			double fTemp38 = (fTemp37 * fTemp36);
			double fTemp39 = ((fTemp31 * faustpower<2>(tan((fConst6 * (fTemp23 / fTemp26))))) * faustpower<2>((1.0 + (fTemp38 * sqrt((fTemp28 / fTemp31))))));
			double fTemp40 = sqrt(((fTemp39 + (fTemp38 * ((2.0 * fTemp35) - (2.0 * fTemp32)))) / fTemp22));
			double fTemp41 = (1.0 + (fTemp38 + fTemp40));
			int iTemp42 = int((fTemp41 == 0.0));
			fRec12[0] = (fSlow9 + (0.99 * fRec12[1]));
			double fTemp43 = exp((0.1151292546497023 * fRec12[0]));
			int iTemp44 = int((fRec12[0] == 0.0));
			double fTemp45 = faustpower<2>(((iTemp44)?1.0:((int((fRec12[0] < 0.0)))?(1.4125375446227544 * fTemp43):(0.7079457843841379 * fTemp43))));
			double fTemp46 = faustpower<2>(fTemp43);
			double fTemp47 = fabs((fTemp46 - fTemp45));
			fRec13[0] = (fSlow10 + (0.99 * fRec13[1]));
			double fTemp48 = pow(2.0,(0.08333333333333333 * (fRec13[0] - 49.0)));
			double fTemp49 = faustpower<2>(fTemp48);
			double fTemp50 = faustpower<2>(((fConst3 * fTemp49) - 9.869604401089358));
			fRec14[0] = (fSlow11 + (0.99 * fRec14[1]));
			double fTemp51 = exp((0.1151292546497023 * fRec14[0]));
			double fTemp52 = (fTemp47 * faustpower<2>(fTemp51));
			double fTemp53 = fabs((fTemp45 - 1.0));
			double fTemp54 = sqrt(((fTemp50 + (fConst4 * (((fTemp46 * fTemp49) * fTemp53) / fTemp52))) / ((fConst4 * ((fTemp49 * fTemp53) / fTemp52)) + fTemp50)));
			double fTemp55 = faustpower<2>(fTemp54);
			double fTemp56 = fabs((fTemp45 - fTemp55));
			double fTemp57 = (fabs((fTemp45 - fTemp54)) - sqrt((fTemp53 * fTemp56)));
			double fTemp58 = fabs((fTemp46 - 1.0));
			double fTemp59 = fabs((fTemp46 - fTemp55));
			double fTemp60 = (fabs((fTemp46 - fTemp54)) - sqrt((fTemp59 * fTemp58)));
			double fTemp61 = sqrt((fTemp59 / fTemp58));
			double fTemp62 = faustpower<2>(tan((fConst5 * fTemp48)));
			double fTemp63 = (fTemp62 * fTemp61);
			double fTemp64 = ((faustpower<2>((1.0 + (fTemp63 * sqrt((fTemp53 / fTemp56))))) * fTemp56) * faustpower<2>(tan((fConst6 * (fTemp48 / fTemp51)))));
			double fTemp65 = sqrt(((fTemp64 + (fTemp63 * ((2.0 * fTemp60) - (2.0 * fTemp57)))) / fTemp47));
			double fTemp66 = (1.0 + (fTemp63 + fTemp65));
			int iTemp67 = int((fTemp66 == 0.0));
			double fTemp68 = (double)input0[i];
			fRec15[0] = (0 - (((fRec15[2] * ((iTemp67)?0.0:(((1 + fTemp63) - fTemp65) / fTemp66))) + (fRec15[1] * ((iTemp67)?0.0:((0 - (2.0 * (1.0 - fTemp63))) / fTemp66)))) - fTemp68));
			double fTemp69 = sqrt((((2.0 * (((fTemp45 * fTemp60) * fTemp62) * fTemp61)) + (fTemp46 * (fTemp64 - (2.0 * (fTemp63 * fTemp57))))) / fTemp47));
			double fTemp70 = (fTemp54 + fTemp63);
			double fTemp71 = ((iTemp44)?fTemp68:((fRec15[0] * ((iTemp67)?0.0:((fTemp69 + fTemp70) / fTemp66))) + ((fRec15[2] * ((iTemp67)?1.0:((fTemp70 - fTemp69) / fTemp66))) + (fRec15[1] * ((iTemp67)?0.0:((0 - (2.0 * (fTemp54 - fTemp63))) / fTemp66))))));
			fRec11[0] = (0 - (((fRec11[2] * ((iTemp42)?0.0:(((1 + fTemp38) - fTemp40) / fTemp41))) + (fRec11[1] * ((iTemp42)?0.0:((0 - (2.0 * (1.0 - fTemp38))) / fTemp41)))) - fTemp71));
			double fTemp72 = sqrt((((2.0 * (((fTemp20 * fTemp35) * fTemp37) * fTemp36)) + (fTemp21 * (fTemp39 - (2.0 * (fTemp38 * fTemp32))))) / fTemp22));
			double fTemp73 = (fTemp29 + fTemp38);
			double fTemp74 = ((iTemp19)?fTemp71:((fRec11[0] * ((iTemp42)?0.0:((fTemp72 + fTemp73) / fTemp41))) + ((fRec11[2] * ((iTemp42)?1.0:((fTemp73 - fTemp72) / fTemp41))) + (fRec11[1] * ((iTemp42)?0.0:((0 - (2.0 * (fTemp29 - fTemp38))) / fTemp41))))));
			fRec7[0] = (0 - (((fRec7[2] * ((iTemp17)?0.0:(((1.0 + (fTemp9 + fTemp14)) - fTemp11) / fTemp15))) + (fRec7[1] * ((iTemp17)?0.0:((0 - (2.0 * ((fTemp9 + (fTemp13 * fTemp16)) - 1.0))) / fTemp15)))) - fTemp74));
			double fTemp75 = (fTemp13 * (1.0 - fTemp9));
			double fTemp76 = ((int((fRec4[0] == 0.0)))?fTemp74:((fRec7[0] * ((iTemp17)?0.0:((fTemp9 * (1.0 + (fTemp12 + fTemp75))) / fTemp15))) + ((fRec7[2] * ((iTemp17)?1.0:((fTemp9 * ((1.0 + (fTemp9 + fTemp75)) - fTemp11)) / fTemp15))) + (fRec7[1] * ((iTemp17)?0.0:(2.0 * ((fTemp9 * ((fTemp9 + (fTemp13 * (0 - fTemp16))) - 1.0)) / fTemp15)))))));
			fRec3[0] = (0 - (((fRec3[2] * ((iTemp8)?0.0:(((1.0 + (fTemp0 + fTemp5)) - fTemp2) / fTemp6))) + (fRec3[1] * ((iTemp8)?0.0:(2.0 * (((fTemp0 + (fTemp4 * (0 - fTemp7))) - 1.0) / fTemp6))))) - fTemp76));
			double fTemp77 = (fTemp4 * (fTemp0 - 1.0));
			double fTemp78 = ((int((fRec0[0] == 0.0)))?fTemp76:((fRec3[0] * ((iTemp8)?0.0:((fTemp0 * (1.0 + (fTemp3 + fTemp77))) / fTemp6))) + ((fRec3[2] * ((iTemp8)?1.0:((fTemp0 * ((1.0 + (fTemp0 + fTemp77)) - fTemp2)) / fTemp6))) + (fRec3[1] * ((iTemp8)?0.0:((((fTemp0 + (fTemp4 * fTemp7)) - 1.0) * (0 - (2.0 * fTemp0))) / fTemp6))))));
			fRec16[0] = (fSlow12 + (0.99 * fRec16[1]));
			double fTemp79 = tan((fConst5 * pow(2.0,(0.08333333333333333 * (fRec16[0] - 49.0)))));
			double fTemp80 = (1.0 / fTemp79);
			double fTemp81 = (1 + ((0.7653668647301795 + fTemp80) / fTemp79));
			double fTemp82 = (1 - (1.0 / faustpower<2>(fTemp79)));
			double fTemp83 = (1 + ((1.8477590650225735 + fTemp80) / fTemp79));
			fRec18[0] = (fTemp78 - ((((1 + ((fTemp80 - 1.8477590650225735) / fTemp79)) * fRec18[2]) + (2 * (fTemp82 * fRec18[1]))) / fTemp83));
			fRec17[0] = (((fRec18[2] + (fRec18[0] + (2 * fRec18[1]))) / fTemp83) - ((((1 + ((fTemp80 - 0.7653668647301795) / fTemp79)) * fRec17[2]) + (2 * (fRec17[1] * fTemp82))) / fTemp81));
			double fTemp84 = ((iSlow13)?((fRec17[2] + (fRec17[0] + (2 * fRec17[1]))) / fTemp81):fTemp78);
			fRec19[0] = (fSlow14 + (0.99 * fRec19[1]));
			double fTemp85 = tan((fConst5 * pow(2.0,(0.08333333333333333 * (fRec19[0] - 49.0)))));
			double fTemp86 = (1.0 / fTemp85);
			double fTemp87 = (1 + ((0.7653668647301795 + fTemp86) / fTemp85));
			double fTemp88 = faustpower<2>(fTemp85);
			double fTemp89 = (1.0 / fTemp88);
			double fTemp90 = (1 - fTemp89);
			double fTemp91 = (1 + ((1.8477590650225735 + fTemp86) / fTemp85));
			fRec21[0] = (fTemp84 - ((((1 + ((fTemp86 - 1.8477590650225735) / fTemp85)) * fRec21[2]) + (2 * (fTemp90 * fRec21[1]))) / fTemp91));
			double fTemp92 = (0 - fTemp89);
			fRec20[0] = (((((fRec21[0] / fTemp88) + (2 * (fRec21[1] * fTemp92))) + (fRec21[2] / fTemp88)) / fTemp91) - ((((1 + ((fTemp86 - 0.7653668647301795) / fTemp85)) * fRec20[2]) + (2 * (fRec20[1] * fTemp90))) / fTemp87));
			fRec22[0] = (fSlow16 + (0.99 * fRec22[1]));
			output0[i] = (FAUSTFLOAT)(exp((0.1151292546497023 * fRec22[0])) * ((iSlow15)?((((fRec20[0] / fTemp88) + (2 * (fRec20[1] * fTemp92))) + (fRec20[2] / fTemp88)) / fTemp87):fTemp84));
			// post processing
			fRec22[1] = fRec22[0];
			fRec20[2] = fRec20[1]; fRec20[1] = fRec20[0];
			fRec21[2] = fRec21[1]; fRec21[1] = fRec21[0];
			fRec19[1] = fRec19[0];
			fRec17[2] = fRec17[1]; fRec17[1] = fRec17[0];
			fRec18[2] = fRec18[1]; fRec18[1] = fRec18[0];
			fRec16[1] = fRec16[0];
			fRec3[2] = fRec3[1]; fRec3[1] = fRec3[0];
			fRec7[2] = fRec7[1]; fRec7[1] = fRec7[0];
			fRec11[2] = fRec11[1]; fRec11[1] = fRec11[0];
			fRec15[2] = fRec15[1]; fRec15[1] = fRec15[0];
			fRec14[1] = fRec14[0];
			fRec13[1] = fRec13[0];
			fRec12[1] = fRec12[0];
			fRec10[1] = fRec10[0];
			fRec9[1] = fRec9[0];
			fRec8[1] = fRec8[0];
			fRec6[1] = fRec6[0];
			fRec5[1] = fRec5[0];
			fRec4[1] = fRec4[0];
			fRec2[1] = fRec2[0];
			fRec1[1] = fRec1[0];
			fRec0[1] = fRec0[0];
		}
	}
};



//----------------------------------------------------------------------------
//  LV2 interface
//----------------------------------------------------------------------------

//#line 391 "lv2.cpp"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* This enables automatic MIDI controller mapping based on the midi:ctrl
   attributes in the Faust source. We have this enabled by default, but you
   may have to disable it if the custom controller mapping gets in the way of
   the automation facilities that the host provides. (But then again if the
   host wants to do its own controller mapping then it probably won't, or at
   least shouldn't, send us the MIDI controllers in the first place.) */
#ifndef FAUST_MIDICC
#define FAUST_MIDICC 1
#endif

#include <lv2/lv2plug.in/ns/lv2core/lv2.h>
#include <lv2/lv2plug.in/ns/ext/dynmanifest/dynmanifest.h>
#if FAUST_MIDICC
#include <lv2/lv2plug.in/ns/ext/event/event-helpers.h>
#include <lv2/lv2plug.in/ns/ext/uri-map/uri-map.h>
#include <lv2/lv2plug.in/ns/ext/urid/urid.h>
#define MIDI_EVENT_URI "http://lv2plug.in/ns/ext/midi#MidiEvent"
#endif

#ifndef URI_PREFIX
#define URI_PREFIX "http://faust-lv2.googlecode.com"
#endif

#ifndef PLUGIN_URI
#define PLUGIN_URI URI_PREFIX "/zameq2"
#endif

/* This allows various manifest data to be generated from the corresponding
   metadata (author, name, description, license) in the Faust source. */
#ifndef FAUST_META
#define FAUST_META 1
#endif

// You can define these for various debugging output items.
//#define DEBUG_META 1 // recognized MIDI controller metadata
//#define DEBUG_MIDI 1 // incoming MIDI messages
//#define DEBUG_MIDICC 1 // controller messages

struct LV2Plugin {
  bool active;		// activation status
  int rate;		// sampling rate
  zameq2 *dsp;		// the dsp
  LV2UI *ui;		// its Faust interface description
  int n_in, n_out;	// number of input and output control ports
  int *ctrls;		// Faust ui elements (indices into ui->elems)
  float **ports;	// corresponding LV2 data
  float *portvals;	// cached port data from the last run
  int *inctrls, *outctrls;	// indices for active and passive controls
  float **inputs, **outputs;	// audio buffers
#if FAUST_MIDICC
  LV2_Event_Buffer* event_port;	// midi input
  std::map<uint8_t,int> ctrlmap; // MIDI controller map
  // Needed host features (the uri-map extension is officially deprecated, but
  // still needed for some if not most hosts at the time of this writing).
  LV2_URI_Map_Feature* uri_map;
  LV2_URID_Map* map;	// the new urid extension
  LV2_URID midi_event;	// midi event uri
  LV2_Event_Feature* event_ref;
#endif

  LV2Plugin() {
    active = false;
    rate = 44100;
    n_in = n_out = 0;
    dsp = NULL;
    ui = NULL;
    ctrls = inctrls = outctrls = NULL;
    ports = inputs = outputs = NULL;
    portvals = NULL;
#if FAUST_MIDICC
    uri_map = NULL; map = NULL;
    midi_event = -1;
    event_ref = NULL;
    event_port = NULL;
#endif
  }
};

static LV2_Handle
instantiate(const LV2_Descriptor*     descriptor,
            double                    rate,
            const char*               bundle_path,
            const LV2_Feature* const* features)
{
  LV2Plugin* plugin = new LV2Plugin;
#if FAUST_MIDICC
  // Scan host features for URID (or URI) map.
  for (int i = 0; features[i]; i++) {
    if (!strcmp(features[i]->URI, LV2_URID_URI "#map")) {
      plugin->map = (LV2_URID_Map*)features[i]->data;
      plugin->midi_event =
	plugin->map->map(plugin->map->handle, MIDI_EVENT_URI);
    } else if (!strcmp(features[i]->URI, LV2_URI_MAP_URI)) {
      plugin->uri_map = (LV2_URI_Map_Feature*)features[i]->data;
      plugin->midi_event =
	plugin->uri_map->uri_to_id(plugin->uri_map->callback_data,
				   LV2_EVENT_URI, MIDI_EVENT_URI);
    } else if (!strcmp(features[i]->URI, LV2_EVENT_URI)) {
      plugin->event_ref = (LV2_Event_Feature *)features[i]->data;
    }
  }
  if (!plugin->map && !plugin->uri_map) {
    fprintf
      (stderr, "%s: host supports neither uri-map or urid:map, giving up\n",
       PLUGIN_URI);
    delete plugin;
    return 0;
  }
#endif
  plugin->rate = rate;
  plugin->dsp = new zameq2();
  plugin->ui = new LV2UI();
  plugin->dsp->init(plugin->rate);
  plugin->dsp->buildUserInterface(plugin->ui);
  // The LV2 ports are numbered as follows: 0..k-1 are the control ports, then
  // come the n audio input ports, finally the m audio output ports (and
  // possibly a MIDI input port).
  int k = plugin->ui->nports, p = 0, q = 0;
  int n = plugin->dsp->getNumInputs(), m = plugin->dsp->getNumOutputs();
  // Allocate tables for the control elements and their LV2 ports.
  plugin->ctrls = (int*)calloc(k, sizeof(int));
  plugin->inctrls = (int*)calloc(k, sizeof(int));
  plugin->outctrls = (int*)calloc(k, sizeof(int));
  plugin->ports = (float**)calloc(k, sizeof(float*));
  plugin->portvals = (float*)calloc(k, sizeof(float));
  assert(k == 0 || (plugin->ctrls && plugin->inctrls && plugin->outctrls &&
		    plugin->ports && plugin->portvals));
  // Scan the Faust UI for active and passive controls which become the
  // input and output control ports of the LV2 plugin, respectively.
  for (int i = 0, j = 0; i < plugin->ui->nelems; i++) {
    switch (plugin->ui->elems[i].type) {
    case UI_T_GROUP: case UI_H_GROUP: case UI_V_GROUP: case UI_END_GROUP:
      // control groups
      break;
    case UI_H_BARGRAPH: case UI_V_BARGRAPH:
      // passive controls (output ports)
      plugin->ctrls[j++] = i;
      plugin->outctrls[q++] = i;
      break;
    default:
      // active controls (input ports)
#if FAUST_MIDICC
      {
	std::map< int, list<strpair> >::iterator it =
	  plugin->ui->metadata.find(i);
	if (it != plugin->ui->metadata.end()) {
	  // Scan for controller mappings.
	  for (std::list<strpair>::iterator jt = it->second.begin();
	       jt != it->second.end(); jt++) {
	    const char *key = jt->first, *val = jt->second;
#if DEBUG_META
	    fprintf(stderr, "ctrl '%s' meta: '%s' -> '%s'\n",
		    plugin->ui[0]->elems[i].label, key, val);
#endif
	    if (strcmp(key, "midi")) continue;
	    unsigned num;
	    if (sscanf(val, "ctrl %u", &num) < 1) continue;
#if 0 // enable this to get feedback about controller assignments
	    fprintf(stderr, "%s: cc %d -> %s\n", PLUGIN_URI, num,
		    plugin->ui->elems[i].label);
#endif
	    plugin->ctrlmap.insert(std::pair<uint8_t,int>(num, p));
	  }
	}
      }
#endif
      plugin->ctrls[j++] = i;
      plugin->inctrls[p++] = i;
      int p = plugin->ui->elems[i].port;
      float val = plugin->ui->elems[i].init;
      plugin->portvals[p] = val;
      break;
    }
  }
  // Realloc the inctrls and outctrls vectors to their appropriate sizes.
  plugin->inctrls = (int*)realloc(plugin->inctrls, p*sizeof(int));
  assert(p == 0 || plugin->inctrls);
  plugin->outctrls = (int*)realloc(plugin->outctrls, q*sizeof(int));
  assert(q == 0 || plugin->outctrls);
  plugin->n_in = p; plugin->n_out = q;
  // Allocate vectors for the audio input and output ports. Like
  // plugin->ports, these will be initialized in the connect_port callback.
  plugin->inputs = (float**)calloc(n, sizeof(float*));
  assert(n == 0 || plugin->inputs);
  plugin->outputs = (float**)calloc(m, sizeof(float*));
  assert(m == 0 || plugin->outputs);
  return (LV2_Handle)plugin;
}

static void
cleanup(LV2_Handle instance)
{
  LV2Plugin* plugin = (LV2Plugin*)instance;
  delete plugin->dsp;
  delete plugin->ui;
  free(plugin->ctrls);
  free(plugin->inctrls);
  free(plugin->outctrls);
  free(plugin->ports);
  free(plugin->portvals);
  free(plugin->inputs);
  free(plugin->outputs);
  delete plugin;
}

static void
connect_port(LV2_Handle instance,
             uint32_t   port,
             void*      data)
{
  LV2Plugin* plugin = (LV2Plugin*)instance;
  int i = port, k = plugin->ui->nports;
  int n = plugin->dsp->getNumInputs(), m = plugin->dsp->getNumOutputs();
  if (i < k)
    plugin->ports[i] = (float*)data;
  else {
    i -= k;
    if (i < n)
      plugin->inputs[i] = (float*)data;
    else {
      i -= n;
      if (i < m)
	plugin->outputs[i] = (float*)data;
#if FAUST_MIDICC
      else if (i == m)
	plugin->event_port = (LV2_Event_Buffer*)data;
#endif
      else
	fprintf(stderr, "%s: bad port number %u\n", PLUGIN_URI, port);
    }
  }
}

#if FAUST_MIDICC
static float ctrlval(const ui_elem_t &el, uint8_t v)
{
  // Translate the given MIDI controller value to the range and stepsize
  // indicated by the Faust control.
  switch (el.type) {
  case UI_BUTTON: case UI_CHECK_BUTTON:
    return (float)(v>=64);
  default:
    /* Continuous controllers. The problem here is that the range 0..127 is
       not symmetric. We'd like to map 64 to the center of the range
       (max-min)/2 and at the same time retain the full control range
       min..max. So let's just pretend that there are 128 controller values
       and map value 127 to the max value anyway. */
    if (v==127)
      return el.max;
    else
      // XXXFIXME: We might want to add proper quantization according to
      // el.step here.
      return el.min+(el.max-el.min)*v/128;
  }
}
#endif

static void
run(LV2_Handle instance, uint32_t n_samples)
{
  LV2Plugin* plugin = (LV2Plugin*)instance;
  int n = plugin->dsp->getNumInputs(), m = plugin->dsp->getNumOutputs();
  AVOIDDENORMALS;
  if (!plugin->active) {
    if (n == m) {
      // copy inputs to outputs
      for (int i = 0; i < m; i++)
	for (unsigned j = 0; j < n_samples; j++)
	  plugin->outputs[i][j] = plugin->inputs[i][j];
    } else {
      // silence
      for (int i = 0; i < m; i++)
	for (unsigned j = 0; j < n_samples; j++)
	  plugin->outputs[i][j] = 0.0f;
    }
    return;
  }
#if FAUST_MIDICC
  if (!plugin->ctrlmap.empty() && plugin->event_port) {
    // Process incoming MIDI events.
    LV2_Event_Iterator i;
    for (lv2_event_begin(&i, plugin->event_port);
	 lv2_event_is_valid(&i);
	 lv2_event_increment(&i)) {
      LV2_Event* ev = lv2_event_get(&i, NULL);
      if (ev->type == 0) {
	if (plugin->event_ref) {
	  plugin->event_ref->lv2_event_unref
	    (plugin->event_ref->callback_data, ev);
	}
      } else if (ev->type == plugin->midi_event) {
	uint8_t *data = (uint8_t*)(ev+1);
#if DEBUG_MIDI
	fprintf(stderr, "midi ev (%u bytes):", ev->size);
	for (unsigned i = 0; i < ev->size; i++)
	  fprintf(stderr, " 0x%0x", data[i]);
	fprintf(stderr, "\n");
#endif
	uint8_t status = data[0] & 0xf0, chan = data[0] & 0x0f;
	if (status == 0xb0) {
	  // interpret all other controller changes according to the MIDI
	  // controller map defined in the Faust plugin itself
	  std::map<uint8_t,int>::iterator it = plugin->ctrlmap.find(data[1]);
	  if (it != plugin->ctrlmap.end()) {
	    // defined MIDI controller
	    int j = plugin->inctrls[it->second];
#if DEBUG_MIDICC
	    fprintf(stderr, "ctrl-change chan %d, ctrl %d, val %d\n", chan+1,
		    data[1], data[2]);
#endif
	    *plugin->ui->elems[j].zone = ctrlval(plugin->ui->elems[j], data[2]);
	  }
	}
      } else {
	fprintf(stderr, "%s: unknown event type %d\n", PLUGIN_URI, ev->type);
      }
    }
  }
#endif
  // Only update the controls if the port value actually changed. This is
  // necessary to preserve the MIDI controller changes (see above). Also note
  // that we do this *after* processing the MIDI controller data so that
  // manual inputs can override these.
  for (int i = 0; i < plugin->n_in; i++) {
    int j = plugin->inctrls[i], k = plugin->ui->elems[j].port;
    float &oldval = plugin->portvals[k], newval = *plugin->ports[k];
    if (newval != oldval)
      *plugin->ui->elems[j].zone = oldval = newval;
  }
  // Let Faust do all the hard work.
  plugin->dsp->compute(n_samples, plugin->inputs, plugin->outputs);
  // Finally grab the passive controls and write them back to the
  // corresponding LV2 ports.
  for (int i = 0; i < plugin->n_out; i++) {
    int j = plugin->outctrls[i], k = plugin->ui->elems[j].port;
    float *z = plugin->ui->elems[j].zone;
    *plugin->ports[k] = *z;
  }
}

static void
activate(LV2_Handle instance)
{
  LV2Plugin* plugin = (LV2Plugin*)instance;
  plugin->dsp->init(plugin->rate);
  for (int i = 0, j = 0; i < plugin->ui->nelems; i++) {
    int p = plugin->ui->elems[i].port;
    float val = plugin->ui->elems[i].init;
    plugin->portvals[p] = val;
  }
  plugin->active = true;
}

static void
deactivate(LV2_Handle instance)
{
  LV2Plugin* plugin = (LV2Plugin*)instance;
  plugin->active = false;
}

const void*
extension_data(const char* uri)
{
  return NULL;
}

static const LV2_Descriptor descriptor = {
  PLUGIN_URI,
  instantiate,
  connect_port,
  activate,
  run,
  deactivate,
  cleanup,
  extension_data
};

extern "C"
LV2_SYMBOL_EXPORT
const LV2_Descriptor*
lv2_descriptor(uint32_t index)
{
  switch (index) {
  case 0:
    return &descriptor;
  default:
    return NULL;
  }
}

//----------------------------------------------------------------------------
//  Dynamic manifest
//----------------------------------------------------------------------------

// NOTE: If your LV2 host doesn't offer this extension then you'll have to
// create a static ttl file with the descriptions of the ports yourself. You
// can also do this by compiling this code to a standalone executable while
// defining the __MAIN__ symbol. Running the executable then prints the
// manifest on stdout.

extern "C"
LV2_SYMBOL_EXPORT
int lv2_dyn_manifest_open(LV2_Dyn_Manifest_Handle *handle,
			  const LV2_Feature *const *features)
{
  LV2Plugin* plugin = new LV2Plugin;
  plugin->dsp = new zameq2();
  plugin->ui = new LV2UI();
  plugin->dsp->init(48000);
  plugin->dsp->buildUserInterface(plugin->ui);
  int k = plugin->ui->nports;
  plugin->ctrls = (int*)calloc(k, sizeof(int));
  assert(k == 0 || plugin->ctrls);
  for (int i = 0, j = 0; i < plugin->ui->nelems; i++) {
    switch (plugin->ui->elems[i].type) {
    case UI_T_GROUP: case UI_H_GROUP: case UI_V_GROUP: case UI_END_GROUP:
      // control groups
      break;
    case UI_H_BARGRAPH: case UI_V_BARGRAPH:
      // passive controls (output ports)
      plugin->ctrls[j++] = i;
      break;
    default:
      // active controls (input ports)
      plugin->ctrls[j++] = i;
      break;
    }
  }
  *handle = (LV2_Dyn_Manifest_Handle)plugin;
  return 0;
}

extern "C"
LV2_SYMBOL_EXPORT
int lv2_dyn_manifest_get_subjects(LV2_Dyn_Manifest_Handle handle,
				  FILE *fp)
{
  fprintf(fp, "@prefix lv2:  <http://lv2plug.in/ns/lv2core#> .\n\
<%s> a lv2:Plugin .\n", PLUGIN_URI);
  return 0;
}

#include <string>
#include <ctype.h>

static string mangle(const string &s)
{
  string t = s;
  size_t n = s.size();
  for (size_t i = 0; i < n; i++)
    if ((i == 0 && !isalpha(t[i]) && t[i] != '_') ||
	(!isalnum(t[i]) && t[i] != '_'))
      t[i] = '_';
  return t;
}

#if FAUST_META
static bool is_xmlstring(const char *s)
{
  // This is just a basic sanity check. The string must not contain any
  // (unescaped) newlines, carriage returns or double quotes.
  return !strchr(s, '\n') && !strchr(s, '\r') && !strchr(s, '"');
}
#endif

extern "C"
LV2_SYMBOL_EXPORT
int lv2_dyn_manifest_get_data(LV2_Dyn_Manifest_Handle handle,
			      FILE *fp,
			      const char *uri)
{
  LV2Plugin* plugin = (LV2Plugin*)handle;
  int k = plugin->ui->nports;
  int n = plugin->dsp->getNumInputs(), m = plugin->dsp->getNumOutputs();
  // Scan the global metadata for plugin name, description, license etc.
  const char *plugin_name = NULL, *plugin_author = NULL, *plugin_descr = NULL,
    *plugin_license = NULL;
#if FAUST_META
  Meta meta;
  plugin->dsp->metadata(&meta);
  for (std::list<strpair>::iterator it = meta.data.begin();
       it != meta.data.end(); it++) {
    const char *key = it->first, *val = it->second;
    if (!val || !is_xmlstring(val)) continue;
    if (!strcmp(key, "name")) {
      if (!plugin_name)
	plugin_name = val;
    } else if (!strcmp(key, "description")) {
      if (!plugin_descr)
	plugin_descr = val;
    } else if (!strcmp(key, "author")) {
      if (!plugin_author)
	plugin_author = val;
    } else if (!strcmp(key, "license")) {
      if (!plugin_license)
	plugin_license = val;
    }
  }
#endif
  if (!plugin_name) plugin_name = "zameq2";
  fprintf(fp, "@prefix doap: <http://usefulinc.com/ns/doap#> .\n\
@prefix foaf: <http://xmlns.com/foaf/0.1/> .\n\
@prefix lv2:  <http://lv2plug.in/ns/lv2core#> .\n\
@prefix epp:  <http://lv2plug.in/ns/ext/port-props#> .\n\
@prefix ev:   <http://lv2plug.in/ns/ext/event#> .\n\
@prefix rdf:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .\n\
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .\n\
@prefix units: <http://lv2plug.in/ns/extensions/units#> .\n\
<%s>\n\
       a lv2:Plugin ;\n\
       doap:name \"%s\" ;\n\
       lv2:binary <zameq2.so> ;\n\
       lv2:optionalFeature epp:supportsStrictBounds ;\n\
       lv2:optionalFeature lv2:hardRtCapable ;\n", PLUGIN_URI, plugin_name);
  if (plugin_author)
    fprintf(fp, "\
       doap:maintainer [ foaf:name \"%s\" ] ;\n", plugin_author);
  if (plugin_descr)
    fprintf(fp, "\
       doap:description \"%s\" ;\n", plugin_descr);
  if (plugin_license)
    fprintf(fp, "\
       doap:license \"%s\" ;\n", plugin_license);
  int idx = 0;
  bool have_midi = false;
  // control ports
  for (int i = 0; i < k; i++, idx++) {
    int j = plugin->ctrls[i];
    assert(idx == plugin->ui->elems[j].port);
    fprintf(fp, "%s [\n", idx==0?"    lv2:port":" ,");
    const char *label = plugin->ui->elems[j].label;
    assert(label);
    string sym = mangle(plugin->ui->elems[j].label);
    switch (plugin->ui->elems[j].type) {
    // active controls (input ports)
    case UI_BUTTON: case UI_CHECK_BUTTON:
    fprintf(fp, "\
	a lv2:ControlPort ;\n\
	a lv2:InputPort ;\n\
	lv2:index %d ;\n\
	lv2:symbol \"%s\" ;\n\
	lv2:name \"%s\" ;\n\
        lv2:portProperty epp:hasStrictBounds ;\n\
        lv2:portProperty lv2:toggled ;\n\
	lv2:default 0.00000 ;\n\
	lv2:minimum 0.00000 ;\n\
	lv2:maximum 1.00000 ;\n", idx, sym.c_str(), label);
      break;
    case UI_NUM_ENTRY: case UI_H_SLIDER: case UI_V_SLIDER:
    fprintf(fp, "\
	a lv2:ControlPort ;\n\
	a lv2:InputPort ;\n\
	lv2:index %d ;\n\
	lv2:symbol \"%s\" ;\n\
	lv2:name \"%s\" ;\n\
	lv2:default %g ;\n\
	lv2:minimum %g ;\n\
	lv2:maximum %g ;\n", idx, sym.c_str(), label,
	    plugin->ui->elems[j].init,
	    plugin->ui->elems[j].min,
	    plugin->ui->elems[j].max);
      break;
    // passive controls (output ports)
    case UI_H_BARGRAPH: case UI_V_BARGRAPH:
    fprintf(fp, "\
	a lv2:ControlPort ;\n\
	a lv2:OutputPort ;\n\
	lv2:index %d ;\n\
	lv2:symbol \"%s\" ;\n\
	lv2:name \"%s\" ;\n\
	lv2:default %g ;\n\
	lv2:minimum %g ;\n\
	lv2:maximum %g ;\n", idx, sym.c_str(), label,
	    plugin->ui->elems[j].min,
	    plugin->ui->elems[j].min,
	    plugin->ui->elems[j].max);
      break;
    default:
      assert(0 && "this can't happen");
      break;
    }
    // Scan for Faust control metadata we understand and add corresponding
    // hints to the LV2 description of the port.
    std::map< int, list<strpair> >::iterator it =
      plugin->ui->metadata.find(j);
    if (it != plugin->ui->metadata.end()) {
      for (std::list<strpair>::iterator jt = it->second.begin();
	   jt != it->second.end(); jt++) {
	const char *key = jt->first, *val = jt->second;
#if FAUST_MIDICC
	unsigned num;
	if (!strcmp(key, "midi") && sscanf(val, "ctrl %u", &num) == 1)
	  have_midi = true;
#endif
	if (!strcmp(key, "unit"))
	  fprintf(fp, "\
	units:unit [\n\
            a            units:Unit ;\n\
            units:name   \"%s\" ;\n\
            units:symbol \"%s\" ;\n\
            units:render \"%%f %s\"\n\
	] ;\n", val, val, val);
	if (strcmp(key, "lv2")) continue;
	if (!strcmp(val, "integer"))
	  fprintf(fp, "\
	lv2:portProperty lv2:integer ;\n");
	else if (!strcmp(val, "hidden") || !strcmp(val, "notOnGUI"))
	  fprintf(fp, "\
	lv2:portProperty epp:notOnGUI ;\n");
	else if (!strncmp(val, "scalepoint", 10) ||
		 !strncmp(val, "scalePoint", 10)) {
	  val += 10;
	  if (!isspace(*val)) continue;
	  char *label = (char*)malloc(strlen(val)+1);
	  float point;
	  int pos;
	  while (sscanf(val, "%s %g%n", label, &point, &pos) == 2) {
	    fprintf(fp, "\
	lv2:scalePoint [ rdfs:label \"%s\"; rdf:value %g ] ;\n",
		    label, point);
	    val += pos;
	  }
	  free(label);
	} else
	  fprintf(stderr, "%s: bad port property '%s:%s'\n", PLUGIN_URI,
		  key, val);
      }
    }
    fprintf(fp, "    ]");
  }
  // audio inputs
  for (int i = 0; i < n; i++, idx++)
    fprintf(fp, "%s [\n\
	a lv2:AudioPort ;\n\
	a lv2:InputPort ;\n\
	lv2:index %d ;\n\
	lv2:symbol \"in%d\" ;\n\
	lv2:name \"in%d\" ;\n\
    ]", idx==0?"    lv2:port":" ,", idx, i, i);
  // audio outputs
  for (int i = 0; i < m; i++, idx++)
    fprintf(fp, "%s [\n\
	a lv2:AudioPort ;\n\
	a lv2:OutputPort ;\n\
	lv2:index %d ;\n\
	lv2:symbol \"out%d\" ;\n\
	lv2:name \"out%d\" ;\n\
    ]", idx==0?"    lv2:port":" ,", idx, i, i);
  if (have_midi) {
    // midi input
    fprintf(fp, "%s [\n\
	a lv2:InputPort ;\n\
	a ev:EventPort ;\n\
	ev:supportsEvent <http://lv2plug.in/ns/ext/midi#MidiEvent> ;\n\
	lv2:index %d ;\n\
	lv2:symbol \"midiin\" ;\n\
	lv2:name \"midiin\"\n\
    ]", idx==0?"    lv2:port":" ,", idx);
    idx++;
  }
  fprintf(fp, "\n.\n");
  return 0;
}

extern "C"
LV2_SYMBOL_EXPORT
void lv2_dyn_manifest_close(LV2_Dyn_Manifest_Handle handle)
{
  LV2Plugin* plugin = (LV2Plugin*)handle;
  delete plugin->dsp;
  delete plugin->ui;
  delete plugin;
}

int main()
{
  LV2_Dyn_Manifest_Handle handle;
  LV2_Feature **features = { NULL };
  int res = lv2_dyn_manifest_open(&handle, features);
  if (res) return res;
  res = lv2_dyn_manifest_get_data(handle, stdout, PLUGIN_URI);
  return res;
}
