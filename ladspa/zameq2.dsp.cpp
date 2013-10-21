//-----------------------------------------------------
// name: "ZamEQ2"
// author: "Damien Zammit"
// copyright: "2013"
// version: "2.1"
// license: "GPLv2"
//
// Code generated with Faust 0.9.62 (http://faust.grame.fr)
//-----------------------------------------------------
/* link with  */
#include <math.h>
#ifndef FAUSTPOWER
#define FAUSTPOWER
#include <cmath>
template <int N> inline float faustpower(float x)          { return powf(x,N); } 
template <int N> inline double faustpower(double x)        { return pow(x,N); }
template <int N> inline int faustpower(int x)              { return faustpower<N/2>(x) * faustpower<N-N/2>(x); } 
template <> 	 inline int faustpower<0>(int x)            { return 1; }
template <> 	 inline int faustpower<1>(int x)            { return x; }
#endif
/************************************************************************

	IMPORTANT NOTE : this file contains two clearly delimited sections :
	the ARCHITECTURE section (in two parts) and the USER section. Each section
	is governed by its own copyright and license. Please check individually
	each section for license and copyright information.
*************************************************************************/

/*******************BEGIN ARCHITECTURE SECTION (part 1/2)****************/

/************************************************************************
    FAUST Architecture File
	Copyright (C) 2003-2011 GRAME, Centre National de Creation Musicale
    ---------------------------------------------------------------------
    This Architecture section is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 3 of
	the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
	along with this program; If not, see <http://www.gnu.org/licenses/>.

	EXCEPTION : As a special exception, you may create a larger work
	that contains this FAUST architecture section and distribute
	that work under terms of your choice, so long as this FAUST
	architecture section is not modified.


 ************************************************************************
 ************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <stack>
#include <string>
#include <iostream>
#include <map>

#include "ladspa.h"
#include "faust/gui/GUI.h"
#include "faust/misc.h"
#include "faust/audio/dsp.h"

#define sym(name) xsym(name)
#define xsym(name) #name

/******************************************************************************
*******************************************************************************

							       VECTOR INTRINSICS

*******************************************************************************
*******************************************************************************/


/********************END ARCHITECTURE SECTION (part 1/2)****************/

/**************************BEGIN USER SECTION **************************/

#ifndef FAUSTFLOAT
#define FAUSTFLOAT float
#endif  

typedef long double quad;

#ifndef FAUSTCLASS 
#define FAUSTCLASS mydsp
#endif

class mydsp : public dsp {
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
	double 	fRec14[3];
	double 	fRec11[3];
	double 	fRec7[3];
	double 	fRec3[3];
	FAUSTFLOAT 	fslider12;
	double 	fRec15[2];
	double 	fRec17[3];
	double 	fRec16[3];
	FAUSTFLOAT 	fcheckbox0;
	FAUSTFLOAT 	fslider13;
	double 	fRec18[2];
	double 	fRec20[3];
	double 	fRec19[3];
	FAUSTFLOAT 	fcheckbox1;
	FAUSTFLOAT 	fslider14;
	double 	fRec21[2];
  public:
	static void metadata(Meta* m) 	{ 
		m->declare("name", "ZamEQ2");
		m->declare("author", "Damien Zammit");
		m->declare("copyright", "2013");
		m->declare("version", "2.1");
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
		fslider1 = 9e+03;
		for (int i=0; i<2; i++) fRec1[i] = 0;
		iConst0 = min(192000, max(1, fSamplingFreq));
		fConst1 = (6.283185307179586 / double(iConst0));
		fslider2 = 1.0;
		for (int i=0; i<2; i++) fRec2[i] = 0;
		fslider3 = 0.0;
		for (int i=0; i<2; i++) fRec4[i] = 0;
		fslider4 = 2.5e+02;
		for (int i=0; i<2; i++) fRec5[i] = 0;
		fslider5 = 1.0;
		for (int i=0; i<2; i++) fRec6[i] = 0;
		fslider6 = 0.0;
		for (int i=0; i<2; i++) fRec8[i] = 0;
		fslider7 = 4e+03;
		for (int i=0; i<2; i++) fRec9[i] = 0;
		iConst2 = faustpower<2>(iConst0);
		fConst3 = (39.47841760435743 / double(iConst2));
		fslider8 = 0.0;
		for (int i=0; i<2; i++) fRec10[i] = 0;
		fConst4 = (9.869604401089358 / double(iConst2));
		fConst5 = (3.141592653589793 / double(iConst0));
		fConst6 = (0.5 / double(iConst0));
		fslider9 = 0.0;
		for (int i=0; i<2; i++) fRec12[i] = 0;
		fslider10 = 0.0;
		for (int i=0; i<2; i++) fRec13[i] = 0;
		fslider11 = 1e+03;
		for (int i=0; i<3; i++) fRec14[i] = 0;
		for (int i=0; i<3; i++) fRec11[i] = 0;
		for (int i=0; i<3; i++) fRec7[i] = 0;
		for (int i=0; i<3; i++) fRec3[i] = 0;
		fslider12 = 1.5e+03;
		for (int i=0; i<2; i++) fRec15[i] = 0;
		for (int i=0; i<3; i++) fRec17[i] = 0;
		for (int i=0; i<3; i++) fRec16[i] = 0;
		fcheckbox0 = 0.0;
		fslider13 = 1.2e+02;
		for (int i=0; i<2; i++) fRec18[i] = 0;
		for (int i=0; i<3; i++) fRec20[i] = 0;
		for (int i=0; i<3; i++) fRec19[i] = 0;
		fcheckbox1 = 0.0;
		fslider14 = 0.0;
		for (int i=0; i<2; i++) fRec21[i] = 0;
	}
	virtual void init(int samplingFreq) {
		classInit(samplingFreq);
		instanceInit(samplingFreq);
	}
	virtual void buildUserInterface(UI* interface) {
		interface->openVerticalBox("ZamEQ2");
		interface->addHorizontalSlider("0 Boost Lowshelf (dB)", &fslider3, 0.0, -2e+01, 2e+01, 0.1);
		interface->addHorizontalSlider("0 Freq Lowshelf (Hz)", &fslider4, 2.5e+02, 1e+01, 2e+04, 1.0);
		interface->addHorizontalSlider("0 Slope Lowshelf", &fslider5, 1.0, 1.0, 1.5, 0.1);
		interface->addHorizontalSlider("1 Boost (dB)", &fslider9, 0.0, -2e+01, 2e+01, 0.1);
		interface->addHorizontalSlider("1 Freq (Hz)", &fslider11, 1e+03, 2e+01, 2e+04, 1.0);
		interface->addHorizontalSlider("1 Q (dB)", &fslider10, 0.0, -18.0, 4e+01, 0.1);
		interface->addHorizontalSlider("2 Boost (dB)", &fslider6, 0.0, -2e+01, 2e+01, 0.1);
		interface->addHorizontalSlider("2 Freq (Hz)", &fslider7, 4e+03, 2e+01, 2e+04, 1.0);
		interface->addHorizontalSlider("2 Q (dB)", &fslider8, 0.0, -18.0, 4e+01, 0.1);
		interface->addHorizontalSlider("3 Boost Highshelf (dB)", &fslider0, 0.0, -2e+01, 2e+01, 0.1);
		interface->addHorizontalSlider("3 Freq Highshelf (Hz)", &fslider1, 9e+03, 1e+01, 2e+04, 1.0);
		interface->addHorizontalSlider("3 Slope Highshelf", &fslider2, 1.0, 1.0, 1.5, 0.1);
		interface->addHorizontalSlider("4 Freq Lowpass", &fslider12, 1.5e+03, 1e+01, 2e+04, 1.0);
		interface->addHorizontalSlider("5 Freq Highpass", &fslider13, 1.2e+02, 1e+01, 2e+04, 1.0);
		interface->addVerticalSlider("6 Master Trim (dB)", &fslider14, 0.0, -12.0, 12.0, 0.1);
		interface->addCheckButton("HP", &fcheckbox1);
		interface->addCheckButton("LP", &fcheckbox0);
		interface->closeBox();
	}
	virtual void compute (int count, FAUSTFLOAT** input, FAUSTFLOAT** output) {
		double 	fSlow0 = (0.010000000000000009 * double(fslider0));
		double 	fSlow1 = (0.010000000000000009 * double(fslider1));
		double 	fSlow2 = (0.010000000000000009 * double(fslider2));
		double 	fSlow3 = (0.010000000000000009 * double(fslider3));
		double 	fSlow4 = (0.010000000000000009 * double(fslider4));
		double 	fSlow5 = (0.010000000000000009 * double(fslider5));
		double 	fSlow6 = (0.010000000000000009 * double(fslider6));
		double 	fSlow7 = (0.010000000000000009 * double(fslider7));
		double 	fSlow8 = (0.010000000000000009 * double(fslider8));
		double 	fSlow9 = (0.010000000000000009 * double(fslider9));
		double 	fSlow10 = (0.010000000000000009 * double(fslider10));
		double 	fSlow11 = double(fslider11);
		double 	fSlow12 = faustpower<2>(fSlow11);
		double 	fSlow13 = (fConst4 * fSlow12);
		double 	fSlow14 = faustpower<2>(((fConst3 * fSlow12) - 9.869604401089358));
		double 	fSlow15 = faustpower<2>(tan((fConst5 * fSlow11)));
		double 	fSlow16 = (fConst6 * fSlow11);
		double 	fSlow17 = (2.0 * fSlow15);
		double 	fSlow18 = (0.010000000000000009 * double(fslider12));
		int 	iSlow19 = int((double(fcheckbox0) == 1.0));
		double 	fSlow20 = (0.010000000000000009 * double(fslider13));
		int 	iSlow21 = int((double(fcheckbox1) == 1.0));
		double 	fSlow22 = (0.010000000000000009 * double(fslider14));
		FAUSTFLOAT* input0 = input[0];
		FAUSTFLOAT* output0 = output[0];
		for (int i=0; i<count; i++) {
			fRec0[0] = (fSlow0 + (0.99 * fRec0[1]));
			double fTemp0 = sqrt(pow(1e+01,(0.05 * fRec0[0])));
			fRec1[0] = (fSlow1 + (0.99 * fRec1[1]));
			double fTemp1 = (fConst1 * fRec1[0]);
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
			double fTemp10 = (fConst1 * fRec5[0]);
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
			double fTemp23 = faustpower<2>(fRec9[0]);
			double fTemp24 = faustpower<2>(((fConst3 * fTemp23) - 9.869604401089358));
			fRec10[0] = (fSlow8 + (0.99 * fRec10[1]));
			double fTemp25 = exp((0.1151292546497023 * fRec10[0]));
			double fTemp26 = (fTemp22 * faustpower<2>(fTemp25));
			double fTemp27 = fabs((fTemp20 - 1.0));
			double fTemp28 = sqrt(((fTemp24 + (fConst4 * (((fTemp23 * fTemp21) * fTemp27) / fTemp26))) / ((fConst4 * ((fTemp23 * fTemp27) / fTemp26)) + fTemp24)));
			double fTemp29 = faustpower<2>(fTemp28);
			double fTemp30 = fabs((fTemp20 - fTemp29));
			double fTemp31 = (fabs((fTemp20 - fTemp28)) - sqrt((fTemp27 * fTemp30)));
			double fTemp32 = fabs((fTemp21 - 1.0));
			double fTemp33 = fabs((fTemp21 - fTemp29));
			double fTemp34 = (fabs((fTemp21 - fTemp28)) - sqrt((fTemp33 * fTemp32)));
			double fTemp35 = sqrt((fTemp33 / fTemp32));
			double fTemp36 = faustpower<2>(tan((fConst5 * fRec9[0])));
			double fTemp37 = (fTemp36 * fTemp35);
			double fTemp38 = ((fTemp30 * faustpower<2>(tan((fConst6 * (fRec9[0] / fTemp25))))) * faustpower<2>((1.0 + (fTemp37 * sqrt((fTemp27 / fTemp30))))));
			double fTemp39 = sqrt(((fTemp38 + (fTemp37 * ((2.0 * fTemp34) - (2.0 * fTemp31)))) / fTemp22));
			double fTemp40 = (1.0 + (fTemp37 + fTemp39));
			int iTemp41 = int((fTemp40 == 0.0));
			fRec12[0] = (fSlow9 + (0.99 * fRec12[1]));
			double fTemp42 = exp((0.1151292546497023 * fRec12[0]));
			int iTemp43 = int((fRec12[0] == 0.0));
			double fTemp44 = faustpower<2>(((iTemp43)?1.0:((int((fRec12[0] < 0.0)))?(1.4125375446227544 * fTemp42):(0.7079457843841379 * fTemp42))));
			double fTemp45 = faustpower<2>(fTemp42);
			double fTemp46 = fabs((fTemp45 - fTemp44));
			fRec13[0] = (fSlow10 + (0.99 * fRec13[1]));
			double fTemp47 = exp((0.1151292546497023 * fRec13[0]));
			double fTemp48 = (fTemp46 * faustpower<2>(fTemp47));
			double fTemp49 = fabs((fTemp44 - 1.0));
			double fTemp50 = sqrt(((fSlow14 + (fSlow13 * ((fTemp45 * fTemp49) / fTemp48))) / (fSlow14 + (fSlow13 * (fTemp49 / fTemp48)))));
			double fTemp51 = faustpower<2>(fTemp50);
			double fTemp52 = fabs((fTemp44 - fTemp51));
			double fTemp53 = (fabs((fTemp44 - fTemp50)) - sqrt((fTemp49 * fTemp52)));
			double fTemp54 = fabs((fTemp45 - 1.0));
			double fTemp55 = fabs((fTemp45 - fTemp51));
			double fTemp56 = (fabs((fTemp45 - fTemp50)) - sqrt((fTemp55 * fTemp54)));
			double fTemp57 = sqrt((fTemp55 / fTemp54));
			double fTemp58 = ((fTemp52 * faustpower<2>(tan((fSlow16 / fTemp47)))) * faustpower<2>((1.0 + (fSlow15 * (fTemp57 * sqrt((fTemp49 / fTemp52)))))));
			double fTemp59 = sqrt(((fTemp58 + (fSlow15 * (fTemp57 * ((2.0 * fTemp56) - (2.0 * fTemp53))))) / fTemp46));
			double fTemp60 = (fSlow15 * fTemp57);
			double fTemp61 = (1.0 + (fTemp60 + fTemp59));
			int iTemp62 = int((fTemp61 == 0.0));
			double fTemp63 = (double)input0[i];
			fRec14[0] = (0 - (((fRec14[2] * ((iTemp62)?0.0:(((1 + fTemp60) - fTemp59) / fTemp61))) + (fRec14[1] * ((iTemp62)?0.0:((0 - (2.0 * (1.0 - fTemp60))) / fTemp61)))) - fTemp63));
			double fTemp64 = sqrt((((fSlow17 * ((fTemp44 * fTemp56) * fTemp57)) + (fTemp45 * (fTemp58 - (fSlow17 * (fTemp57 * fTemp53))))) / fTemp46));
			double fTemp65 = (fTemp50 + fTemp60);
			double fTemp66 = ((iTemp43)?fTemp63:((fRec14[0] * ((iTemp62)?0.0:((fTemp64 + fTemp65) / fTemp61))) + ((fRec14[2] * ((iTemp62)?1.0:((fTemp65 - fTemp64) / fTemp61))) + (fRec14[1] * ((iTemp62)?0.0:((0 - (2.0 * (fTemp50 - fTemp60))) / fTemp61))))));
			fRec11[0] = (0 - (((fRec11[2] * ((iTemp41)?0.0:(((1 + fTemp37) - fTemp39) / fTemp40))) + (fRec11[1] * ((iTemp41)?0.0:((0 - (2.0 * (1.0 - fTemp37))) / fTemp40)))) - fTemp66));
			double fTemp67 = sqrt((((2.0 * (((fTemp20 * fTemp34) * fTemp36) * fTemp35)) + (fTemp21 * (fTemp38 - (2.0 * (fTemp37 * fTemp31))))) / fTemp22));
			double fTemp68 = (fTemp28 + fTemp37);
			double fTemp69 = ((iTemp19)?fTemp66:((fRec11[0] * ((iTemp41)?0.0:((fTemp67 + fTemp68) / fTemp40))) + ((fRec11[2] * ((iTemp41)?1.0:((fTemp68 - fTemp67) / fTemp40))) + (fRec11[1] * ((iTemp41)?0.0:((0 - (2.0 * (fTemp28 - fTemp37))) / fTemp40))))));
			fRec7[0] = (0 - (((fRec7[2] * ((iTemp17)?0.0:(((1.0 + (fTemp9 + fTemp14)) - fTemp11) / fTemp15))) + (fRec7[1] * ((iTemp17)?0.0:((0 - (2.0 * ((fTemp9 + (fTemp13 * fTemp16)) - 1.0))) / fTemp15)))) - fTemp69));
			double fTemp70 = (fTemp13 * (1.0 - fTemp9));
			double fTemp71 = ((int((fRec4[0] == 0.0)))?fTemp69:((fRec7[0] * ((iTemp17)?0.0:((fTemp9 * (1.0 + (fTemp12 + fTemp70))) / fTemp15))) + ((fRec7[2] * ((iTemp17)?1.0:((fTemp9 * ((1.0 + (fTemp9 + fTemp70)) - fTemp11)) / fTemp15))) + (fRec7[1] * ((iTemp17)?0.0:(2.0 * ((fTemp9 * ((fTemp9 + (fTemp13 * (0 - fTemp16))) - 1.0)) / fTemp15)))))));
			fRec3[0] = (0 - (((fRec3[2] * ((iTemp8)?0.0:(((1.0 + (fTemp0 + fTemp5)) - fTemp2) / fTemp6))) + (fRec3[1] * ((iTemp8)?0.0:(2.0 * (((fTemp0 + (fTemp4 * (0 - fTemp7))) - 1.0) / fTemp6))))) - fTemp71));
			double fTemp72 = (fTemp4 * (fTemp0 - 1.0));
			double fTemp73 = ((int((fRec0[0] == 0.0)))?fTemp71:((fRec3[0] * ((iTemp8)?0.0:((fTemp0 * (1.0 + (fTemp3 + fTemp72))) / fTemp6))) + ((fRec3[2] * ((iTemp8)?1.0:((fTemp0 * ((1.0 + (fTemp0 + fTemp72)) - fTemp2)) / fTemp6))) + (fRec3[1] * ((iTemp8)?0.0:((((fTemp0 + (fTemp4 * fTemp7)) - 1.0) * (0 - (2.0 * fTemp0))) / fTemp6))))));
			fRec15[0] = (fSlow18 + (0.99 * fRec15[1]));
			double fTemp74 = tan((fConst5 * fRec15[0]));
			double fTemp75 = (1.0 / fTemp74);
			double fTemp76 = (1 + ((0.7653668647301795 + fTemp75) / fTemp74));
			double fTemp77 = (1 - (1.0 / faustpower<2>(fTemp74)));
			double fTemp78 = (1 + ((1.8477590650225735 + fTemp75) / fTemp74));
			fRec17[0] = (fTemp73 - ((((1 + ((fTemp75 - 1.8477590650225735) / fTemp74)) * fRec17[2]) + (2 * (fTemp77 * fRec17[1]))) / fTemp78));
			fRec16[0] = (((fRec17[2] + (fRec17[0] + (2 * fRec17[1]))) / fTemp78) - ((((1 + ((fTemp75 - 0.7653668647301795) / fTemp74)) * fRec16[2]) + (2 * (fRec16[1] * fTemp77))) / fTemp76));
			double fTemp79 = ((iSlow19)?((fRec16[2] + (fRec16[0] + (2 * fRec16[1]))) / fTemp76):fTemp73);
			fRec18[0] = (fSlow20 + (0.99 * fRec18[1]));
			double fTemp80 = tan((fConst5 * fRec18[0]));
			double fTemp81 = (1.0 / fTemp80);
			double fTemp82 = (1 + ((0.7653668647301795 + fTemp81) / fTemp80));
			double fTemp83 = faustpower<2>(fTemp80);
			double fTemp84 = (1.0 / fTemp83);
			double fTemp85 = (1 - fTemp84);
			double fTemp86 = (1 + ((1.8477590650225735 + fTemp81) / fTemp80));
			fRec20[0] = (fTemp79 - ((((1 + ((fTemp81 - 1.8477590650225735) / fTemp80)) * fRec20[2]) + (2 * (fTemp85 * fRec20[1]))) / fTemp86));
			double fTemp87 = (0 - fTemp84);
			fRec19[0] = (((((fRec20[0] / fTemp83) + (2 * (fRec20[1] * fTemp87))) + (fRec20[2] / fTemp83)) / fTemp86) - ((((1 + ((fTemp81 - 0.7653668647301795) / fTemp80)) * fRec19[2]) + (2 * (fRec19[1] * fTemp85))) / fTemp82));
			fRec21[0] = (fSlow22 + (0.99 * fRec21[1]));
			output0[i] = (FAUSTFLOAT)(exp((0.1151292546497023 * fRec21[0])) * ((iSlow21)?((((fRec19[0] / fTemp83) + (2 * (fRec19[1] * fTemp87))) + (fRec19[2] / fTemp83)) / fTemp82):fTemp79));
			// post processing
			fRec21[1] = fRec21[0];
			fRec19[2] = fRec19[1]; fRec19[1] = fRec19[0];
			fRec20[2] = fRec20[1]; fRec20[1] = fRec20[0];
			fRec18[1] = fRec18[0];
			fRec16[2] = fRec16[1]; fRec16[1] = fRec16[0];
			fRec17[2] = fRec17[1]; fRec17[1] = fRec17[0];
			fRec15[1] = fRec15[0];
			fRec3[2] = fRec3[1]; fRec3[1] = fRec3[0];
			fRec7[2] = fRec7[1]; fRec7[1] = fRec7[0];
			fRec11[2] = fRec11[1]; fRec11[1] = fRec11[0];
			fRec14[2] = fRec14[1]; fRec14[1] = fRec14[0];
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



/***************************END USER SECTION ***************************/

/*******************BEGIN ARCHITECTURE SECTION (part 2/2)***************/

//-----------------------------------portCollector--------------------------------------
//
// portCollector is passed to the buildUserInterface method of a dsp object
// in order to build a description of its inputs, outputs and control ports.
// This description is used to fill a LADSPA_Descriptor
//
//--------------------------------------------------------------------------------------

//--------------------------------useful constants--------------------------------------

#define MAXPORT 1024
static const int ICONTROL 	= LADSPA_PORT_INPUT|LADSPA_PORT_CONTROL;
static const int OCONTROL 	= LADSPA_PORT_OUTPUT|LADSPA_PORT_CONTROL;
static const int RANGE 		= LADSPA_PORT_INPUT|LADSPA_PORT_CONTROL;

static const char* inames[] = {
					"input00", "input01", "input02", "input03", "input04",
					"input05", "input06", "input07", "input08", "input09",
					"input10", "input11", "input12", "input13", "input14",
					"input15", "input16", "input17", "input18", "input19",
					"input20", "input21", "input22", "input23", "input24",
					"input25", "input26", "input27", "input28", "input29",
					"input30", "input31", "input32", "input33", "input34",
					"input35", "input36", "input37", "input38", "input39"
};

static const char* onames[] = {
					"output00", "output01", "output02", "output03", "output04",
					"output05", "output06", "output07", "output08", "output09",
					"output10", "output11", "output12", "output13", "output14",
					"output15", "output16", "output17", "output18", "output19",
					"output20", "output21", "output22", "output23", "output24",
					"output25", "output26", "output27", "output28", "output29",
					"output30", "output31", "output32", "output33", "output34",
					"output35", "output36", "output37", "output38", "output39"
};

class portCollector : public UI
{
 private:

	//--------------------------------------------------------------------------------------

	const int				fInsCount;					// number of audio input ports
	const int				fOutsCount;					// number of audio output ports
	int						fCtrlCount;					// number of control ports

	LADSPA_PortDescriptor 	fPortDescs[MAXPORT];		// table of port descriptors to be used in a LADSPA_Descriptor
	const char* 			fPortNames[MAXPORT];		// table of port names to be used in a LADSPA_Descriptor
	LADSPA_PortRangeHint 	fPortHints[MAXPORT];		// table of port hints to be used in a LADSPA_Descriptor

    std::string					fPluginName;				// toplevel prefix used as plugin name
    std::stack<std::string>			fPrefix;					// current prefix for controls name


	//--------------------------------------------------------------------------------------
    std::string simplify(const std::string& src)
	{
		int		i=0;
		int		level=2;
        std::string	dst;

		while (src[i] ) {

			switch (level) {

				case 0 :
				case 1 :
				case 2 :
					// Skip the begin of the label "--foo-"
					// until 3 '-' have been read
					if (src[i]=='-') { level++; }
					break;

				case 3 :
					// copy the content, but skip non alphnum
					// and content in parenthesis
					switch (src[i]) {
						case '(' :
						case '[' :
							level++;
							break;

						case '-' :
							dst += '-';
							break;

						default :
							if (isalnum(src[i])) {
								dst+= tolower(src[i]);
							}

					}
					break;

				default :
					// here we are inside parenthesis and
					// we skip the content until we are back to
					// level 3
					switch (src[i]) {

						case '(' :
						case '[' :
							level++;
							break;

						case ')' :
						case ']' :
							level--;
							break;

						default :
							break;
					}

			}
			i++;
		}
		return (dst.size() > 0) ? dst :src;
	}

	void addPortDescr(int type, const char* label, int hint, float min=0.0, float max=0.0)
	{
        std::string fullname = simplify(fPrefix.top() + "-" + label);
		char * str = strdup(fullname.c_str());

		fPortDescs[fInsCount + fOutsCount + fCtrlCount] = type;
		fPortNames[fInsCount + fOutsCount + fCtrlCount] = str;
		fPortHints[fInsCount + fOutsCount + fCtrlCount].HintDescriptor = hint;
		fPortHints[fInsCount + fOutsCount + fCtrlCount].LowerBound = min;
		fPortHints[fInsCount + fOutsCount + fCtrlCount].UpperBound = max;
		fCtrlCount++;
	}

	void openAnyBox(const char* label)
	{
		if (fPrefix.size() == 0) {
			// top level label is used as plugin name
			fPluginName = label;
			fPrefix.push(label);

		} else {
            std::string s;
			if (label && label[0]) {
				s = fPrefix.top() + "-" + label;
			} else {
				s = fPrefix.top();
			}
			fPrefix.push(s);
		}
	}

 public:

	//--------------------------------Collect the audio ports-------------------------------

	portCollector(int ins, int outs) : UI(), fInsCount(ins), fOutsCount(outs), fCtrlCount(0)
	{
		for (int i = 0; i < ins; i++) {
			fPortDescs[i] = LADSPA_PORT_INPUT | LADSPA_PORT_AUDIO;
			fPortNames[i] = inames[i];
			fPortHints[i].HintDescriptor = 0;
		}
		for (int j = 0; j < outs; j++) {
			fPortDescs[ins + j] = LADSPA_PORT_OUTPUT | LADSPA_PORT_AUDIO;
			fPortNames[ins + j] = onames[j];
			fPortHints[ins + j].HintDescriptor = 0;
		}
	};

	virtual ~portCollector() {}

	//------------------------------Collect the control ports-------------------------------

	virtual void addButton(const char* label, float* zone) {
		addPortDescr(ICONTROL, label, LADSPA_HINT_TOGGLED);
	}

	virtual void addToggleButton(const char* label, float* zone) {
		addPortDescr(ICONTROL, label, LADSPA_HINT_TOGGLED);
	}

	virtual void addCheckButton(const char* label, float* zone) {
		addPortDescr(ICONTROL, label, LADSPA_HINT_TOGGLED);
	}

	virtual void addVerticalSlider(const char* label, float* zone, float init, float min, float max, float step) {
		addPortDescr(ICONTROL, label, LADSPA_HINT_BOUNDED_BELOW | LADSPA_HINT_BOUNDED_ABOVE, min, max);
	}

	virtual void addHorizontalSlider(const char* label, float* zone, float init, float min, float max, float step) {
		(max < 18000.0) ? addPortDescr(ICONTROL, label, LADSPA_HINT_BOUNDED_BELOW | LADSPA_HINT_BOUNDED_ABOVE, min, max) : addPortDescr(ICONTROL, label, LADSPA_HINT_BOUNDED_BELOW | LADSPA_HINT_BOUNDED_ABOVE | LADSPA_HINT_LOGARITHMIC, min, max);
	}

	virtual void addNumEntry(const char* label, float* zone, float init, float min, float max, float step) {
		addPortDescr(ICONTROL, label, LADSPA_HINT_BOUNDED_BELOW | LADSPA_HINT_BOUNDED_ABOVE, min, max);
	}

	// -- passive widgets

	virtual void addNumDisplay(const char* label, float* zone, int precision) {
		addPortDescr(OCONTROL, label, 0, -10000, +10000);
	}
	virtual void addTextDisplay(const char* label, float* zone, const char* names[], float min, float max) {
		addPortDescr(OCONTROL, label, LADSPA_HINT_BOUNDED_BELOW | LADSPA_HINT_BOUNDED_ABOVE, min, max);
	}
	virtual void addHorizontalBargraph(const char* label, float* zone, float min, float max) {
		addPortDescr(OCONTROL, label, LADSPA_HINT_BOUNDED_BELOW | LADSPA_HINT_BOUNDED_ABOVE, min, max);
	}
	virtual void addVerticalBargraph(const char* label, float* zone, float min, float max){
		addPortDescr(OCONTROL, label, LADSPA_HINT_BOUNDED_BELOW | LADSPA_HINT_BOUNDED_ABOVE, min, max);
	}

	virtual void openFrameBox(const char* label)		{ openAnyBox(label); }
	virtual void openTabBox(const char* label)		{ openAnyBox(label); }
	virtual void openHorizontalBox(const char* label)	{ openAnyBox(label); }
	virtual void openVerticalBox(const char* label)	{ openAnyBox(label); }

	virtual void closeBox() 					{ fPrefix.pop(); }

	virtual void show() {}
	virtual void run() 	{}

	//---------------------------------Fill the LADSPA descriptor---------------------------

	// generate an ID from a plugin name
	int makeID (const char* s) {
		int h = 0;
		for (int i = 0; s[i]; i++) {
			h = (h << 3) + (s[i] & 7);
		}
		return 1+h%1000;
	}

	// fill a ladspa descriptor with the information collected on ports
	void fillPortDescription (LADSPA_Descriptor * descriptor) {
		const char* name = sym(mydsp);
		descriptor->PortCount 			= fCtrlCount+fInsCount+fOutsCount;
		descriptor->PortDescriptors 	= fPortDescs;
		descriptor->PortNames 			= fPortNames;
		descriptor->PortRangeHints 		= fPortHints;

		descriptor->Label = strdup(name);
		descriptor->UniqueID = makeID(name);
//		descriptor->Label = strdup(fPluginName.c_str());
//		descriptor->UniqueID = makeID(fPluginName.c_str());
		descriptor->Properties = LADSPA_PROPERTY_HARD_RT_CAPABLE;
		descriptor->Name = "ZamEQ2";
//		descriptor->Name = strdup(fPluginName.c_str());
		descriptor->Maker = "Damien Zammit";
		descriptor->Copyright = "GPL";
	}
};

//--------------------------------------portData----------------------------------------
//
// portData : a user interface used to associate the data buffers and the ports
//
//--------------------------------------------------------------------------------------

class portData : public UI
{

 private:

	//--------------------------------------------------------------------------------------

	const int				fInsCount;					// number of audio input ports
	const int				fOutsCount;					// number of audio output ports
	int						fCtrlCount;					// number of control ports

	float* 					fPortZone[MAXPORT];			//
	float* 					fPortData[MAXPORT];

	//--------------------------------------------------------------------------------------

	void addZone(float* zone)
	{
		fPortZone[fInsCount + fOutsCount + fCtrlCount] = zone;
		fCtrlCount++;
	}

 public:

	//--------------------------------Collect the audio ports-------------------------------

	portData(int ins, int outs) : UI(), fInsCount(ins), fOutsCount(outs), fCtrlCount(0) {};
	virtual ~portData() {}

	//------------------------------Collect the control zones-------------------------------

	virtual void addButton(const char* label, float* zone) 			{ addZone(zone); }
	virtual void addToggleButton(const char* label, float* zone)  	{ addZone(zone); }
	virtual void addCheckButton(const char* label, float* zone)  		{ addZone(zone); }

	virtual void addVerticalSlider(const char* label, float* zone, float init, float min, float max, float step) 		{ addZone(zone); }
	virtual void addHorizontalSlider(const char* label, float* zone, float init, float min, float max, float step) 	{ addZone(zone); }
	virtual void addNumEntry(const char* label, float* zone, float init, float min, float max, float step)  			{ addZone(zone); }

	// -- passive widgets

	virtual void addNumDisplay(const char* label, float* zone, int precision) 						{ addZone(zone); }
	virtual void addTextDisplay(const char* label, float* zone, const char* names[], float min, float max) 	{ addZone(zone); }
	virtual void addHorizontalBargraph(const char* label, float* zone, float min, float max) 			{ addZone(zone); }
	virtual void addVerticalBargraph(const char* label, float* zone, float min, float max)			{ addZone(zone); }

	virtual void openFrameBox(const char* label)		{ }
	virtual void openTabBox(const char* label)		{ }
	virtual void openHorizontalBox(const char* label)	{ }
	virtual void openVerticalBox(const char* label)	{ }
	virtual void closeBox() 					{ }

	virtual void show() {}
	virtual void run() 	{}

	//---------------------------------interaction with LADSPA------------------------------

	void setPortData (unsigned long port, LADSPA_Data* data) {
		fPortData[port] = data;
	}

	void updateCtrlZones() {
		for (int i = fInsCount+fOutsCount; i < fInsCount+fOutsCount+fCtrlCount; i++)	*fPortZone[i] = *fPortData[i];
	}

	float** getInputs() {
		return &fPortData[0];
	}

	float** getOutputs() {
		return &fPortData[fInsCount];
	}
};

//--------------------------------Faust-LADSPA plugin-----------------------------------
//
// Plugin structure, callbacks and LADSPA_descriptor(i) entry point
//
//--------------------------------------------------------------------------------------

LADSPA_Descriptor* 	gDescriptor = 0;

struct PLUGIN
{
	unsigned long	fSampleRate;
	portData*		fPortData;
	dsp*			fDsp;

	PLUGIN(unsigned long r, portData* d, dsp* p) : fSampleRate(r), fPortData(d), fDsp(p) {}
};

LADSPA_Handle instantiate_method (const struct _LADSPA_Descriptor * Descriptor, unsigned long SampleRate)
{
	dsp*		p = new mydsp();
	portData* 	d = new portData(p->getNumInputs(), p->getNumOutputs());

	p->buildUserInterface(d);
	return new PLUGIN (SampleRate, d, p);
}

void connect_method (LADSPA_Handle Instance, unsigned long Port, LADSPA_Data * DataLocation)
{
	PLUGIN* p = (PLUGIN*) Instance;
	p->fPortData->setPortData(Port, DataLocation);
}

void activate_method (LADSPA_Handle Instance)
{
	PLUGIN* p = (PLUGIN*) Instance;
	p->fDsp->init(p->fSampleRate);
}

void run_method (LADSPA_Handle Instance, unsigned long SampleCount)
{
	PLUGIN* p = (PLUGIN*) Instance;
	p->fPortData->updateCtrlZones();
	AVOIDDENORMALS;
	p->fDsp->compute(SampleCount, p->fPortData->getInputs(), p->fPortData->getOutputs());
}

void deactivate_method (LADSPA_Handle Instance)
{}

void cleanup_method (LADSPA_Handle Instance)
{
	PLUGIN* p = (PLUGIN*) Instance;
	delete p->fPortData;
	delete p->fDsp;
	delete p;
}

//--------------------------------------------------------------------------------------

void init_descriptor(LADSPA_Descriptor* descriptor)
{
	descriptor->UniqueID = 123456;
	descriptor->Label = "none";
	descriptor->Properties = LADSPA_PROPERTY_HARD_RT_CAPABLE;
	descriptor->Name = "none";
	descriptor->Maker = "Yann Orlarey";
	descriptor->Copyright = "GPL";

	descriptor->ImplementationData = 0;

	// description des methods
	descriptor->instantiate = instantiate_method;
	descriptor->connect_port = connect_method;
	descriptor->activate = activate_method;
	descriptor->run = run_method;
	descriptor->run_adding = 0;
	descriptor->set_run_adding_gain = 0;
	descriptor->deactivate = deactivate_method;
	descriptor->cleanup = cleanup_method;
}

//--------------------------------------------------------------------------------------

const LADSPA_Descriptor * ladspa_descriptor(unsigned long Index)
{
    if (Index == 0) {
		if (gDescriptor == 0)
		{
			// allocate temporaries dsp and portCollector to build the plugin description
			mydsp* p = new mydsp();
			if (p) {
				portCollector*	c=new portCollector(p->getNumInputs(), p->getNumOutputs());
				p->buildUserInterface(c);
				gDescriptor = new LADSPA_Descriptor;
				init_descriptor(gDescriptor);
				c->fillPortDescription(gDescriptor);
				delete p;
			} else {
				printf("Memory Error : unable to allocate the dsp object\n");
			}
		}
		return gDescriptor;
	} else {
		return NULL;
	}
}

/********************END ARCHITECTURE SECTION (part 2/2)****************/


