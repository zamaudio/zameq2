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
		interface->openVerticalBox("ZamEQ2");
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

	string					fPluginName;				// toplevel prefix used as plugin name
	stack<string>			fPrefix;					// current prefix for controls name


	//--------------------------------------------------------------------------------------
	string simplify(const string& src)
	{
		int		i=0;
		int		level=2;
		string	dst;

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
		string fullname = simplify(fPrefix.top() + "-" + label);
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
			string s;
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
		addPortDescr(ICONTROL, label, LADSPA_HINT_BOUNDED_BELOW | LADSPA_HINT_BOUNDED_ABOVE, min, max);
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

		descriptor->Label = "ZamEQ2";
		descriptor->UniqueID = makeID(name);
//		descriptor->Label = strdup(fPluginName.c_str());
//		descriptor->UniqueID = makeID(fPluginName.c_str());
		descriptor->Properties = LADSPA_PROPERTY_HARD_RT_CAPABLE;
		descriptor->Name = "ZamEQ2";
//		descriptor->Name = strdup(fPluginName.c_str());
		descriptor->Maker = "Damien Zammit";
		descriptor->Copyright = "2013";
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


