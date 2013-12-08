#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "lv2/lv2plug.in/ns/lv2core/lv2.h"

#define ZAMEQ2_URI "http://zamaudio.com/lv2/zameq2"

typedef enum {
	ZAMEQ2_INPUT  = 0,
	ZAMEQ2_OUTPUT = 1,

	ZAMEQ2_BOOSTDB1 = 2,
	ZAMEQ2_Q1 = 3,
	ZAMEQ2_FREQ1 = 4,
	
	ZAMEQ2_BOOSTDB2 = 5,
	ZAMEQ2_Q2 = 6,
	ZAMEQ2_FREQ2 = 7,
	
	ZAMEQ2_BOOSTDBL = 8,
	ZAMEQ2_SLOPEDBL = 9,
	ZAMEQ2_FREQL = 10,

	ZAMEQ2_BOOSTDBH = 11,
	ZAMEQ2_SLOPEDBH = 12,
	ZAMEQ2_FREQH = 13
} PortIndex;

typedef struct {
	float* input;
	float* output;

	float* boostdb1;
	float* q1;
	float* freq1;

	float* boostdb2;
	float* q2;
	float* freq2;

	float* boostdbl;
	float* slopedbl;
	float* freql;

	float* boostdbh;
	float* slopedbh;
	float* freqh;

	double x1,x2,y1,y2;
	double x1a,x2a,y1a,y2a;
	double zln1,zln2,zld1,zld2;
	double zhn1,zhn2,zhd1,zhd2;
	double a0x,a1x,a2x,b0x,b1x,b2x,gainx;
	double a0y,a1y,a2y,b0y,b1y,b2y,gainy;
	double Bl[3];
	double Al[3];
	double Bh[3];
	double Ah[3];
	float srate;
} ZamEQ2;



static LV2_Handle
instantiate(const LV2_Descriptor*     descriptor,
            double                    rate,
            const char*               bundle_path,
            const LV2_Feature* const* features)
{
	int i;
	ZamEQ2* zameq2 = (ZamEQ2*)malloc(sizeof(ZamEQ2));
	zameq2->srate = rate;
	zameq2->x1=zameq2->x2=zameq2->y1=zameq2->y2=0.f;
	zameq2->x1a=zameq2->x2a=zameq2->y1a=zameq2->y2a=0.f;
	zameq2->a0x=zameq2->a1x=zameq2->a2x=zameq2->b0x=zameq2->b1x=zameq2->b2x=zameq2->gainx=0.f;
	zameq2->a0y=zameq2->a1y=zameq2->a2y=zameq2->b0y=zameq2->b1y=zameq2->b2y=zameq2->gainy=0.f;
	zameq2->zln1=zameq2->zln2=zameq2->zld1=zameq2->zld2=0.f;
	zameq2->zhn1=zameq2->zhn2=zameq2->zhd1=zameq2->zhd2=0.f;

	for (i = 0; i < 3; ++i) {
		zameq2->Bl[i] = zameq2->Al[i] = zameq2->Bh[i] = zameq2->Ah[i] = 0.f;
		zameq2->Bl[i] = zameq2->Al[i] = zameq2->Bh[i] = zameq2->Ah[i] = 0.f;
		zameq2->Bl[i] = zameq2->Al[i] = zameq2->Bh[i] = zameq2->Ah[i] = 0.f;
	}

	return (LV2_Handle)zameq2;
}

static void
connect_port(LV2_Handle instance,
             uint32_t   port,
             void*      data)
{
	ZamEQ2* zameq2 = (ZamEQ2*)instance;

	switch ((PortIndex)port) {
	case ZAMEQ2_INPUT:
		zameq2->input = (float*)data;
		break;
	case ZAMEQ2_OUTPUT:
		zameq2->output = (float*)data;
		break;
	case ZAMEQ2_BOOSTDB1:
		zameq2->boostdb1 = (float*)data;
		break;
	case ZAMEQ2_Q1:
		zameq2->q1 = (float*)data;
		break;
	case ZAMEQ2_FREQ1:
		zameq2->freq1 = (float*)data;
		break;
	case ZAMEQ2_BOOSTDB2:
		zameq2->boostdb2 = (float*)data;
		break;
	case ZAMEQ2_Q2:
		zameq2->q2 = (float*)data;
		break;
	case ZAMEQ2_FREQ2:
		zameq2->freq2 = (float*)data;
		break;
	case ZAMEQ2_BOOSTDBL:
		zameq2->boostdbl = (float*)data;
		break;
	case ZAMEQ2_SLOPEDBL:
		zameq2->slopedbl = (float*)data;
		break;
	case ZAMEQ2_FREQL:
		zameq2->freql = (float*)data;
		break;
	case ZAMEQ2_BOOSTDBH:
		zameq2->boostdbh = (float*)data;
		break;
	case ZAMEQ2_SLOPEDBH:
		zameq2->slopedbh = (float*)data;
		break;
	case ZAMEQ2_FREQH:
		zameq2->freqh = (float*)data;
		break;
	}
}

// Force already-denormal float value to zero
static inline void 
sanitize_denormal(double& value) {
    if (!isnormal(value)) {
        value = 0.f;
    }
}

static inline int 
sign(double x) {
        return (x >= 0.f ? 1 : -1);
}

static inline double 
from_dB(double gdb) {
        return (exp(gdb/20.f*log(10.f)));
}

static inline double
to_dB(double g) {
        return (20.f*log10(g));
}

static void
activate(LV2_Handle instance)
{
}

static void
peq(double G0, double G, double GB, double w0, double Dw,
        double *a0, double *a1, double *a2, double *b0, double *b1, double *b2, double *gn) {

        double F,G00,F00,num,den,G1,G01,G11,F01,F11,W2,Dww,C,D,B,A;
        F = fabs(G*G - GB*GB);
        G00 = fabs(G*G - G0*G0);
        F00 = fabs(GB*GB - G0*G0);
        num = G0*G0 * (w0*w0 - M_PI*M_PI)*(w0*w0 - M_PI*M_PI)
                + G*G * F00 * M_PI*M_PI * Dw*Dw / F;
        den = (w0*w0 - M_PI*M_PI)*(w0*w0 - M_PI*M_PI)
                + F00 * M_PI*M_PI * Dw*Dw / F;
        G1 = sqrt(num/den);
        G01 = fabs(G*G - G0*G1);
        G11 = fabs(G*G - G1*G1);
        F01 = fabs(GB*GB - G0*G1);
        F11 = fabs(GB*GB - G1*G1);
        W2 = sqrt(G11 / G00) * tan(w0/2.f)*tan(w0/2.f);
        Dww = (1.f + sqrt(F00 / F11) * W2) * tan(Dw/2.f);
        C = F11 * Dww*Dww - 2.f * W2 * (F01 - sqrt(F00 * F11));
        D = 2.f * W2 * (G01 - sqrt(G00 * G11));
        A = sqrt((C + D) / F);
        B = sqrt((G*G * C + GB*GB * D) / F);
        *gn = G1;
        *b0 = (G1 + G0*W2 + B) / (1.f + W2 + A);
        *b1 = -2.f*(G1 - G0*W2) / (1.f + W2 + A);
        *b2 = (G1 - B + G0*W2) / (1.f + W2 + A);
        *a0 = 1.f;
        *a1 = -2.f*(1.f - W2) / (1.f + W2 + A);
        *a2 = (1 + W2 - A) / (1.f + W2 + A);

        sanitize_denormal(*b1);
        sanitize_denormal(*b2);
        sanitize_denormal(*a0);
        sanitize_denormal(*a1);
        sanitize_denormal(*a2);
        sanitize_denormal(*gn);
        if (!isnormal(*b0)) { *b0 = 1.f; }
}

static bool
lowshelfeq(double G0, double G, double GB, double w0, double Dw, double q,
		double B[], double A[]) {
 	double alpha,b0,b1,b2,a0,a1,a2;
	G = powf(10.f,G/20.f); 
	double AA  = sqrt(G);
	
	alpha = sin(w0)/2.f * sqrt( (AA + 1.f/AA)*(1.f/q - 1.f) + 2.f );
	b0 =    AA*( (AA+1.f) - (AA-1.f)*cos(w0) + 2.f*sqrt(AA)*alpha );
        b1 =  2.f*AA*( (AA-1.f) - (AA+1.f)*cos(w0)                   );
        b2 =    AA*( (AA+1.f) - (AA-1.f)*cos(w0) - 2.f*sqrt(AA)*alpha );
        a0 =        (AA+1.f) + (AA-1.f)*cos(w0) + 2.f*sqrt(AA)*alpha;
        a1 =   -2.f*( (AA-1.f) + (AA+1.f)*cos(w0)                   );
        a2 =        (AA+1.f) + (AA-1.f)*cos(w0) - 2.f*sqrt(AA)*alpha;
	
	B[0] = b0/a0;
        B[1] = b1/a0;
        B[2] = b2/a0;
        A[0] = 1.f;
        A[1] = a1/a0;
        A[2] = a2/a0;

	return true;
}

static bool
highshelfeq(double G0, double G, double GB, double w0, double Dw, double q,
		double B[], double A[]) {
        double alpha,b0,b1,b2,a0,a1,a2;
        G = powf(10.f,G/20.f);
        double AA  = sqrt(G);

        alpha = sin(w0)/2.f * sqrt( (AA + 1.f/AA)*(1.f/q - 1.f) + 2.f );
        b0 =    AA*( (AA+1.f) + (AA-1.f)*cos(w0) + 2.f*sqrt(AA)*alpha );
        b1 =  -2.f*AA*( (AA-1.f) + (AA+1.f)*cos(w0)                   );
        b2 =    AA*( (AA+1.f) + (AA-1.f)*cos(w0) - 2.f*sqrt(AA)*alpha );
        a0 =        (AA+1.f) - (AA-1.f)*cos(w0) + 2.f*sqrt(AA)*alpha;
        a1 =   2.f*( (AA-1.f) - (AA+1.f)*cos(w0)                   );
        a2 =        (AA+1.f) - (AA-1.f)*cos(w0) - 2.f*sqrt(AA)*alpha;

        B[0] = b0/a0;
        B[1] = b1/a0;
        B[2] = b2/a0;
        A[0] = 1.f;
        A[1] = a1/a0;
        A[2] = a2/a0;

        return true;
}

static void
run(LV2_Handle instance, uint32_t n_samples)
{
	ZamEQ2* zameq2 = (ZamEQ2*)instance;

	const float* const input  = zameq2->input;
	float* const       output = zameq2->output;

	const double        boostdb1 = *(zameq2->boostdb1);
	const double        q1 = *(zameq2->q1);
	const double        freq1 = *(zameq2->freq1);
	
	const double        boostdb2 = *(zameq2->boostdb2);
	const double        q2 = *(zameq2->q2);
	const double        freq2 = *(zameq2->freq2);
	
	const double        boostdbl = *(zameq2->boostdbl);
	const double        slopedbl = *(zameq2->slopedbl);
	const double        freql = *(zameq2->freql);

	const double        boostdbh = *(zameq2->boostdbh);
	const double        slopedbh = *(zameq2->slopedbh);
	const double        freqh = *(zameq2->freqh);

	double dcgain = 1.f;
	
	double boost1 = from_dB(boostdb1);
  	double fc1 = freq1 / zameq2->srate;
	double w01 = fc1*2.f*M_PI;
	double bwgain1 = (boostdb1 == 0.f) ? 1.f : (boostdb1 < 0.f) ? boost1*from_dB(3.f) : boost1*from_dB(-3.f);
	double bw1 = fc1 / q1;

	double boost2 = from_dB(boostdb2);
  	double fc2 = freq2 / zameq2->srate;
	double w02 = fc2*2.f*M_PI;
	double bwgain2 = (boostdb2 == 0.f) ? 1.f : (boostdb2 < 0.f) ? boost2*from_dB(3.f) : boost2*from_dB(-3.f);
	double bw2 = fc2 / q2;

	double boostl = from_dB(boostdbl);
	double All = sqrt(boostl);
	double bwl = 2.f*M_PI*freql/ zameq2->srate;
	double bwgaindbl = to_dB(All);
	
	double boosth = from_dB(boostdbh);
	double Ahh = sqrt(boosth);
	double bwh = 2.f*M_PI*freqh/ zameq2->srate;
	double bwgaindbh = to_dB(Ahh);

	peq(dcgain,boost1,bwgain1,w01,bw1,&zameq2->a0x,&zameq2->a1x,&zameq2->a2x,&zameq2->b0x,&zameq2->b1x,&zameq2->b2x,&zameq2->gainx);
	peq(dcgain,boost2,bwgain2,w02,bw2,&zameq2->a0y,&zameq2->a1y,&zameq2->a2y,&zameq2->b0y,&zameq2->b1y,&zameq2->b2y,&zameq2->gainy);
	lowshelfeq(0.f,boostdbl,bwgaindbl,2.f*M_PI*freql/zameq2->srate,bwl,slopedbl,zameq2->Bl,zameq2->Al);
	highshelfeq(0.f,boostdbh,bwgaindbh,2.f*M_PI*freqh/zameq2->srate,bwh,slopedbh,zameq2->Bh,zameq2->Ah);

	for (uint32_t pos = 0; pos < n_samples; pos++) {
		double tmp,tmpl, tmph;
		double in = input[pos];
		sanitize_denormal(zameq2->x1);
		sanitize_denormal(zameq2->x2);
		sanitize_denormal(zameq2->y1);
		sanitize_denormal(zameq2->y2);
		sanitize_denormal(zameq2->x1a);
		sanitize_denormal(zameq2->x2a);
		sanitize_denormal(zameq2->y1a);
		sanitize_denormal(zameq2->y2a);
		sanitize_denormal(zameq2->zln1);
		sanitize_denormal(zameq2->zln2);
		sanitize_denormal(zameq2->zld1);
		sanitize_denormal(zameq2->zld2);
		sanitize_denormal(zameq2->zhn1);
		sanitize_denormal(zameq2->zhn2);
		sanitize_denormal(zameq2->zhd1);
		sanitize_denormal(zameq2->zhd2);
		sanitize_denormal(in);

		//lowshelf
		tmpl = in * zameq2->Bl[0] + 
			zameq2->zln1 * zameq2->Bl[1] +
			zameq2->zln2 * zameq2->Bl[2] -
			zameq2->zld1 * zameq2->Al[1] -
			zameq2->zld2 * zameq2->Al[2];
		zameq2->zln2 = zameq2->zln1;
		zameq2->zld2 = zameq2->zld1;
		zameq2->zln1 = in;
		zameq2->zld1 = tmpl;
		
		//highshelf
		tmph = tmpl * zameq2->Bh[0] + 
			zameq2->zhn1 * zameq2->Bh[1] +
			zameq2->zhn2 * zameq2->Bh[2] -
			zameq2->zhd1 * zameq2->Ah[1] -
			zameq2->zhd2 * zameq2->Ah[2];
		zameq2->zhn2 = zameq2->zhn1;
		zameq2->zhd2 = zameq2->zhd1;
		zameq2->zhn1 = tmpl;
		zameq2->zhd1 = tmph;
		
		//parametric1
		tmp = tmph * zameq2->b0x + zameq2->x1 * zameq2->b1x + zameq2->x2 * zameq2->b2x - zameq2->y1 * zameq2->a1x - zameq2->y2 * zameq2->a2x;
		zameq2->x2 = zameq2->x1;
		zameq2->y2 = zameq2->y1;
		zameq2->x1 = tmph;
		zameq2->y1 = tmp;
		
		//parametric2
		output[pos] = tmp * zameq2->b0y + zameq2->x1a * zameq2->b1y + zameq2->x2a * zameq2->b2y - zameq2->y1a * zameq2->a1y - zameq2->y2a * zameq2->a2y;
		zameq2->x2a = zameq2->x1a;
		zameq2->y2a = zameq2->y1a;
		zameq2->x1a = tmp;
		zameq2->y1a = output[pos];
	}
}

static void
deactivate(LV2_Handle instance)
{
}

static void
cleanup(LV2_Handle instance)
{
	free(instance);
}

const void*
extension_data(const char* uri)
{
	return NULL;
}

static const LV2_Descriptor descriptor = {
	ZAMEQ2_URI,
	instantiate,
	connect_port,
	activate,
	run,
	deactivate,
	cleanup,
	extension_data
};

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
