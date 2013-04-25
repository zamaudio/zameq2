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
	ZAMEQ2_FREQL = 10
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

	float x1,x2,y1,y2;
	float x1a,x2a,y1a,y2a;
	float zln1a,zln2a,zld1a,zld2a;
	float zln1b,zln2b,zld1b,zld2b;
	float a0x,a1x,a2x,b0x,b1x,b2x,gainx;
	float a0y,a1y,a2y,b0y,b1y,b2y,gainy;
	float B[3][5];
	float A[3][5];
	float Bh[3][5];
	float Ah[3][5];
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
	zameq2->zln1a=zameq2->zln2a=zameq2->zld1a=zameq2->zld2a=0.f;
	zameq2->zln1b=zameq2->zln2b=zameq2->zld1b=zameq2->zld2b=0.f;	
	zameq2->a0x=zameq2->a1x=zameq2->a2x=zameq2->b0x=zameq2->b1x=zameq2->b2x=zameq2->gainx=0.f;
	zameq2->a0y=zameq2->a1y=zameq2->a2y=zameq2->b0y=zameq2->b1y=zameq2->b2y=zameq2->gainy=0.f;
	for (i = 0; i < 5; ++i) {
		zameq2->B[0][i] = zameq2->A[0][i] = zameq2->Bh[0][i] = zameq2->Ah[0][i] = 0.f;
		zameq2->B[1][i] = zameq2->A[1][i] = zameq2->Bh[1][i] = zameq2->Ah[1][i] = 0.f;
		zameq2->B[2][i] = zameq2->A[2][i] = zameq2->Bh[2][i] = zameq2->Ah[2][i] = 0.f;
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
	}
}

// Works on little-endian machines only
static inline bool 
is_nan(float& value ) {
    if (((*(uint32_t *) &value) & 0x7fffffff) > 0x7f800000) {
      return true;
    }
    return false;
}

// Force already-denormal float value to zero
static inline void 
sanitize_denormal(float& value) {
    if (is_nan(value)) {
        value = 0.f;
    }
}

static inline int 
sign(float x) {
        return (x >= 0.f ? 1 : -1);
}

static inline float 
from_dB(float gdb) {
        return (exp(gdb/20.f*log(10.f)));
}

static inline float
to_dB(float g) {
        return (20.f*log10(g));
}

static float
arcsinh (float x) {
	return (log(2.f*x+sqrt(x*x+1.f)));
}

static void
activate(LV2_Handle instance)
{
}

static void
peq(float G0, float G, float GB, float w0, float Dw,
        float *a0, float *a1, float *a2, float *b0, float *b1, float *b2, float *gn) {

        float F,G00,F00,num,den,G1,G01,G11,F01,F11,W2,Dww,C,D,B,A;
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
        if (is_nan(*b0)) { *b0 = 1.f; }
}

static bool
shelfeq(float G0, float G, float GB, float w0, float Dw, float q,
		float B[][5], float A[][5], float Bh[][5], float Ah[][5]) {
	float r,L,c0,WB,e,g,g0,b,D,phi,si,b0h,b1h,b2h,a1h,a2h,b0,b1,b2,b3,b4,a0,a1,a2,a3,a4,alpha;
	int N = 2;
	int i;

	G0 = powf(10.f,G0/20.f);
	G = powf(10.f,G/20.f); 
	GB = powf(10.f,GB/20.f); 

	r = N % 2; L = (N-r)/2.f;

	Bh[0][0] = 1.f;
	Bh[0][1] = 0.f;
	Bh[0][2] = 0.f;

	Ah[0][0] = 1.f;
	Ah[0][1] = 0.f;
	Ah[0][2] = 0.f;
		
	B[0][0] = 1.f;
	B[0][1] = 0.f;
	B[0][2] = 0.f;
	B[0][3] = 0.f;
	B[0][4] = 0.f;
	
	A[0][0] = 1.f;
	A[0][1] = 0.f;
	A[0][2] = 0.f;
	A[0][3] = 0.f;
	A[0][4] = 0.f;

	//if (G==G0) 
	//	return false; 

	c0 = cos(w0); 

	if (w0==0.f) 
		c0 = 1.f;

	if (w0==M_PI/2.f)
		c0 = 0.f;

	if (w0==M_PI)
		c0 = -1.f;

	//WB = tan(Dw/2.f);
	WB = tan(Dw/(q*q));
	e = sqrt((G*G - GB*GB)/(GB*GB - G0*G0)); 
	g = powf(G,1.f/N); 
	g0 = powf(G0,1.f/N); 

	//Butterworth
	b = WB / (powf(e,1.f/N));
  
	if (r==1) {                                       
		D = b + 1.f;
		Bh[0][0] = (g*b+g0)/D;
		Bh[0][1] = (g*b-g0)/D;
		Bh[0][2] = 0.f;

    		Ah[0][0] = 1.f;
    		Ah[0][1] = (b-1.f)/D;
    		Ah[0][2] = 0.f;

    		B[0][0] = (g0+g*b)/D;
    		B[0][1] = (-2.f*g0*c0)/D;
    		B[0][2] = (g0-g*b)/D;
    		B[0][3] = 0.f; 
    		B[0][4] = 0.f;

    		A[0][0] = 1.f;
    		A[0][1] = (-2.f*c0)/D;
    		A[0][2] = (1.f-b)/D;
    		A[0][3] = 0.f; 
    		A[0][4] = 0.f;
	}    

	for (i = 0; i <= L; ++i) {
		phi = (2.f*i-1.f)*M_PI/(2.f*N);                                
		si = sin(phi);
		D = b*b + 2.f*b*si + 1.f;
		b0h = (g*g*b*b + 2.f*g0*g*b*si + g0*g0)/D;        
		b1h = 2.f*(g*g*b*b - g0*g0)/D;
		b2h = (g*g*b*b - 2.f*g0*g*b*si + g0*g0)/D;
		a1h = 2.f*(b*b - 1.f)/D;
		a2h = (b*b - 2.f*b*si + 1.f)/D;
		
	        Bh[1+i][0] = b0h;
	        Bh[1+i][1] = b1h;
	        Bh[1+i][2] = b2h;

	        Ah[1+i][0] = 1.f;
		Ah[1+i][1] = a1h;
		Ah[1+i][2] = a2h;

		b0 = (g*g*b*b + g0*g0 + 2.f*g*g0*si*b)/D;
		b1 = -4.f*c0*(g0*g0 + g*g0*si*b)/D;
		b2 = 2.f*((1.f+2.f*c0*c0)*g0*g0 - g*g*b*b)/D;
		b3 = -4.f*c0*(g0*g0 - g*g0*si*b)/D;
		b4 = (g*g*b*b + g0*g0 - 2.f*g*g0*si*b)/D;
		a1 = -4.f*c0*(1.f + si*b)/D;
		a2 = 2.f*(1.f+2.f*c0*c0 - b*b)/D;
		a3 = -4.f*c0*(1.f - si*b)/D;
		a4 = (b*b - 2.f*si*b + 1.f)/D;

		B[1+i][0] = b0;
		B[1+i][1] = b1;
		B[1+i][2] = b2;
		B[1+i][3] = b3;
		B[1+i][4] = b4;
		A[1+i][0] = 1.f;
		A[1+i][1] = a1;
		A[1+i][2] = a2;
		A[1+i][3] = a3;
		A[1+i][4] = a4;
	}

	// LP or HP shelving filter
	if (c0==1.f || c0==-1.f)  { 
		for (i = 1; i <= L+1; ++i) {
			B[i][0] = Bh[i][0];
        	        B[i][1] = Bh[i][1]*c0;
	                B[i][2] = Bh[i][2];
	                B[i][3] = 0.f;
	                B[i][4] = 0.f;
	                A[i][0] = Ah[i][0];
	                A[i][1] = Ah[i][1]*c0;
	                A[i][2] = Ah[i][2];
	                A[i][3] = 0.f;
	                A[i][4] = 0.f;

			Bh[i][3] = 0.f;
			Bh[i][4] = 0.f;
			Ah[i][3] = 0.f;
			Ah[i][4] = 0.f;
			sanitize_denormal(B[i][1]);
			sanitize_denormal(B[i][2]);
			sanitize_denormal(A[i][1]);
			sanitize_denormal(A[i][2]);
			if (is_nan(B[i][0])) { B[i][0] = 1.f; }

		}
	}

 
	float AA  = sqrt(G);
	
	alpha = sin(w0)/2.f * sqrt( (AA + 1.f/AA)*(1.f/q - 1.f) + 2.f );
	b0 =    AA*( (AA+1.f) - (AA-1.f)*cos(w0) + 2.f*sqrt(AA)*alpha );
        b1 =  2.f*AA*( (AA-1.f) - (AA+1.f)*cos(w0)                   );
        b2 =    AA*( (AA+1.f) - (AA-1.f)*cos(w0) - 2.f*sqrt(AA)*alpha );
        a0 =        (AA+1.f) + (AA-1.f)*cos(w0) + 2.f*sqrt(AA)*alpha;
        a1 =   -2.f*( (AA-1.f) + (AA+1.f)*cos(w0)                   );
        a2 =        (AA+1.f) + (AA-1.f)*cos(w0) - 2.f*sqrt(AA)*alpha;
	
	B[2][0] = b0/a0;
        B[2][1] = b1/a0;
        B[2][2] = b2/a0;
        A[2][0] = 1.f;
        A[2][1] = a1/a0;
        A[2][2] = a2/a0;

	return true;
}

static void
run(LV2_Handle instance, uint32_t n_samples)
{
	ZamEQ2* zameq2 = (ZamEQ2*)instance;

	const float* const input  = zameq2->input;
	float* const       output = zameq2->output;

	const float        boostdb1 = *(zameq2->boostdb1);
	const float        q1 = *(zameq2->q1);
	const float        freq1 = *(zameq2->freq1);
	
	const float        boostdb2 = *(zameq2->boostdb2);
	const float        q2 = *(zameq2->q2);
	const float        freq2 = *(zameq2->freq2);
	
	const float        boostdbl = *(zameq2->boostdbl);
	const float        slopedbl = *(zameq2->slopedbl);
	const float        freql = *(zameq2->freql);

	float dcgain = 1.f;
	
	float boost1 = from_dB(boostdb1);
  	float fc1 = freq1 / zameq2->srate;
	float w01 = fc1*2.f*M_PI;
	float bwgain1 = (boostdb1 == 0.f) ? 1.f : (boostdb1 < 0.f) ? boost1*from_dB(3.f) : boost1*from_dB(-3.f);
	float bw1 = fc1 / q1;

	float boost2 = from_dB(boostdb2);
  	float fc2 = freq2 / zameq2->srate;
	float w02 = fc2*2.f*M_PI;
	float bwgain2 = (boostdb2 == 0.f) ? 1.f : (boostdb2 < 0.f) ? boost2*from_dB(3.f) : boost2*from_dB(-3.f);
	float bw2 = fc2 / q2;

	float boostl = from_dB(boostdbl);
	float fcl = zameq2->srate/4.f;
	float w0l = 0.f;
	float Al = sqrt(boostl);
	//float ql = 1.f / (sqrt((Al+1.f/Al)*(1.f/slopel-1.f)+2.f)); 
	
	float bwl = 2.f*M_PI*freql/ zameq2->srate;

	//bwl = 2.f/log(2.f)*arcsinh(1.f/(2.f*ql));
	

	float bwgaindbl = to_dB(Al);
	//float bwgaindbl = to_dB((boostl-dcgain)*exp(-1.f/ql)+dcgain);
	
	peq(dcgain,boost1,bwgain1,w01,bw1,&zameq2->a0x,&zameq2->a1x,&zameq2->a2x,&zameq2->b0x,&zameq2->b1x,&zameq2->b2x,&zameq2->gainx);
	peq(dcgain,boost2,bwgain2,w02,bw2,&zameq2->a0y,&zameq2->a1y,&zameq2->a2y,&zameq2->b0y,&zameq2->b1y,&zameq2->b2y,&zameq2->gainy);
	shelfeq(0.f,boostdbl,bwgaindbl,2.f*M_PI*freql/zameq2->srate,bwl,slopedbl,zameq2->B,zameq2->A,zameq2->Bh,zameq2->Ah);

	//printf("Gdb=%f fc=%f A=%f Q=%f bw=%f GBdb=%f\n",boostdbl,fcl, Al, ql, bwl, bwgaindbl);
	printf("B=[%f, %f, %f; %f, %f, %f]\n",zameq2->B[1][0],zameq2->B[1][1],zameq2->B[1][2],zameq2->B[2][0],zameq2->B[2][1],zameq2->B[2][2]);
	printf("A=[%f, %f, %f; %f, %f, %f]\n",zameq2->A[1][0],zameq2->A[1][1],zameq2->A[1][2],zameq2->A[2][0],zameq2->A[2][1],zameq2->A[2][2]);

		

	for (uint32_t pos = 0; pos < n_samples; pos++) {
		float tmp,tmpla, tmplb;
		float in = input[pos];
		sanitize_denormal(zameq2->x1);
		sanitize_denormal(zameq2->x2);
		sanitize_denormal(zameq2->y1);
		sanitize_denormal(zameq2->y2);
		sanitize_denormal(zameq2->x1a);
		sanitize_denormal(zameq2->x2a);
		sanitize_denormal(zameq2->y1a);
		sanitize_denormal(zameq2->y2a);
		sanitize_denormal(zameq2->zln1a);
		sanitize_denormal(zameq2->zln2a);
		sanitize_denormal(zameq2->zld1a);
		sanitize_denormal(zameq2->zld2a);
		sanitize_denormal(zameq2->zln1b);
		sanitize_denormal(zameq2->zln2b);
		sanitize_denormal(zameq2->zld1b);
		sanitize_denormal(zameq2->zld2b);
		sanitize_denormal(in);

		//lowshelf a
		tmpla = in * zameq2->B[2][0] + 
			zameq2->zln1a * zameq2->B[2][1] +
			zameq2->zln2a * zameq2->B[2][2] -
			zameq2->zld1a * zameq2->A[2][1] -
			zameq2->zld2a * zameq2->A[2][2];
		zameq2->zln2a = zameq2->zln1a;
		zameq2->zld2a = zameq2->zld1a;
		zameq2->zln1a = in;
		zameq2->zld1a = tmpla;
		
		//lowshelf b
/*		tmplb = tmpla * zameq2->B[2][0] + 
			zameq2->zln1b * zameq2->B[2][1] +
			zameq2->zln2b * zameq2->B[2][2] -
			zameq2->zld1b * zameq2->A[2][1] -
			zameq2->zld2b * zameq2->A[2][2];
		zameq2->zln2b = zameq2->zln1b;
		zameq2->zld2b = zameq2->zld1b;
		zameq2->zln1b = tmpla;
		zameq2->zld1b = tmplb;
*/		
		//parametric1
		tmp = tmpla * zameq2->b0x + zameq2->x1 * zameq2->b1x + zameq2->x2 * zameq2->b2x - zameq2->y1 * zameq2->a1x - zameq2->y2 * zameq2->a2x;
		zameq2->x2 = zameq2->x1;
		zameq2->y2 = zameq2->y1;
		zameq2->x1 = tmpla;
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
