#include <lv2plugin.hpp>
#include <stdint.h>
#include <stdio.h>
#include <cmath>

using namespace LV2;


class ZamEQ2 : public Plugin<ZamEQ2> {
public:
  
  ZamEQ2(double rate)
    : Plugin<ZamEQ2>(5) {
    
    srate = rate;
    x1=x2=y1=y2=0.f;
    a0x=a1x=a2x=b0x=b1x=b2x=0.f;
    a0y=a1y=a2y=b0y=b1y=b2y=0.f;
    gainx=gainy=0.f;
  }

  void peq(float G0, float G, float GB, float w0, float Dw,
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

  // Works on little-endian machines only
  inline bool is_nan(float& value ) {
    if (((*(uint32_t *) &value) & 0x7fffffff) > 0x7f800000) {
      return true;
    }
    return false;
  }

  // Force already-denormal float value to zero
  inline void sanitize_denormal(float& value) {
    if (is_nan(value)) {
        value = 0.f;
    }
  }

  inline int sign(float x) {
  	return (x >= 0.f ? 1 : -1);
  }

  inline float from_dB(float gdb) {
  	return (exp(gdb/20.f*log(10.f)));
  }

  inline float to_dB(float g) {
  	return (20.f*log10(g));
  }

  float x1,x2,y1,y2;
  float a0x,a1x,a2x,b0x,b1x,b2x;
  float a0y,a1y,a2y,b0y,b1y,b2y;
  float gainx,gainy;
  float srate;

 void run(uint32_t nframes) {
  
  float dcgain1 = 1.f; 
  //float boost1 = exp(*p(0)/20.f*log(10.f));
  //float q1 = *p(1);
  //float fc1 = *p(2)/srate;
  
  //12dB boost at 1000Hz Q=5
  float boost1 = exp(12.f/20.f*log(10.f));
  float q1 = 5.f;
  float fc1 = 1000.f/srate;

  float w01 = fc1*2.f*M_PI;
  
  float bwgain1 = (*p(0) == 0.f) ? 1.f : (*p(0) < 0.f) ? boost1*from_dB(3.f) : boost1*from_dB(-3.f);
  float bw1 = fc1 / q1;
  
 /* 
  float dcgain2 = *p(5);
  float boost2 = *p(6);
  float bwgain2 = *p(7);
  float fc2 = *p(8)*2.f*M_PI / srate;
  float bw2 = *p(9)*2.f*M_PI / srate;
 */
  
  uint32_t i = 0;

  peq(dcgain1,boost1,bwgain1,w01,bw1,&a0x,&a1x,&a2x,&b0x,&b1x,&b2x,&gainx);
  //printf("B = [%f, %f, %f]\n",b0x, b1x, b2x);
  //printf("A = [%f, %f, %f]\n\n",a0x, a1x, a2x);


  //peq(dcgain2,boost2,bwgain2,fc2,bw2,&a0y,&a1y,&a2y,&b0y,&b1y,&b2y,&gainy);
  for (i = 0; i < nframes; ++i) {
    sanitize_denormal(x1);
    sanitize_denormal(x2);
    sanitize_denormal(y1);
    sanitize_denormal(y2);
    p(4)[i] = (p(3)[i] * b0x + x1 * b1x + x2 * b2x - y1 * a1x - y2 * a2x);
    sanitize_denormal(p(4)[i]);
    x2 = x1;
    y2 = y1;
    x1 = p(3)[i];
    y1 = p(4)[i];
  }
 }  

};

static int _ = ZamEQ2::register_class("http://zamaudio.com/lv2/zameq2");
