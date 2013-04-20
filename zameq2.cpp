#include <lv2plugin.hpp>
#include <stdint.h>
#include <stdio.h>
#include <cmath>

#define PEAKING		0
#define LOWSHELF	1
#define HIGHSHELF	2

using namespace LV2;

class ZamEQ2 : public Plugin<ZamEQ2> {
public:
  
  ZamEQ2(double rate)
    : Plugin<ZamEQ2>(6) {
    
    srate = rate;
    y1=y2=0.f;
    z[12]=z[11]=z[10]=z[9]=z[8]=z[7]=z[6]=z[5]=z[4]=z[3]=z[2]=z[1]=0.f;
    x[5]=x[4]=x[3]=x[2]=x[1]=0.f;
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
  
  float boostdb;
  float f0;
  float Q;
  int type;

  float a1,a2,b1,b2,c,d,dt,srate;
  float BB[12];
  float CC[12];
  float z[13],y1,y2,x[6];
  float K;
  float w0;
 
  float A11(float t) {
	if (Q >= 0.5) {
		return (exp(-a1*w0*t/2.f)*cos(d*w0*t/2.f)-(a1/d)*exp(-a1*w0*t/2.f)*sin(d*w0*t/2.f));
	} else {
		return (exp(-a1*w0*t/2.f)*cosh(d*w0*t/2.f)-(a1/d)*exp(-a1*w0*t/2.f)*sinh(d*w0*t/2.f));
	}
  }

  float A12(float t) {
	if (Q >= 0.5) {
		return ((2.f/d)*exp(-a1*w0*t/2.f)*sin(d*w0*t/2.f));
	} else {
		return ((2.f/d)*exp(-a1*w0*t/2.f)*sinh(d*w0*t/2.f));
	}
  }

  float A21(float t) {
	if (Q >= 0.5) {
		return (-(2.f*a2/d)*exp(-a1*w0*t/2.f)*sin(d*w0*t/2.f));
	} else {
		return (-(2.f*a2/d)*exp(-a1*w0*t/2.f)*sinh(d*w0*t/2.f));
	}
  }

  float A22(float t) {
	if (Q >= 0.5) {
		return (exp(-a1*w0*t/2.f)*cos(d*w0*t/2.f)+(a1/d)*exp(-a1*w0*t/2.f)*sin(d*w0*dt/2.f));
	} else {
		return (exp(-a1*w0*t/2.f)*cosh(d*w0*t/2.f)+(a1/d)*exp(-a1*w0*t/2.f)*sinh(d*w0*t/2.f));
	}
  }

  float sinc(float t) {
	return ((t==0.f) ? 1.f : sin(M_PI*t)/(M_PI*t));
  }

  float simpson1(int j, float a, float b, float n) {
        float ret;
	double h = (b - a) / n;
        float s;
        int i;
        s = Bf(a,j) + Bf(b,j);

        for (i = 1; i < n; i += 2) {
                s += 4.f * Bf(a + i * h, j);
        }

        for (i = 2; i < n-1; i += 2) {
                s += 2.f * Bf(a + i * h, j);
        }

        ret = s * h / 3.f;
//	sanitize_denormal(ret);
	return (ret);
  }

  float simpson2(int j, float a, float b, float n) {
        float ret;
	double h = (b - a) / n;
        float s;
        int i;
        s = Cf(a,j) + Cf(b,j);

        for (i = 1; i < n; i += 2) {
                s += 4.f * Cf(a + i * h, j);
        }

        for (i = 2; i < n-1; i += 2) {
                s += 2.f * Cf(a + i * h, j);
        }

        ret = s * h / 3.f;
//	sanitize_denormal(ret);
	return (ret);
  }


  float Bf(float s, int j) {
        return ((A11(dt-s)*b1*w0 + A12(dt-s)*b2*w0)*sinc((s+j*dt)/dt)*(0.54f+0.46f*cos(M_PI*(s+j*dt)/(5.f*dt))));
  }

  float Cf(float s, int j) {
        return ((A21(dt-s)*b1*w0 + A22(dt-s)*b2*w0)*sinc((s+j*dt)/dt)*(0.54f+0.46f*cos(M_PI*(s+j*dt)/(5.f*dt))));
  }

 
  void run(uint32_t nframes) {
  int ii;
  uint32_t i;

  boostdb = *p(0);
  f0 = *p(1);
  Q = *p(2) ;
  type = *p(3);

  K = exp(boostdb/40.f*log(10.f));
  w0 = 2.f*M_PI*f0;
  dt = 1.f / srate;

    switch (type) {
        case PEAKING:
                a1 = 1.f / (Q*K);
                a2 = 1.f;
                b1 = (K - 1.f) / (K*Q);
                b2 = 0.f;
                c = 1.f;
                break;

        case LOWSHELF:
                a1 = 1.f / (sqrt(K)*Q);
                a2 = 1.f / K;
                b1 = (sqrt(K) - 1.f) / (sqrt(K)*Q);
                b2 = (K - 1.f) / K;
                c = 1.f;
                break;

        case HIGHSHELF:
                a1 = sqrt(K) / Q;
                a2 = K;
                b1 = (sqrt(K) - K*K*sqrt(K)) / Q;
                b2 = K - K*K*K;
                c = K*K;
                break;
        default:
		printf("wtf type=%d\n",type);
    }

    d = sqrt(std::abs(4.f*a2-a1*a1));
	

    for (ii = 0; ii <= 10; ++ii) {
	BB[ii] = simpson1(ii-5, 0.f, dt, 10.f);
	CC[ii] = simpson2(ii-5, 0.f, dt, 10.f);
//	printf("%d:%f:%f\n",ii,BB[ii],CC[ii]);
    }


    float aa22 = A22(dt);
    float aa12 = A12(dt);

    for (i = 0; i < nframes; ++i) {
      int k;
      sanitize_denormal(p(4)[i]);
      sanitize_denormal(z[11]);
      sanitize_denormal(z[10]);
      sanitize_denormal(z[9]);
      sanitize_denormal(z[8]);
      sanitize_denormal(z[7]);
      sanitize_denormal(z[6]);
      sanitize_denormal(z[5]);
      sanitize_denormal(z[4]);
      sanitize_denormal(z[3]);
      sanitize_denormal(z[2]);
      sanitize_denormal(z[1]);
      for (k = 0; k <= 10; ++k) {
 //       sanitize_denormal(BB[k]);
 //       sanitize_denormal(CC[k]);
      }
      sanitize_denormal(y2);
      sanitize_denormal(y1);


      p(5)[i] = -BB[10]*aa22-CC[10]*aa12*z[12] + 
	(-BB[9]*aa22 - CC[9]*aa12 + BB[10])*z[11] + 
	(-BB[8]*aa22 - CC[8]*aa12 + BB[9])*z[10] + 
	(-BB[7]*aa22 - CC[7]*aa12 + BB[8])*z[9] + 
	(-BB[6]*aa22 - CC[6]*aa12 + BB[7])*z[8] + 
	(-BB[5]*aa22 - CC[5]*aa12 + BB[6] + 
		c*(A11(dt)*A22(dt)-A21(dt)*A12(dt)))*z[7] + 
	(-BB[4]*aa22 - CC[4]*aa12 + BB[5] - 
		c*(A11(dt)+A22(dt)))*z[6] + 
	(-BB[3]*aa22 - CC[3]*aa12 + BB[4] +
		c)*z[5] + 
	(-BB[2]*aa22 - CC[2]*aa12 + BB[3])*z[4] + 
	(-BB[1]*aa22 - CC[1]*aa12 + BB[2])*z[3] + 
	(-BB[0]*aa22 - CC[0]*aa12 + BB[1])*z[2] + 
	(BB[0]*z[1]) -
	p(4)[i]*1.f - (A11(dt)+A22(dt))*y1 - (A11(dt)*A22(dt)-A12(dt)*A21(dt))*y2;

      if(p(5)[i] > 1.f) { p(5)[i] = 1.f; }
      if(p(5)[i] < -1.f) { p(5)[i] = -1.f; }

      z[12] = z[11];
      z[11] = z[10];
      z[10] = z[9];
      z[9] = z[8];
      z[8] = z[7];
      z[7] = z[6];
      z[6] = z[5];
      z[5] = z[4];
      z[4] = z[3];
      z[3] = z[2];
      z[2] = z[1];
      z[1] = p(4)[i];

      y2 = y1;
      y1 = p(5)[i];
    }
  } 

};

static int _ = ZamEQ2::register_class("http://zamaudio.com/lv2/zameq2");
