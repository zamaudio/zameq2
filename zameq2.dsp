declare name "ZamEQ2";
declare author "Damien Zammit";
declare copyright "2013";
declare version "1.00";
declare license "GPLv2";

import("math.lib");
import("filter.lib");
import("music.lib");

sign(x) = ((x >= 0.0) - 0.5)*2.0;
from_dB(gdb) = (exp(gdb/20.0*log(10.0)));
to_dB(g) = (20.0*log10(g));

peq(G0, G, GB, w0, Dw) = TF2(bb0,bb1,bb2,aa1,aa2)
	with {
		F = abs(G*G - GB*GB);
		G00 = abs(G*G - G0*G0);
		F00 = abs(GB*GB - G0*G0);
		num = G0*G0 * (w0*w0 - PI*PI)*(w0*w0 - PI*PI)
			+ G*G * F00 * PI*PI * Dw*Dw / F;
		den = (w0*w0 - PI*PI)*(w0*w0 - PI*PI)
			+ F00 * PI*PI * Dw*Dw / F;
		G1 = sqrt(num/den);
		G01 = abs(G*G - G0*G1);
		G11 = abs(G*G - G1*G1);
		F01 = abs(GB*GB - G0*G1);
		F11 = abs(GB*GB - G1*G1);
		W2 = sqrt(G11 / G00) * tan(w0/2.0)*tan(w0/2.0);
		Dww = (1.0 + sqrt(F00 / F11) * W2) * tan(Dw/2.0);
		C = F11 * Dww*Dww - 2.0 * W2 * (F01 - sqrt(F00 * F11));
		D = 2.0 * W2 * (G01 - sqrt(G00 * G11));
		A = sqrt((C + D) / F);
		B = sqrt((G*G * C + GB*GB * D) / F);
		
		aa = 1.0 + W2 + A;
		b0 = (G1 + G0*W2 + B);
		b1 = -2.0*(G1 - G0*W2);
		b2 = (G1 - B + G0*W2);

		a1 = -2.0*(1.0 - W2);
		a2 = (1 + W2 - A);
		
		bb0 = if((aa==0.0),0.0,b0/aa);
		bb1 = if((aa==0.0),0.0,b1/aa);
		bb2 = if((aa==0.0),1.0,b2/aa);
		
		aa1 = if((aa==0.0),0.0,a1/aa);
		aa2 = if((aa==0.0),0.0,a2/aa);
	};

lowshelfeq(G0, GG, GB, w0, Dw, q) = TF2(bb0,bb1,bb2,aa1,aa2)
	with {
		G = pow(10.0, GG/20.0); 
		AA  = sqrt(G);
		
		alpha = sin(w0)/2.0 * sqrt( (AA + 1.0/AA)*(1.0/q - 1.0) + 2.0 );
		b0 =    AA*( (AA+1.0) + (-AA+1.0)*cos(w0) + 2.0*sqrt(AA)*alpha );
		b1 =  2.0*AA*( (AA-1.0) + (-AA-1.0)*cos(w0)                   );
		b2 =    AA*( (AA+1.0) + (-AA+1.0)*cos(w0) - 2.0*sqrt(AA)*alpha );
		a0 =        (AA+1.0) + (AA-1.0)*cos(w0) + 2.0*sqrt(AA)*alpha;
		a1 =   -2.0*( (AA-1.0) + (AA+1.0)*cos(w0)                   );
		a2 =        (AA+1.0) + (AA-1.0)*cos(w0) - 2.0*sqrt(AA)*alpha;
		
		bb0 = if((a0==0.0),0.0,b0/a0);
		bb1 = if((a0==0.0),0.0,b1/a0);
		bb2 = if((a0==0.0),1.0,b2/a0);
		
		aa1 = if((a0==0.0),0.0,a1/a0);
		aa2 = if((a0==0.0),0.0,a2/a0);
	};

highshelfeq(G0, GG, GB, w0, Dw, q) = TF2(bb0,bb1,bb2,aa1,aa2)
	with {
		G = pow(10.0, GG/20.0); 
		AA  = sqrt(G);

		alpha = sin(w0)/2.0 * sqrt( (AA + 1.0/AA)*(1.0/q - 1.0) + 2.0 );
		b0 =    AA*( (AA+1.0) + (AA-1.0)*cos(w0) + 2.0*sqrt(AA)*alpha );
		b1 =  -2.0*AA*( (AA-1.0) + (AA+1.0)*cos(w0)                   );
		b2 =    AA*( (AA+1.0) + (AA-1.0)*cos(w0) - 2.0*sqrt(AA)*alpha );
		a0 =        (AA+1.0) + (-AA+1.0)*cos(w0) + 2.0*sqrt(AA)*alpha;
		a1 =   2.0*( (AA-1.0) + (-AA-1.0)*cos(w0)                   );
		a2 =        (AA+1.0) + (-AA+1.0)*cos(w0) - 2.0*sqrt(AA)*alpha;
		
		bb0 = if((a0==0.0),0.0,b0/a0);
		bb1 = if((a0==0.0),0.0,b1/a0);
		bb2 = if((a0==0.0),1.0,b2/a0);
		
		aa1 = if((a0==0.0),0.0,a1/a0);
		aa2 = if((a0==0.0),0.0,a2/a0);
	};

pianokey2hz(x) = 440.0*pow(2.0, (x-49.0)/12); // piano key 49 = A440

boostdb1 = hslider("1 Boost (dB)", 0.0, -20.0, 20.0, 0.1) : smooth(0.99);
q1 = hslider("1 Q (dB)", 0.0, -10.0, 40.0, 0.1) : smooth(0.99) : from_dB;
//freq1 = hslider("1 Freq (Hz) [unit:PK]", 1000.0, 20.0, 20000.0, 1.0);
freq1 = hslider("1 Freq [unit:PK]", 64, 1, 100, 0.01) : smooth(0.99) : pianokey2hz;

boostdb2 = hslider("2 Boost (dB)", 0.0, -20.0, 20.0, 0.1) : smooth(0.99);
q2 = hslider("2 Q (dB)", 0.0, -10.0, 40.0, 0.1) : smooth(0.99) : from_dB;
//freq2 = hslider("2 Freq (Hz)", 4000.0, 20.0, 20000.0, 1.0) : smooth(0.99);
freq2 = hslider("2 Freq [unit:PK]", 88, 1, 100, 0.01) : smooth(0.99) : pianokey2hz;

boostdbl = hslider("0 Boost Lowshelf (dB)", 0.0, -20.0, 20.0, 0.1) : smooth(0.99);
slopedbl = hslider("0 Slope Lowshelf", 1.0, 1.0, 1.5, 0.1) : smooth(0.99);
//freql = hslider("0 Freq Lowshelf (Hz)", 250.0, 10.0, 20000.0, 1.0) : smooth(0.99);
freql = hslider("0 Freq Lowshelf [unit:PK]", 40, 1, 100, 0.01) : smooth(0.99) : pianokey2hz;

boostdbh = hslider("3 Boost Highshelf (dB)", 0.0, -20.0, 20.0, 0.1) : smooth(0.99);
slopedbh = hslider("3 Slope Highshelf", 1.0, 1.0, 1.5, 0.1) : smooth(0.99);
//freqh = hslider("3 Freq Highshelf (Hz)", 9000.0, 10.0, 20000.0, 1.0) : smooth(0.99);
freqh = hslider("3 Freq Highshelf [unit:PK]", 100, 1, 100, 0.01) : smooth(0.99) : pianokey2hz;

freqlp = hslider("4 Freq Lowpass [unit:PK]", 70, 1, 100, 0.01) : smooth(0.99) : pianokey2hz;
bypasslp = checkbox("LP");
freqhp = hslider("5 Freq Highpass [unit:PK]", 20, 1, 100, 0.01) : smooth(0.99) : pianokey2hz;
bypasshp = checkbox("HP");

mainlevel = vslider("6 Master Trim (dB)", 0.0, -12.0, 12.0, 0.1) : smooth(0.99) : from_dB;

peakingeq1(x) = if((boostdb1 == 0.0), x, x : peq(dcgain,boost1,bwgain1,w01,bw1))
	with {
                dcgain = 1.0;
                boost1 = from_dB(boostdb1);
                fc1 = freq1 / SR;
                w01 = fc1*2.0*PI;
                bwgain1 = if((boostdb1 == 0.0), 1.0, (if((boostdb1 < 0.0), boost1*from_dB(3.0), boost1*from_dB(-3.0))));
                bw1 = fc1 / q1;
	};

peakingeq2(x) = if((boostdb2 == 0.0), x, x : peq(dcgain,boost2,bwgain2,w02,bw2))
	with {
		dcgain = 1.0;
		boost2 = from_dB(boostdb2);
                fc2 = freq2 / SR;
                w02 = fc2*2.0*PI;
                bwgain2 = if((boostdb2 == 0.0), 1.0, (if((boostdb2 < 0.0), boost2*from_dB(3.0), boost2*from_dB(-3.0))));
                bw2 = fc2 / q2;
	};
		
lowshelf1(x) = if((boostdbl == 0.0), x, x : lowshelfeq(0.0,boostdbl,bwgaindbl,bwl,bwl,slopedbl)) 
	with {
		boostl = from_dB(boostdbl);
		All = sqrt(boostl);
		bwl = 2.0*PI*freql/SR;
		bwgaindbl = to_dB(All);
	};

highshelf1(x) = if((boostdbh == 0.0), x, x : highshelfeq(0.0,boostdbh,bwgaindbh,bwh,bwh,slopedbh))
	with {
		boosth = from_dB(boostdbh);
		Ahh = sqrt(boosth);
		bwh = 2.0*PI*freqh/SR;
		bwgaindbh = to_dB(Ahh);
	};
	
lowpass1(x) = if((bypasslp == 1.0), x : lowpass(4,freqlp), x);

highpass1(x) = if((bypasshp == 1.0), x : highpass(4,freqhp), x);

process = *(peakingeq1 : peakingeq2 : lowshelf1 : highshelf1 : lowpass1 : highpass1, mainlevel);
