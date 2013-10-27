zameq2
======

ZamEQ2 - Equalizer plugin

Yes, this is yet another equalizer plugin.
The difference with this one is that you get 2 parametric sections using
some of the latest Orfanidis butterworth filters.
These filters are supposed to behave more like analog EQs close to the
nyquist frequency.  Also, you get a shelving circuit on each end and a HP/LP
section which can be toggled.

Feel free to leave comments on my blog as all feedback is appreciated.

http://www.zamaudio.com/?p=870


Install instructions:

Faust >=0.9.62 and LADSPA SDK is required to compile 
this plugin. Build scripts are provided for LADSPA.

	make && sudo make install

