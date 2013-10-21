all: ladspa/zameq2.so

ladspa/zameq2.so: ladspa/zameq2.dsp.cpp
	./compileladspa zameq2.dsp

ladspa/zameq2.dsp.cpp:
	./genladspa zameq2.dsp

install:
	mkdir -p /usr/local/lib/ladspa
	cp -a ladspa/zameq2.so /usr/local/lib/ladspa
	
uninstall:
	rm -f /usr/local/lib/ladspa/zameq2.so

clean:
	rm -f ladspa/zameq2.so
