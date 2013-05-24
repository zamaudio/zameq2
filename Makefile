all: lv2/zameq2.lv2/zameq2.ttl ladspa/zameq2.so

ladspa: ladspa/zameq2.so

lv2: lv2/zameq2.lv2/zameq2.ttl

lv2/zameq2.lv2/zameq2.ttl: lv2/zameq2.cpp
	./compilelv2 zameq2.dsp

lv2/zameq2.cpp:
	./genlv2 zameq2.dsp	

ladspa/zameq2.so: ladspa/zameq2.dsp.cpp
	./compileladspa zameq2.dsp

ladspa/zameq2.dsp.cpp:
	./genladspa zameq2.dsp

install:
	mkdir -p /usr/local/lib/lv2
	mkdir -p /usr/local/lib/ladspa
	cp -a lv2/zameq2.lv2 /usr/local/lib/lv2
	cp -a ladspa/zameq2.so /usr/local/lib/ladspa
	
uninstall:
	rm -fr /usr/local/lib/lv2/zameq2.lv2
	rm -f /usr/local/lib/ladspa/zameq2.so

clean:
	rm -fr lv2/zameq2.lv2
	rm -f ladspa/zameq2.so

distclean:
	rm -fr lv2/*
	rm -fr ladspa/*
