all: rms flipG192 mld align

rms: rms.c tinywavein_c.h
	gcc -Wall -O2 -lm $< -o $@

flipG192: flipG192.c
	gcc -Wall -O2 -lm $< -o $@

mld:
	make -C mld_src/build config=release
	cp mld_src/bin/mld .

align: align.c tinywavein_c.h tinywaveout_c.h
	gcc -Wall -O2 -lm $< -o $@

clean:
	rm -rf rms flipG192 mld align mld_src/bin/ mld_src/build/obj
