#******************************************************************************
#                        ETSI TS 103 634 V1.6.1                               *
#              Low Complexity Communication Codec Plus (LC3plus)              *
#                                                                             *
# Copyright licence is solely granted through ETSI Intellectual Property      *
# Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
# estoppel or otherwise.                                                      *
#*****************************************************************************/

all: rms flipG192 mld align

rms: rms.c tinywavein_c.h
	gcc -Wall -O2 $< -lm -o $@

flipG192: flipG192.c
	gcc -Wall -O2 $< -lm -o $@

mld:
	make -C mld_src/build config=release
	cp mld_src/bin/mld .

align: align.c tinywavein_c.h tinywaveout_c.h
	gcc -Wall -O2 $< -lm -o $@

clean:
	rm -rf rms flipG192 mld align mld_src/bin/ mld_src/build/obj
