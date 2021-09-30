 ================================================================
 LC3plus Precision Requirement script for High Resolution V1.0.0
 ================================================================

 /******************************************************************************
 *                        ETSI TS 103 634 V1.3.1                               *
 *              Low Complexity Communication Codec Plus (LC3plus)              *
 *                                                                             *
 * Copyright licence is solely granted through ETSI Intellectual Property      *
 * Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
 * estoppel or otherwise.                                                      *
 ******************************************************************************/

The script precision_requirements.py performs THD+N and SNR measurements according
to [1] and produces results to check if the requirements on precision for the high
resolution mode are fulfilled for a vendor implementation of LC3plus. 

Pre-requisites
==============

 - python3
 - python matplotlib module
 - SoX V14.4.2 (http://sox.sourceforge.net)
 - Unix based operating system, or Cygwin if using Windows

Usage:
======

python3 precision_requirements.py LC3plus_precision.cfg

Before execution:
 - Build LC3plus floating point reference executable.
 - set paths to encoder and decoder executables under test in LC3plus_precision.cfg.

The script requires a configuration file which contains paths to executables, settings and
operating points to be tested. 



Usage of configuration file
============================

The configuration file is separated in a [globals] and [tests] section by square brackets. 
Within each section, variables can be set by 'variable=value'.

In the [globals] section, general parameters have to be set like paths to binaries and 
thresholds (see LC3plus_precision.cfg).

(reference|test)_(encoder|decoder) : path to binary
example:
reference_encoder = ../src/floating_point/LC3plus -E -q -hrmode -formatG192 -frame_ms {frame_ms} {infile} {bitstream} {bitrate}
reference_decoder = ../src/floating_point/LC3plus -D -q -formatG192 -bps {bps} {bitstream} {outfile}
For test_(encoder|decoder) the strings in curly brackets have to stay the same, order and flags can change.

(enc|dec)_(thd|snr)_threshold      : threshold for worst case value out of all measured frequency points
(enc|dec)_(thd|snr)_threshold_1k   : threshold for (thd|snr) at 1kHz sine frequency

For 'encdec' mode the more critical threshold out of set enc_..._threshold, dec_..._threshold is selected.

num_cores       : specifies number of used processing cores for parallelization
plot            : (1 | 0) generate plots for every operating point
freq_grid       : (linlog | onethirdoctave) used frequency grid
freq_grid_steps : number of steps for frequency grid
bps             : (24 | 32) bits per sample of wave file

In the [tests] section operating points can be defined in 'config =' by Mode, Frame Size, Sampling rate 
and Bitrate as described in LC3plus_precision.cfg. Mode must be (encode|decode|encdec). Every new line 
of operating points has to be indented.

Results
=======

After all operating points are processed at all frequencies, results are saved in LC3plus_precision.html. 
The html file contains a table showing the measured THD and SNR values. If the value has passed the 
corresponding threshold criteria the cell is presented blue, otherwise red.

If 'plot = 1' was set in the config file, a click on the THD/SNR value in the table opens a plot showing
all measured frequencies of the operating point.

The plots are saved in png_LC3plus_precision, the generated wave files in wav_LC3plus_precision.

References
==========
[1] ETSI TS 103 634 Chapter 7.3.5.4	Precision tests (THD+N)
