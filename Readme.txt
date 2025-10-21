/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                        Software Version V1.8.0ETSI                          *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

Description
-----------
This package contains the fixed-point and floating-point source code of the
LC3plus codec as electronic attachment to ETSI TS 103 634.

The software uses fix-point arithmetic utilizing the ITU-T STL2009 including
the latest updates introduced by 3GPP and/or IEEE 754 floating-point arithmetic.

The following structure outlines the content of this package.

-/conformance         : conformance script with example configuration file and
                        Readme
-/fec_control_unit    : example fec control unit script and Readme
-/src/fixed_point     : fixed-point source files, makefile and Visual Studio
                        solution
-/src/floating_point  : floating-point source files, makefile and Visual Studio
                        solution
-/testvec             : testvector package with script. Verification that LC3plus
                        build on your system operates as expected by comparing to
                        pre-calculated MD5 hashes
-/tools               : LC3plus helper tools and examples

Please refer to the respective Readme files for more information.


Features
--------
    - Supported sampling rates: 8 kHz, 16 kHz, 24 kHz, 32 kHz, 44.1 kHz,
                                48 kHz, 96 kHz
    - Multichannel support by multi-mono coding
    - Support of audio sample depth: 16 bits and 24 bits
    - Frame duration of 10 ms, 7.5 ms, 5 ms, 2.5 ms and 1.25 ms
    - Supported bit rates as outlined in Table 5.1 and Table 5.2 of TS 103 634
    - Packet loss concealment as defined in ETSI TS 103 634

Changelog
---------
  Latest non-bitexact changes are encapsulated in defines listed in defines.h
  for review.
  
    - V1.8.0ETSI 2025-10-21 (ETSI TS 103 634 V1.6.1)
       - General
            - Updated source headers to version V1.6.1
  
    - V1.8.0ETSI 2025-09-12 (ETSI TS 103 634 V1.5.6)
       - General
            - Finalize support for 1.25 ms frame duration for both fixed and floating point
            - For hrmode, 1.25 ms is not supported

    - V1.7.6ETSI 2025-06-02 (ETSI TS 103 634 V1.5.5)
       - General
            - Improved support for 1.25 ms frame duration

    - V1.7.5ETSI 2025-03-10 (ETSI TS 103 634 V1.5.3)
       - General
            - Added support for 1.25 ms frame duration
  
    - V1.7.4ETSI 2024-05-01 (ETSI TS 103 634 V1.5.1)
       - General
            - Accepted defines CR8, CR9, CR10 and CR11
  
    - V1.7.3ETSI 2024-04-01 (ETSI TS 103 634 V1.4.5)
       - General
            - Implemented CR1 on ETSI TS 103 634 V1.4.4 (changes marked with CR11)

    - V1.7.2ETSI 2024-03-13 (ETSI TS 103 634 V1.4.4)
       - General
            - Implemented CR1 on ETSI TS 103 634 V1.4.3 (changes marked with CR10)

    - V1.7.1ETSI 2024-02-09 (ETSI TS 103 634 V1.4.4)
       - General
            - Implemented CR1 on ETSI TS 103 634 V1.4.2 (changes marked with CR9)
  
    - V1.7.0ETSI 2023-09-15 (ETSI TS 103 634 V1.4.2)
       - General
            - Implemented CR1, CR2 and CR3 on ETSI TS 103 634 V1.4.1 (changes marked with CR8)
  
    - V1.6.9ETSI 2023-03-01 (ETSI TS 103 634 V1.4.1)
       - General
            - Updated source headers to version V1.4.1
            
    - V1.6.8ETSI 2023-01-20 (ETSI TS 103 634 V1.3.4)
       - General
            - Accepted defines CR6 and CR7
  
    - V1.6.7ETSI 2022-11-18 (ETSI TS 103 634 V1.3.3)
       - General
            - Updated software version and ETSI TS version
            - Implemented CR1 on TS 103 634 V1.3.2 (changes marked with defines CR7)

    - V1.6.6ETSI 2022-06-22 (ETSI TS 103 634 V1.3.2)
       - General
            - Updated software version and ETSI TS version
            - Implemented CR1 on TS 103 634 V1.3.1 (changes marked with defines CR6)

    - V1.6.4ETSI 2021-10-01 (ETSI TS 103 634 V1.3.1)
       - General
            - Updated software version and ETSI TS version
  
    - V1.6.3ETSI 2021-08-01 (ETSI TS 103 634 V1.2.5)
       - General
            - Implemented CR1 on TS 103 634 V1.2.4 (changes marked with defines CR6)
            - Accepted CR1 to CR5
  
    - V1.6.2ETSI 2021-06-23 (ETSI TS 103 634 V1.2.4)
       - General
            - Implemented CR1 on TS 103 634 V1.2.3 (changes marked with defines CR5)
            - Major Features
                 - Added high resolution audio support and improved precision coding for regular modes
                   in fixed-point code (wrapped in the define ENABLE_HR_MODE)
                 - Added support for fallback bitrates for the high resolution mode
                   for critical conditions

    - V1.6.0ETSI 2021-05-25 (ETSI TS 103 634 V1.2.3)
       - General
            - Implemented CR1 on TS 103 634 V1.2.2 (changes marked with defines CR4)

    - V1.4.12ETSI 2021-03-03 (ETSI TS 103 634 V1.2.2)
       - General
            - Implemented CR1 on TS 103 634 V1.2.1 (changes marked with defines CR3)
            - Major Features
                 - Added concealment and forward error protection for floating point code
                 - Limited bit rate to 320 kbps for regular modes
    - V1.4.10ETSI 2020-07-23 (ETSI TS 103 634 V1.2.1)
       - General
            - Accepted CR2 defines
    - V1.4.9E 2020-03-23 (ETSI TS 103 634 V1.1.3)
        - General
            - Fix PLC fade-out for frame lengths < 10ms, implemented
              under define CR2_B_NONBE_PLC3_ENHANCEMENTS
        - Channel Coder
            - added epok_flags according to CR2_G, implemented
              under define CR2_G_EPOK_FLAGS
            - fixed incorrect value of error_report for
              certain cases with bfi=1 according to CR2_J, implemented
              under define CR2_J_NONBE_ERROR_REPORT_FIX
    - V1.4.6E 2020-01-22 (ETSI TS 103 634 V1.1.2)
        - General
            - Renamed NONBE defines according to corresponding change in CR

    - V1.4.5E 2020-01-10 (ETSI TS 103 634 V1.1.2)
        - General
            - Added Channel Coder Converter
            - Several Packet Loss Concealment simplifications
            - Aligned LTPF normcorr window

    - V1.4.2 2019-07-25 (ETSI TS 103 634 V1.1.1)
        - General
            - Updated copyright header

    - V1.4.1 2019-06-26 (ETSI TS 103 634 V0.1.1)
        - General
            - Added several sanity checks for corrupt bitstreams
            - Fix delay compensation (dc == 2)
            - Changed setup struct size defines to actual maximum values


Building
--------
    Unix platforms:
        - Go to src/fixed_point or src/floating_point folder
        - Call "make"
        - Executable path and name "./LC3plus"

    Windows platforms:
        - Go to src/fixed_point/msvc or src/floating_point/msvc folder
        - Open up solution file LC3plus.sln and build it
        - Standard config executable path and name ".\Win32\Release\LC3plus.exe"
        - The solution is optimized for Visual Studio 2017


The floating-point source code introduces a floating-point data type and
corresponding functions via macros labeled LC3_XYZ. They use single
precision floating point format to represent data and perform mathematical
operations. 


Usage
-----
    The following example commands explain the usage of the LC3plus binary. A
    complete list is available by calling ./LC3plus -h.

    To call encoder+decoder at the same time
        ./LC3plus INPUT.wav OUTPUT.wav BITRATE

    To call encoder only
        ./LC3plus -E INPUT.wav OUTPUT.bin BITRATE

    To call decoder only
        ./LC3plus -D INPUT.bin OUTPUT.wav

    To specify output bits per sample
        ./LC3plus -bps NUM INPUT.wav OUTPUT.wav BITRATE
    where NUM is either 16 (default) or 24.

    To specify bitrate switching file instead of fixed bitrate
        ./LC3plus -swf FILE INPUT.wav OUTPUT.wav BITRATE
    where FILE is a binary file containing the bitrate as a
    sequence of 64-bit values.

    To specify delay compensation
        ./LC3plus -dc DC INPUT.wav OUTPUT.wav BITRATE
    where DC is 0, 1 (default) or 2 meaning that:
        0: Don't use delay compensation
        1: Compensate delay in decoder (default)
        2: Split delay equally in encoder and decoder

    To disable frame counter (quiet mode)
        ./LC3plus -q INPUT.wav OUTPUT.wav BITRATE

    To activate verbose mode (print switching commands)
        ./LC3plus -v INPUT.wav OUTPUT.wav BITRATE

    To use the G192 bitstream format
        ./LC3plus -E -formatG192 INPUT.wav OUTPUT.g192 BITRATE
    Note that an additional file OUTPUT.cfg will be created containing
    the decoder information. To specify the configuration file,
    the flag -cfgG192 FILE can be used where FILE is the path to the
    configuration file.
    Note that the same flags (-formatG192 and -cfgG192) shall be used for decoding,
    if the G192 format was chosen for encoding.

    To call decoder with frame loss simulations
        ./LC3plus -D -epf patternfile INPUT.bin OUTPUT.wav
        where the patternfile is a binary file containing a sequence of
        16-bit values, non-zero values indicating a frame loss

    To write error detection pattern (from arithmetic decoder) into a 16-bit binary file
        ./LC3plus -D -edf patternfile INPUT.bin OUTPUT.wav
        where the patternfile is a binary file containing a sequence of
        16-bit values, non-zero values indicating a detected frame loss

    The high-resolution mode is mandatory for usage with 96 kHz sampling frequency 
    and can also be used with 48 kHz sampling frequency. 
    To enable the high-resolution mode call

        ./LC3plus -hrmode INPUT.wav OUTPUT.wav BITRATE

    To activate error protection mode call

        ./LC3plus -epmode NUM INPUT.wav OUTPUT.wav BITRATE

        where NUM is the error protection mode using the LC3plus channel coder.
        NUM must be in the range of 0 - 4 where:
              - 0 = error protection disabled
              - 1 = minimum error protection, bit error detection only
              - 2 = moderate error protection, e.g. up to 3 bits correctable for
                    a 40 bytes frame
              - 3 = strong error protection, e.g. up to 9 bits correctable for
                    a 40 bytes frame
              - 4 = maximum error protection, e.g. up to 18 bits correctable for
                    a 40 bytes frame
              -NUM can also be a binary switching file holding the values for the
                    current epmode.

    The actual number of correctable bits depends on the bit error distribution
    in a frame where the codes are strongest when bit errors accumulate (e.g.
    burst errors) and weakest when bit errors are uniformly distributed.

    Error protection can be in encoder only mode:
         ./LC3plus -E -epmode NUM INPUT.wav OUTPUT.bin BITRATE

    The decoder needs the information, if error protection is enabled
         ./LC3plus -D -epmode NUM INPUT.bin OUTPUT.wav
         with NUM one of the following: 0 = off, 1 = on. The channel-decoder
         identifies the error protection mode by itself

    The parameter bitrate, bandwidth and epmode also allow the usage of switching
    files instead of a fixed number.
    If the switching file is shorter than the
    input it is looped from the beginning. The switching files contain data stored
    as little-endian 64-bit values.

    For the bitrate, each value represents a frame's bitrate.
    For the bandwidth, each value represents a frame's cutoff bandwidth in Hz.
    For error protection, each value represents a frame's epmode multiplied with
    100, i.e. epmode = value / 100.
    Switching files can be created with 'gen_rate_profile' available from
    ITU-T G191.

    For supported bitrate configurations refer to ETSI TS 103 634 Table 5.1.

    When the number of bytes per frame is smaller than or equal to 160, left
    and right channel are coded jointly by the LC3plus channel coder. Above this
    threshold, left and right channel are channel encoded separately. The
    cross-over bitrates are given in the table below.

    Frame Duration [ms]|  Max. Bitrate for Joint Channel Coding [kbps]
    ===================================================================
    1.25               |  1024
    -------------------------------------------------------------------
    2.5                |  512
    -------------------------------------------------------------------
    5                  |  256
    -------------------------------------------------------------------
    7.5                |  170.667
    -------------------------------------------------------------------
    10                 |  128
    -------------------------------------------------------------------

    The software package also includes a Channel Coder Converter (CCC) which converts unprotected
    and protected LC3plus payloads within each other.

    Building:
    On Unix:
        - make ccConvert

    On Windows:
        - Go to src/fixed_point/msvc_ccc
        - Open up solution file ccConvert.sln and build it
        - Standard config executable path and name ".\Win32\Release\ccConvert.exe"

    Usage:
    The CCC can be used as follows:
    Pack mode (convert to protected):     ./ccConvert -pack gross_bytes ep_mode in.lc3 out.lc3
    Unpack mode (convert to unprotected): ./ccConvert -unpack in.lc3 out.lc3

    where gross_bytes is the gross number of bytes available and ep_mode the error protection
    mode which is one of 0,1,2,3,4.
    For G.192 format, the filename extension should be ".g192".


Binary File Format
------------------
    Note: The binary file format is intended only for testing purposes. It is
    not standardized and may change in the future.
    The file starts with a config header structure as follows.

    Field            Bytes   Content
    ------------------------------------------------------------------
    file_id           2       file identifier, value 0xcc1c
    header_size       2       total config header size in bytes
    samplingrate      2       sampling frequency / 100
    bitrate           2       bitrate / 100
    channels          2       number of channels
    frame_ms          2       Frame duration in ms * 100
    epmode            2       error protection, 0: off, non-zero: on
    signal_len        2       input signal length in samples
    signal_len_red    2       signal_len >> 16
    hrmode            2       high resolution mode, 0: off, 1: on

    All fields are stored in little-endian byte order. The config header could
    be extended in the future so readers should seek to header_size to skip any
    unknown fields.
    The header is immediately followed by a series of coded audio frames, where
    each frame consists of a two-byte frame length information and the current
    coded frame.
    For multichannel audio, the coded audio frames for each channel are concatenated and the 
    length refers to the sum of the frame lengths.

