/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                        Software Version V1.4.10ETSI                         *
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

-/src/fixed_point     : fixed-point source files, makefile and Visual Studio
                        solution
-/src/floating_point  : floating-point source files, makefile and Visual Studio
                        solution
-/conformance         : conformance script with example configuration file and
                        Readme
-/testvec             : testvector package with script. Verification that LC3plus
                        build on your system operates as expected by comparing to
                        pre-calculated MD5 hashes
-/tools               : LC3plus helper tools and examples

Please refer to the respective Readme files for more information.


Features
--------
    - Supported sampling rates: 8 kHz, 16 kHz, 24 kHz, 32 kHz, 44.1 kHz,
                                      48 kHz, 96 kHz
    - Supported bit rates: 20 bytes ... 400 bytes per frame per channel
    - Multichannel support by multi-mono coding
    - Support of audio sample depth: 16 bits and 24 bits
    - Frame duration of 10 ms, 5 ms and 2.5 ms
    - Supported bit rates:
        - frame duration 10 ms:  16 kbps ... 320 kbps per channel
        - frame duration 5 ms:   32 kbps ... 640 kbps per channel
        - frame duration 2.5 ms: 64 kbps ... 1280 kbps per channel
    - High resolution audio mode for 48 kHz and 96 kHz and higher SNR.
    -- Supported bit rates for 48 kHz
        - frame duration 10 ms:  124.8 ... 500 kbps per channel
        - frame duration 5 ms:   148.8 ... 600 kbps per channel
        - frame duration 2.5 ms: 172.8 ... 672 kbps per channel
    -- Supported bit rates for 96 kHz
        - frame duration 10 ms:  149.6 ... 500 kbps per channel
        - frame duration 5 ms:   174.4 ... 600 kbps per channe
        - frame duration 2.5 ms: 198.4 ... 672 kbps per channel

    - Packet loss concealment as defined in ETSI TS 103 634

The fixed-point source code currently has limited support for the
following functions:
    - High resolution mode

The floating-point source code currently has limited support for the
following functions:
    - Packet loss concealment
    - Error protection (LC3plus channel coder)
    - Partial concealment


Changelog
---------
  Latest non-bitexact changes are encapsulated in defines listed in defines.h
  for review.

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
corresponding functions via macros labeled LC3_XYZ. By default they use single
precision floating point format to represent data and perform mathematical
operations. Double precision floating point format can be acitivated by enabling
the define LC3_DOUBLE_PRECISION.


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


    The high-resolution mode is only available in the floating-point source code.
    It is mandatory for usage with 96 kHz sampling frequency and can also be used
    with 48 kHz sampling frequency. To enable the high-resolution mode call

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

    If error protection is enabled, the compressed LC3plus frame is reduced
    depending on the robustness according to the following table:

    Frame Duration| Min. Bitrate  |  Max. Bitrate
    Channel Mode  | [kbps]        |  [kbps]
    ==============================================
    2.5 ms        |  128          |  960
    Mono          |               |
    ----------------------------------------------
    5 ms          |  64           |  480
    Mono          |               |
    ----------------------------------------------
    10 ms         |  32           |  240
    Mono          |               |
    ----------------------------------------------
    2.5 ms        |  256          |  1920
    Stereo        |               |
    ----------------------------------------------
    5 ms          |  128          |  960
    Stereo        |               |
    ----------------------------------------------
    10 ms         |  64           |  480
    Stereo        |               |

    While the number of bytes per frame is smaller than or equal to 160, left
    and right channel are coded jointly by the LC3plus channel coder. Above this
    threshold, left and right channel are channel encoded separately. The
    cross-over bitrates are given in the table below.

    Frame Duration [ms]|  Max. Bitrate for Joint Channel Coding [kbps]
    ===================================================================
    2.5                |  512
    -------------------------------------------------------------------
    5                  |  256
    -------------------------------------------------------------------
    10                 |  128
    -------------------------------------------------------------------

    The minimum and maximum allowed bitrate per channel depends on the frame duration
    since the minimum and maximum frame data size of 20 bytes ... 400 bytes is kept
    independently of the frame duration. Therefore, the bitrate per channel is
    restricted according to the table below.

    Frame Duration [ms]|  Min. Bitrate [kbps]  |  Max. Bitrate [kbps]
    ===================================================================
    2.5                |  64                   |  1280
    -------------------------------------------------------------------
    5                  |  32                   |  640
    -------------------------------------------------------------------
    10                 |  16                   |  320
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

    Field         Bytes   Content
    ------------------------------------------------------------
    file_id       2       file identifier, value 0xcc1c
    header_size   2       total config header size in bytes
    samplingrate  2       sampling frequency / 100
    bitrate       2       bitrate / 100
    channels      2       number of channels
    frame_ms      2       Frame duration in ms * 100
    epmode        2       error protection, 0: off, non-zero: on
    signal_len    4       input signal length in samples

    Note: the floating-point implementation adds an additional field
    hrmode        2       high resolution mode on/off

    All fields are stored in little-endian byte order. The config header could
    be extended in the future so readers should seek to header_size to skip any
    unknown fields.
    The header is immediately followed by a series of coded audio frames, where
    each frame consists of a two-byte frame length information and the current
    coded frame.

