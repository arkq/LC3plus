LC3plus Lossless Mode Testvector Verification
=============================================

This package lets you verify a third-party LC3plus lossless mode
implementation against the ETSI reference. It covers both residual
prioritization rules (b_relative = 0 and 1).

Files
-----
  testvecCheck_lossless.pl  - Perl driver script
  md5_bin_fx_lossless.txt   - Reference MD5 hashes of encoded streams
  md5_dec_fx_lossless.txt   - Reference MD5 hashes of decoded WAVs
  Readme_lossless.txt       - This file

Lossless mode is implemented in the fixed-point reference only;
the floating-point build does not contain a lossless path.

Test matrix
-----------
For each sampling rate fs in { 44 100, 48 000, 96 000 } Hz with a
16 bit-per-sample mono input from input/, four configurations are
verified:

  1. VBR lossless, b_relative = 0   (-lossless -rel_prio 0)
  2. VBR lossless, b_relative = 1   (-lossless -rel_prio 1)
  3. CBR lossless, b_relative = 0   (-lossless -padding -rel_prio 0 <bitrate>)
  4. CBR lossless, b_relative = 1   (-lossless -padding -rel_prio 1 <bitrate>)

Total: 12 test cases.

Configurations 1+2 exercise the bitstream construction (verifies that
b_relative is written/read correctly and that the rest of the lossless
encoder is bit-exact regardless of which prioritization rule is signaled).

Configurations 3+4 exercise the prioritization algorithm in
compute_resbits_priority(), because the CBR bitrate forces priority
decisions on frames that do not fit. Bitrates per fs are chosen above
the HR fallback rate (Table 5.2 of TS 103 634) so the encoder stays in
lossless mode most of the time but has to drop some LSB planes on
demanding frames.

How to run
----------
1. Build the fixed-point executable:
     cd ../src/fixed_point && cmake -B build -S . && cmake --build build -j

2. Verify your implementation against the reference hashes:
     cd ../../testvec
     ./testvecCheck_lossless.pl -fixed -log

   The script encodes and decodes each test configuration with your
   build and compares the resulting MD5 hashes against the reference
   values in md5_bin_fx_lossless.txt and md5_dec_fx_lossless.txt.
   It exits with code 0 on full pass, non-zero on any mismatch.
   The log file lc3plus_lossless_testvectors_report_<timestamp>.txt
   contains per-case reference vs. observed hashes.
