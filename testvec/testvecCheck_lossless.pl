#!/usr/bin/perl -w

#******************************************************************************
#                        ETSI TS 103 634 V1.7.1                               *
#              Low Complexity Communication Codec Plus (LC3plus)              *
#                                                                             *
# Copyright licence is solely granted through ETSI Intellectual Property      *
# Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
# estoppel or otherwise.                                                      *
#*****************************************************************************/

# ============================================================
# LC3plus ETSI Lossless Testvectors script V0.0.1
# ============================================================
#
# Verifies bit-exactness of the LC3plus lossless mode encoder/decoder
# for all supported sampling rates and both prioritization rules
# (b_relative = 0 and b_relative = 1, set via the -rel_prio CLI flag).
#
# Covers:
#   - VBR lossless (no bitrate restriction) — exercises bitstream
#     construction with b_relative both written and read correctly.
#   - CBR lossless with padding — exercises the prioritization
#     iteration in compute_resbits_priority() for both rules.
#
# Usage analogous to testvecCheck.pl:
#   ./testvecCheck_lossless.pl -fixed [-create] [-clean] [-log] [-verbose]
#
# Note: 24-bit testing requires dedicated 24-bit input WAVs with
# low-level noise (-100 dBFS) — see Spec clause 7.3.8. Not run here
# until such WAVs are provided in input/.
##################################################
my $EXE_TST_FX = './../src/fixed_point/LC3plus_HR';
my $md5_bin_fx = './md5_bin_fx_lossless.txt';
my $md5_dec_fx = './md5_dec_fx_lossless.txt';
my $MY_MD5  = 'md5sum';
##################################################

use strict;
use File::Basename;
use File::Path 'rmtree';
use POSIX;

my $VERSION = 'V0.0.1';
my $timestamp = strftime "%m_%d_%Y_%H_%M_%S", localtime;
my $tmp_folder = "lc3plus_lossless_testvectors".$timestamp;
my $inputFile = "./input/thetest";
my $output_folder_stream_tst = $tmp_folder."/bitstream_tst";
my $output_folder_decoded_tst = $tmp_folder."/decoded_tst";
my $report = "lc3plus_lossless_testvectors_report_".$timestamp.".txt";
my $fh;
my $quiet = '>/dev/null 2>&1';
my $EXE_TST;
my $md5_bin;
my $md5_dec;

my $testvectors_fail = 0;

# -------- Lossless test matrix --------
# Lossless mode supports fs in {44100, 48000, 96000, 192000} per Tab. 5.2a.
# 192 kHz is included only if input/thetest192.wav is provided.
# Each config is tested with -rel_prio 0 and -rel_prio 1.
# Two modes per fs:
#   (a) VBR  (no bitrate argument): variable, encoder chooses
#   (b) CBR  (-padding + bitrate):  constrained, exercises the priority iteration
# Bitrates for CBR are chosen above the HR fallback rate (Tab. 5.2) to keep
# lossless mode active most of the time but force priority decisions on
# difficult frames.

# fs => [list of CBR bitrates to test]; VBR is always tested.
my %FS_BR = (
    44100 => [400000],     # ~ ~50 kBytes/s, well above fallback rate at 44.1 kHz
    48000 => [500000],     # bit above fallback at 48 kHz
    96000 => [1000000],    # bit above fallback at 96 kHz
);

# Optionally enable when input/thetest192.wav exists:
if (-e $inputFile."192.wav")
{
    $FS_BR{192000} = [2500000];
}

my @REL_PRIO = (0, 1);

# Set MD5 command according to OS
getOS();

# Args
my ($fixed, $float, $create, $test, $clean, $log) = (0, 0, 0, 0, 0, 0);
getArgs(\$EXE_TST, \$md5_bin, \$md5_dec, \$fixed, \$float, \$create, \$clean, \$log, \$quiet, \@ARGV);

$test = 1;
if ($create) { $test = 0; }

checkExe($EXE_TST, $md5_bin, $md5_dec);
checkMD5($MY_MD5, $quiet);
createDirs($tmp_folder, $output_folder_stream_tst, $output_folder_decoded_tst);
checkInputs($inputFile);

if ($log)
{
    open($fh, '>', $report) or die "Could not open file '$report' $!";
}

my $md5stream;
my $md5decoded;
if ($create)
{
    open($md5stream, '>', $md5_bin) or die "Could not open file '$md5_bin' for write: $!";
    open($md5decoded, '>', $md5_dec) or die "Could not open file '$md5_dec' for write: $!";
}
else
{
    open($md5stream, '<', $md5_bin) or die "Could not open file '$md5_bin' for read: $!";
    open($md5decoded, '<', $md5_dec) or die "Could not open file '$md5_dec' for read: $!";
}

print("Lossless testvectors script started...\n");

foreach my $sr (sort { $a <=> $b } keys %FS_BR)
{
    my $input = "$inputFile".(floor($sr/1000)).".wav";
    next unless -e $input;
    my $base = basename($input);
    $base =~ s/\.[^.]+$//;

    # --- (a) VBR lossless: -lossless [no bitrate] -rel_prio X ---
    foreach my $rp (@REL_PRIO)
    {
        my $tag = $base."_lossless_vbr_rp${rp}";
        my $stream  = "$output_folder_stream_tst/$tag.lc3plus";
        my $decoded = "$output_folder_decoded_tst/$tag.wav";
        my $cmd_enc = "$EXE_TST -E -lossless -rel_prio $rp $input $stream $quiet";
        my $cmd_dec = "$EXE_TST -D $stream $decoded $quiet";
        run_pair($cmd_enc, $cmd_dec, $stream, $decoded, $md5stream, $md5decoded, $log, $fh, \$testvectors_fail, $create);
    }

    # --- (b) CBR lossless with padding for each configured bitrate ---
    foreach my $br (@{$FS_BR{$sr}})
    {
        foreach my $rp (@REL_PRIO)
        {
            my $tag = $base."_lossless_cbr${br}_rp${rp}";
            my $stream  = "$output_folder_stream_tst/$tag.lc3plus";
            my $decoded = "$output_folder_decoded_tst/$tag.wav";
            my $cmd_enc = "$EXE_TST -E -lossless -padding -rel_prio $rp $input $stream $br $quiet";
            my $cmd_dec = "$EXE_TST -D $stream $decoded $quiet";
            run_pair($cmd_enc, $cmd_dec, $stream, $decoded, $md5stream, $md5decoded, $log, $fh, \$testvectors_fail, $create);
        }
    }
}

close $md5stream;
close $md5decoded;

if ($clean) { cleanup($tmp_folder); }

if ($test)
{
    my $result = $testvectors_fail ? "NOT passed" : "passed";
    print("\nLossless testvector check $result!\n");
    if ($log) { print $fh "\nLossless testvector check $result!\n"; close $fh; }
}

if ($log) { print("See logfile: $report\n"); }
exit($testvectors_fail ? 1 : 0);


# ----------------- helpers -----------------

sub run_pair
{
    my ($cmd_enc, $cmd_dec, $stream, $decoded, $md5s, $md5d, $log, $fh, $fail_ref, $create) = @_;

    system($cmd_enc);
    if ($log) { print $fh "$cmd_enc\n"; }
    system($cmd_dec);
    if ($log) { print $fh "$cmd_dec\n"; }

    if ($create)
    {
        record_md5($stream,  $md5s);
        record_md5($decoded, $md5d);
    }
    else
    {
        compare_md5($stream,  $md5s, $log, $fh, $fail_ref);
        compare_md5($decoded, $md5d, $log, $fh, $fail_ref);
    }
}

sub record_md5
{
    my ($file, $fh_out) = @_;
    my $md5 = qx($MY_MD5 $file);
    my ($hash) = $md5 =~ /([a-f0-9]{32})/i;
    print $fh_out basename($file).":".$hash."\n";
}

sub compare_md5
{
    my ($file, $fh_in, $log, $fh_log, $fail_ref) = @_;
    my $line = <$fh_in>;
    chomp($line) if defined $line;
    my @kv = split(/:/, $line // "");
    my $base = basename($file);
    if (!defined $kv[0] || $kv[0] ne $base)
    {
        print("Error: hash list out of sync (expected $base, got '".($kv[0]//"")."')\n");
        ${$fail_ref} = 1;
        return;
    }
    my $md5 = qx($MY_MD5 $file);
    my ($hash) = $md5 =~ /([a-f0-9]{32})/i;
    if ($log)
    {
        print $fh_log "Check $base: ref=$kv[1] got=$hash\n";
    }
    if (!defined $hash || $kv[1] ne $hash)
    {
        ${$fail_ref} = 1;
        print("MISMATCH: $base (ref=$kv[1] got=".($hash//"<none>").")\n");
    }
}

sub cleanup { my ($dir) = @_; rmtree($dir); }

sub checkMD5
{
    my ($cmd, $quiet) = @_;
    my $ret = system("$cmd $0 $quiet");
    if ($ret != 0)
    {
        print("Error: cannot find md5 command: $cmd\n");
        exit(1);
    }
}

sub checkInputs
{
    my ($base) = @_;
    my @needed = (44, 48, 96);
    foreach my $sr (@needed)
    {
        my $f = $base.$sr.".wav";
        if (! -e $f)
        {
            print("Cannot find input file $f. Please provide it in ./input/\n");
            exit(1);
        }
    }
}

sub createDirs
{
    my ($tmp, $a, $b) = @_;
    mkdir($tmp); mkdir($a); mkdir($b);
}

sub getArgs
{
    my ($EXE_TST, $md5_bin, $md5_dec, $fixed, $float, $create, $clean, $log, $quiet, $args) = @_;
    my @arg = @{$args};
    my $help = 0;
    foreach my $a (@arg)
    {
        if    ($a eq '-create')  { ${$create} = 1; }
        elsif ($a eq '-fixed')   { ${$fixed}  = 1; }
        elsif ($a eq '-float')   { ${$float}  = 1; }
        elsif ($a eq '-clean')   { ${$clean}  = 1; }
        elsif ($a eq '-log')     { ${$log}    = 1; }
        elsif ($a eq '-verbose') { ${$quiet} = ''; }
        elsif ($a =~ /^-+h(elp)?$/i) { $help = 1; }
        else { print("Unknown parameter: $a\n"); printUsage(); }
    }
    if ($help) { printUsage(); }
    if (${$float})
    {
        print("Lossless mode is only available in the fixed-point reference. Use -fixed.\n"); exit(1);
    }
    if (!${$fixed})
    {
        print("Please select -fixed (lossless is fixed-point only).\n"); printUsage();
    }
    ${$EXE_TST} = $EXE_TST_FX; ${$md5_bin} = $md5_bin_fx; ${$md5_dec} = $md5_dec_fx;
}

sub checkExe
{
    my ($exe, $a, $b) = @_;
    if (! -e $exe) { print("ERROR: cannot find executable $exe\n"); exit(1); }
    if (! -e $a)   { print("ERROR: cannot find MD5 reference: $a\n"); exit(1); }
    if (! -e $b)   { print("ERROR: cannot find MD5 reference: $b\n"); exit(1); }
}

sub getOS
{
    my $os = $^O;
    if ($os eq "darwin") { $MY_MD5 = "md5"; }
    else                 { $MY_MD5 = "md5sum"; }
}

sub printUsage
{
    print("\nLC3plus ETSI Lossless Testvectors script $VERSION\n");
    print("Verifies bit-exactness of lossless mode for both b_relative values.\n");
    print("Usage: testvecCheck_lossless.pl -fixed [-create] [-clean] [-log] [-verbose]\n");
    print("  -fixed   : test fixed-point executable (lossless is fixed-point only)\n");
    print("  -create  : regenerate reference MD5 lists (do NOT use for verification)\n");
    print("  -clean   : remove temporary directory after run\n");
    print("  -log     : write detailed log to file\n");
    print("  -verbose : show executable output\n");
    exit(0);
}
