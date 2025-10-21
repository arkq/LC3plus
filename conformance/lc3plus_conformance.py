#!/usr/bin/env python3

# Version 1.2.6

#******************************************************************************
#                        ETSI TS 103 634 V1.6.1                               *
#              Low Complexity Communication Codec Plus (LC3plus)              *
#                                                                             *
# Copyright licence is solely granted through ETSI Intellectual Property      *
# Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
# estoppel or otherwise.                                                      *
#*****************************************************************************/


import argparse
import collections
import concurrent.futures
import configparser
import datetime
import filecmp
import hashlib
import io
import itertools
import logging
import math
import os
import pathlib
import re
import shlex
import shutil
import struct
import subprocess
import sys
import uuid
import wave
import zipfile
import locale
try:
    import numpy
except ImportError:
    sys.exit('Numpy missing! Try running "pip3 install numpy".')
import time


# INTERNALS
DEBUG_SETTING = 0

# convenience
Executor = concurrent.futures.ThreadPoolExecutor
TestEnv  = collections.namedtuple('TestEnv', 'profile test config test_dir workers opid opid_log wav_only cmp_only')
WavInfo  = collections.namedtuple('WavInfo', 'ch rate frames')

# constants
MAX_DELAY = 322
MAX_SAMPLES_PER_FRAME = 480
SAMPLERATES = [8000, 16000, 24000, 32000, 44100, 48000, 96000]
FRAME_SIZES = [1.25, 2.5, 5.0, 7.5, 10.0]
CHANNELS = [1, 2]
SQAM_URL    = 'https://qc.ebu.io/testmaterial/523/1/download/'
SQAM_SHA256 = '7d6fcd0fc42354637291792534b61bf129612f221f8efef97b62e8942a8686aa'
SOX_URL     = 'https://sourceforge.net/projects/sox/files/sox/14.4.2/sox-14.4.2-win32.zip'
SOX_SHA256  = '8072cc147cf1a3b3713b8b97d6844bb9389e211ab9e1101e432193fad6ae6662'
SOX_EXE     = pathlib.Path('SoX/sox-14.4.2/sox.exe')
INF         = float('inf')

PC_BER_FILES = {
  1.25  : "etc/pc_ber_2.5ms.dat",
  2.5  : "etc/pc_ber_2.5ms.dat",
  5.0  : "etc/pc_ber_5ms.dat",
  7.5  : "etc/pc_ber_10ms.dat",
  10.0 : "etc/pc_ber_10ms.dat"
}

PLC_BURST_FILES = {
  1.25  : "etc/plc_bfer_1.25ms_eid.dat",
  2.5  : "etc/plc_bfer_2.5ms_eid.dat",
  5.0  : "etc/plc_bfer_5ms_eid.dat",
  7.5  : "etc/plc_bfer_7.5ms_eid.dat",
  10.0 : "etc/plc_bfer_10ms_eid.dat"
}

PLC_RANDOM_FILES = {
  1.25  : "etc/plc_fer_1.25ms_eid.dat",
  2.5  : "etc/plc_fer_2.5ms_eid.dat",
  5.0  : "etc/plc_fer_5ms_eid.dat",
  7.5  : "etc/plc_fer_7.5ms_eid.dat",
  10.0 : "etc/plc_fer_10ms_eid.dat"
}

MIN_HRMODE_BYTES = {# (frame_size, sampling_rate): min_bytes_per_frame
    (10,48000): 156,
    (10,96000): 187,
    (7.5,48000): 117,
    (7.5,96000): 140,
    (5,48000): 93,
    (5,96000): 109,
    (2.5,48000): 54,
    (2.5,96000): 62,
    (1.25,48000): 54,
    (1.25,96000): 62,
}

# test items
ITEM_DIR = pathlib.Path('test_items')
ITEMS = {  # start, frag, SQAM name
    'ABBA'                : ( 7,  8, '69.flac'),
    'Castanets'           : ( 0,  8, '27.flac'),
    'Eddie_Rabbitt'       : ( 0,  8, '70.flac'),
    'Female_Speech_German': ( 0,  8, '53.flac'),
    'Glockenspiel'        : ( 0, 10, '35.flac'),
    'Piano_Schubert'      : ( 0,  8, '60.flac'),
    'Violoncello'         : ( 0, 10, '10.flac'),
    'Harpsichord'         : (39,  9, '40.flac'),
    'Male_Speech_English' : ( 0,  8, '50.flac'),
}
ITEMS_PLC = ['ABBA', 'Castanets', 'Female_Speech_German', 'Harpsichord' , 'Male_Speech_English']
ITEM_LOW_PASS   = 'White_Noise_LP20'
ITEM_BAND_LIMIT = 'Female_Speech_German'
ITEM_LFE        = 'Sine_Sweep'

# sampling rate, band widths, bytes / frame
BAND_LIMITS = {
    48000: ([4000, 8000, 12000, 16000], 115),
    32000: ([4000, 8000, 12000], 80),
    24000: ([4000, 8000], 60),
    16000: ([4000], 40),
}
BAND_WIDTHS = {
    48000: [4000, 8000, 12000, 16000, 20000],
    32000: [4000, 8000, 12000, 16000],
    24000: [4000, 8000, 12000],
    16000: [4000, 8000],
}

# config default values
TESTS = [
    'sqam',
    'band_limiting',
    'low_pass',
    'bitrate_switching',
    'bandwidth_switching',
    'plc',
    'plc_burst',
    'pc',
    'ep_correctable',
    'ep_non_correctable',
    'ep_mode_switching',
    'ep_combined',
    'ep_combined_nc',
    'lfe',
]
TEST_MODES = ['encode', 'encdec', 'decode']
DEFAULTS_GLOBAL = {
    'option_bandwidth' : '',
    'option_ep_debug'  : '',
    'option_ep_mode'   : '',
    'option_plc_mode'  : '',
    'peaq_bin'         : '',
    'peaq_odg_regex'   : '',
    'ref_encoder'      : '',
    'ref_decoder'      : '',
    'tst_encoder'      : '',
    'tst_decoder'      : '',
    'config_bitrate_unit' : 'bps',
    'config_bitrate_format' : 'list',
    'delay'                 : '0',
    'hrmode'                : '0',
}
DEFAULTS_TEST = {'configs': []}
for test in TESTS:
    DEFAULTS_TEST['test_' + test] = False
for test, mode in itertools.product(TESTS, TEST_MODES):
    DEFAULTS_TEST['{}_{}_eng_threshold'.format(test, mode)] = [70,90]
    DEFAULTS_TEST['{}_{}_mld_threshold'.format(test, mode)] = 4
    DEFAULTS_TEST['{}_{}_odg_threshold'.format(test, mode)] = [0.06, 0.08, 0.12]
    DEFAULTS_TEST['{}_{}_rms_threshold'.format(test, mode)] = 14
    DEFAULTS_TEST['{}_{}_mad_threshold'.format(test, mode)] = 0.00148
    DEFAULTS_TEST['{}_{}_metric'.format(test, mode)]        = 'rms'
METRIC_DEFAULTS = {
    'low_pass_encode_metric'           : 'eng',
    'low_pass_encdec_metric'           : 'eng',
    'plc_decode_metric'                : 'mld',
    'plc_burst_decode_metric'          : 'mld',
    'pc_decode_metric'                 : 'mld',
    'ep_non_correctable_decode_metric' : 'mld',
    'ep_non_correctable_encdec_metric' : 'mld',
    'ep_non_correctable_encode_metric' : 'mld',
    'ep_combined_nc_decode_metric'     : 'mld',
    'ep_combined_nc_encdec_metric'     : 'mld',
    'ep_combined_nc_encode_metric'     : 'mld',
    'ep_combined_decode_metric'        : 'odg',
    'ep_combined_encode_metric'        : 'rms',
    'ep_combined_encdec_metric'        : 'mld',
    'ep_combined_encode_metric'        : 'mld',
    'ep_correctable_encode_metric'     : 'odg',
    'ep_correctable_decode_metric'     : 'rms',
    'ep_mode_switching_encode_metric'  : 'odg',
    'ep_mode_switching_decode_metric'  : 'rms',
    'lfe_encode_metric'                : 'mld',
    'sqam_encode_metric'               : 'odg',
    'sqam_encdec_metric'               : 'odg',
    'band_limiting_encode_metric'      : 'odg',
    'band_limiting_encdec_metric'      : 'odg',
    'bitrate_switching_encode_metric'  : 'odg',
    'bitrate_switching_encdec_metric'  : 'odg',
    'bandwidth_switching_encode_metric': 'odg',
    'bandwidth_switching_encdec_metric': 'odg',
    'ep_correctable_encode_metric'     : 'odg',
    'ep_correctable_encdec_metric'     : 'odg',
    'ep_mode_switching_encode_metric'  : 'odg',
    'ep_mode_switching_encdec_metric'  : 'odg',
}
DEFAULTS_TEST.update(METRIC_DEFAULTS)

# html output stuff
LABEL = {
    'sqam'               : 'SQAM',
    'band_limiting'      : 'Band Limitation',
    'low_pass'           : 'Low Pass',
    'bitrate_switching'  : 'Bitrate Switching',
    'bandwidth_switching': 'Bandwidth Switching',
    'plc'                : 'Packet Loss Concealment - 10% Random Loss',
    'plc_burst'          : 'Packet Loss Concealment - 10% Burst Loss',
    'pc'                 : 'Partial Concealment',
    'ep_correctable'     : 'Channel Coder for Correctable Frames',
    'ep_non_correctable' : 'Channel Coder for Non-Correctable Frames',
    'ep_mode_switching'  : 'Error protection mode switching',
    'ep_combined'        : 'Combined Channel Coder for Correctable Frames',
    'ep_combined_nc'     : 'Combined Channel Coder for Non-Correctable Frames',
    'lfe'                : 'Low Frequency Effects',
}
HEADER_ALL = ['Mode', 'Item', 'Frame Size', 'Samplerate', 'Bitrate']
HEADER_EP  = ['EP Mode']
HEADER_EPD = ['BFI', 'EPMR', 'ER']
HEADER_EPD_BFI = ['BFI']
HEADER_BL  = ['Bandwidth']
HEADER_METRIC = {
    'rms': ['Max. Abs. Diff', 'RMS [dB]', 'RMS Reached [bits]'],
    'odg': ['ODG<sub>ref</sub>', '&#916<sub>ODG</sub>'],
    'mld': ['MLD'],
    'eng': ['E<sub>diff</sub> [dB]'],
    None : []
}
HTML_HEAD = ('<!DOCTYPE html><head><meta charset="UTF-8"><title>{title} Report</title>'
             '<style>{style}</style></head>\n<body><h2>Conformance for "{title}" {state}!</h2>\n')
HTML_DIV = '<div><h3>{label} - {percent}%</h3>\n'
STYLE = ('body {font-family:sans-serif; color:#f8f8f2; background-color:#272822; font-size:80%} div {border:1px solid '
         '#8f908a; border-radius:4px; overflow:hidden; display:table; margin-left:30px; margin-bottom:30px} h2 {text-a'
         'lign:left; margin-left:30px} h3 {text-align:left; margin:4px} table {border-spacing:0px; width:100%} th {pad'
         'ding:4px; border-top:1px solid #8f908a} td {padding:4px} tr:nth-child(even) {background-color:rgba(255,255,2'
         '55,0.1)} td.pass {background-color:rgba(0,192,255,0.4)} td.fail {background-color:rgba(255,0,0,0.4)} td.warn'
         '{background-color:rgba(214,137,16,0.4)} a:link {color: white;text-decoration: none;} a:visited {color: '
         'white} a:hover {color: white;font-weight:bold;} a:active {color: rgb(155, 155, 155)}')

# Constants for bitstream editing
G192_SYNC_WORD        = int('0xCC1C', 16)
G192_GOOD_FRAME       = int('0x6B21', 16)
G192_BAD_FRAME        = int('0x6B20', 16)
G192_ONE              = int('0x0081', 16)
G192_ZERO             = int('0x007F', 16)
G192_REDUNDANCY_FRAME = int('0x6B22', 16)

# List of conditions for silencing frames to ensure stable conformance testing
# when different LTPF behaviour is observed between fl and fx implementations
LTPF_BITSTREAM_CONDITIONS = [
    {'item': 'Male_Speech_English',  "channels": 2, 'sampling_rate': 16000, 'frame_size': 1.25, 'bitrate': 512000, 'ep_mode': 4,'critical_frames': list(range(5069, 5080))},
    {'item': 'Male_Speech_English',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': 256000, 'ep_mode': 3,'critical_frames': list(range(5069, 5080))},
    {'item': 'Male_Speech_English',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': 256000, 'ep_mode': '1-4','critical_frames': list(range(2238, 2249))},
    {'item': 'Male_Speech_English',  "channels": 2, 'sampling_rate': 32000, 'frame_size': 1.25, 'bitrate': 512000, 'ep_mode': 4,'critical_frames': list(range(4477, 4490))},
    {'item': 'Female_Speech_German',  "channels": 1, 'sampling_rate': 16000, 'frame_size': 1.25, 'bitrate': 128000, 'ep_mode': 0,'critical_frames': list(range(4179, 4190))},
    {'item': 'Female_Speech_German',  "channels": 1, 'sampling_rate': 16000, 'frame_size': 1.25, 'bitrate': '128000-128000', 'ep_mode': 0,'critical_frames': list(range(4179, 4190))},
    {'item': 'Female_Speech_German',  "channels": 1, 'sampling_rate': 16000, 'frame_size': 1.25, 'bitrate': 256000, 'ep_mode': 4,'critical_frames': list(range(4179, 4190))},
    {'item': 'Female_Speech_German',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': 128000, 'ep_mode': 0,'critical_frames': list(range(2034, 2045))},
    {'item': 'Female_Speech_German',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': '128000-128000', 'ep_mode': 0,'critical_frames': list(range(2034, 2045))},
    {'item': 'Female_Speech_German',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': 256000, 'ep_mode': 4,'critical_frames': list(range(2034, 2045))},
    {'item': 'Female_Speech_German',  "channels": 1, 'sampling_rate': 48000, 'frame_size': 1.25, 'bitrate': 256000, 'ep_mode': '1-4','critical_frames': list(range(5690, 5701))},
    {'item': 'Female_Speech_German',  "channels": 1, 'sampling_rate': 48000, 'frame_size': 1.25, 'bitrate': 256000, 'ep_mode': 3,'critical_frames': list(range(812, 923))},
    {'item': 'Female_Speech_German',  "channels": 2, 'sampling_rate': 32000, 'frame_size': 1.25, 'bitrate': 512000, 'ep_mode': 4,'critical_frames': list(range(730, 741))},
    {'item': 'Female_Speech_German',  "channels": 2, 'sampling_rate': 32000, 'frame_size': 1.25, 'bitrate': 512000, 'ep_mode': 3,'critical_frames': list(range(4067, 4078))},
    {'item': 'Female_Speech_German',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': 256000, 'ep_mode': '1-4','critical_frames': list(range(1098, 1109))},
    {'item': 'Female_Speech_German',  "channels": 1, 'sampling_rate': 16000, 'frame_size': 1.25, 'bitrate': 256000, 'ep_mode': '1-4','critical_frames': list(range(3004, 3015))},
    {'item': 'Harpsichord',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': 128000, 'ep_mode': 0,'critical_frames': list(range(3786, 4082))},
    {'item': 'Harpsichord',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': '128000-128000', 'ep_mode': 0,'critical_frames': list(range(3786, 4082))},
    {'item': 'Harpsichord',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': 256000, 'ep_mode': 4,'critical_frames': list(range(3786, 4082))},
    {'item': 'Piano_Schubert',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': '128000-128000', 'ep_mode': 0,'critical_frames': list(range(2815, 2826)) + list(range(4662, 4673))},
    {'item': 'Piano_Schubert',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': 128000, 'ep_mode': 0,'critical_frames': list(range(2815, 2826))+ list(range(4662, 4673))},
    {'item': 'Piano_Schubert',  "channels": 1, 'sampling_rate': 48000, 'frame_size': 1.25, 'bitrate': 128000, 'ep_mode': 0,'critical_frames': list(range(2793, 2817))},
    {'item': 'Glockenspiel',  "channels": 1, 'sampling_rate': 16000, 'frame_size': 1.25, 'bitrate': 128000, 'ep_mode': 0,'critical_frames': list(range(6487, 6523))},
    {'item': 'Violoncello',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': 128000, 'ep_mode': 0,'critical_frames': list(range(3106, 3117))},
    {'item': 'Violoncello',  "channels": 1, 'sampling_rate': 24000, 'frame_size': 1.25, 'bitrate': '128000-128000', 'ep_mode': 0,'critical_frames': list(range(3106, 3117))}
]

# convenience wrapper for os.makedirs
def makedirs(path):
    os.makedirs(str(path), exist_ok=True)
    return path


# convenience wrapper for shutil.rmtree
def removedir(path):
    shutil.rmtree(str(path), ignore_errors=True)


# Run command and return output. cmd can be string or list. Commands with .exe suffix are automatically
# called with wine unless wine=False. Set unicode=False to get binary output. Set hard_fail=False to
# to ignore nonzero return codes.
def call(cmd, wine=True, unicode=True, hard_fail=True, log_output=True, opid='', opid_log={}, wine_error_counter=0):
    if isinstance(cmd, str):
        cmd = [x for x in shlex.split(cmd) if x]
    if sys.platform != 'cygwin' and wine and cmd[0].lower().endswith('.exe'):
        cmd = ['wine'] + cmd
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=unicode)
    out = p.communicate()[0] or (b'', '')[unicode]
    quoted_cmd = ' '.join(map(shlex.quote, cmd))
    logging.debug( '[{}] '.format(opid) + quoted_cmd)
    if opid:
        if opid_log.__contains__(opid) == False:
            opid_log[opid] = ''
        opid_log[opid] += '[CMD] ' + quoted_cmd + '\n\n' + out + '\n\n'
    if unicode and log_output:
        logging.debug(out) 
        # check if wine error occured and repeat call 3 times
        if 'client error' in out or 'Connection reset by peer' in out or 'wine: chdir to' in out or 'wine: created' in out:
            if wine_error_counter < 3:
                print('WINE ERROR OCCURRED ' + str(wine_error_counter+1) + ' times')
                print('[CMD]: ' + quoted_cmd + '\n[OUT]: ' + out )
                return call(cmd, wine=wine, unicode=unicode, hard_fail=hard_fail, log_output=log_output, opid=opid, opid_log=opid_log, wine_error_counter=wine_error_counter+1)

    if hard_fail and p.returncode != 0:
        raise OSError(quoted_cmd + ' failed!')
    return out


# return url as bytes object, validate against hash
def download(url, sha256=None):
    try:
        buf = call('curl --silent -L "{}"'.format(url), unicode=False)
    except OSError:
        sys.exit('Failed to download {}!'.format(url))
    if sha256 and hashlib.sha256(buf).hexdigest() != sha256:
        sys.exit('Failed to validate hash for {}!'.format(url))
    return buf


def download_sox():
    if not is_file(SOX_EXE):
        print('Downloading SoX ...')
        buf = download(SOX_URL, SOX_SHA256)
        zipfile.ZipFile(io.BytesIO(buf)).extractall(str(SOX_EXE.parent.parent))
        if sys.platform == 'cygwin':
            call('chmod -R +x "{}"'.format(SOX_EXE.parent))


def exe_exists(exe, wine=False):
    try:
        out = call(exe, wine=wine, hard_fail=False)
    except OSError:
        return False
    return not (wine and out.startswith('wine: ')) # detect wine: cannot find


def check_system(args, globvars):
    if sys.platform == 'win32':
        sys.exit('This script must run under cygwin')
    tools = ['curl', 'gcc', 'make']
    if sys.platform != 'cygwin' and sys.platform != 'darwin':
        tools += ['wine']
    if args.system_sox:
        tools += ['sox']
    for tool in tools:
        if not exe_exists(tool):
            sys.exit('Failed to find {} executable'.format(tool))
    if not args.wav_only: # no peaq binary needed since only wave files are created 
        if globvars['peaq_bin'] and not exe_exists(globvars['peaq_bin'], wine=True):
            sys.exit('Failed to find PEAQ executable. Please adjust config file')
    if not exe_exists(globvars['ref_encoder'], wine=True):
        sys.exit('Failed to find LC3plus ref_encoder executable ('+globvars['ref_encoder']+'). Please adjust config file.')
    if not exe_exists(globvars['ref_decoder'], wine=True):
        sys.exit('Failed to find LC3plus ref_decoder executable ('+globvars['ref_decoder']+'). Please adjust config file.')
    if not args.cmp_only: # no test coder needed since files are provided
        if not exe_exists(globvars['tst_encoder'], wine=True):
            sys.exit('Failed to find LC3plus tst_encoder executable ('+globvars['tst_encoder']+'). Please adjust config file.')
        if not exe_exists(globvars['tst_decoder'], wine=True):
            sys.exit('Failed to find LC3plus tst_decoder executable ('+globvars['tst_decoder']+'). Please adjust config file.')

# search s for expr and return first match or exit
def regex_search(expr, s):
    if not re.search(expr, s):
        sys.exit('No match for regular expression "{}"!'.format(expr))
    return re.search(expr, s).group(1)


# convert byte objects to signed int16
def byte_to_float(b, frames, channels):
    return struct.unpack("%ih" % (frames * channels), b)


def build_tools():
    call('make -C tools')


# return info about wav file
def wav_info(path):
    wav = wave.open(str(path))
    return WavInfo(wav.getnchannels(), wav.getframerate(), wav.getnframes())


# call sox with args in repeatable mode, lazy skips execution if output already exists
def sox(*args, lazy=False, opid='',opid_log={}):
    wavs = [x for x in map(str, args) if x.endswith('.wav')]
    if not (lazy and is_file(wavs[-1])): # last .wav is assumed to be output
        call('{} -R {}'.format(SOX_EXE, ' '.join(map(str, args))), opid=opid, opid_log=opid_log)


# convert file to given bits per sample (bps) 
def convert_bps(input, output, bps, lazy=False):
    sox(input, '-b', bps, output, lazy=lazy)


def trim(input, output, start, end, nChannels, lazy=False):
    if not (lazy and is_file(output)):
        tmp = output.parent / uuid_file('trim_', '.wav')
        if nChannels == 1:
            sox(input, tmp, 'trim', start, end, 'remix -')
        else:
            sox(input, tmp, 'trim', start, end)
        wi = wav_info(tmp)
        sox(tmp, output, 'fade', 0.5, wi.frames / wi.rate, 0.7)
        tmp.unlink()


# resample wav using sox
def resample(input, output, rate, lazy=False, opid='',opid_log={}):
    sox(input, output, 'rate -vs', rate, lazy=lazy, opid=opid, opid_log=opid_log)


# apply lowpass filter using sox
def low_pass(input, output, fc, lazy=False, opid='', opid_log={}):
    sox(input, output, 'sinc -{}'.format(fc), lazy=lazy, opid=opid, opid_log=opid_log)

# generate rate switching file with unique name, returns path
# byteswise step down and upwards starting in the middle
def generate_switching_file2(env, fs, *values):
    path   = env.test_dir / uuid_file('swf_', '.dat')
    fps = 2000 / fs
    minBR  = min(values) / fps / 8
    maxBR  = max(values) / fps / 8
    up   = numpy.linspace(minBR, maxBR, num=int(fps)).round()
    down = numpy.linspace(maxBR, minBR, num=int(fps)).round()
    p=numpy.repeat(numpy.concatenate((down,up)),20)
    p=numpy.delete(p,[x for x in range(0, int(fps/2))])
    p=p*8*fps
    p=p.astype('int32') 
    p.tofile(str(path))
    return path, min(values)

# get switching file, returns path
def get_switching_file(env, *values):
    layers_sorted = sorted(values)
    layers = '_'.join(map(str, layers_sorted))
    switching_dir = pathlib.Path('switching_files')

    # the files have been generated and stored using tools/gen-rate-profile.exe -layers {} {} 10 0 15000 123456789
    base_name = f"swf_layers_{layers}_frame_10_min_0_max_15000_seed_123456789.dat"
    switching_file_path = switching_dir / base_name

    # Check if switching_file exists
    if not switching_file_path.exists():
        raise FileNotFoundError(f"Switching file {switching_file_path} not found in '{pathlib.Path.cwd()}'.")

    return switching_file_path


# compares binary equality of files
def compare_bin(file1, file2, opid='', opid_log={}):
    logging.debug('compare_bin: %s %s', file1, file2)
    opid_log[opid] += 'compare_bin: {} {}\n\n'.format(file1, file2)
    return filecmp.cmp(str(file1), str(file2))


# generate unique file name with extension
def uuid_file(prefix='', suffix=''):
    return prefix + str(uuid.uuid4()) + suffix


# like str() but with special case for float
def fstr(x):
    return '{:.3g}'.format(x) if type(x) == float else str(x)


# like str() but with special case for list
def lstr(x):
    return '-'.join(map(str, x)) if type(x) in (list, tuple) else str(x)


# returns true if path is a file
def is_file(path):
    return os.path.isfile(str(path))


# calculate bitrate from bytes per frame
def get_bitrate(bytes_per_frame, frame_size):
    return int(min(bytes_per_frame * 8000 / frame_size, 320000))


# apply func to list of argumets,
def thread_executor(func, args, workers):
    list(ThreadPoolExecutor(workers).map(lambda x: func(*x), args)) # list() to collect futures


def prepare_items(workers):
    sqam_dir = pathlib.Path('SQAM')
    item_dir = makedirs(ITEM_DIR)
    if not sqam_dir.exists():
        print('Downloading test items ...')
        buf = download(SQAM_URL, SQAM_SHA256)
        zipfile.ZipFile(io.BytesIO(buf)).extractall(str(sqam_dir))

    print('Preparing test items ...')
    # trim items
    with Executor(workers) as ex:
        for item, (st, fr, flac) in ITEMS.items():
            for ch in CHANNELS:
                infile  = sqam_dir / flac
                outfile = item_dir / (item + '_{}ch.wav').format(ch)
                ex.submit(trim, infile, outfile, st, fr, ch, lazy=True)
    # resample items
    with Executor(workers) as ex:
        for item, sr in itertools.product(ITEMS, SAMPLERATES):
            for ch in CHANNELS:
                infile  = item_dir / (item + '_{}ch.wav').format(ch)
                outfile = item_dir / '{}_{}_{}ch.wav'.format(item, sr, ch)
                ex.submit(resample, infile, outfile, sr, lazy=True)
    with Executor(workers) as ex:
        # 20 kHz lowpass
        for item, sr in itertools.product(ITEMS, SAMPLERATES):
            if sr >= 44100:
                for ch in CHANNELS:
                    infile  = item_dir / '{}_{}_{}ch.wav'.format(item, sr, ch)
                    outfile = item_dir / '{}_{}_{}ch_lp20.wav'.format(item, sr, ch)
                    ex.submit(low_pass, infile, outfile, 20000, lazy=True)
        # band limit
        for sr, (bws, _) in BAND_LIMITS.items():
            for bw in bws:
                for ch in CHANNELS:
                    infile  = item_dir / '{}_{}_{}ch.wav'.format(ITEM_BAND_LIMIT, sr, ch)
                    outfile = item_dir / '{}_{}_{}ch_bw{}.wav'.format(ITEM_BAND_LIMIT, sr, ch, bw)
                    ex.submit(low_pass, infile, outfile, bw, lazy=True)
        # LP20 item with 4 seconds of white noise above 20kHz
        outfile = item_dir / (ITEM_LOW_PASS + '_48000_1ch.wav')
        synth   = 'synth 4 white fir etc/hp_fir_coef.txt'
        ex.submit(sox, '-n -r 48000 -c 1 -b 16', outfile, synth, lazy=True)
    with Executor(workers) as ex:
        # LFE item with 4 seconds of sinesweep from 10-250Hz
        for sr in SAMPLERATES:
            outfile = item_dir / '{}_{}_1ch.wav'.format(ITEM_LFE, str(sr))
            synth   = 'synth 4 sine 10:250 gain -0.5'
            ex.submit(sox, '-n -r',sr ,'-c 1 -b 16 -e signed', outfile, synth, lazy=True)
            if sr >= 44100:
                outfile2 = item_dir / '{}_{}_1ch_lp20.wav'.format(ITEM_LFE, sr)
                while not os.path.exists(outfile):
                    time.sleep(0.01)
                ex.submit(low_pass, outfile, outfile2, 20000, lazy=True)

#    # create 1ch items
#    with Executor(workers) as ex:
#        for path in os.listdir(str(item_dir)):
#            if '_2ch' in path:
#                infile  = item_dir / path
#                outfile = item_dir / path.replace('_2ch', '_1ch')
#                ex.submit(sox, infile, outfile, 'remix -', lazy=True)
    # create 24bit items
    with Executor(workers) as ex:
        for path in os.listdir(str(item_dir)):
            if 'ch.wav' in path:
                infile = item_dir / path
                outfile = item_dir / path.replace('ch.wav','ch_24b.wav')
                ex.submit(convert_bps, infile, outfile, '24', lazy=True)


def parse_config(path):
    def strip_comment(line):
        return line.split('#', 1)[0].strip()

    def split_list(line):
        return [x.strip() for x in strip_comment(line).split(',')]

    def parse_conf_line(line, test_config):
        mode, fs, sr, br = split_list(line)
        fs, sr = float(fs), int(sr)
        if '[' in br:
            br = br[1:-1]
            br = list(br.split(':'))
            br = [int(i) for i in br]
        elif ':' in br:
            br_start, br_step, br_stop = map(int, br.split(':'))
            br = list(range(br_start, br_stop + 1, br_step))
        else:
            br = [int(br)]

        # convert formats
        ## iteration to list
        if (test_config['config_bitrate_format'] == 'iteration'):
            br = list(range(br[0],br[2],br[1]))
        ## convert bytes per frame to bps
        if (test_config['config_bitrate_unit'] == 'bytes'):
            if sr == 44100:
                fs2 = fs*48000/44100
            else:
                fs2 = fs
            br = [math.ceil(x*8*1000/fs) for x in br] 

        if fs not in FRAME_SIZES:
            sys.exit('Unsupported frame size: {}!'.format(line))
        if sr not in SAMPLERATES:
            sys.exit('Unsupported sampling rate: {}!'.format(line))
#        if min(br) < 16000 or max(br) > 672000:
#            sys.exit('Invalid bitrate: {}!'.format(line))

        return mode, fs, sr, br

    def parse_bool(val):
        if val not in ('0', '1'):
            raise ValueError
        return val == '1'

    if not is_file(path):
        sys.exit('No such file: ' + path)

    globals_required = ['enabled_tests', 'ref_encoder', 'ref_decoder', 'tst_encoder', 'tst_decoder']
    globals_all      = list(DEFAULTS_GLOBAL) + globals_required
    test_keys = ['test_' + t for t in TESTS]
    globels   = DEFAULTS_GLOBAL.copy()
    configs   = {}

    try:
        parser = configparser.ConfigParser()
        parser.read(path)
        # parse global section
        for key in parser['globals']:
            globels[key] = strip_comment(parser['globals'][key])
            if key not in globals_all:
                sys.exit('Unknown key "{}" in config'.format(key))
        globels['enabled_tests'] = split_list(parser['globals']['enabled_tests'])
        # trigger KeyError for required keys
        map(lambda key: globels[key], globals_required)
        # parse test sections
        for profile in globels['enabled_tests']:
            configs[profile] = {**globels, **DEFAULTS_TEST}
            test_config = {
                'config_bitrate_unit' : DEFAULTS_GLOBAL['config_bitrate_unit'],
                'config_bitrate_format' : DEFAULTS_GLOBAL['config_bitrate_format'],
            }
            for key in parser[profile]:
                try:
                    val = strip_comment(parser[profile][key])
                    if key in test_keys:
                        configs[profile][key] = parse_bool(val)
                    elif key in test_config:
                        test_config[key] = val
                    elif key == 'configs':
                        lines = parser[profile][key].splitlines()
                        configs[profile][key] = []
                        for l in lines:
                            if parse_conf_line(l, test_config) not in configs[profile][key]:
                                configs[profile][key].append(parse_conf_line(l, test_config))
                    elif key.endswith('_threshold') and key in DEFAULTS_TEST:
                        if args.bit_exact:
                            continue
                        configs[profile][key] = eval(val)
                    elif key.endswith('_metric') and key in DEFAULTS_TEST:
                        if val not in ('rms', 'odg', 'mld', 'eng'):
                            raise ValueError
                        configs[profile][key] = val
                    elif key == 'bitflip_seed':
                        configs[profile][key] = val
                    else:
                        sys.exit('Unknown key "{}" in config'.format(key))
                except ValueError:
                    sys.exit('Invalid value in config: {} = {}'.format(key, parser[profile][key]))
    except KeyError as e:
        sys.exit('Missing "{}" in config'.format(e.args[0]))
    except configparser.DuplicateOptionError as e:
        sys.exit('Duplicate key "{}" in config'.format(e.args[1]))

    return globels, configs


# splits up files into channels, yields channel files
def split_channels(env, *wav_files):
    channels = call('{} --i -c {}'.format(SOX_EXE, wav_files[0]))
    if channels == '1':
        yield wav_files
    else:
        for ch in range(1, int(channels) + 1):
            channel_files = []
            for wav_f in wav_files:
                split_file = wav_f.with_suffix('.ch{}.wav'.format(ch))
                sox(wav_f, split_file, 'remix', ch, opid=env.opid, opid_log=env.opid_log)
                channel_files.append(split_file)
            yield channel_files


def run_rms(env, ref, tst, threshold):
    rms, diff, bits = -INF, 0, 24
    for ref_chan_x, tst_chan_x in split_channels(env, ref, tst):
        out = call('tools/rms {} {} {}'.format(ref_chan_x, tst_chan_x, threshold), opid=env.opid, opid_log=env.opid_log)
        diff_samp = int(regex_search(r'different samples\s+: (\d+)', out))
        if diff_samp != 0:
            rms  = max(rms, float(regex_search(r'Overall RMS value\s+: (\S+) dB ---', out)))
            diff = max(diff, float(regex_search(r'Maximum difference\s+: (\S+) ---', out)))
            bits = min(bits, int(regex_search(r'RMS criteria\s+: (\d+) bit', out)))
    return rms, diff, bits


def run_peaq(env, reference, test):
    odg = 5
    for ref_chan_x, tst_chan_x in split_channels(env, reference, test):
        ref48 = ref_chan_x.with_suffix('.48k.16b.wav')
        tst48 = tst_chan_x.with_suffix('.48k.16b.wav')
        # resample to 48kHz and convert to 16b
        call(f'{SOX_EXE} -R {ref_chan_x} -b 16 {ref48} rate -vs 48000',opid=env.opid, opid_log=env.opid_log)
        call(f'{SOX_EXE} -R {tst_chan_x} -b 16 {tst48} rate -vs 48000',opid=env.opid, opid_log=env.opid_log)
        out = call(env.config['peaq_bin'].format(reference=ref48, test=tst48), opid=env.opid, opid_log=env.opid_log)
        odg = min(odg, float(regex_search(env.config['peaq_odg_regex'], out)))
    return odg


def run_mld(env, reference, test):
    mld = 0
    for ref_chan_x, tst_chan_x in split_channels(env, reference, test):
        tst48 = ref_chan_x.with_suffix('.48k.wav')
        ref48 = tst_chan_x.with_suffix('.48k.wav')
        resample(ref_chan_x, ref48, 48000, opid=env.opid, opid_log=env.opid_log)
        resample(tst_chan_x, tst48, 48000, opid=env.opid, opid_log=env.opid_log)
        out = call('tools/mld -d {} {}'.format(ref48, tst48), opid=env.opid, opid_log=env.opid_log)
        mld = max(mld, float(regex_search(r'maximum loudness difference:\s*(\S+)', out)))
    return mld


# calculate energy of wav
def energy(env, test):
    logging.debug('[{}] '.format(env.opid, opid_log=env.opid_log) + 'energy: %s', str(test))
    bps = int(call('{} --i -b {}'.format(SOX_EXE, test)))
    if bps != 16:
        raise ValueError('energy test not implemented for {} bit'.format(bps))
    with wave.open(str(test), 'rb') as tst_wf:
        bytes_tst = tst_wf.readframes(tst_wf.getnframes())
        tst = byte_to_float(bytes_tst, tst_wf.getnframes(), tst_wf.getnchannels())
        eng = sum(numpy.square(tst))
        eng = 10 * math.log10(eng) if eng != 0 else -INF
        return eng


def get_payload_size(epmode, pcmode, slotbytes):
    n_codewords = math.ceil(2*slotbytes/15)
    payloadsize = slotbytes

    if epmode > 0:
        # RS redundancy
        if epmode == 1:
            payloadsize = payloadsize - 1
        else:
            payloadsize = payloadsize - n_codewords*(epmode - 1)

        # CRC1
        if slotbytes == 40:
            if epmode == 1:
                payloadsize = payloadsize - 3
            else:
                payloadsize = payloadsize - 2
        else:
            payloadsize = payloadsize - 3

        # CRC2
        if pcmode > 0 and epmode > 2 and slotbytes >= 80:
            payloadsize = payloadsize - 2

    return payloadsize


def bps_to_bpf(bitrate, frame_size):
    return math.floor(float(bitrate)/8 * frame_size/1000)


def bpf_to_bps(bytes, frame_size):
    return math.ceil(bytes * 8 * 1000 / frame_size)


def get_net_bytes(bitrate, frame_size, epmode):
    gross_bytes = bps_to_bpf(bitrate, frame_size)
    net_bytes = get_payload_size(epmode, 1, gross_bytes)
    return net_bytes


# calculate threshold depending on ODG of reference and for hrmode on net bytes per frame
# higher ref ODG correlate with higher delta ODG
def get_delta_threshold(hr_info, sr, odg_ref, delta_odg_thresh):
    if hr_info['hrmode']:
        # for ep_mode switching calculate the net bytes with ep mode 4 (highest)
        if isinstance(hr_info['ep_mode'], pathlib.PosixPath):
            hr_info['ep_mode'] = 4
        real_bytes = get_net_bytes(hr_info['bitrate'], hr_info['frame_size'], hr_info['ep_mode'])
        if real_bytes < MIN_HRMODE_BYTES[(hr_info['frame_size'], sr)]:
            hr_info['hrmode'] = 0
            return get_delta_threshold(hr_info, sr, odg_ref, [0.06, 0.08, 0.12])
    if sr == 8000:
        return 0.15
    if type(delta_odg_thresh).__name__ != 'list':
        return delta_odg_thresh
    elif odg_ref > -1.0:
        return delta_odg_thresh[0]
    elif odg_ref > -2.3: 
        return delta_odg_thresh[1]
    else:
        return delta_odg_thresh[2]


def get_eng_threshold(frame_size, engergy_thresholds):
    if frame_size < 2.5:
        return engergy_thresholds[1]
    else:
        return engergy_thresholds[0]


# compare output wavs by metric rms, odg, mld, eng
def compare_wav(env, hr_info, mode, sr, infile, file_ref, file_tst):
    mkey   = '{}_{}_metric'.format(env.test, mode)
    metric = env.config[mkey]
    tkey   = '{}_{}_{}_threshold'.format(env.test, mode, metric)
    thresh = env.config[tkey]

    ref_al = file_ref.with_suffix('.aligned.wav')
    tst_al = file_tst.with_suffix('.aligned.wav')
    infile_copy = infile.with_suffix('.copy_for_{}t.wav'.format(file_ref.stem))
    call('cp {} {}'.format(infile, infile_copy), opid=env.opid, opid_log=env.opid_log)
    
    infile_len = int(call('{} --i -s {}'.format(SOX_EXE, infile_copy)))
    call('tools/align {} {} {} {}'.format(file_ref, ref_al, infile_len, env.config['delay']), opid=env.opid, opid_log=env.opid_log)
    call('tools/align {} {} {} {}'.format(file_tst, tst_al, infile_len, env.config['delay']), opid=env.opid, opid_log=env.opid_log)

    if metric == 'rms':
        rms, diff, bits = run_rms(env, ref_al, tst_al, thresh)
        rms_thr  = 20 * math.log10(2 ** (-thresh + 1) / 12 ** 0.5)
        tkey     = '{}_{}_mad_threshold'.format(env.test, mode)
        diff_thr = env.config[tkey]
        ok_rms   = rms <= rms_thr
        ok_diff  = diff <= diff_thr
        ok_bits  = bits >= thresh
        ok       = ok_rms and ok_diff
        values   = [(diff,  ('fail', 'pass')[ok_diff], diff_thr),
                    (rms, ('fail', 'pass')[ok_rms], rms_thr),
                    (bits, ('warn', 'none')[ok_bits], thresh)]
    if metric == 'odg':
        odg_ref  = run_peaq(env, infile_copy, ref_al)
        odg_tst  = run_peaq(env, infile_copy, tst_al)
        odg_diff = abs(odg_ref - odg_tst)
        thresh = get_delta_threshold(hr_info, sr, odg_ref, thresh)
        ok       = odg_diff <= thresh
        values   = [(odg_ref, '', None),
                    (odg_diff, ('fail', 'pass')[ok], thresh)]
    if metric == 'mld':
        mld    = run_mld(env, ref_al, tst_al)
        ok     = mld <= thresh
        values = [(mld, ('fail', 'pass')[ok], thresh)]
    if metric == 'eng':
        d_eng  = energy(env, file_tst)
        thresh = get_eng_threshold(hr_info['frame_size'], thresh)
        ok     = d_eng <= thresh
        values = [(d_eng, ('fail', 'pass')[ok], thresh)]

    return ok, values

def check_bfi_file(file, opid='', opid_log={}):
    logging.debug('check_bfi_file: %s', file)
    opid_log[opid] += 'check_bfi_file: {}\n\n'.format(file)
    with open(file, 'rb') as file:
        while True:
            chunk = file.read(4)
            if not chunk:
                break

            index = chunk.find(b'\x01')
            if index != -1:
                return False
            
            index = chunk.find(b'\x02')
            if index != -1:
                return False
            
    return True

# compare outout files of ep debug flag
def compare_errors(env, file_ref, file_tst):
    ok_all, values = True, []
    if env.test == 'sqam':
        ok = check_bfi_file(str(file_ref) + '.bfi', opid=env.opid, opid_log=env.opid_log)
        ok_all = ok_all and ok
        ok = check_bfi_file(str(file_tst) + '.bfi', opid=env.opid, opid_log=env.opid_log)
        ok_all = ok_all and ok
        values += [('bad', 'ok')[ok_all]]
    else:
        for ext in ['.bfi', '.epmr', '.error_report']:
            ok = compare_bin(str(file_ref) + ext, str(file_tst) + ext, opid=env.opid, opid_log=env.opid_log)
            values += [(('bad', 'ok')[ok], ('fail', 'pass')[ok], None)]
            ok_all = ok_all and ok
            
    return ok_all, values


# ensure inputs exist and nothig is overwritten
def check_io_files(input, output):
    if not is_file(input):
        raise FileNotFoundError(input)
    if is_file(output):
        raise FileExistsError(output)


def encode_reference(env, input, output, frame_size, bitrate, bandwidth=None, ep_mode=0, lfe=0):
    check_io_files(input, output)
    cmd = env.config['ref_encoder']
    opt = []
    if lfe:
        opt += ['-lfe 1']
    if bandwidth:
        opt += ['-bandwidth', bandwidth]
    if ep_mode:
        opt += ['-epmode', ep_mode]
    options = ' '.join(map(str, opt))
    call(cmd.format(input=input, output=output, frame_size=frame_size, bitrate=bitrate, options=options), opid=env.opid, opid_log=env.opid_log)


def encode_test(env, input, output, frame_size, bitrate, bandwidth=None, ep_mode=0,lfe=0):
    check_io_files(input, output)
    cmd = env.config['tst_encoder']
    opt = []
    if lfe:
        opt += ['-lfe 1']
    if bandwidth:
        opt += [env.config['option_bandwidth'].format(arg=bandwidth)]
    if ep_mode:
        opt += [env.config['option_ep_mode'].format(arg=ep_mode)]
    options = ' '.join(opt)
    call(cmd.format(input=input, output=output, frame_size=frame_size, bitrate=bitrate, options=options), opid=env.opid, opid_log=env.opid_log)


def decode_reference(env, input, output, error_file=None):
    check_io_files(input, output)
    cmd = env.config['ref_decoder']
    opt = []
    if error_file:
        opt += ['-ep_dbg', error_file]
    options = ' '.join(map(str, opt))
    call(cmd.format(input=input, output=output, options=options), opid=env.opid, opid_log=env.opid_log)


def decode_test(env, input, output, error_file=None):
    check_io_files(input, output)
    cmd = env.config['tst_decoder']
    opt = []
    if error_file:
        opt += [env.config['option_ep_debug'].format(arg=error_file)]
    options = ' '.join(map(str, opt))
    call(cmd.format(input=input, output=output, options=options), opid=env.opid, opid_log=env.opid_log)


def apply_error_pattern(env, input, output, mode, pattern):
    assert mode in ('fer', 'ber', 'flip')
    check_io_files(input, output)
    if mode == 'fer':
        if sys.platform == 'darwin':
            cmd = 'tools/eid-xor_mac -vbr -bs g192 -ep byte -fer {} {} {}'
        else: 
            cmd = 'tools/eid-xor_linux -vbr -bs g192 -ep byte -fer {} {} {}'
        call(cmd.format(input, pattern, output), opid=env.opid, opid_log=env.opid_log)
    if mode == 'ber':
        if sys.platform == 'darwin':
            cmd = 'tools/eid-xor_mac -vbr -bs g192 -ep byte -ber {} {} {}'
        else:
            cmd = 'tools/eid-xor_linux -vbr -bs g192 -ep byte -ber {} {} {}'
        call(cmd.format(input, pattern, output), opid=env.opid, opid_log=env.opid_log)
    if mode == 'flip':
        if 'bitflip_seed' in env.config.keys():
            seed = env.config['bitflip_seed']
        else: 
            seed = '1911'
        cmd = 'tools/flipG192 {} {} {} {} ' + seed + ' 0'
        flips, frames = pattern
        call(cmd.format(input, output, flips, frames), opid=env.opid, opid_log=env.opid_log)
    # copy the config file of g192 bitstreams
    if is_file(str(input) + '.cfg'):
        call('cp {} {}'.format(str(input) + '.cfg', str(output) + '.cfg'), opid=env.opid, opid_log=env.opid_log)

def read_bit(bytes, bp, mask):
    if bytes[bp] & mask:
        bit = 1
    else:
        bit = 0

    if mask == 128:
        mask = 1
        bp -= 1
    else:
        mask <<= 1

    return bit, bp, mask

def read_uint(bytes, bp, mask, nbits):
    val, bp, mask = read_bit(bytes, bp, mask)

    for i in range(1, nbits):
        bit, bp, mask = read_bit(bytes, bp, mask)
        val += bit << i

    return val, bp, mask

def getLastNzBits(N):
    minBits = int(numpy.ceil(numpy.log2(N >> 1)))

    if ((1 << minBits) - (N >> 1)) < 2:
        minBits += 1
    
    return minBits

def decode_and_edit_bitstream(config, data, bitstream, nbytes, bytes_read):
    # bw_cutoff_bits is set in FillDecSetup() in setup_dec_lc3.c as decoder->BW_cutoff_bits = BW_cutoff_bits_all[decoder->fs_idx];
    # decoder->fs_idx = ( decoder->fs ) == 96000 ? 5 : ( decoder->fs ) / 10000
    BW_cutoff_bits_all =  [ 0, 1, 2, 2, 3, 0 ] # Defined in constans.c
    fs_idx = 5 if config['samplingrate'] == 96000 else int(config['samplingrate']/10000)
    bw_cutoff_bits = BW_cutoff_bits_all[fs_idx]

    # N = decoder->yLen which is set in setup_dec_lc3.c in set_dec_frame_params()
    frame_length = int(numpy.ceil( config['samplingrate'] * 10 / 1000 ))
    if config['hrmode'] == 1:
        N = frame_length
    else:
        N = min( 400, frame_length ) # 400 is MAX_BW
    if config['frame_ms'] == 7.5:
        N = int((N * 3) / 4)
    else:
        N = int(N / int(10/config['frame_ms']))
    if config['frame_ms']  == 1.25:
        bw_cutoff_bits = 0

    # Start decoding the side information
    bp = 0
    bp_side = nbytes - 1
    mask_side = 1

    # Bandwidth
    if bw_cutoff_bits > 0:
        bw_cutoff_idx, bp_side, mask_side = read_uint( bw_cutoff_bits, data, bp_side, mask_side, bw_cutoff_bits )
    else:
        bw_cutoff_idx = fs_idx

    # Last non-zero tuple
    lastnz, bp_side, mask_side = read_uint( data, bp_side, mask_side, getLastNzBits( N ) )
    lastnz = ( lastnz + 1 ) * 2

    # LSB mode bit
    lsbMode, bp_side, mask_side = read_bit( data, bp_side, mask_side )

    # Global Gain
    gg_idx, bp_side, mask_side = read_uint( data, bp_side, mask_side, 8 )

    # TNS activation flag
    if config['frame_ms']  == 1.25:
        tns_numfilters = 0
    else:
        if bw_cutoff_idx < 3 or config['frame_ms']  == 1.25:
            tns_numfilters = 1
        else:
            tns_numfilters = 2
    tns_order = numpy.zeros(tns_numfilters)
    for i in range(tns_numfilters):
        tns_order[i], bp_side, mask_side = read_bit( data, bp_side, mask_side )

    # LTPF activation flag
    pitch_present, bp_side, mask_side = read_bit( data, bp_side, mask_side )

    # Skipping 38 bits in SNS data
    if (8 - int(numpy.log2(mask_side))) > 1:
        bp_side -= 1
    bp_side -= (38 - (8 - int(numpy.log2(mask_side)))) // 8

    mask_side = 1 << (38 - (8 - int(numpy.log2(mask_side)))) % 8

    # LTPF pitch data
    if pitch_present:
        mask_side_b4 = mask_side
        bp_side_b4 = bp_side
        
        # Read ltpf flag by indexing the bitstream
        ltpf_idx_in_bitstream = bytes_read-16*(nbytes-1-(bp_side-1)) + int(numpy.log2(mask_side)*2)
        x, = struct.unpack('<H', bitstream[ltpf_idx_in_bitstream : ltpf_idx_in_bitstream + 2])
        ltpf_active_v2 = 1 if x == G192_ONE else 0 

        # Read ltpf flag as in the C code
        ltpf_active, bp_side, mask_side = read_uint( data, bp_side, mask_side, 1 )

        assert( ltpf_active == ltpf_active_v2 )

        # If LTPF is active, deactivate it by modifying the bitsream.
        if ltpf_active:
            bitstream[ltpf_idx_in_bitstream : ltpf_idx_in_bitstream + 2] = struct.pack('<H', G192_ZERO)
        else:
            bitstream[ltpf_idx_in_bitstream : ltpf_idx_in_bitstream + 2] = struct.pack('<H', G192_ONE)

        pitch_index, bp_side, mask_side = read_uint( data, bp_side, mask_side, 9 )
    else:
        ltpf_active = 0
        pitch_index = 0

    return bitstream

def read_cfg_header(g192_cfg_file):
    with open(g192_cfg_file, 'rb') as f_g192:
            data = f_g192.read()

    data = struct.unpack('<HhhhhhhHhh', data[:20])

    config = {}
    config['file_id']           = data[0]
    config['header_size']       = data[1]
    config['samplingrate']      = data[2] * 100
    config['bitrate']           = data[3] * 100
    config['channels']          = data[4]
    config['frame_ms']          = data[5] / 100
    config['epmode']            = data[6]
    config['signal_len']        = data[7]
    config['signal_len_red']    = data[8] << 16
    config['hrmode']            = data[9]

    assert( hex(config['file_id']) == '0xcc1c' )

    return config

def pack_unpack_bitstream_epmode(g192_file_in, g192_file_out, fs, br, ep_mode, task='unpack'):
    try:
        print(os.getcwd())
        ccConvert = '../src/fixed_point/ccConvert'
        if task == 'unpack':
            cmd = [ccConvert, '-unpack', g192_file_in, g192_file_out]
        elif task == 'pack':
            gross_bytes = int((br / 8) / (1000 / fs))
            cmd = [ccConvert, '-pack', str(gross_bytes), str(ep_mode), g192_file_in, g192_file_out]
        call(cmd)
    except:
        raise Exception

# Toggles LTPF activations flag
def edit_ltpf_bitstream(g192_file, frames, fs, br, ep_mode, opid='', opid_log={}):
    logging.debug( '[{}] '.format(opid) + f'edit_ltpf_bitstream: {g192_file}')
    if opid:
        if opid_log.__contains__(opid) == False:
            opid_log[opid] = ''
        opid_log[opid] += f'[CMD] # edit_ltpf_bitstream: {g192_file}\n\n'
    g192_cfg_file = g192_file + '.cfg'

    config = read_cfg_header(g192_cfg_file)

    problematic_frames = numpy.array(frames)
    problematic_frames -= 20

    bytes_read = 0

    # unpack bitstream in case of epmode
    if ep_mode > 0:
        g192_file_unpacked = g192_file.replace('.g192', '_unpacked.g192') #'out_unpacked.g192'
        pack_unpack_bitstream_epmode(g192_file, g192_file_unpacked, fs, br, ep_mode, task='unpack')
        with open(g192_file_unpacked, 'rb') as f_g192:
            bitstream = bytearray(f_g192.read())
    else:
        with open(g192_file, 'rb') as f_g192:
            bitstream = bytearray(f_g192.read())


    with open(g192_file, 'rb') as f_g192:
        bitstream = bytearray(f_g192.read())

    for x, frame in enumerate(problematic_frames):
        for f in range(frame+1) if x==0 else range(problematic_frames[x-1]+1, frame+1):
            frame_indicator, length = struct.unpack('<HH', bitstream[bytes_read:bytes_read+4])
            bytes_read += 4

            assert( frame_indicator == G192_GOOD_FRAME )

            nbytes = length // 8
            data = numpy.zeros(nbytes, dtype=numpy.uint8)

            for i in range(nbytes):
                byte = numpy.uint8(0)
                
                for b in range(8):
                    x, = struct.unpack('<H', bitstream[bytes_read:bytes_read+2])
                    if x == G192_ONE:
                        byte |= numpy.uint8(1) << b
                    bytes_read += 2
                data[i] = byte

            if f == frame:
                bitstream = decode_and_edit_bitstream(config, data, bitstream, nbytes, bytes_read)

    if ep_mode > 0:
        # if epmode, pack the bitstream back into protected mode
        # Write the modified bitstream
        g192_file_unpacked_edit = g192_file.replace('.g192', '_unpacked.edit.g192')
        call(f'cp {g192_file_unpacked}.cfg {g192_file_unpacked_edit}.cfg')
        with open(g192_file_unpacked_edit, 'wb') as f_g192:
            f_g192.write(bitstream)
        #g192_file_new_packed = 'out_modified_packed.g192'
        pack_unpack_bitstream_epmode(g192_file_unpacked_edit, g192_file, fs, br, ep_mode, task='pack')
    else:
        # Write the modified bitstream
        g192_file_new = g192_file
        with open(g192_file_new, 'wb') as f_g192:
            f_g192.write(bitstream) 

# rename switcing files
def swf_bitrate_to_label (tmp):
    return '{}-to-{}'.format(tmp[0], tmp[-1])


# create file names for test
def make_files(env, files, *args):
    tmp = args
    if type (tmp[4]).__name__ == 'list':
        tmp = list(tmp)
        tmp[4] = swf_bitrate_to_label(tmp[4])
        tmp = tuple(tmp)
    protoyp = '_'.join(map(lstr, tmp)) + '_'
    return tuple(env.test_dir / (protoyp + f) for f in files)


# permutate test configs
def sqam_configs(config, items, ch=1, lp20=True, modes=None):
    for mode, fs, sr, brs in config['configs']:
        if modes and mode not in modes:
            continue
        for item, br in itertools.product(items, brs):
            if config['hrmode'] == '1':
                suffix = '_24b'
            else:
                suffix = '_lp20' if lp20 and sr in (44100, 48000) else ''
            infile = ITEM_DIR / '{}_{}_{}ch{}.wav'.format(item, sr, ch, suffix)
            yield mode, item, fs, sr, br, infile


# apply test func to list of tests, multithreadded
def test_executor(env, func, tests):
    ex = Executor(env.workers)
    return list(ex.map(lambda args: func(*args), tests))


def check_ltpf_conditions(fs,sr,br,epmode,item_path):
    for condition in LTPF_BITSTREAM_CONDITIONS:
        if '.dat' in str(epmode):
            epmode_tmp = '1-4'
        else:
            epmode_tmp = epmode
        if (fs == condition['frame_size'] and sr == condition['sampling_rate'] and br == condition['bitrate'] and epmode_tmp == condition['ep_mode'] and condition['item'] in item_path and f"_{str(condition['channels'])}ch_" in item_path):
                return condition['critical_frames']
    return []

def silence_samples(input_wav, output_wav, start_sample, len_silence):
    with wave.open(str(input_wav), 'rb') as wav_in:
        params = wav_in.getparams()  # Get parameters of the wav file
        frames = wav_in.readframes(params.nframes)  # Read all frames

    # Convert frames to numpy array + Create a writable copy
    audio_data = numpy.frombuffer(frames, dtype=numpy.int16).copy() 

    start_index = start_sample * params.nchannels
    # Silence the critical samples
    audio_data[start_index:start_index + len_silence * params.nchannels] = 0

    # Write the modified audio to the output WAV file
    tmp_path = output_wav.parent / f"{output_wav.name}.{uuid.uuid4().hex}.tmp"
    with wave.open(str(tmp_path), 'wb') as wav_out:
        wav_out.setparams(params)
        wav_out.writeframes(audio_data.tobytes())
    os.replace(tmp_path, output_wav)

# process a single test item
# performs encoding, erroring, decoding and evaluation
# returns tuple of bool, list (pass condition, metric values)
# mode: encode/decode/encdec, fs: frame size, sr: sampling rate, br: bitrate
def process_item(env, mode, item, fs, sr, br, infile, bandwidth=None, ep_mode=0, error_mode=None, error_pattern=None, lfe=0):
    bw_name = 'swf' if is_file(bandwidth) else bandwidth
    ep_name = '1-4' if is_file(ep_mode) else ep_mode
    fmt     = '  {} {:20} {:3g} ms {:5} Hz {:>6} bit/s ep:{}'
    br_label = br
    if type(br).__name__ == 'list' or type(br).__name__ == 'tuple':
        br_label = swf_bitrate_to_label(br)
    opid = '{}__{}__{}__{}__{}__{}__{}__{}'.format(env[0],env[1],mode,item,fs,sr,br_label,ep_name) 
    env = env._replace(opid = opid) # operation point ID to log commands associated with a operationpoint (=unique configuration)
    wav_only = env[7]['active']
    wav_only_dir = env[7]['path']
    cmp_only = env[8]['active']
    cmp_only_dir = env[7]['path']
    print(fmt.format(mode, item, fs, sr, lstr(br), ep_name))
    if isinstance(br, int):
        br_param = br
    elif isinstance(br, tuple):
        br_param = lstr(br)

    file_names = ['r.g192', 't.g192', 're.g192', 'te.g192', 'r.wav', 't.wav', 'rd', 'td']
    file_tuple = make_files(env, file_names, mode, infile.stem, fs, sr, br, bw_name, ep_name)
    ref_bin, tst_bin, ref_err, tst_err, ref_wav, tst_wav, ref_dbg, tst_dbg = file_tuple
    # evaluate channel coder output only for decode_ep_* tests
    if not (env.test.startswith('ep_') and mode == 'decode'):
        err_ok, err_val, ref_dbg, tst_dbg = True, [], None, None

    try:
        # generate bitrate switching file if needed
        sw_min_bitrate = 0
        if type(br) in (list, tuple):
            #br = generate_switching_file(env, *br)
            br, sw_min_bitrate = generate_switching_file2(env, fs, *br)
        # encode
        encode_reference(env, infile, ref_bin, fs, br, bandwidth=bandwidth, ep_mode=ep_mode, lfe=lfe)
        # for ltpf with adaptive gain, edit the bitstream
        critical_frames = check_ltpf_conditions(fs,sr,br_param,ep_mode,ref_bin.stem)
        if 0: # len(critical_frames) > 0:
            edit_ltpf_bitstream(str(ref_bin), critical_frames, fs, br, ep_mode, opid=env.opid, opid_log=env.opid_log)
        if not cmp_only:
            if mode in ('encode', 'encdec'):
                encode_test(env, infile, tst_bin, fs, br, bandwidth=bandwidth, ep_mode=ep_mode, lfe=lfe)
                critical_frames = check_ltpf_conditions(fs,sr,br_param,ep_mode,tst_bin.stem)
                if 0: # len(critical_frames) > 0:
                    edit_ltpf_bitstream(str(tst_bin), critical_frames, fs, br, ep_mode, opid=env.opid, opid_log=env.opid_log)
            # apply errors
            if error_mode:
                apply_error_pattern(env, ref_bin, ref_err, error_mode, error_pattern)
                ref_bin = ref_err
                if mode in ('encode', 'encdec'):
                    apply_error_pattern(env, tst_bin, tst_err, error_mode, error_pattern)
                    tst_bin = tst_err
        # decode
        decode_reference(env, ref_bin, ref_wav, error_file=ref_dbg)
        if not cmp_only:
            if mode == 'encode':
                decode_reference(env, tst_bin, tst_wav, error_file=tst_dbg)
            if mode == 'encdec':
                decode_test(env, tst_bin, tst_wav, error_file=tst_dbg)
            if mode == 'decode':
                decode_test(env, ref_bin, tst_wav, error_file=tst_dbg)
        # compare outputs
        if wav_only:
            ok, err_ok, val, err_val = True, True, [], []
            dest_path = pathlib.Path(str(tst_wav).replace(tst_wav.parts[0], wav_only_dir))
            call('mkdir -p ' + str(dest_path.parent))
            call('cp ' + str(tst_wav) + ' ' + str(dest_path))
            if tst_dbg:
                call('cp ' + str(tst_dbg) + '.bfi ' + str(dest_path.parent))
                call('cp ' + str(tst_dbg) + '.epmr ' + str(dest_path.parent))
                call('cp ' + str(tst_dbg) + '.error_report ' + str(dest_path.parent))
        else:
            if cmp_only:
                tst_wav = pathlib.Path(str(tst_wav).replace(tst_wav.parts[0], cmp_only_dir))
            hr_info = {'frame_size':fs, 'ep_mode':ep_mode, 'hrmode':int(env.config['hrmode']), 'bitrate':sw_min_bitrate if sw_min_bitrate else br}
            if critical_frames:
                for i, crit_frame in enumerate(critical_frames):
                    frame_length = int(numpy.ceil( sr * 10 / 1000 ))
                    N = int(frame_length / int(10/fs))
                    start_sample = crit_frame * N
                    new_infile = re.sub(r'(?i)^(.*?)(?:\.edit(\d+))?(\.wav)$', lambda m: f"{m.group(1)}.edit{int(m.group(2) or 0)+1}{m.group(3)}", str(infile))
                    new_infile = pathlib.Path(new_infile)
                    silence_samples(infile, new_infile, start_sample, N)
                    infile = new_infile
                    new_ref_wav = re.sub(r'(?i)^(.*?)(?:\.edit(\d+))?(\.wav)$', lambda m: f"{m.group(1)}.edit{int(m.group(2) or 0)+1}{m.group(3)}", str(ref_wav))
                    new_ref_wav = pathlib.Path(new_ref_wav)
                    silence_samples(ref_wav, new_ref_wav, start_sample, N)
                    ref_wav = new_ref_wav
                    new_tst_wav = re.sub(r'(?i)^(.*?)(?:\.edit(\d+))?(\.wav)$', lambda m: f"{m.group(1)}.edit{int(m.group(2) or 0)+1}{m.group(3)}", str(tst_wav))
                    new_tst_wav = pathlib.Path(new_tst_wav)
                    silence_samples(tst_wav, new_tst_wav, start_sample, N)
                    tst_wav = new_tst_wav
            ok, val = compare_wav(env, hr_info, mode, sr, infile, ref_wav, tst_wav)
            if ref_dbg and tst_dbg:
                if cmp_only:
                    dest_path = pathlib.Path(str(tst_wav).replace(tst_wav.parts[0], cmp_only_dir))
                    tst_dbg = str(dest_path.parent) + '/' + str(tst_dbg.parts[-1])
                if (mode == 'decode' or mode == 'encdec') and (env.test == 'sqam' or env.test.startswith('ep_')):
                    err_ok, err_val = compare_errors(env, ref_dbg, tst_dbg)
        return ok and err_ok, err_val + val, opid

    except (OSError, FileNotFoundError, FileExistsError, KeyError) as e:
        logging.error('[%s] process_item: %s: %s', opid, type(e).__name__, str(e))
        return False, [], opid


def test_sqam(env):
    print('Testing SQAM ...')
    def func(mode, item, fs, sr, br, infile):
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile)
        return [ok, mode, item, fs, sr, br] + values + [opid]

    tests = sqam_configs(env.config, ITEMS)
    return test_executor(env, func, tests)


def test_band_limiting(env):
    print('Testing band limitation ...')
    def func(mode, item, fs, sr, br, bw):
        infile = ITEM_DIR / '{}_{}_1ch_bw{}.wav'.format(item, sr, bw)
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile)
        return [ok, mode, item, fs, sr, br, bw] + values + [opid]

    tests = set()
    for mode, fs, sr, _ in env.config['configs']:
        if sr >= 16000:
            bw_limits, frame_bytes = BAND_LIMITS[sr]
            br = get_bitrate(frame_bytes, fs)
            for bw in bw_limits:
                tests.add((mode, ITEM_BAND_LIMIT, fs, sr, br, bw))
    return test_executor(env, func, list(tests))


def test_low_pass(env):
    print('Testing low pass ...')
    def func(mode, item, fs, sr, br, infile):
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile)
        return [ok, mode, item, fs, sr, br] + values + [opid]

    items = [ITEM_LOW_PASS]
    modes = ['encode', 'encdec']
    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, items , modes=modes, lp20=False):
        if sr >= 44100:
            tests.append((mode, item, fs, sr, br, infile))
    return test_executor(env, func, tests)


def test_lfe(env):
    print('Testing LFE mode ...')
    def func(mode, item, fs, sr, br, infile):
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile, lfe=1)
        return [ok, mode, item, fs, sr, br] + values + [opid]

    items = [ITEM_LFE]
    modes = ['encode']
    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, items , modes=modes):
        tests.append((mode, item, fs, sr, br, infile))
    return test_executor(env, func, tests)


def test_bitrate_switching(env):
    print('Testing bitrate switching ...')
    def func(mode, item, fs, sr, bitrates):
        if env.config['hrmode'] == '1':
            br     = (min(bitrates), max(bitrates)) # take all steps
            infile = ITEM_DIR / '{}_{}_1ch{}.wav'.format(item, sr, '_24b')
        else:
            br     = (int(160000 / fs), max(bitrates))
            infile = ITEM_DIR / '{}_{}_1ch.wav'.format(item, sr)
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile)
        if type(br).__name__ == 'list':
            br = swf_bitrate_to_label(br)
        return [ok, mode, item, fs, sr, br] + values + [opid]

    tests = []
    for mode, fs, sr, bitrates in env.config['configs']:
        for item in ITEMS:
            tests.append((mode, item, fs, sr, bitrates))
    return test_executor(env, func, tests)


def test_bandwidth_switching(env):
    print('Testing bandwidth switching ...')
    def func(mode, item, fs, sr, br, infile):
        bwf = get_switching_file(env, *BAND_WIDTHS[sr])
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile, bwf)
        return [ok, mode, item, fs, sr, br] + values + [opid]

    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, ITEMS, modes=['encode','encdec']):
        if sr >= 16000:
            tests.append((mode, item, fs, sr, br, infile))
    return test_executor(env, func, tests)


def test_plc(env):
    print('Testing packet loss concealment ...')
    def func(mode, item, fs, sr, br, infile):
        pattern = PLC_RANDOM_FILES[fs]
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile, None, 0, 'fer', pattern)
        return [ok, mode, item, fs, sr, br] + values + [opid]

    tests = sqam_configs(env.config, ITEMS_PLC, modes=['decode'])
    return test_executor(env, func, tests)


def test_plc_burst(env):
    print('Testing packet loss concealment ...')
    def func(mode, item, fs, sr, br, infile):
        pattern = PLC_BURST_FILES[fs]
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile, None, 0, 'fer', pattern)
        return [ok, mode, item, fs, sr, br] + values + [opid]

    tests = sqam_configs(env.config, ITEMS_PLC, modes=['decode'])
    return test_executor(env, func, tests)


def test_pc(env):
    print('Testing partial concealment ...')
    def func(mode, item, fs, sr, br, infile):
        pattern = PC_BER_FILES[fs]
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile, None, 4, 'ber', pattern)
        return [ok, mode, item, fs, sr, br] + values + [opid]

    tests = sqam_configs(env.config, ITEMS_PLC, modes=['decode'])
    return test_executor(env, func, tests)


def test_ep_correctable(env):
    print('Testing channel coder for correctable frames ...')
    def func(mode, item, fs, sr, br, infile, ep_mode):
        pattern = 'etc/ep_{}bit_{}ms_{}epm_{}cor_{}cmb.dat'.format(int(br), fs, int(ep_mode),1,0) # correctable 1, combined 0)
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile, None, ep_mode, 'fer', pattern)
        return [ok, mode, item, fs, sr, br, ep_mode] + values + [opid]

    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, ITEMS_PLC, modes=['encode','decode']):
        for ep_mode in [1, 2, 3, 4]:
            tests.append((mode, item, fs, sr, br, infile, ep_mode))
    return test_executor(env, func, tests)


def test_ep_non_correctable(env):
    print('Testing channel coder for non-correctable frames ...')
    def func(mode, item, fs, sr, br, infile, ep_mode):
        pattern = 'etc/ep_{}bit_{}ms_{}epm_{}cor_{}cmb.dat'.format(int(br), fs, int(ep_mode),0,0) # correctable 0, combined 0)
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile, None, ep_mode, 'fer', pattern)
        return [ok, mode, item, fs, sr, br, ep_mode] + values + [opid]

    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, ITEMS_PLC, modes=['decode']):
        for ep_mode in [1, 2, 3, 4]:
            tests.append((mode, item, fs, sr, br, infile, ep_mode))
    return test_executor(env, func, tests)


def test_ep_mode_switching(env):
    print('Testing ep-mode switching ...')
    ep_mode = get_switching_file(env, 100, 200, 300, 400)
    def func(mode, item, fs, sr, br, infile):
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile, None, ep_mode, None, None)
        return [ok, mode, item, fs, sr, br, '1-4'] + values + [opid]

    tests = sqam_configs(env.config, ITEMS_PLC, modes=['encode','decode'])
    return test_executor(env, func, tests)


def test_ep_combined(env):
    print('Testing combined channel coder for correctable frames ...')
    def func(mode, item, fs, sr, br, infile, ep_mode):
        pattern = 'etc/ep_{}bit_{}ms_{}epm_{}cor_{}cmb.dat'.format(int(br), fs, int(ep_mode),1,1) # correctable 1, combined 1)
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile, None, ep_mode, 'fer', pattern)
        return [ok, mode, item, fs, sr, br, ep_mode] + values + [opid]

    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, ITEMS_PLC, ch=2, modes=['encode','decode']):
        bytes_per_frame = (br / 800) / (10 / fs)
        # Combined channel coder works for bitrates up to 128 kbps; 160 = 128000 / 800
        if bytes_per_frame <= 160:
            for ep_mode in [1, 2, 3, 4]:
                #print('pattern: {},{},{},{}'.format(br,fs,ep_mode, '1','1'))
                tests.append((mode, item, fs, sr, br, infile, ep_mode))

    return test_executor(env, func, tests)


def test_ep_combined_nc(env):
    print('Testing combined channel coder for non-correctable frames ...')
    def func(mode, item, fs, sr, br, infile, ep_mode):
        pattern = 'etc/ep_{}bit_{}ms_{}epm_{}cor_{}cmb.dat'.format(int(br), fs, int(ep_mode),0,1) # correctable 0, combined 1)
        ok, values, opid = process_item(env, mode, item, fs, sr, br, infile, None, ep_mode, 'fer', pattern)
        return [ok, mode, item, fs, sr, br, ep_mode] + values + [opid]

    tests = []
    for mode, item, fs, sr, br, infile in sqam_configs(env.config, ITEMS_PLC, ch=2, modes=['decode']):
        bytes_per_frame = (br / 800) / (10 / fs)
        if bytes_per_frame <= 160:
            for ep_mode in [1, 2, 3, 4]:
                #print('pattern: {},{},{},{}'.format(br,fs,ep_mode, '0','1'))
                tests.append((mode, item, fs, sr, br, infile, ep_mode))
    return test_executor(env, func, tests)


def pass_ratio(results):
    num_passed = sum(ok for ok, *_ in results)
    return num_passed / len(results) if results else 1.0


def profile_passed(all_results):
    flat_results = itertools.chain(*all_results.values())
    return all(ok for ok, *_ in flat_results)


def gen_td(value):
    if type(value) in (tuple, list) and len(value) == 3 and value[1] != '-to-':
        value, clazz, thresh = value
        clazz  = ' class={}'.format(clazz) if clazz else ''
        thresh = ' ({})'.format(fstr(thresh)) if thresh != None else ''
        return '<td{}>{}{}</td>'.format(clazz, fstr(value), thresh)
    else:
        return '<td>{}</td>'.format(lstr(value))


def gen_table(test, mode, config, results):
    mkey   = '{}_{}_metric'.format(test, mode)
    metric = config[mkey]
    header = HEADER_ALL.copy()
    if test == 'band_limiting':
        header += HEADER_BL
    if test == 'sqam' and (mode == 'decode' or mode == 'encdec'):
        header += HEADER_EPD_BFI
    if test.startswith('ep_'):
        header += HEADER_EP
        if mode == 'decode':
            header += HEADER_EPD
    header += HEADER_METRIC[metric]
    buf = '<table>\n<tr>'
    buf += ''.join('<th>{}</th>'.format(x) for x in header)
    buf += '</tr>\n'
    for values in results:
        link = str(pathlib.Path(config['log_dir']) / (values[-1] + '.txt'))
        buf += '<tr><td><a href="{}">{}</a></td>'.format(link, values[1])
        buf += ''.join(map(gen_td, values[2:-1]))
        buf += '</tr>\n'
    return buf + '</table>\n'


def gen_div(test, config, results):
    percent = math.floor(100 * pass_ratio(results))
    buf     = HTML_DIV.format(label=LABEL[test], percent=percent)
    for mode in TEST_MODES:
        mode_results = [r for r in results if r[1] == mode]
        if mode_results:
            buf += gen_table(test, mode, config, mode_results)
    return buf + '</div>\n'


def save_html(profile, config, all_results, html_file, wav_only):
    ok    = profile_passed(all_results)
    state = ('failed', 'passed')[ok]
    if wav_only:
        if state == 'passed':
            state = 'All files generated!'
        else:
            state = 'File generation failed!'     
    buf = HTML_HEAD.format(title=profile, style=STYLE, state=state)
    for test in TESTS:
        if test in all_results:
            buf += gen_div(test, config, all_results[test])
    buf += '</body>\n'

    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(buf)


def save_op_cmds(opid_log, log_dir):
    for opid in opid_log.keys():
        opid_file = str(pathlib.Path(log_dir) / (opid + '.txt'))
        with open(opid_file, 'w') as fid:
            fid.write(opid_log[opid])


def main(args):
    args.workers = min(max(args.workers, 1), os.cpu_count())
    time_stamp   = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
    log_file     = '{}_{}.log'.format(args.config.replace('.cfg',''), time_stamp)
    stem_config = str(pathlib.Path(args.config).stem)
    log_dir      = makedirs('{}_{}_log'.format(stem_config, time_stamp))
    opid_log     = {}
    log_handlers = [logging.FileHandler(log_file)]
    if args.verbose:
        log_handlers += [logging.StreamHandler(sys.stdout)]
    logging.basicConfig(level=logging.DEBUG, handlers=log_handlers)
    work_dir = makedirs(pathlib.Path('lc3plus_conformance_' + time_stamp))

    logging.debug(' '.join([str(arg) for arg in sys.argv]))
    if args.test_items_only:
        if not args.system_sox:
            download_sox()
        prepare_items(args.workers)
        print('File preparation complete. All input test files stored in "test_items"')
        sys.exit(0)

    try:
        all_passed = True
        globels, configs = parse_config(args.config)
        profiles = globels['enabled_tests']
        # remove these point for remote debug
        if DEBUG_SETTING == 0:
            check_system(args, globels)
            build_tools()
            if not args.system_sox:
                download_sox()
            prepare_items(args.workers)

        for profile in profiles:
            print('Running tests for "{}" ...'.format(profile))
            config      = configs[profile]
            all_results = {}
            wav_only_dict = {'active':args.wav_only, 'path':'test_files_' + stem_config}
            cmp_only_dict = {'active':args.cmp_only, 'path':'test_files_' + stem_config}
            for test in TESTS:
                test_test = 'test_' + test
                if config[test_test]:
                    test_dir    = makedirs(work_dir / profile / test)
                    test_env    = TestEnv(profile, test, config, test_dir, args.workers,'', opid_log, wav_only_dict, cmp_only_dict)
                    test_func   = globals()[test_test]
                    test_result = test_func(test_env)
                    if not test_result:
                        print('{} in "{}" is enabled with no suitable configuration!'.format(test, profile))
                        all_passed = False
                    all_results[test] = test_result
                    if not args.keep_files:
                        removedir(work_dir / profile / test)

            if all_results:
                all_passed = all_passed and profile_passed(all_results)
                html_file  = '{}_{}.html'.format(profile, time_stamp)
                print('Saving results ...')
                config['log_dir'] = log_dir
                save_html(profile, config, all_results, html_file, args.wav_only)
                save_op_cmds(opid_log, log_dir)
            else:
                print('No tests in "{}" were enabled!'.format(profile))

        print('\nLogfile:', log_file)
        if args.wav_only:
            if all_passed:
                print('\nAll test files created successfully and saved in test_files_{}\n'.format(stem_config))
            else:
                print('\nFile generation process failed! See {}\n'.format(log_file))
        else:
            print('Results:', '         \n'.join('%s_%s.html' % (p, time_stamp) for p in profiles))
            print('\nConformance test', 'passed.' if all_passed else 'failed!', '\n')
        sys.exit(0 if all_passed else 2)
    except KeyboardInterrupt:
        print('\rExiting. Please wait while workers shut down ...')
    finally:
        if not args.keep_files:
            removedir(work_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='LC3plus conformance tool - checks if your version of the LC3plus cod'
                                                 'ec is conform to the reference provided by Fraunhofer & Ericsson.')
    parser.add_argument('-keep', action='store_true', dest='keep_files', help="Don't delete workdir at end of test")
    parser.add_argument('-system_sox', action='store_true', help='Use system sox')
    parser.add_argument('-bit_exact', action='store_true', help='reset thresholds to test for bit exactness')
    parser.add_argument('-v', action='store_true', dest='verbose', help='Activate verbose output')
    parser.add_argument('-w', dest='workers', type=int, default=os.cpu_count(), help='Number of worker threads')
    parser.add_argument('config', help='Conformance config file')
    parser.add_argument('-wav_only', action='store_true', help='encode wav files whithout applying peaq comparison')
    parser.add_argument('-cmp_only', action='store_true', help='use wav files from directory for peaq comparison')
    parser.add_argument('-test_items_only', action='store_true', help='only prepare test_items folder without running conformance. Can be called with dummy.cfg')
    args = parser.parse_args()

    if args.system_sox:
        SOX_EXE = 'sox'
    if args.bit_exact:
        for test, mode in itertools.product(TESTS, TEST_MODES):
            DEFAULTS_TEST['{}_{}_mld_threshold'.format(test, mode)] = 0
            DEFAULTS_TEST['{}_{}_odg_threshold'.format(test, mode)] = [0,0,0]
            DEFAULTS_TEST['{}_{}_rms_threshold'.format(test, mode)] = 16
            DEFAULTS_TEST['{}_{}_mad_threshold'.format(test, mode)] = 0
    if args.cmp_only and args.wav_only:
        print('Error: -wav_only and -cmp_only can only be used exclusively!')
        exit(1)

    main(args)
