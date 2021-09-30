# ================================================================
# LC3plus Precision Requirement script for High Resolution V1.0.0
# ================================================================

# /******************************************************************************
# *                        ETSI TS 103 634 V1.3.1                               *
# *              Low Complexity Communication Codec Plus (LC3plus)              *
# *                                                                             *
# * Copyright licence is solely granted through ETSI Intellectual Property      *
# * Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
# * estoppel or otherwise.                                                      *
# ******************************************************************************/

import sys
import shlex
import subprocess
import wave
import logging
from thd_plus_n import thd_plus_n, generate_sinusoid, snr
import math as m
import sys
import os
import shutil
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import configparser
import multiprocessing
from multiprocessing import Process, Queue
import datetime

rdict = {"freq":0, "thd": 0, "snr": 0}

HTML_HEADER = """<!DOCTYPE html><head><meta charset="UTF-8"><title>{cfg_file_name} Report</title>
<style>body {{font-family:sans-serif; color:#f8f8f2; background-color:#272822; font-size:80%}} div {{border:1px solid #8f908a; border-radius:4px; overflow:hidden; display:table; margin-left:30px; margin-bottom:30px}} h2 {{text-align:left; margin-left:30px}} h3 {{text-align:left; margin:4px}} table {{border-spacing:0px; width:100%}} th {{padding:4px; border-top:1px solid #8f908a}} td {{padding:4px}} tr:nth-child(even) {{background-color:rgba(255,255,255,0.1)}} td.pass {{background-color:rgba(0,192,255,0.4)}} td.fail {{background-color:rgba(255,0,0,0.4)}} td.warn{{background-color:rgba(214,137,16,0.4)}} a:link {{color: white;text-decoration: none;}} a:visited {{color: white}} a:hover {{color: white;font-weight:bold;}} a:active {{color: rgb(155, 155, 155)}}</style></head>
<body><h2>Precision requirements for {cfg_file_name} TEST_STATUS</h2>
<div><h3>THD/SNR Test - worst case and 1kHz values</h3>
<table>
<tr><th>Mode</th><th>Frame Size</th><th>Samplerate</th><th>Bitrate</th><th>SNR [dB]</th><th>SNR 1kHz [dB]</th><th>THD [dB]</th><th>THD 1kHz [dB]</th></tr>"""
HTML_ROW = '<tr><td>{mode}</td><td>{frame_ms}</td><td>{fs}</td><td>{bitrate}</td><td class="{snr_results_status}">{snr} ({snr_threshold})</td><td class="{snr_results_status}">{snr_1k} ({snr_threshold_1k})</td><td class="{thd_result_status}">{thd} ({thd_threshold})</td><td class="{thd_result_status}">{thd_1k} ({thd_threshold_1k})</td></tr>'
HTML_ROW_PNG = '<tr><td>{mode}</td><td>{frame_ms}</td><td>{fs}</td><td>{bitrate}</td><td class="{snr_results_status}"><a href="{png_path}">{snr} ({snr_threshold})</a></td><td class="{snr_results_status_1k}"><a href="{png_path}">{snr_1k} ({snr_threshold_1k})</a></td><td class="{thd_result_status}"><a href="{png_path}">{thd} ({thd_threshold})</a></td><td class="{thd_result_status_1k}"><a href="{png_path}">{thd_1k} ({thd_threshold_1k})</a></td></tr>'
GLOBAL_KEYS = ['reference_encoder','reference_decoder','test_encoder','test_decoder','enc_thd_threshold','enc_snr_threshold','enc_thd_threshold_1k','enc_snr_threshold_1k','dec_thd_threshold','dec_snr_threshold','dec_thd_threshold_1k','dec_snr_threshold_1k','num_cores','plot',]
CONFIGS = ['test_mode','frame_ms','fs','bitrate']

# Run command and return output. cmd can be string or list. Commands with .exe suffix are automatically
# called with wine unless wine=False. Set unicode=False to get binary output. Set hard_fail=False to
# to ignore nonzero return codes.
def call(cmd, wine=True, unicode=True, hard_fail=True, log_output=True):
    if isinstance(cmd, str):
        cmd = [x for x in shlex.split(cmd) if x]
    if sys.platform != 'cygwin' and wine and cmd[0].lower().endswith('.exe'):
        cmd = ['wine'] + cmd
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=unicode)
    out = p.communicate()[0] or (b'', '')[unicode]
    quoted_cmd = ' '.join(map(shlex.quote, cmd))
    logging.debug(quoted_cmd)
    if unicode and log_output:
        logging.debug(out)
    if hard_fail and p.returncode != 0:
        raise OSError(quoted_cmd + ' failed!\n' + out )
        exit(1)
    return out


def bytes_to_float_24(b):
    maxval = 1<<23
    period = 1<<24
    scale = 1/((1<<23)-1)
    
    ints = [scale * float((((int(b[3*n]) | int(b[3*n+1]) << 8 | int(b[3*n+2]) << 16) + maxval) % period) - maxval)
            for n in range(len(b)//3)] 

    return ints


def float_to_bytes_24(x):
    mask = 255
    scale = (1<<23)-1
    out = bytes()
    for samp in x:
        val = int(round(scale*samp))
        out += bytes([val & mask, (val >> 8) & mask, (val >> 16) & mask])
    return out


def get_op_string(op):
    ops = {}
    for key in op:
        ops[key] = str(op[key])
    return '_'.join([ops['test_mode'],ops['frame_ms'],ops['fs'],ops['bitrate']])


def measure_thd_plus_n(call_string_enc, call_string_dec, freq, op):
    y = generate_sinusoid(freq, op['fs'], 1, -3)
    infile = op['wav_folder'] + '/' + "sin_{}Hz.wav".format(freq)
    operation_point = get_op_string(op)
    bitstream = op['wav_folder'] + '/' +"sin_{}Hz_{}.g192".format(freq, operation_point)
    outfile = op['wav_folder'] + '/' +"sin_{}Hz_{}_dec.wav".format(freq, operation_point)
    outfile_raw = op['wav_folder'] + '/' +"sin_{}Hz_{}_dec.raw".format(freq, operation_point)
    outfile2 = op['wav_folder'] + '/' +"sin_{}Hz_{}_dec_preproc.wav".format(freq, operation_point)
    outfile2_raw = op['wav_folder'] + '/' +"sin_{}Hz_{}_dec_preproc.raw".format(freq, operation_point)
    sinctaps = 10000
    
    # 24 bit PCM
    with wave.open(infile, 'wb') as wf:
        wf.setparams((1,3,op['fs'],len(y),'NONE', 'not compressed'))
        b = float_to_bytes_24(y)
        wf.writeframes(b)

    if op['bps'] == 32:
        tmp = '%s.tmp.wav' % (infile)
        sox_cmd = 'sox %s -e floating-point -b 32 %s' % (infile, tmp)
        log = call(sox_cmd)
        log = call('mv %s %s' % (tmp, infile))

    log = call(call_string_enc.format(infile=infile, bitstream=bitstream))
    log = call(call_string_dec.format(bitstream=bitstream, outfile=outfile))

    sox_cmd = 'sox {infile} {outfile} '.format(infile = outfile, outfile=outfile2)

    sox_log = call(sox_cmd)

    if op['bps'] == 32:
        tmp = '%s.tmp.wav' % (outfile2)
        log = call('sox %s -b 24 -e signed-integer %s' % (outfile2, tmp))
        shutil.move(tmp, outfile2)

        tmp = '%s.tmp.wav' % (outfile)
        log = call('sox %s -b 24 -e signed-integer %s' % (outfile, tmp))
        shutil.move(tmp, outfile)

    log = call('sox %s %s' % (outfile2, outfile2_raw))
    log = call('sox %s %s' % (outfile, outfile_raw))

    fid = open(outfile2_raw, 'rb')
    b = fid.read()
    y_out2 = bytes_to_float_24(b)

    fid = open(outfile_raw, 'rb')
    b = fid.read()
    y_out = bytes_to_float_24(b)

    t = thd_plus_n(y_out2[sinctaps:-sinctaps], y_out[sinctaps:-sinctaps], freq, op['fs'])
    s = snr(y[4800:-4800], y_out[4800:-4800])
                      
    return freq, t, s
    

def usage(msg):
    print("error: {}".format(msg))
    print("usage: python3 precision_requirements.py LC3plus_precision.cfg")
    sys.exit(1)


def result_to_html(html, op_result, op, png_path):
    snr_results_status, thd_result_status  = 'pass', 'pass'
    snr_results_status_1k, thd_result_status_1k  = 'pass', 'pass'

    if op_result['snr'] < op_result['snr_threshold']: 
        snr_results_status = 'fail'
    if op_result['thd'] > op_result['thd_threshold']:
        thd_result_status = 'fail'
    if op_result['snr_1k'] < op_result['snr_threshold_1k']: 
        snr_results_status_1k = 'fail'
    if op_result['thd_1k'] > op_result['thd_threshold_1k']:
        thd_result_status_1k = 'fail'
    if png_path:
        new_row = HTML_ROW_PNG.format(mode=op['test_mode'], frame_ms=op['frame_ms'], fs=op['fs'], \
                              bitrate=op['bitrate'], snr_results_status=snr_results_status, snr_results_status_1k=snr_results_status_1k,\
                              snr=round(op_result['snr'],2), snr_1k=round(op_result['snr_1k'],2), snr_threshold=op_result['snr_threshold'], snr_threshold_1k=op_result['snr_threshold_1k'],\
                              thd_result_status=thd_result_status, thd_result_status_1k=thd_result_status_1k, thd=round(op_result['thd'],2), thd_1k=round(op_result['thd_1k'],2),\
                              thd_threshold=op_result['thd_threshold'], thd_threshold_1k=op_result['thd_threshold_1k'], png_path=png_path)
    else:
        new_row = HTML_ROW.format(mode=op['test_mode'], frame_ms=op['frame_ms'], fs=op['fs'], \
                                bitrate=op['bitrate'], snr_results_status=snr_results_status, snr_results_status_1k=snr_results_status_1k, \
                                snr=round(op_result['snr'],2), snr_1k=round(op_result['snr_1k'],2), snr_threshold=op_result['snr_threshold'], snr_threshold_1k=op_result['snr_threshold_1k'], \
                                thd_result_status=thd_result_status, thd_result_status_1k=thd_result_status_1k, thd=round(op_result['thd'],2), thd_1k=round(op_result['thd_1k'],2),\
                                thd_threshold=op_result['thd_threshold'], thd_threshold_1k=op_result['thd_threshold_1k'])
    op_status = 'passed'
    if 'fail' in [snr_results_status, thd_result_status, snr_results_status_1k, thd_result_status_1k]:
        op_status = 'failed'

    return (html + new_row, op_status)


def measure_thd_plus_n_parallel(call_string_enc, call_string_dec, freq, op):
    freq, thd, snr = measure_thd_plus_n(call_string_enc, call_string_dec, freq, op)
    rdict["thd"] = thd
    rdict["snr"] = snr
    rdict["freq"] = freq
    q.put(rdict.copy())


def plot(op, freqs, tvec, svec, png_folder):
    title = 'LC3plus High Resolution: {}, {}ms, {}kHz, {:.1f}kbps'.format(op['test_mode'], op['frame_ms'], op['fs']//1000, op['bitrate']/1000)
    filename = 'thd+n_{}_{}ms_{}kHz_{:.1f}kbps'.format(op['test_mode'], op['frame_ms'], op['fs']//1000, op['bitrate']/1000)
    filename += '.png'

    fig, ax = plt.subplots()
    fig.suptitle(title)
    ax.set_xscale('log')
    ax.plot(freqs, tvec, label='THD + N')
    ax.plot(freqs, [-s for s in svec],'.' ,label = '-SNR')
    ax.legend()
    plt.xlabel('Hz')
    plt.ylabel('dB')
    if not os.path.isdir(png_folder):
        os.makedirs(png_folder)
    plt.savefig(png_folder+'/'+filename)
    print('data plotted to {}'.format(png_folder+'/'+filename))
    return png_folder+'/'+filename


def parse_globals(config):
    missing_keys = set(GLOBAL_KEYS) - set(config['globals'].keys())
    if bool(missing_keys):
        print('Global key "{}" missing.'.format(list(missing_keys)[0]))
        exit()

    globals = {}
    for key in GLOBAL_KEYS:
        globals[key] = config['globals'][key]

    if not os.path.isfile(globals['reference_encoder'].split(' ')[0]):
        print('reference_encoder "{}" not found.'.format(globals['reference_encoder']))
        exit()
    if not os.path.isfile(globals['reference_decoder'].split(' ')[0]):
        print('reference_decoder "{}" not found.'.format(globals['reference_decoder']))
        exit()
    if not os.path.isfile(globals['test_encoder'].split(' ')[0]):
        print('test_encoder "{}" not found.'.format(globals['test_encoder']))
        exit()
    if not os.path.isfile(globals['test_decoder'].split(' ')[0]):
        print('test_decoder "{}" not found.'.format(globals['test_decoder']))
        exit()

    if not 'freq_grid' in config['globals'].keys(): 
        globals['freq_grid'] = 'linlog'
    else:
        globals['freq_grid'] = config['globals']['freq_grid']

    if not 'freq_grid_steps' in config['globals'].keys():
        globals['freq_grid_steps'] = 200
    else:
        globals['freq_grid_steps'] = config['globals']['freq_grid_steps'].format()

    if globals['freq_grid'] == 'linlog':
        print('Using frequency grid with linear spacing in log domain with {} points'.format(globals['freq_grid_steps']))
    elif globals['freq_grid'] == 'onethirdoctave':
        print('Using frequency grid of one-third octave spacing')
    else:
        print('Frequency grid {} unkown, use linlog or onethirdoctave'.format(globals['freq_grid']))
        exit()

    if not 'bps' in config['globals'].keys(): 
        globals['bps'] = 24
    else:
        globals['bps'] = int(config['globals']['bps'])

    if int(globals['num_cores']) > multiprocessing.cpu_count():
        globals['num_cores'] = multiprocessing.cpu_count()
    else:
        globals['num_cores'] = int(globals['num_cores'])
    print('Using {} processing cores'.format(globals['num_cores']))

    return globals


def parse_op(line):
    op = {} # operation point
    try:
        op['test_mode'] = line.split(',')[0].strip()
        op['frame_ms'] = float(line.split(',')[1].strip())
        op['fs'] = int(line.split(',')[2].strip())
        op['bitrate'] = int(line.split(',')[3].strip())
    except IndexError:
        print('Not enough test entries in config.')

    if op['test_mode'] not in ['encode', 'decode', 'encdec']:
        print('Mode must be (encode|decode|encdec)!')
        exit()
    return op


def print_op_status(op):
    op_print = 'Running... '
    for key in CONFIGS:
        op_print += '{}: {}, '.format(key,op[key])
    print(op_print)


def select_threshold(test_mode, xxx_threshold):
    key_enc = 'enc_' + xxx_threshold
    key_dec = 'dec_' + xxx_threshold
    if test_mode == 'encode':
        return float(globals[key_enc])
    if test_mode == 'decode':
        return float(globals[key_dec])
    if test_mode == 'encdec': # select the more critical of enc_threshold/dec_threshold
        if xxx_threshold[:3] == 'snr':
            return max(float(globals[key_enc]), float(globals[key_dec]))
        elif xxx_threshold[:3] == 'thd':
            return min(float(globals[key_enc]), float(globals[key_dec]))


q = Queue()
if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage("not enough input arguments")

    if sys.argv[1] == '-h':
        usage("Printing help")

    conf_path = sys.argv[1]
    if not os.path.isfile(conf_path):
        usage('could not find {}'.format(conf_path))
    
    # logging
    time_stamp   = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
    log_file     = 'precision_requirements_{}.log'.format(time_stamp)
    log_handlers = [logging.FileHandler(log_file)]
    logging.basicConfig(level=logging.DEBUG, handlers=log_handlers)
    logging.debug(' '.join([str(arg) for arg in sys.argv]))

    # parse [globals] config
    config = configparser.ConfigParser()
    config.read(conf_path)

    globals = parse_globals(config)

    cfg_file_name = conf_path.split('/')[-1].rsplit('.',1)[0]
    html = HTML_HEADER.format(cfg_file_name=conf_path.split('/')[-1])
    wav_folder = 'wav_' + cfg_file_name
    index_1kHz = 0
    
    if os.path.isdir(wav_folder):
        print('removing wave folder from previous run: {}'.format(wav_folder))
        shutil.rmtree(wav_folder)
    os.makedirs(wav_folder)
    # parse [tests] config
    tests = config['tests']['configs']
    test_status = 'passed'
    for line in tests.split('\n'):
        op = parse_op(line)
        print_op_status(op)
       
        if globals['freq_grid'] == 'linlog':
            # n point log space
            nSteps = int (globals['freq_grid_steps'])
            stepSize = (m.log10(op['fs']/2) - 1) / nSteps
            freqs = [10**(n*stepSize+1) for n in range(nSteps)]
            index_1k = m.floor ((m.log10(1000)-1)/stepSize)
            if freqs[index_1k] != 1000:
                freqs.insert (index_1k+1, 1000)
                index_1k = index_1k + 1

        elif globals['freq_grid'] == 'onethirdoctave':
            # one-third octave steps
            max_steps = int (m.log10(op['fs']/2)*10)
            freqs = [10 ** (0.1 * ((n+12))) for n in range(max_steps-12+1)]
            index_1k = int (m.log10(1000)*10) - 12
        else:
            exit ('unknow frequency grid')
        
        op['wav_folder'] = wav_folder
        op['bps'] = globals['bps']

        if op['test_mode'] == 'encode':
            encoder_exe = globals['test_encoder']
            decoder_exe = globals['reference_decoder']
        if op['test_mode'] == 'decode':
            encoder_exe = globals['reference_encoder']
            decoder_exe = globals['test_decoder']
        if op['test_mode'] == 'encdec':
            encoder_exe = globals['test_encoder']
            decoder_exe = globals['test_decoder']

        call_string_enc = encoder_exe.format(frame_ms = op['frame_ms'], infile = '{infile}', bitstream = '{bitstream}',bitrate = op['bitrate'])
        call_string_dec = decoder_exe.format(bps = globals['bps'], bitstream = '{bitstream}', outfile = '{outfile}')

        rdictlist = []
        plist = []
        cores_in_use = 0

        try:
            for freq in freqs:
                p = Process(target = measure_thd_plus_n_parallel, args=(call_string_enc, call_string_dec, freq, op))
                p.start()
                plist.append(p)
                cores_in_use = cores_in_use + 1
                if (cores_in_use == globals['num_cores'] or freq == freqs[-1]):
                    for p in plist:
                        p.join()
                        rdict = q.get()
                        rdictlist.append(rdict.copy())
                        cores_in_use = cores_in_use - 1
                        print("freq: {:9.3f}  SNR: {:6.2f}  THD: {:6.2f}".format(rdict["freq"], rdict["snr"], rdict["thd"]))
                            
                    plist = []
        except KeyboardInterrupt:
            print('Received keyboard interrupt, stopping')
        
        rdictlist_sorted = sorted(rdictlist, key = lambda k: k["freq"])
        freqs = [d["freq"] for d in rdictlist_sorted]
        tvec = [d["thd"] for d in rdictlist_sorted]
        svec = [d["snr"] for d in rdictlist_sorted]

        op_result = {}
        op_result['snr'] = min(svec)
        op_result['thd'] = max(tvec)
        op_result['snr_1k'] = svec[index_1k]
        op_result['thd_1k'] = tvec[index_1k]

        op_result['snr_threshold'] = select_threshold(op['test_mode'], 'snr_threshold')
        op_result['thd_threshold'] = select_threshold(op['test_mode'], 'thd_threshold')
        op_result['snr_threshold_1k'] = select_threshold(op['test_mode'], 'snr_threshold_1k')
        op_result['thd_threshold_1k'] = select_threshold(op['test_mode'], 'thd_threshold_1k')

        png_path = ''    
        if globals['plot']:
            # plot figure
            png_path = plot(op, freqs, tvec, svec, 'png_'+cfg_file_name)
        
        html, op_status = result_to_html(html, op_result, op, png_path)
        if op_status == 'failed':
            test_status = 'failed'
    html = html.replace('TEST_STATUS', test_status)
    html += '</table></div></body>'
    html_file_name = cfg_file_name+'.html'
    with open(html_file_name,'w') as fid:
        fid.writelines(html)
    print('results saved in {}'.format(html_file_name))
    #shutil.rmtree(wav_folder)
