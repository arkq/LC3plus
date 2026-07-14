#******************************************************************************
#                        ETSI TS 103 634 V1.7.1                               *
#              Low Complexity Communication Codec Plus (LC3plus)              *
#                                                                             *
# Copyright licence is solely granted through ETSI Intellectual Property      *
# Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
# estoppel or otherwise.                                                      *
#*****************************************************************************/

import argparse 
import numpy as np
from scipy.io import wavfile
from itertools import pairwise
import sys


def check_lossless_samples(Nchannels, WavIn, DecWavOut):
    
    fs, WavInSamples = wavfile.read(WavIn)
    fs, DecWavOutSamples = wavfile.read(DecWavOut)
    
    nSamplesNotLossless = np.sum( np.abs(WavInSamples - DecWavOutSamples) > 0 )
    
    if Nchannels == 2:
        losslessSamples = 2*len(WavInSamples) - nSamplesNotLossless
        ratio =  losslessSamples / (2*len(WavInSamples))
    else:
        losslessSamples = len(WavInSamples) - nSamplesNotLossless
        ratio =  losslessSamples / len(WavInSamples)  
    
    if ratio == 1:
        print(f'lossless samples:\t{ratio * 100}%')
    else:
        print(f'lossless samples:\t{(ratio * 100):.10f}%')
    print(f'non-lossless samples:\t{nSamplesNotLossless}')


def parse_file_header(file):

    file_identifier = file.read(2)                              # file identifier, value 0xcc1c
    header_size = int.from_bytes(file.read(2),"little")         # total config header size in bytes
    samplerate = int.from_bytes(file.read(2),"little") * 100    # sample rate 
    bitrate = int.from_bytes(file.read(2),"little") * 100       # bitrate 
    Nchannels = int.from_bytes(file.read(2),"little")           # Nchannels
    frame_ms = int.from_bytes(file.read(2),"little") / 100      # frame duration in ms
    epmode = int.from_bytes(file.read(2),"little")              # eps mode
    signal_len =  int.from_bytes(file.read(2),"little")         # signal len
    signal_len_2 = int.from_bytes(file.read(2),"little")        # signal len 2
        
    if header_size > 18:
            hr_mode = int.from_bytes(file.read(2),"little")         #hrmode
        
    if header_size > 20:      
            wav_index =  int.from_bytes(file.read(2),"little")      #bits per sample / wavindex 
            if (wav_index == 2):
                wavformat = 24
            elif(wav_index == 1):
                wavformat = 16
            else:
                wavformat = 16
    else:
        wavformat = 16

    return  samplerate, bitrate, Nchannels,frame_ms, wavformat


def read_bitstream(bitstream, samplerate, bitrate, frame_ms, Nchannels, padding, g192):
        
        payload_counts = []
        bytesFull      = []
        paddingFull    = []


        if samplerate == 44100:
            frame_ms_scaled = frame_ms * (480 / 441)
        else:
            frame_ms_scaled = frame_ms

        while True:

            if g192:
                #read SNYC_WORD
                SNYC_WORD = int.from_bytes(bitstream.read(2),"little")
                if SNYC_WORD == 0:
                    break
                if SNYC_WORD != 0x6b20 and SNYC_WORD != 0x6b21:
                    print()
                    print("ERROR: Snyc Word of g192 frame could not be interpreted")
                    sys.exit()
                #read frame length 
                nBytes = int.from_bytes(bitstream.read(2),"little")                

            else:                  
                #read frame length 
                nBytes = int.from_bytes(bitstream.read(2),"little")
                          
            if not nBytes:
                    break
                
            #skip the rest 
            if g192:
                bits = bitstream.read(int(nBytes * 2))[::-1]
            else:
                bits = bitstream.read(nBytes)

            
            if padding:   
                bytesPerFrame = int(bitrate / 8000 * frame_ms_scaled)
               
                #if Nchannels == 2:
                #    bytesPerFrame = bytesPerFrame / 2
                
                bytesFull.append(bytesPerFrame)                       
                sumPadding = 0
                byte = 0
                if g192:
                        for i in range(0,len(bits),16):
                            #pick always 2 bytes
                            FrameBytes =  bytes(reversed(bits[i:i+16]))
                            for j in range(0,len(FrameBytes),2):
                                TwoBytes = FrameBytes[j:j+2]
                                if TwoBytes.hex() == '8100':
                                    byte |= 1 << int(j / 2)
                            if byte == 0x07:
                                sumPadding += 1
                                byte = 0
                            else:
                                break
                else:
                    for bit in reversed(bits):
                        if bit == 0x07: 
                            sumPadding += 1
                        else:
                            break  
                                                   
                #nBytes = nBytes - sumPadding
                if g192:
                    payload_counts.append(8 * (int(nBytes / 8) - sumPadding))
                    paddingFull.append(sumPadding)
                else:
                    payload_counts.append(8*(nBytes - sumPadding))
                    paddingFull.append(sumPadding)                   
                
            else:
                if g192:
                    payload_counts.append( nBytes )
                else: 
                    payload_counts.append( 8 * nBytes)
                                
        if Nchannels == 2:
            combined = [payload_counts[i] + payload_counts[i + 1] for i in range(0, len(payload_counts), 2)]
            payload_counts = combined

        return payload_counts, frame_ms_scaled, paddingFull, bytesFull


def analyse_bitstream(path, outfile, padding, g192, g192_configfile):
    
    payload_counts  = []
    if padding:
        paddingFull = []
        bytesFull   = []
    
    if g192:
        if not g192_configfile:
            if len(path) == 2:
                if '_L' in path[0]:
                    g192_configfile = path[0].replace("_L","")  + '.cfg'
                elif '_R' in path[0]:
                    g192_configfile = path[0].replace("_R","")  + '.cfg'
            else:
                g192_configfile = path[0] + '.cfg' 
        with open(g192_configfile, 'rb') as file:
            samplerate, bitrate, Nchannels, frame_ms, wavformat = parse_file_header(file)
        
        if len(path) == 2:
            with open(path[0], 'rb') as file:
                payload_counts_L, frame_ms_scaled, paddingFull_L, bytesFull_L = read_bitstream(file, samplerate, bitrate, frame_ms, 1, padding, g192)
            file.close()
            with open(path[1], 'rb') as file:
                payload_counts_R, frame_ms_scaled, paddingFull_R, bytesFull_R = read_bitstream(file, samplerate, bitrate, frame_ms, 1, padding, g192)
            file.close()

            #combine channels
            payload_counts = payload_counts_L + payload_counts_R
            paddingFull = paddingFull_L + paddingFull_R
            bytesFull = bytesFull_L + bytesFull_R
        
        else:
            with open(path[0], 'rb') as file:
                payload_counts, frame_ms_scaled, paddingFull, bytesFull = read_bitstream(file, samplerate, bitrate, frame_ms, Nchannels, padding, g192)
    else:
        with open(path[0], 'rb') as file:
            samplerate, bitrate, Nchannels, frame_ms, wavformat = parse_file_header(file)
            payload_counts, frame_ms_scaled, paddingFull, bytesFull = read_bitstream(file, samplerate,bitrate,frame_ms,Nchannels,padding,0)

    max_payload = (1000 / frame_ms_scaled) * max(payload_counts) if payload_counts else 0
    avg_payload = (1000 / frame_ms_scaled) * sum(payload_counts) / len(payload_counts) if payload_counts else 0
          
    print()
    print(f'nChannels:\t\t{Nchannels}')
    print(f'frame duration:\t\t{frame_ms} ms')
    print(f'samplerate\t\t{samplerate} Hz')
    print(f'nFrames:\t\t{len(payload_counts)}')     

    print(f'Average bitrate:\t{avg_payload:.3f} bps')
    print(f'Maximum bitrate:\t{max_payload:.3f} bps')
    print(f'Average compression:\t{(avg_payload / (Nchannels * samplerate * wavformat)):.10}')
        
    if padding:
        
        if Nchannels == 2:
             sumBytesFull = sum(bytesFull) / 2
        else:
             sumBytesFull =  sum(bytesFull)

        print(f'Padded Bytes:\t\t{(100 * ( sum(paddingFull) / sumBytesFull ) ):.3f} %')
        print(f'Padding Offset:\t\t{int(abs(sum(paddingFull) + int(sum(payload_counts) / 8) - sumBytesFull ))}')
            
        
    if outfile != None:
        payload_bytes = np.array(payload_counts) / 8
        np.savetxt(outfile, payload_bytes, fmt='%d')
            
    return Nchannels
        

if __name__ == "__main__":
    
    print ('###############################################################################')
    print ('#                           ETSI TS 103 634 V1.7.1                            #')
    print ('#                         LC3plus bitstream analyser                          #')
    print ('# Copyright licence is solely granted through ETSI Intellectual Property      #')
    print ('# Rights Policy, 3rd April 2019. No patent licence is granted by implication, #')
    print ('# estoppel or otherwise.                                                      #')                                                                      #')
    print ('###############################################################################')
    
    parser = argparse.ArgumentParser(prog='bitstream_analyser.py',
                                     description='Program to extract audio frames and get number of lossless samples')

    parser.add_argument('--input','-i',
                        type=str,
                        nargs='*',
                        help='Path to LC3plus bitstream',
                        required=True
                        )
    
    parser.add_argument('--g192Format','-g192',
                        help='Bitstream format is interpreted as g192',
                        required=False,
                        action='store_true'
                        )
    
    parser.add_argument('--g192ConfigFile','-g192cfg',
                        type=str,
                        help='Config file of g192 Bitstream',
                        required=False
                        )
    
    parser.add_argument('--WavIn','-wi',
                        type=str,
                        help='Path to reference wavfile',
                        required=False
                        )
    
    parser.add_argument('--DecWavOut','-dwi',
                        type=str,
                        help='Path to decoded output file',
                        required=False
                        )
    parser.add_argument('--FrameLenInfoFile','-fl',
                        type=str,
                        help='Path to csv file containing all frame lengths in bytes',
                        required=False,
                        default=None
                        )
    parser.add_argument('--padding','-pd',
                        help='Switch to padding mode',
                        required=False,
                        action='store_true'
                        )
    

    args = parser.parse_args()
    if args.g192Format:
        print(f'Analysing g192 bitstream!\n')
    if len(args.input) == 2 and not args.g192Format:
        print(f'Two bitstreams given but g192 format not set!')
        sys.exit(-1)
    if len(args.input) > 2 or len(args.input) < 1:
        print(f'Wrong number of input files ({args.input}). One input file allowed or two input files for g192 stereo analysis')
        sys.exit(-1)
    if len(args.input) == 1:
        pass #args.input = args.input[0]
    Nchannels = analyse_bitstream(args.input, args.FrameLenInfoFile, args.padding, args.g192Format,args.g192ConfigFile)
    if args.WavIn and args.DecWavOut:
        check_lossless_samples(Nchannels,args.WavIn,args.DecWavOut)
            
        
    print()

