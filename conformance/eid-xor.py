#!/usr/bin/env python3
"""
eid-xor.py - Version 1.2 (Python port)

Program Description:
~~~~~~~~~~~~~~~~~~~~

This example program performs a "logical" exclusive-or operation of
a file (representing an encoded speech bitstream) with another file
(representing a bit error pattern).

The file containing an encoded speech bitstream can be in a compact
binary format, in the G.192 serial bitstream format (which uses
16-bit softbits), or in the byte-oriented G.192 format.

The file containing the error pattern will be in one of three
possible formats: G.192 16-bit softbit format (without synchronism
header for bit errors), byte-oriented version of the G.192 format,
and compact, hard-bit binary (bit) mode.

Original C version by:
    Simao Ferraz de Campos Neto
    Comsat Laboratories

Python port maintains exact compatibility with the C version output.
"""

import sys
import struct
import argparse
from enum import IntEnum
from typing import BinaryIO, Callable, List, Optional, Tuple

# Constants from softbit.h
G192_ZERO = 0x007F
G192_ONE = 0x0081
G192_FER = 0x6B20  # Frame erasure
G192_SYNC = 0x6B21  # Good frame (sync header)

# Buffer length
EID_BUFFER_LENGTH = 256


class Format(IntEnum):
    G192 = 0
    BYTE = 1
    COMPACT = 2


class ErrorType(IntEnum):
    BER = 0  # Bit Error Rate
    FER = 1  # Frame Erasure Rate


def format_str(fmt: int) -> str:
    """Return string representation of format."""
    names = ["g192", "byte", "compact"]
    if 0 <= fmt < len(names):
        return names[fmt]
    return "unknown"


def type_str(t: int) -> str:
    """Return string representation of error type."""
    return "FER" if t == ErrorType.FER else "BER"


def eid_xor(a: int, b: int) -> int:
    """Perform XOR operation on two samples."""
    return G192_ONE if (a ^ b) else G192_ZERO


def insert_errors(a: List[int], b: List[int], n: int) -> Tuple[List[int], int]:
    """
    Insert errors by XOR-ing the input data arrays.
    Returns the disturbed data and count of disturbed bits.
    """
    c = []
    disturbed = 0
    for i in range(n):
        bit = eid_xor(a[i], b[i])
        if bit != a[i]:
            disturbed += 1
        c.append(bit)
    return c, disturbed


def read_g192(f: BinaryIO, count: int) -> List[int]:
    """Read G.192 format data (16-bit little-endian)."""
    data = f.read(count * 2)
    if not data:
        return []
    # Unpack as unsigned shorts (little-endian)
    n = len(data) // 2
    return list(struct.unpack(f'<{n}H', data[:n*2]))


def read_byte(f: BinaryIO, count: int) -> List[int]:
    """Read byte format data and convert to G192 equivalents.
    
    Byte format uses lower byte of G192 values:
    - 0x7F -> G192_ZERO (0x007F)
    - 0x81 -> G192_ONE (0x0081)
    - 0x21 -> G192_SYNC (0x6B21)
    - 0x20 -> G192_FER (0x6B20)
    """
    data = f.read(count)
    if not data:
        return []
    
    result = []
    for b in data:
        if b == 0x7F:
            result.append(G192_ZERO)
        elif b == 0x81:
            result.append(G192_ONE)
        elif b == 0x21:
            result.append(G192_SYNC)
        elif b == 0x20:
            result.append(G192_FER)
        else:
            # For other values, try to infer based on the pattern
            # FER indicators have upper nibble 0x20, softbits don't
            if (b & 0xF0) == 0x20:
                result.append(0x6B00 | b)  # FER-style indicator
            else:
                result.append(b)  # Keep as-is (softbit style)
    return result


def read_bit_ber(f: BinaryIO, count: int) -> List[int]:
    """Read compact bit format for BER (bit errors)."""
    result = []
    bytes_needed = (count + 7) // 8
    data = f.read(bytes_needed)
    if not data:
        return []
    
    for byte in data:
        for bit_pos in range(8):
            if len(result) >= count:
                break
            # LSB first
            bit = (byte >> bit_pos) & 1
            result.append(G192_ONE if bit else G192_ZERO)
    
    return result[:count]


def read_bit_fer(f: BinaryIO, count: int) -> List[int]:
    """Read compact bit format for FER (frame erasures)."""
    result = []
    bytes_needed = (count + 7) // 8
    data = f.read(bytes_needed)
    if not data:
        return []
    
    for byte in data:
        for bit_pos in range(8):
            if len(result) >= count:
                break
            # LSB first
            bit = (byte >> bit_pos) & 1
            result.append(G192_FER if bit else G192_SYNC)
    
    return result[:count]


def save_g192(f: BinaryIO, data: List[int]) -> int:
    """Save G.192 format data (16-bit little-endian)."""
    packed = struct.pack(f'<{len(data)}H', *data)
    f.write(packed)
    return len(data)


def save_byte(f: BinaryIO, data: List[int]) -> int:
    """Save byte format data, converting from G192 equivalents.
    
    Converts G192 values to their byte equivalents:
    - G192_ZERO (0x007F) -> 0x7F
    - G192_ONE (0x0081) -> 0x81
    - G192_SYNC (0x6B21) -> 0x21
    - G192_FER (0x6B20) -> 0x20
    """
    result = []
    for d in data:
        if d == G192_ZERO:
            result.append(0x7F)
        elif d == G192_ONE:
            result.append(0x81)
        elif d == G192_SYNC:
            result.append(0x21)
        elif d == G192_FER:
            result.append(0x20)
        elif (d & 0xFF00) == 0x6B00:
            result.append(d & 0xFF)  # FER-style indicator
        else:
            result.append(d & 0xFF)  # Keep lower byte
    f.write(bytes(result))
    return len(data)


def save_bit(f: BinaryIO, data: List[int]) -> int:
    """Save compact bit format data."""
    result = []
    for i in range(0, len(data), 8):
        byte = 0
        for bit_pos in range(8):
            if i + bit_pos < len(data):
                # Check if it's a '1' (either G192_ONE or error/erasure)
                val = data[i + bit_pos]
                if val == G192_ONE or val == G192_FER:
                    byte |= (1 << bit_pos)
        result.append(byte)
    f.write(bytes(result))
    return len(data)


def check_eid_format(f: BinaryIO, filename: str) -> Tuple[int, int]:
    """
    Check the format and type of an EID file.
    Returns (format, type).
    
    This mimics the behavior in the C version's check_eid_format from softbit.c
    """
    # Save position
    pos = f.tell()
    
    # Read enough data to analyze - need to check multiple samples
    header = f.read(10)
    if len(header) < 2:
        f.seek(pos)
        return Format.COMPACT, ErrorType.BER
    
    # Check for G.192 format (16-bit) by looking at multiple samples
    if len(header) >= 4:
        word1 = struct.unpack('<H', header[:2])[0]
        word2 = struct.unpack('<H', header[2:4])[0]
        
        # Check for sync header (FER type indicators) - must be consistent
        if (word1 & 0xFFF0) == 0x6B20:
            f.seek(pos)
            return Format.G192, ErrorType.FER
        
        # Check for softbit values (BER type indicators)
        # Both first samples should be valid softbits
        if (word1 == G192_ZERO or word1 == G192_ONE) and \
           (word2 == G192_ZERO or word2 == G192_ONE):
            f.seek(pos)
            return Format.G192, ErrorType.BER
    
    # Check for byte format by examining multiple bytes
    if len(header) >= 2:
        byte1 = header[0]
        byte2 = header[1]
        
        # Check for byte-mode sync header (more strict check)
        # Sync byte should be 0x21 (good) or 0x20 (bad)
        if byte1 == 0x21 or byte1 == 0x20:
            f.seek(pos)
            return Format.BYTE, ErrorType.FER
        
        # Check for byte-mode softbit values - need multiple matches
        if (byte1 == 0x7F or byte1 == 0x81) and \
           (byte2 == 0x7F or byte2 == 0x81):
            f.seek(pos)
            return Format.BYTE, ErrorType.BER
    
    # Default to compact format
    f.seek(pos)
    return Format.COMPACT, ErrorType.BER


def display_usage(level: int = 0):
    """Display program usage."""
    print("eid-xor.py - Version 1.2 (Python port)\n")
    
    if level:
        print("""Program Description:

This example program performs a "logical" exclusive-or operation of
a file (representing an encoded speech bitstream) with another file
(representing a bit error pattern.

The file containing an encoded speech bitstream can be in a compact
binary format, in the G.192 serial bitstream format (which uses
16-bit softbits), or in the byte-oriented G.192 format.

The file containing the error pattern will be in one of three
possible formats: G.192 16-bit softbit format (without synchronism
header for bit errors), byte-oriented version of the G.192 format,
and compact, hard-bit binary (bit) mode. These are described in the
following.

The headerless G.192 serial bitstream format is as described in
G.192, with the exceptions listed below. The main feature is that
the softbits and frame erasure indicators are right-aligned at
16-bit word boundaries (unsigned short): 
'0'=0x007F and '1'=0x0081, and good/bad frame = 0x6B21/0x6B20

In the byte-oriented softbit serial bitstream, only the lower byte
of the softbits defined in G.192 are used. Hence:
'0'=0x7F and '1'=0x81, and good/bad frame = 0x21/0x20

In the compact (bit) mode, only hard bits are saved. Each byte will
have information about eight bits or frames. The LBbs will refer to
bits or frames that occur first in time. Here, '1' means that a bit
is in error or that a frame should be erased, and a '0', otherwise.

Conventions:
~~~~~~~~~~~~

Bitstreams can be disturbed in two ways: by bit errors, or by frame
erasures. The STL EID supports three basic modes: random/bit errors
(labeled BER), simple frame erasure (labeled FER), and Bellcore
model burst frame erasure (labeled BFER). Here are some conventions
that apply to the particular formats for each of these three EID
operating modes.

BER: bitstream generated by this program are composed of bits 1/0,
     *without* synchronism headers or any other frame delimitation
     (i.e., only bits affecting the payload are present). Frame
     boundaries are defined by the user's application only. The
     following applies:
     G.192 mode: file will contain either 0x007F (no disturbance) or
                0x0081 (bit error)
     Byte mode:  file will contain either 0x7F (no disturbance) or
                0x81 (bit error)
     Compact mode: each bit in the file will indicate whether a
                disturbance occurred (bit 1) or not (bit 0).
                Lower order bits apply to bits occurring first
                in time.

FER/BFER: bitstream generate by this program is composed only by
     the indication of whether a frame should be erased or not. No
     payload is present. The following applies:
     G.192 mode: file will contain either 0x6B21 (no disturbance) or
                0x6B20 (frame erasure)
     Byte mode:  file will contain either 0x21 (no disturbance) or
                0x20 (frame erasure)
     Compact mode: each bit in the file will indicate whether a frame
                erasure occurred (bit 1) or not (bit 0). Lower order
                bits apply to bits occurring first in time.
""")
    else:
        print("Program to insert bit errors and frame erasures in bitstream")
        print("files using a previously generated error pattern. Three formats")
        print("are acceptable: g192, byte, and (compact) bit.\n")
    
    print("""Usage:
eid-xor.py [Options] in_bs err_pat_bs out_bs
Where:
 in_bs ...... input encoded speech bitstream file
 err_pat .... error pattern bitstream file
 out_bs ..... disturbed encoded speech bitstream file    

Options:
 -frame # ... Set the frame size to # (for headerless G.192
              bitstreams or for compact binary files).
 -bs mode ... Mode for bitstream (g192, byte, or bit)
 -ep mode ... Mode for error pattern (g192, byte, or bit)
 -ber ....... Error pattern is a bit error pattern (needed for bit format)
 -fer ....... Error pattern is a frame erasure pattern (for bit format)
 -vbr ....... Enables variable bit rate operation (different frame sizes)
 -q ......... Quiet operation
 -? ......... Displays this message
 -help ...... Displays a complete help message""")
    sys.exit(-128 & 0xFF)


def parse_format(s: str) -> int:
    """Parse format string to Format enum value."""
    s = s.lower()
    if 'g192' in s:
        return Format.G192
    elif 'byte' in s:
        return Format.BYTE
    elif 'bit' in s or 'compact' in s:
        return Format.COMPACT
    return -1


def main():
    # Parse arguments manually to match C behavior
    args = sys.argv[1:]
    
    # Default values
    ep_type = ErrorType.BER
    bs_format = Format.G192
    ep_format = Format.G192
    fr_len = 0
    vbr = False
    quiet = False
    start_frame = 1
    
    # Parse options
    positional = []
    i = 0
    while i < len(args):
        arg = args[i]
        if arg == '-start':
            start_frame = int(args[i + 1])
            i += 2
        elif arg == '-frame':
            fr_len = int(args[i + 1])
            i += 2
        elif arg == '-bs':
            fmt = parse_format(args[i + 1])
            if fmt < 0:
                sys.stderr.write("Invalid BS format type. Aborted\n")
                sys.exit(5)
            bs_format = fmt
            i += 2
        elif arg == '-ep':
            fmt = parse_format(args[i + 1])
            if fmt < 0:
                sys.stderr.write("Invalid error pattern format type. Aborted\n")
                sys.exit(5)
            ep_format = fmt
            i += 2
        elif arg in ['-ber', '-BER']:
            ep_type = ErrorType.BER
            i += 1
        elif arg in ['-fer', '-FER', '-bfer', '-BFER']:
            ep_type = ErrorType.FER
            i += 1
        elif arg == '-vbr':
            vbr = True
            i += 1
        elif arg == '-q':
            quiet = True
            i += 1
        elif arg == '-?':
            display_usage(0)
        elif arg == '-help' or arg == '--help':
            display_usage(1)
        elif arg.startswith('-'):
            sys.stderr.write(f'ERROR! Invalid option "{arg}" in command line\n\n')
            display_usage(0)
        else:
            positional.append(arg)
            i += 1
    
    if len(positional) < 3:
        if len(positional) == 0:
            display_usage(0)
        # Interactive mode - prompt for missing files
        if len(positional) < 1:
            positional.append(input("_Input bit stream file ..................: "))
        if len(positional) < 2:
            positional.append(input("_Error pattern file .....................: "))
        if len(positional) < 3:
            positional.append(input("_Output bit stream file .................: "))
    
    ibs_file = positional[0]
    ep_file = positional[1]
    obs_file = positional[2]
    
    # Print file names like the C tool does
    if not quiet:
        sys.stderr.write(f"_Input bit stream file ..................: {ibs_file}\n")
        sys.stderr.write(f"_Error pattern file .....................: {ep_file}\n")
        sys.stderr.write(f"_Output bit stream file .................: {obs_file}\n")
    
    # Starting frame is from 0 to number_of_frames-1
    start_frame -= 1
    
    # Open files
    try:
        Fibs = open(ibs_file, 'rb')
    except IOError:
        sys.stderr.write("Could not open input bitstream file\n")
        sys.exit(1)
    
    try:
        Fep = open(ep_file, 'rb')
    except IOError:
        sys.stderr.write("Could not open error pattern file\n")
        sys.exit(1)
    
    try:
        Fobs = open(obs_file, 'wb')
    except IOError:
        sys.stderr.write("Could not create output file\n")
        sys.exit(1)
    
    # Check consistency - inspect input bitstream format
    detected_format, tmp_type = check_eid_format(Fibs, ibs_file)
    
    if detected_format != bs_format:
        sys.stderr.write(f"*** Switching bitstream format from {format_str(bs_format)} to {format_str(detected_format)} ***\n")
        bs_format = detected_format
    
    # Check for sync header in bitstream
    sync_header = False
    if tmp_type == ErrorType.FER:
        if bs_format == Format.G192:
            # Check for G.192 sync header
            tmp = read_g192(Fibs, 2)
            if len(tmp) >= 2:
                frame_len_1 = tmp[1]
                # Seek to next expected header position
                Fibs.seek(tmp[1] * 2, 1)  # SEEK_CUR
                tmp2 = read_g192(Fibs, 2)
                if len(tmp2) >= 2 and (tmp2[0] & 0xFFF0) == 0x6B20:
                    fr_len = frame_len_1
                    sync_header = True
                    if frame_len_1 != tmp2[1]:
                        vbr = True
            Fibs.seek(0)
        elif bs_format == Format.BYTE:
            # Check for byte-wise sync header
            data = Fibs.read(2)
            if len(data) >= 2:
                frame_len_1 = data[1]
                Fibs.seek(data[1], 1)
                data2 = Fibs.read(2)
                if len(data2) >= 2 and (data2[0] & 0xF0) == 0x20:
                    fr_len = frame_len_1
                    sync_header = True
                    if frame_len_1 != data2[1]:
                        vbr = True
            Fibs.seek(0)
    
    # Default frame length
    blk = EID_BUFFER_LENGTH
    if fr_len == 0:
        fr_len = blk
    
    # Check error pattern format
    detected_ep_format, detected_ep_type = check_eid_format(Fep, ep_file)
    
    if detected_ep_format != ep_format:
        sys.stderr.write(f"*** Switching error pattern format from {format_str(ep_format)} to {format_str(detected_ep_format)} ***\n")
        ep_format = detected_ep_format
    
    if detected_ep_type != ep_type:
        if ep_format == Format.COMPACT:
            sys.stderr.write(f"*** Cannot infer error pattern type. Using {type_str(ep_type)} ***\n")
        else:
            sys.stderr.write(f"*** Switching error pattern type from {type_str(ep_type)} to {type_str(detected_ep_type)} ***\n")
            ep_type = detected_ep_type
    
    # VBR operation checks
    if vbr and (bs_format == Format.COMPACT or not sync_header):
        vbr = False
        if bs_format == Format.COMPACT:
            sys.stderr.write("VBR operation disabled for compact bitstreams!\n")
        else:
            sys.stderr.write("VBR operation disabled for headerless bitstreams!\n")
    
    # Output format check
    if bs_format == Format.COMPACT and ep_type == ErrorType.FER:
        obs_format = Format.G192
    else:
        obs_format = bs_format
    
    # Reset sync for compact BS
    if bs_format == Format.COMPACT and sync_header:
        sys.stderr.write("*** Disabling SYNC header for compact bitstream ***\n")
        sync_header = False
    
    # Select I/O functions
    if bs_format == Format.BYTE:
        read_data = read_byte
    elif bs_format == Format.G192:
        read_data = read_g192
    else:
        read_data = read_bit_ber
    
    if ep_format == Format.BYTE:
        read_patt = read_byte
    elif ep_format == Format.G192:
        read_patt = read_g192
    elif ep_type == ErrorType.BER:
        read_patt = read_bit_ber
    else:
        read_patt = read_bit_fer
    
    if obs_format == Format.BYTE:
        save_data = save_byte
    elif obs_format == Format.G192:
        save_data = save_g192
    else:
        save_data = save_bit
    
    # Sample size
    ibs_sample_len = 1 if bs_format == Format.BYTE else (2 if bs_format == Format.G192 else 0)
    
    # VBR: scan for max frame size
    if vbr:
        max_fr_len = fr_len
        while True:
            Fibs.seek(ibs_sample_len, 1)
            offset_data = read_data(Fibs, 1)
            if not offset_data:
                break
            offset = offset_data[0]
            if offset > max_fr_len:
                max_fr_len = offset
            Fibs.seek(offset * ibs_sample_len, 1)
        Fibs.seek(0)
        fr_len = max_fr_len
    
    # Frame lengths
    bs_len = fr_len + 2 if sync_header else fr_len
    ep_len = fr_len
    
    ori_bs_len = bs_len
    ori_fr_len = fr_len
    
    # Prepare erased frame template
    erased_frame = [0] * bs_len
    if sync_header:
        erased_frame[0] = G192_FER
        erased_frame[1] = fr_len
    
    # Counters
    disturbed = 0.0
    processed = 0.0
    wraps = 0
    
    # Main processing
    if ep_type == ErrorType.FER:
        k = 0
        ep_buffer = []
        ep_true_len = 0
        
        while True:
            # Read one frame from BS
            if vbr:
                bs = read_data(Fibs, 2)
                if len(bs) != 2:
                    break
                fr_len = bs[1]
                bs_len = fr_len + 2 if sync_header else fr_len
                if fr_len != 0:
                    payload = read_data(Fibs, fr_len)
                    bs.extend(payload)
            else:
                bs = read_data(Fibs, bs_len)
            
            if not bs:
                break
            
            items = len(bs)
            
            if items < bs_len:
                if sync_header:
                    sys.stderr.write("*** File size for this bitstream file is not multiple  ***\n")
                    sys.stderr.write("*** of the given frame length. Check that the correct  ***\n")
                    sys.stderr.write("*** frame size was used (is this a variable-frame size ***\n")
                    sys.stderr.write("*** file?) and that the bitstream is not corrupted.***\n")
                    sys.exit(9)
                else:
                    sys.stderr.write("*** File size for this HEADERLESS bitstream is not ***\n")
                    sys.stderr.write("*** multiple of the given frame length. Check that ***\n")
                    sys.stderr.write("*** the correct frame size was selected & that the ***\n")
                    sys.stderr.write("*** bitstream file is not corrupted.***\n")
                    bs_len = fr_len = items
            
            # Read erasure flags
            while k == 0:
                ep_buffer = read_patt(Fep, ep_len)
                ep_true_len = k = len(ep_buffer)
                if k <= 0:
                    Fep.seek(0)
                    wraps += 1
            
            processed += 1
            
            # Save original or erased frame
            ep_idx = ep_true_len - k
            if ep_buffer[ep_idx] == G192_FER:
                if vbr:
                    erased_frame[1] = fr_len
                # Adjust erased frame size if needed
                ef = erased_frame[:bs_len] if len(erased_frame) >= bs_len else erased_frame + [0] * (bs_len - len(erased_frame))
                save_data(Fobs, ef)
                disturbed += 1
            else:
                save_data(Fobs, bs)
            
            k -= 1
    
    else:  # BER
        while True:
            # Read one frame from BS
            if vbr:
                bs = read_data(Fibs, 2)
                if len(bs) != 2:
                    break
                fr_len = bs[1]
                bs_len = fr_len + 2 if sync_header else fr_len
                if fr_len != 0:
                    payload = read_data(Fibs, fr_len)
                    bs.extend(payload)
            else:
                bs = read_data(Fibs, bs_len)
            
            if not bs:
                break
            
            items = len(bs)
            
            if items < bs_len:
                if sync_header:
                    sys.stderr.write("*** File size for this bitstream file is not multiple  ***\n")
                    sys.stderr.write("*** of the given frame length. Check that the correct  ***\n")
                    sys.stderr.write("*** frame size was used (is this a variable-frame size ***\n")
                    sys.stderr.write("*** file?) and that the bitstream is not corrupted.***\n")
                    sys.exit(9)
                else:
                    bs_len = fr_len = items
            
            # Read error pattern
            ep = read_patt(Fep, ep_len)
            
            if len(ep) < ep_len:
                if len(ep) < 0:
                    sys.stderr.write(f"Error reading {ep_file}\n")
                    sys.exit(7)
                k = ep_len - len(ep)
                Fep.seek(0)
                ep.extend(read_patt(Fep, k))
                wraps += 1
            
            # Extract payload
            if sync_header:
                header = bs[:2]
                payload = bs[2:]
            else:
                header = []
                payload = bs
            
            # Insert errors
            payload_len = min(len(payload), len(ep), fr_len)
            new_payload, err_count = insert_errors(payload[:payload_len], ep[:payload_len], payload_len)
            
            # Update counters
            disturbed += err_count
            processed += fr_len
            
            # Reconstruct and save frame
            output = header + new_payload
            if len(output) < bs_len:
                output.extend(payload[payload_len:])
            
            save_data(Fobs, output)
    
    # Restore frame lengths
    bs_len = ori_bs_len
    fr_len = ori_fr_len
    
    # Print summary
    header_str = "(G.192 header) " if sync_header else "(headerless) .."
    sys.stderr.write(f"# Bitstream format {header_str}...... : {format_str(bs_format)}\n")
    if bs_format != obs_format:
        sys.stderr.write(f"# Out bitstream format {header_str}.. : {format_str(obs_format)}\n")
    
    ep_type_str = "(frame erasure) " if ep_type == ErrorType.FER else "(bit error) ...."
    sys.stderr.write(f"# Pattern format {ep_type_str}....... : {format_str(ep_format)}\n")
    sys.stderr.write(f"# Error pattern files wrapped ...........: {wraps} times\n")
    sys.stderr.write(f"# Frame size ............................: {fr_len}\n")
    
    unit = "frames " if ep_type == ErrorType.FER else "bits .."
    sys.stderr.write(f"# Processed {unit}..................... : {processed:.0f} \n")
    sys.stderr.write(f"# Distorted {unit}..................... : {disturbed:.0f} \n")
    
    rate_label = "Frame erasure rate" if ep_type == ErrorType.FER else "Bit error rate ..."
    if processed > 0:
        rate = 100.0 * disturbed / processed
    else:
        rate = 0.0
    sys.stderr.write(f"# {rate_label}.....................: {rate:f} %\n")
    
    # Close files
    Fibs.close()
    Fep.close()
    Fobs.close()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())