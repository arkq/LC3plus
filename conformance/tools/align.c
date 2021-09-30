/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "tinywavein_c.h"
#include "tinywaveout_c.h"


#define MAX_CHANNEL_NUMBER 64

static void printUsage(void);
void WriteWav(WAVEFILEOUT* out_file, int *sample_buf, int nSamples, int bps);


int main(int ac, char * av[])
{
    char     *in_filename = NULL, *out_filename = NULL;
    int       samplerate = 0, channels = 0, input_len = 0, bps = 0;
    WAVEFILEIN  *in_file;
    WAVEFILEOUT *out_file;
    int sample_buf[MAX_CHANNEL_NUMBER];
    int zero_buf[MAX_CHANNEL_NUMBER];
    int nSamplesRead = 0, ref_len = 0, delay = 0, count = 0, end_1, end_2;


    if ( ac < 4 )
    {
        printUsage();
    }
    
    in_filename = av[1];
    out_filename = av[2];
    ref_len = atoi(av[3]);
    delay = atoi(av[4]);
    
    if ( delay < 0 )
    {
        printf("delay must be a positiv value!\n");
        exit(1);
    }
    
    in_file = OpenWav(in_filename, &samplerate, &channels, &input_len, &bps);
    out_file = CreateWav(out_filename, samplerate, channels, bps);

    // set parameter for writing output file
    if (ref_len > input_len - delay){
        end_1 = input_len;
    }
    else if (ref_len <= input_len - delay){
        end_1 = ref_len + delay;
    }
    end_2 = ref_len + delay;

    for (int i = 0; i < channels; i++) { zero_buf[i] = 0; }

    while(1)
    {
        if (count < input_len){
            ReadWavShort(in_file, sample_buf, channels, &nSamplesRead); //read only one sample per iteration
        }
        if (delay <= count && count < end_1){
            WriteWav(out_file, sample_buf, channels, bps);
        }
        if (end_1 <= count && count < end_2){
            WriteWav(out_file, zero_buf, channels, bps);
        }        
        count ++;
        if (count > end_2){
            break;
        }
    }
    CloseWavIn(in_file);
    CloseWav(out_file);
    printf("Produced aligned file: %s\n", out_filename);
    printf("cut of: %d samples at the beginning\n", delay);
    printf("cut of: %d samples at the end\n", input_len - end_1);
    int padded = 0;
    if (end_2 - input_len > 0 ){ padded = end_2 - input_len;}
    printf("padded: %d zeros at the end\n", padded);
    return 0;
}

void printUsage(void)
{
    printf("   Usage: ./align infile.wav outfile.wav ref_len delay \n\n");
    printf("   infile.wav  : coded test file [16b/24b,mono/stereo] with known delay, but possibly different lenght \n");
    printf("   outfile.wav : aligned file, will be trimmed/padded with 0s to maintain ref_len\n");
    printf("   ref_len     : length of the reference file to be aligned with\n");
    printf("   delay       : known delay/lag between infile.wav and reference file\n\n");

    printf("   Procedure:\n");
    printf("   Cut of 'delay' samples of infile.wav at the beginning and fill up with zeros or trim overlap at the end to maintain length 'ref_len'.\n");
    exit(1);
}

void WriteWav(WAVEFILEOUT* out_file, int *sample_buf, int nSamples, int bps)
{
    if ( bps == 16 ){
        short sample_buf_short[nSamples];
        for (int i = 0; i < nSamples; i++ ){
            sample_buf_short[i] = (short) sample_buf[i];
        }
        WriteWavShort(out_file, sample_buf_short, nSamples);
        //WriteWavShort(out_file, (short *) sample_buf, nSamples);
    }
    else if (bps == 24 ){
        WriteWavLong(out_file, sample_buf, nSamples);
    }
    else {
        printf("  Only 16 or 24 bits per samples supported!");
        exit(1);
    }
}
