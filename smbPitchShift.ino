/****************************************************************************

  NAME: smbPitchShift.cpp
  VERSION: 1.2
  HOME URL: http://www.dspdimension.com
  KNOWN BUGS: none

  SYNOPSIS: Routine for doing pitch shifting while maintaining
  duration using the Short Time Fourier Transform.

  DESCRIPTION: The routine takes a pitchShift factor value which is between 0.5
  (one octave down) and 2. (one octave up). A value of exactly 1 does not change
  the pitch. numSampsToProcess tells the routine how many samples in indata[0...
  numSampsToProcess-1] should be pitch shifted and moved to outdata[0 ...
  numSampsToProcess-1]. The two buffers can be identical (ie. it can process the
  data in-place). fftFrameSize defines the FFT frame size used for the
  processing. Typical values are 1024, 2048 and 4096. It may be any value <=
  MAX_FRAME_LENGTH but it MUST be a power of 2. osamp is the STFT
  oversampling factor which also determines the overlap between adjacent STFT
  frames. It should at least be 4 for moderate scaling ratios. A value of 32 is
  recommended for best quality. sampleRate takes the sample rate for the signal
  in unit Hz, ie. 44100 for 44.1 kHz audio. The data passed to the routine in
  indata[] should be in the range [-1.0, 1.0), which is also the output range
  for the data, make sure you scale the data accordingly (for 16bit signed integers
  you would have to divide (and multiply) by 32768).

  COPYRIGHT 1999-2009 Stephan M. Bernsee <smb [AT] dspdimension [DOT] com>

              The Wide Open License (WOL)

  Permission to use, copy, modify, distribute and sell this software and its
  documentation for any purpose is hereby granted without fee, provided that
  the above copyright notice and this license appear in all source copies.
  THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF
  ANY KIND. See http://www.dspguru.com/wol.htm for more information.

*****************************************************************************/
#include <arm_math.h>
#define MAX_FRAME_LENGTH FFT_SIZE

static float gLastPhasetmp = 0;
static int stepShift = 0;
static float gInFIFO[MAX_FRAME_LENGTH];
//static float gOutFIFO[MAX_FRAME_LENGTH];
static float gFFTworksp[2 * MAX_FRAME_LENGTH];
static float gLastPhase[MAX_FRAME_LENGTH];
static float gSumPhase[MAX_FRAME_LENGTH];
static float gOutputAccum[2 * MAX_FRAME_LENGTH];
static float gAnaFreq[MAX_FRAME_LENGTH];
static float gAnaMagn[MAX_FRAME_LENGTH];
static float gSynFreq[MAX_FRAME_LENGTH];
static float gSynMagn[MAX_FRAME_LENGTH];
static long gRover = false, gInit = false;

const long fftFrameSize2 = FFT_SIZE / 2;
arm_cfft_radix4_instance_f32 fft_inst;

void initSmb() {
  memset(gInFIFO,      0, MAX_FRAME_LENGTH * sizeof(float));
  memset(gFFTworksp,   0, 2 * MAX_FRAME_LENGTH * sizeof(float));
  memset(gLastPhase,   0, (MAX_FRAME_LENGTH) *sizeof(float));
  memset(gSumPhase,    0, (MAX_FRAME_LENGTH) *sizeof(float));
  memset(gOutputAccum, 0, 2 * MAX_FRAME_LENGTH * sizeof(float));
  memset(gAnaFreq,     0, (MAX_FRAME_LENGTH) * sizeof(float));
  memset(gAnaMagn,     0, (MAX_FRAME_LENGTH) * sizeof(float));
  gInit = true;
}
// -----------------------------------------------------------------------------------------------------------------
void smbPitchShift(const float pitchShift, const long numSampsToProcess, const long osamp, const float sampleRate, float *indata, float *outdata) {

  float magn, phase, tmp, window, real, imag;
  long i, k, qpd, index;
  // set up some handy variables
  const long stepSize = FFT_SIZE / osamp;
  const float freqPerBin = sampleRate / (float)FFT_SIZE;
  const float expct = 2.*M_PI * (float)stepSize / (float)FFT_SIZE;
  const float inFifoLatency = FFT_SIZE - stepSize;

  stepShift = 0;

  if (gRover == false) gRover = inFifoLatency;

  for (int b = 0; b < inFifoLatency; b++) {
    gInFIFO[b] = indata[b];
  }

  // main processing loop
  for ( i = inFifoLatency; i < numSampsToProcess; i++) {

    // As long as we have not yet collected enough data just read in
    gInFIFO[gRover] = indata[i];
    gRover++;

    // now we have enough data for processing
    if (gRover >= FFT_SIZE) {

      gRover = inFifoLatency;

      // do windowing and re,im interleave
      for (k = 0; k < FFT_SIZE; k++) {
        //window = -.5 * cosf(2.*M_PI * (float)k / (float)FFT_SIZE) + .5;
        gFFTworksp[2 * k] = gInFIFO[k] * win[k];
        gFFTworksp[2 * k + 1] = 0.;
      }

      // ***************** ANALYSIS *******************

      // do transform
      arm_cfft_radix4_init_f32(&fft_inst, 256, 0, 1);
      arm_cfft_radix4_f32(&fft_inst, gFFTworksp);
      
      // this is the analysis step
      for (k = 0; k <= fftFrameSize2; k++) {

        // de-interlace FFT buffer
        real = gFFTworksp[2 * k];
        imag = gFFTworksp[2 * k + 1];

        // compute magnitude and phase
        magn = 2.*sqrtf(real * real + imag * imag);
        phase = atan2f(imag, real);

        // compute phase difference
        tmp = phase - gLastPhase[k];
        gLastPhase[k] = phase;

        // subtract expected phase difference
        tmp -= (float)k * expct;

        // map delta phase into +/- Pi interval
        qpd = tmp / M_PI;
        if (qpd >= 0) qpd += qpd & 1;
        else qpd -= qpd & 1;
        tmp -= M_PI * (float)qpd;

        // get deviation from bin frequency from the +/- Pi interval
        tmp = osamp * tmp / (2.*M_PI);

        // compute the k-th partials' true frequency
        tmp = (float)k * freqPerBin + tmp * freqPerBin;

        // store magnitude and true frequency in analysis arrays
        gAnaMagn[k] = magn;
        gAnaFreq[k] = tmp;
      }

      // ***************** PROCESSING *******************
      // this does the actual pitch shifting
      memset(gSynMagn, 0, fftFrameSize2 * sizeof(float));
      memset(gSynFreq, 0, fftFrameSize2 * sizeof(float));
      for (k = 0; k <= fftFrameSize2; k++) {
        index = k * pitchShift;

        if (index <= fftFrameSize2) {
          gSynMagn[index] += gAnaMagn[k];
          gSynFreq[index] = gAnaFreq[k] * pitchShift;
        }
      }

      // ***************** SYNTHESIS *******************
      // this is the synthesis step

      for (k = 0; k <= fftFrameSize2; k++) {

        // get magnitude and true frequency from synthesis arrays
        magn = gSynMagn[k];
        tmp = gSynFreq[k];

        // subtract bin mid frequency
        tmp -= (float)k * freqPerBin;

        // get bin deviation from freq deviation
        tmp /= freqPerBin;

        // take osamp into account
        tmp = 2.*M_PI * tmp / osamp;

        // add the overlap phase advance back in
        tmp += (float)k * expct;

        // accumulate delta phase to get bin phase
        gSumPhase[k] += tmp;
        phase = gSumPhase[k];

        // get real and imag part and re-interleave
        gFFTworksp[2 * k] = magn * arm_cos_f32(phase);//cosf(phase);
        gFFTworksp[2 * k + 1] = magn * arm_sin_f32(phase);//;//sinf(phase);
      }

      // zero negative frequencies
      for (k = FFT_SIZE + 2; k < 2 * FFT_SIZE; k++) gFFTworksp[k] = 0.;

      // do inverse transform
      arm_cfft_radix4_init_f32(&fft_inst, 256, 1, 1);
      arm_cfft_radix4_f32(&fft_inst, gFFTworksp);

      // do windowing and add to output accumulator
      for (k = 0; k < FFT_SIZE; k++) {
        //window = -.5 * cosf(2.*M_PI * (float)k / (float)FFT_SIZE) + .5;
        //window = -.5 * cosf(2.*M_PI * (float)k / (float)FFT_SIZE) + .5;
        gOutputAccum[k] += 2.*win[k] * (gFFTworksp[2 * k]*256) / (fftFrameSize2 * osamp);
      }

      // shift output data buffer
      for (k = 0; k < stepSize; k++) {
        outdata[k + (stepShift * stepSize)] = gOutputAccum[k];
      }
      stepShift++;

      // shift accumulator
      memmove(gOutputAccum, gOutputAccum + stepSize, FFT_SIZE * sizeof(float));

      // move input FIFO
      for (k = 0; k < inFifoLatency; k++) gInFIFO[k] = gInFIFO[k + stepSize];

    }
  }
}
