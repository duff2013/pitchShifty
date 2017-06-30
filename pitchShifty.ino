// S. M. Bernsee pitch shift algorithm for teensy
/*******************************************
 *           editable items
 ******************************************/
#define FFT_SIZE    256
#define OVER_SAMPLE 16
// input waveform
const float frequency  = 440;
const float sampleRate = 44100;
const float amplitude  = 1; 
/******************************************/
// phase accumulator
float phi         = 0;
// phase increment per sample
const float delta = 2 * PI * frequency / sampleRate;

float input[FFT_SIZE * 8];
float output[FFT_SIZE * 8];

void setup() {
  // semitones to shift 12 is one octive
  int semitones = 12;
  // length of the input buffer
  int bufferLengthFrames = FFT_SIZE * 8;
  // shift value for algorithm
  float pitchShift = powf(2., semitones / 12.);
  
  while (!Serial);
  delay(120);
  
  for (int i = 0; i < FFT_SIZE * 8; i++) {
    input[i] = (amplitude * sinf(phi));
    // increment phase accumulator
    phi += delta;
  }
  /*
   * smbPitchShift params:
   * 1: "pitchShift" - > semitones to shift up
   * 2: "bufferLengthFrames" -> number of samples in input buffer must be larger than FFT_SIZE
   * 3: "FFT_SIZE" -> size of the FFT, needs to be a power of 2
   * 4: "OVER_SAMPLE" -> fifo buffer overlap factor, more the better but slower, has to be divisable by FFT_SIZE
   * 5: "sampleRate" -> sample rate for sin generation
   * 6: "input" -> input buffer
   * 7: "output" -> output buffer
   */
  smbPitchShift(pitchShift, bufferLengthFrames, FFT_SIZE, OVER_SAMPLE, sampleRate, input, output);
  
  for (int i = 0; i < bufferLengthFrames - FFT_SIZE; i++) {
    Serial.printf("%f,%f\n", input[i], output[i]);
    delay(10);
  }
}

void loop() {
  // put your main code here, to run repeatedly:

}
