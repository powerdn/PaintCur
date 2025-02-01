function Hd = mylowpass2
%MYLOWPASS2 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.3 and DSP System Toolbox 9.5.
% Generated on: 14-Nov-2020 01:27:27

% Butterworth Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are in Hz.
Fs = 400e3;  % Sampling Frequency

Fpass = 50;          % Passband Frequency
Fstop = 90;          % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 60;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);

% [EOF]
