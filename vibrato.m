%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%vibrato.m
%Program author: Akihiro Inui
%Program details: Read audio file and apply simple vibrate
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
close all; clear all;  clc;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 1. Preamble and variable declaration
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Read audio file 
[x,Fs] = audioread('filename');
 
% Sample rate
Ts = 1/Fs;

% Frequency for sine wave
f0 = 2;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 2. Apply Vibrato to y
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Length of signal in smaples
n = 0:length(x);

% Create sine wave for vibrato
vib = exp(-2*pi*f0*n/Fs).*sin(2*pi*f0*n/Fs); 

% Initialise output
y = zeros(1,length(x));

% Multiple the sinusoid to y
for i = 1:length(x)
y(i)= x(i)*vib(i);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 3. Listen and visualise the output signal
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Normalise output
MaxAmp = max(abs(y));
y = y/MaxAmp;

% Listen to the output signal
soundsc(y,Fs)

% Create a time vector
t = (0:Ts:max(size(x))/Fs - Ts)';

% Plot output signal in frequency domain
plot(t,y)
title('Vibrato output')
xlabel('Time (s)')
ylabel('Amplitude')
