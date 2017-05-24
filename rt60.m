%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%rt60.m
%Program author: Akihiro Inui
%Program details: Calculate RT 60
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear all; close all; clc;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 1. Preamble and variable declaration
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Read audio file
[x,Fs] = audioread('filename');

% Reverberation time(T1=RT30,T2=RT60)
T1 = -5;
T2 = -35;
T3 = -60;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 2. Calculate Reverberation time
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Calculate sum of input signal and flip
y = flipud(cumsum(flipud(x .^ 2)));

% Create time vector
t = ([0:length(y)-1]/Fs)';

% Find a maximum value of y in Db expression
m = max(10*log10(y));

% Set the start point to zero
z = 10*log10(y)-m;

% Sum the signal by -5, -35 dB from 0dB
ind = [sum(z >= T1)  sum(z >= T2)];

% Linear Approximation
p = polyfit(z(ind(1):ind(2)), t(ind(1):ind(2)), 1);

% Derive RT60
rt = polyval(p, T3);