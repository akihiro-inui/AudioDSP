%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%am.m
%Program author: Akihiro Inui
%Program details: Implement Amplitude Modulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear all; close all; clc;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 1. Preamble and variable declaration
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Read audio file
[x,Fs] = audioread('filename');

% Normalise the input signal
xn = x/max(abs(x(:))); 

% Sample period
Ts = 1/Fs;

% Create a time vector
t = [0:Ts:max(size(x))/Fs - Ts]';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 2. Amplitude Modulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Amplitude modulation
A=1;     % Amplitude
fc=0.75; % Carrier frequency
a= A*cos(2*pi*fc*t);
y=(1+a).*xn;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 3. Visualise signals
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figfont = 16;
figure
plot(t,a,t,xn,t,y)
xlabel('Time (s)','FontSize',figfont)
ylabel('Amplitude (arbitrary)','FontSize',figfont)
title(['Amplitude modulation'],'FontSize',figfont)
l1 = legend('Modulator','Original','Amplitude modulated');
set(l1,'FontSize',figfont-2)
set(gca,'FontSize',figfont-2)
grid on

% Play back the unmodulated and modulated sounds
%soundsc(x,Fs)
%soundsc(y,Fs)