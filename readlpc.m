%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%readlpcc.m
%Program author: Akihiro Inui
%Extract LPC from audio data
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Output is LPC Mean value and standard deviation
% Input values are Frame Rate: F, Label:label, Number of Coefficient:lpcccoeff
% This file has to be set in a directory has sound files with lpcc.m

function [Cmean,Cstd] = readlpc(F,label,lpcccoeff)

% Read directory into D
D = dir;

% Number of audio files
Nf = length(D)-4;

% Create a empty vector and Matrix to input LPCC
l = zeros(1,lpcccoeff+1);
Cmean = zeros(Nf,lpcccoeff+1);
Mdsum = zeros(1,lpcccoeff);
Cstd = zeros(length(D)-4,lpcccoeff);

% Main loop
countr = 0;

% sampling frequency
Fs = 44100;

for n = 3:length(D)-2
    names = getfield(D(n),'name');          % Read names of files
    [x,Fs] = audioread(names);
    lp = lpcc(x,Fs,F,lpcccoeff);            % Extract LPC
    l = lp(2:lpcccoeff+1,:);                % Extract the first row 
    l = l';
    [r,i] = size(lp);
    Mean = mean(l,'omitnan');               % Calculate mean vlaues for each frame(Ignore silence)
    
    % Calculate standard deviation
    for q =1:i
        Cdelta(:,q) = lp(2:lpcccoeff+1,q)-Mean';         % Departures from mean values
    end
    Cd = Cdelta.*Cdelta;
    
    for p =1:lpcccoeff
    Mdsum(1,p) = sum(Cd(p,:),'omitnan');                % Sum the all squared departures
    end
    
    Cstd(n-2,1:lpcccoeff) = sqrt((1/(r-1))*Mdsum(1,:)); % Compute standard deviation and store
    Cmean(n-2,1:lpcccoeff) = Mean(1,:);                 % Store mean values of each frame
    countr = countr+1;
end

% Label data 
Cmean(:,lpcccoeff+1) = label;
Cstd(:,lpcccoeff+1) = label;

end