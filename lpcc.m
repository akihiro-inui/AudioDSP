%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%lpcc.m
%Program author: Akihiro Inui
%read linear predictive coefficients
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% This file has to be in a directory has sound files with readlpc.m
% Output is lpc

function C = lpcc(input,samplingRate,frameRate,lpcccoeff)


% Avoid error
if (nargin < 2) samplingRate = 16000; end
if (nargin < 3) frameRate = 100; end

% Make a window
windowSize = 256;
hamWindow = 0.54 - 0.46*cos(2*pi*(0:windowSize-1)/windowSize);

% Apply preemphasis filter to the input signal
preEmphasized = filter([1 -.97], 1, input);

% Calculate the number of columns
windowStep = samplingRate/frameRate;
cols = fix((length(input)-windowSize)/windowStep); 

% Set empty matrix to store entire output signal
C = zeros(lpcccoeff+1, cols);

% Initialise variables
a = 1;
b = a + windowStep;

% Calculate lpc coefficients to each segment
for start=0:cols-1
    if b >= length(preEmphasized)   % If the frame is out of maximum length
       b = length(preEmphasized);   % Set the last frame to the end point
       a = b-windowStep;
    end
    seg = preEmphasized(a:b);       % Take out signals from each segment
    if seg(:) == 0
    c = zeros(1,lpcccoeff+1);
    else
    c = lpc(seg,lpcccoeff);         % Calculate lpc
    end
    C(:,start+1) = c';              % Store lpc into matrix        
    a = b+1;                        % Update variables
    b = a + windowStep;
end



    