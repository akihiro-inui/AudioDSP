%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Karplus_Strong.m
%Program author: Akihiro Inui
%Program details: Implement Karplus-Strong Algorithm
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 1. Preamble and variable declaration
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Sample rate
Fs = 44100;
% Sample period
Ts = 1/Fs;
% Fundamental frequency
f0 = 110;
% Length of delay line
N = round(Fs/f0 - 0.5);
% Loss factor
rho = 0.95;
% Dynamics filter coefficient
R = 0.95;
% duration of simulation(sec)
tEnd = 3;
% length of output vector
M = tEnd*Fs;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 2. Create input vector
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create white noise input
v = -1 + (1+1)*rand(M,1);

% Create filtered vector u
u = zeros(N,1);

% Initialise output
u(1) = (1-R)*v(1);

% Calculate dynamics filter
for n = 2:N
        u(n) = (1-R)*v(n) + R*u(n-1);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 3. Create input/output vector
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Create input vector
x = zeros(N,1);
% Create output vector
y = zeros(N+M,1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 4. Calculate Karplus
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The equation of Karplus-Strong is                   %
% y(n) = x(n)+(rho/2)*(y(n-N)+y(n-(N+1)))             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implement Karplus-Strong Algorithm
% for n < N, x(n)=u(n) and for n >= N, x(n)= 0

% n < N
for n = 1:N-1
        x(n) = u(n);
        y(n) = x(n);
end

% n = N, N+1
y(N) = (rho/2)*y(1);
y(N+1) = y(N);

%n >= N
for n = N+2:N+M
    y(n) = rho*(y(n-N)+y(n-(N+1)))/2;
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 5. Listen to the output signal
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Normalise output
MaxAmp = max(abs(y));
y = y/MaxAmp;

% Sound
soundsc(y,Fs)
