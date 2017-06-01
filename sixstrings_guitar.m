%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%sixstrings_guitar.m
%Program author: Akihiro Inui
%Simulate six strings guitar
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear all; close all; clc;

%%%%% flags

plot_on = 0;                % in-loop plotting on (1) or off (0)
itype = 2;                  % type of input: 1: struck, 2: plucked
bctype = 1;                 % boundary condition type: 1: simply supported, 2: clamped
outtype = 1;                % output type: 1: string displacement, 2: string velocity

%%%%% parameters
% tension (N)
T = [123.88 135.915 147.19 145.32 120.48 105.407];       
% string radius (m)
r = [0.00045 0.00042 0.00040 0.00034 0.00041 0.00030];  
% Young's modulus (Pa)
E = [1.92e11 1.92e11 1.92e11 1.92e11 1.92e11 1.92e11]; 
% density (kg/m^3)
rho = [17600 12500 8200 6250 2250 2050]; 
% T60 (s)
T60 = [4 4 4 4 4 4];  
% length (m)
L = [0.647 0.647 0.647 0.647 0.647 0.647];  
% I/O
SR = 44100;                 % sample rate (Hz)
k = 1/SR;                   % time step
Tf = 3.0;                  % duration of simulation (s)

xi = 0.1;                   % coordinate of excitation (normalised, 0-1)
famp = 1;                   % peak amplitude of excitation (N)
dur = 0.001;                % duration of excitation (s)
exc_st = 0.1;               % start time of excitation (s)

xo = 0.6;                   % coordinate of output (normalised, 0-1)
window_dur = 0.01;          % duration of fade-out window (s).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Error check for flags
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Check only 1 or 2 are selected
% Return error if other values were selected
if (itype == 1 || itype == 2) == 0
           error('itype has to be 1 or 2');
end

if (bctype == 1 || bctype == 2) == 0
           error('bctype has to be 1 or 2');
end

if (outtype == 1 || outtype == 2) == 0
           error('outtype has to be 1 or 2')
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Error check for the durartion
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if dur > Tf || dur < 0
    error('Stability condition violated');
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Error check for parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (plot_on == 1 || plot_on == 0) == 0
    error('The value has to be 0 or 1');
end

% Create empty vectors for error check
coeff = zeros(1,6); herror = zeros(1,6);
short = zeros(1,6); long = zeros(1,6);

% Error check loop
for e = 1:6
    if dur<=0 || E(e)<=0 || exc_st<=0 || famp<=0 || L(e)<=0 || r(e)<=0 || rho(e)<=0 ...
    || SR<=0 || T(e)<=0 || T60(e)<=0 || Tf<=0 || window_dur<=0 || xi<=0 || xo<=0
        error('Variables have to be positive number');
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Error check for N
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % coefficient for calculation


    coeff(e) = T(e)/(rho(e)*pi*(r(e)^2)*(SR^2));

    % herror is another calculation method of hmin
    herror(e) = sqrt(0.5*((coeff(e)) + sqrt((coeff(e)^2)+(4*E(e)*r(e)^2)/(rho(e) * SR^2))));  

    if L(e) > floor(10000*herror(e))
        error('Stability condition violated'); 
    end


    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Error check for locations
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Calculate short and far distance from 0 side
    short(e) = L(e)/(floor(L(e)/herror(e)));
    long(e) = L(e)*(1 - 1/(floor(L(e)/herror(e))));


    if xi < short(e) || xo < short(e)
            error('Stability condition violated');
    end

    if xi > long(e) || xo > long(e)
            error('Stability condition violated');
    end

end


%%%%% derived parameters

A  = pi.*r.^2;                 % string cross-sectional area
I = 0.25.*pi.*r.^4;            % string moment of intertia

c = sqrt(T./(rho.*A));        % wave speed
K = sqrt(E.*I./(rho.*A));      % stiffness constant
sig = 6*log(10)./T60;        % loss parameter

%%%%% grid


hmin = sqrt(0.5.*(c.^2.*k.^2+sqrt(c.^4.*k.^4+16.*K.^2.*k.^2)));    % minimal grid spacing

N = floor(L./hmin);          % number of segments (N+1 is number of grid points)
h = L./N;                    % adjusted grid spacing

lambda = c.*k./h;             % Courant number
mu = K.*k./h.^2;               % numerical stiffness constant


%%%%% I/O

Nf = floor(SR*Tf);          % number of time steps



% Ceiling instead of floor so that decimal part will not be ignored
li = ceil(xi.*N./L);         % grid index of excitation
lo = ceil(xo.*N./L);         % grid index of output



% create force signal
% Initialise input force signal
f = zeros(Nf,1);               
durint = floor(dur*SR);                 % duration of force signal, in samples
exc_st_int = floor(exc_st*SR);          % start time index for excitation


% Use switch to adopt itype
switch itype 
    case 1  % struck case
             for n=exc_st_int:exc_st_int+durint-1
             f(n)=famp*0.5*(1-cos(2*pi*(n-exc_st_int)/durint));
             end
    case 2  % plucked case
             for n=exc_st_int:exc_st_int+durint-1
             f(n) = famp*0.5*(1-cos(pi*(n-exc_st_int)/durint));
             end
end

%%%%% scheme coefficients

% interior

fac = 1./(1+sig.*k);
a2 = -fac.*mu.^2;
a1 = fac.*(lambda.^2+4*mu.^2);
a0 = fac.*(2-2*lambda.^2-6*mu.^2);
b0 = fac.*(1-sig.*k);

% boundary

if(bctype==1)
    % simply supported
    a00 = fac.*(2-2.*lambda.^2-5.*mu.^2);
end
if(bctype==2)
    % clamped
    a00 = fac.*(2-2.*lambda.^2-6.*mu.^2);
end


% input

d0 = fac.*(k.^2./(h.*rho.*A));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E string%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialise scheme variables
N1 = N(1); li1 = li(1); lo1 = lo(1); 

% allocate boundary conditions 
a01 = a0(1); a001 = a00(1); b01 = b0(1); 
a11 = a1(1); a21 = a2(1); d01 = d0(1);

u2 = zeros(N1+1,1);              % state
u1 = u2;                        % state
u = u2;                         % state

y1 = zeros(Nf,1);                % output

%%%%% main loop

for n=1:Nf
    % interior update
    u(3:N1-1) = a01*u1(3:N1-1)+a11*(u1(2:N1-2)+u1(4:N1))+a21*(u1(1:N1-3)+u1(5:N1+1))-b01*u2(3:N1-1); 
    % boundary updates 
    u(2) = a001*u1(2)+a11*u1(3)+a21*u1(4)-b01*u2(2);
    u(N1) = a001*u1(N1)+a11*u1(N1-1)+a21*u1(N1-2)-b01*u2(N1);
    % send in input
    u(li) = u(li)+d01*f(n);   
    % read output
    if(outtype==1)
        y1(n) = u(lo1);
    end  
    if(outtype==2)
        y1(n) = (u(lo1)-u1(lo1))/k;    %Take derivative
    end
    % plot
    if(plot_on==1)
        % draw current state
        figure(1)
        plot([0:N1]'*h, u, 'k');
        axis([0 L(1) -0.005 0.005])
        drawnow
    end
    % shift state
    u2 = u1;
    u1 = u;
end

%%% Windowing

% Size of window
windowSize = 2*window_dur*SR;
% Make a hann window
hanWindow = 0.5 - 0.5*cos(2*pi*(0:windowSize-1)/windowSize);

% Trancate half
han = hanWindow(windowSize/2:windowSize);

% location of output signal to be windowed
st = (Tf-window_dur)*SR; ed = Tf*SR;

% Apply the half window to output signal
ywindowed1 = y1(st:ed).*han';

% combine with original y
out1 = [y1(1:st-1);ywindowed1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% A string%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialise scheme variables
N2 = N(2); li2 = li(2); lo2 = lo(2); 

% allocate boundary conditions 
a02 = a0(2); a002 = a00(2); b02 = b0(2); 
a12 = a1(2); a22 = a2(2); d02 = d0(2);

u2 = zeros(N2+1,1);              % state
u1 = u2;                        % state
u = u2;                         % state

y2 = zeros(Nf,1);                % output

%%%%% main loop

for n=1:Nf
    % interior update
    u(3:N2-1) = a02*u1(3:N2-1)+a12*(u1(2:N2-2)+u1(4:N2))+a22*(u1(1:N2-3)+u1(5:N2+1))-b02*u2(3:N2-1); 
    % boundary updates 
    u(2) = a002*u1(2)+a12*u1(3)+a22*u1(4)-b02*u2(2);
    u(N2) = a002*u1(N2)+a12*u1(N2-1)+a22*u1(N2-2)-b02*u2(N2);
    % send in input
    u(li2) = u(li2)+d02*f(n);   
    % read output
    if(outtype==1)
        y2(n) = u(lo2);
    end  
    if(outtype==2)
        y2(n) = (u(lo2)-u1(lo2))/k;    %Take derivative
    end
    % plot
    if(plot_on==1)
        % draw current state
        figure(1)
        plot([0:N2]'*h, u, 'k');
        axis([0 L(2) -0.005 0.005])
        drawnow
    end
    % shift state
    u2 = u1;
    u1 = u;
end

%%% Windowing

% Size of window
windowSize = 2*window_dur*SR;
% Make a hann window
hanWindow = 0.5 - 0.5*cos(2*pi*(0:windowSize-1)/windowSize);

% Trancate half
han = hanWindow(windowSize/2:windowSize);

% location of output signal to be windowed
st = (Tf-window_dur)*SR; ed = Tf*SR;

% Apply the half window to output signal
ywindowed2 = y2(st:ed).*han';

% combine with original y
out2 = [y2(1:st-1);ywindowed2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D string%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialise scheme variables
N3 = N(3); li3 = li(3); lo3 = lo(3); 

% allocate boundary conditions 
a03 = a0(3); a003 = a00(3); b03 = b0(3); 
a13 = a1(3); a23 = a2(3); d03 = d0(3);

u2 = zeros(N3+1,1);              % state
u1 = u2;                        % state
u = u2;                         % state

y3 = zeros(Nf,1);                % output

%%%%% main loop

for n=1:Nf
    % interior update
    u(3:N3-1) = a03*u1(3:N3-1)+a13*(u1(2:N3-2)+u1(4:N3))+a23*(u1(1:N3-3)+u1(5:N3+1))-b03*u2(3:N3-1); 
    % boundary updates 
    u(2) = a003*u1(2)+a13*u1(3)+a23*u1(4)-b03*u2(2);
    u(N3) = a003*u1(N3)+a13*u1(N3-1)+a23*u1(N3-2)-b03*u2(N3);
    % send in input
    u(li3) = u(li3)+d03*f(n);   
    % read output
    if(outtype==1)
        y3(n) = u(lo3);
    end  
    if(outtype==2)
        y3(n) = (u(lo3)-u1(lo3))/k;    %Take derivative
    end
    % plot
    if(plot_on==1)
        % draw current state
        figure(1)
        plot([0:N3]'*h, u, 'k');
        axis([0 L(3) -0.005 0.005])
        drawnow
    end
    % shift state
    u2 = u1;
    u1 = u;
end

%%% Windowing

% Size of window
windowSize = 2*window_dur*SR;
% Make a hann window
hanWindow = 0.5 - 0.5*cos(2*pi*(0:windowSize-1)/windowSize);

% Trancate half
han = hanWindow(windowSize/2:windowSize);

% location of output signal to be windowed
st = (Tf-window_dur)*SR; ed = Tf*SR;

% Apply the half window to output signal
ywindowed3 = y3(st:ed).*han';

% combine with original y
out3 = [y3(1:st-1);ywindowed3];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G string%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialise scheme variables
N4 = N(4); li4 = li(4); lo4 = lo(4); 

% allocate boundary conditions 
a04 = a0(4); a004 = a00(4); b04 = b0(4); 
a14 = a1(4); a24 = a2(4); d04 = d0(4);

u2 = zeros(N4+1,1);              % state
u1 = u2;                        % state
u = u2;                         % state

y4 = zeros(Nf,1);                % output

%%%%% main loop

for n=1:Nf
    % interior update
    u(3:N4-1) = a04*u1(3:N4-1)+a14*(u1(2:N4-2)+u1(4:N4))+a24*(u1(1:N4-3)+u1(5:N4+1))-b04*u2(3:N4-1); 
    % boundary updates 
    u(2) = a004*u1(2)+a14*u1(3)+a24*u1(4)-b04*u2(2);
    u(N4) = a004*u1(N4)+a14*u1(N4-1)+a24*u1(N4-2)-b04*u2(N4);
    % send in input
    u(li4) = u(li4)+d04*f(n);   
    % read output
    if(outtype==1)
        y4(n) = u(lo4);
    end  
    if(outtype==2)
        y4(n) = (u(lo4)-u1(lo4))/k;    %Take derivative
    end
    % plot
    if(plot_on==1)
        % draw current state
        figure(1)
        plot([0:N4]'*h, u, 'k');
        axis([0 L(4) -0.005 0.005])
        drawnow
    end
    % shift state
    u2 = u1;
    u1 = u;
end

%%% Windowing

% Size of window
windowSize = 2*window_dur*SR;
% Make a hann window
hanWindow = 0.5 - 0.5*cos(2*pi*(0:windowSize-1)/windowSize);

% Trancate half
han = hanWindow(windowSize/2:windowSize);

% location of output signal to be windowed
st = (Tf-window_dur)*SR; ed = Tf*SR;

% Apply the half window to output signal
ywindowed4 = y4(st:ed).*han';

% combine with original y
out4 = [y4(1:st-1);ywindowed4];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% B string%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialise scheme variables
N5 = N(5); li5 = li(5); lo5 = lo(5); 

% allocate boundary conditions 
a05 = a0(5); a005 = a00(5); b05 = b0(5); 
a15 = a1(5); a25 = a2(5); d05 = d0(5);

u2 = zeros(N5+1,1);              % state
u1 = u2;                        % state
u = u2;                         % state

y5 = zeros(Nf,1);                % output

%%%%% main loop

for n=1:Nf
    % interior update
    u(3:N5-1) = a05*u1(3:N5-1)+a15*(u1(2:N5-2)+u1(4:N5))+a25*(u1(1:N5-3)+u1(5:N5+1))-b05*u2(3:N5-1); 
    % boundary updates 
    u(2) = a005*u1(2)+a15*u1(3)+a25*u1(4)-b05*u2(2);
    u(N5) = a005*u1(N5)+a15*u1(N5-1)+a25*u1(N5-2)-b05*u2(N5);
    % send in input
    u(li5) = u(li5)+d05*f(n);   
    % read output
    if(outtype==1)
        y5(n) = u(lo5);
    end  
    if(outtype==2)
        y5(n) = (u(lo5)-u1(lo5))/k;    %Take derivative
    end
    % plot
    if(plot_on==1)
        % draw current state
        figure(1)
        plot([0:N5]'*h, u, 'k');
        axis([0 L(5) -0.005 0.005])
        drawnow
    end
    % shift state
    u2 = u1;
    u1 = u;
end

%%% Windowing

% Size of window
windowSize = 2*window_dur*SR;
% Make a hann window
hanWindow = 0.5 - 0.5*cos(2*pi*(0:windowSize-1)/windowSize);

% Trancate half
han = hanWindow(windowSize/2:windowSize);

% location of output signal to be windowed
st = (Tf-window_dur)*SR; ed = Tf*SR;

% Apply the half window to output signal
ywindowed5 = y5(st:ed).*han';

% combine with original y
out5 = [y5(1:st-1);ywindowed5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% High E string%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialise scheme variables
N6 = N(6); li6 = li(6); lo6 = lo(6); 

% allocate boundary conditions 
a06 = a0(6); a006 = a00(6); b06 = b0(6); 
a16 = a1(6); a26 = a2(6); d06 = d0(6);

u2 = zeros(N6+1,1);              % state
u1 = u2;                        % state
u = u2;                         % state

y6 = zeros(Nf,1);                % output

%%%%% main loop

for n=1:Nf
    % interior update
    u(3:N6-1) = a06*u1(3:N6-1)+a16*(u1(2:N6-2)+u1(4:N6))+a26*(u1(1:N6-3)+u1(5:N6+1))-b06*u2(3:N6-1); 
    % boundary updates 
    u(2) = a006*u1(2)+a16*u1(3)+a26*u1(4)-b06*u2(2);
    u(N6) = a006*u1(N6)+a16*u1(N6-1)+a26*u1(N6-2)-b06*u2(N6);
    % send in input
    u(li6) = u(li6)+d06*f(n);   
    % read output
    if(outtype==1)
        y6(n) = u(lo6);
    end  
    if(outtype==2)
        y6(n) = (u(lo6)-u1(lo6))/k;    %Take derivative
    end
    % plot
    if(plot_on==1)
        % draw current state
        figure(1)
        plot([0:N6]'*h, u, 'k');
        axis([0 L(6) -0.005 0.005])
        drawnow
    end
    % shift state
    u2 = u1;
    u1 = u;
end

%%% Windowing

% Size of window
windowSize = 2*window_dur*SR;
% Make a hann window
hanWindow = 0.5 - 0.5*cos(2*pi*(0:windowSize-1)/windowSize);

% Trancate half
han = hanWindow(windowSize/2:windowSize);

% location of output signal to be windowed
st = (Tf-window_dur)*SR; ed = Tf*SR;

% Apply the half window to output signal
ywindowed6 = y6(st:ed).*han';

% combine with original y
out6 = [y6(1:st-1);ywindowed6];



% Normalise output signals
Max1 = max(abs(out1)); Max2 = max(abs(out2)); Max3 = max(abs(out3));
Max4 = max(abs(out4)); Max5 = max(abs(out5)); Max6 = max(abs(out6));

out1 = out1/Max1; out2 = out2/Max2; out3 = out3/Max3;
out4 = out4/Max4; out5 = out5/Max5; out6 = out6/Max6;

%%%%% play sound

out=out1+out2+out3+out4+out5+out6;

% Play sound

soundsc(out,SR)