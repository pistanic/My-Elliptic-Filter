function [H] = myElliptic( k, Ap, Aa)
clc

% MYELLIPTIC Calulates the Transfer function and loss charictaristics of 
% a digital approximation of the Elliptic filter using the method described 
% in: Digital Signal Processing - Signals, Systems and Filters 
% by: Andreas Antoniou 
% 
% Normalized transfer function: 
% Hn(s) = H0/D0(s)*PI[i=1 -> r]{(s^2 + a0[i])/(s^2 + b1[i]*s + b0[i0])}
%
% Where r = |(n-1)/2 for odd n
%           |(n/2)   for even n
% and 
%   D0(s) = |s + sig0 for odd n   
%           |1        for even n 
%
%
%
% Plots Produced by this function: 
% - Normalized Transferfunction plot (freqz)
% - Loss charictaristics assoicated with obtined transferfunction. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               USAGE
%
% Selectivity factor - k is defined as the passband edge over stopband 
% edge. K = Wp/Wa
% Maximum pass band attenuation : Ap 
% Minimum stop band attenuation : Aa 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k1 = sqrt(1-k^2); % (Eq: 10.86)

q0 = 0.5* (1-sqrt(k1) ) / (1+sqrt(k1)); % (Eq: 10.87)

a = [1 5 9 13];

b = [1 2 15 150];

q = (q0.^a)*b'; % (Eq: 10.88)

D = ( 10^(0.1*Aa) -1 ) / ( 10^(0.1*Ap) -1 ); % (Eq: 10.89)

nf = log10(16*D)/(log10(1/q)); % Calulates exact number for n 

% following line rounds exact number for n, nf to integer
n = round((nf+0.5)); % (Eq: 10.90)

AA = 10*log10( ( (10^(0.1*Ap)-1) / (16*q^n) ) + 1 );

lam = 1/2/n*log( (10^(0.05*Ap) +1 ) / ( 10^(0.05*Aa) -1 ) ); %(Eq: 10.91) 
 
m = 10;

sigN = (-1).^(0:m).*q.^( (0:m).*((0:m)+1) ).*sinh( (2*(0:m)+1)*lam );

sigD = (-1).^(1:m).*q.^((1:m).^2).*cosh(2*(1:m)*lam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INCORRECT CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sig0 = abs( (2*(q^(1/4))*sum(sigN)) / (1+ 2*sum(sigD)) );  %(Eq: 10.92)
sig02 = sig0^2; 
W = sqrt( ( 1+ (k*sig02) ) * ( 1+ (sig02/k) ) ); %(Eq: 10.93)

if(mod(n,2) == 0)
    r = n/2;
    u = (1:r) - 0.5;
else
    r = (n-1)/2;
    u = 1:r;
end
OmegaN = @(u) (-1).^(0:m).*q.^((0:m).*( (0:m)+1 )).*sin((2*(0:m)+1)*pi*u/n);
OmegaD = @(u) (-1).^(1:m).*q.^((1:m).^2).*cos(2*(1:m)*pi*u/n);
Omega = zeros(size(u));
for i = 1:length(Omega)
    Omega(i) = ( 2*q^(1/4) )*sum( OmegaN(u(i) ) )/( 1+ 2*sum(OmegaD(u(i)) ) );
end

Omega2 = power(Omega,2);

% This Loop Calulates the Coeffictions for the transferfunction
for t = 1:r
    V(t) = sqrt( (1- (k*Omega2(t)) )*( 1- (Omega2(t)/k) ) ); % (Eq: 10.95)
    a0(t) = 1/Omega2(t); % (Eq: 10.96) 
end
b0 = ( (sig0*V).^2 + (Omega*W).^2) ./ (1+(sig0*Omega).^2).^2;
b1 = 2*sig0*V ./ (1+ (sig0*Omega).^2);

if(mod(n,2) == 0)
    H0 = 10^(-0.05*Ap)*prod(b0./a0);
    D0 = 1; 
else 
    H0 = sig0 *prod(b0 ./ a0); % (Eq: 10.99)
    D0 = [1 sig0]; 
end 

Num = [1 0 a0(1)];
Den = [1 b1(1) b0(1)];

% Calulate Transfer Function 
if(length(a0) > 1)  
    for t = 2:length(a0 - 1)
        Num = conv(Num, [1, 0, a0(t)]);
        Den = conv(Den, [1, b1(t), b0(t)]);
    end 

Num = H0*Num;
Den = conv(Den,D0);
HnS = tf(Num, Den);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DEBUGGING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%format long 
%q0
%q
%k1 
%D
% sig0 = 0.623480181667140
%Omega
%Omega2
%deltaN;
%deltaD;
%a0;
%b0;
%b1;
%H0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

option = bodeoptions;
option.Grid = 'on';
option.PhaseVisible = 'off'; 
option.FreqUnits = 'Hz';
figure(1) 
axL = gca; 
set(gcf,'Position',[154 129 866 537]);
bodeplot(1/HnS,option);
axL.YLim = [-10 100];
title(['Elliptic Loss Charictaristics']);
% Transfer Function in S
HnS 
% Actual Stopband Loss
fprintf('Actual Stopband loss = ');
AA  
fprintf('* * * * * * * * * * Additional OutPuts * * * * * * * * * *')
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Z Domain Analysis 

%%%%%%%%%%%%%%%%%%%%%%%%%% Bilinear Transform %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[numZ,denZ] = bilinear(Num,Den, 0.5);
fprintf('Z-Domain Transfer Function');
HnZ = tf(numZ,denZ,1024);
HnZ % Transfer Function in Z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%         - Z-Domain Polt of myElliptinc plot        %%%%%%%%%%%%%%%%%%%
%%%%%%%         - Z-Domain Comparison Plot vs ellip function %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
figure(2);  
N = 1024; % Number of points to evaluate at
upp = pi; % Evaluate only up to fs/2
% Create the vector of angular frequencies at one more point.
% After that remove the last element (Nyquist frequency)
w = linspace(0, pi, N+1); 
w(end) = [];
ze = exp(-1j*w); % Pre-compute exponent
H = polyval(numZ, ze)./polyval(denZ, ze); % Evaluate transfer function and take the amplitude
Ha = abs(H);
Hdb_  = 20*log10(Ha) % Convert to dB scale
wn_   = w/pi;
% Plot and set axis limits
plot(wn_, Hdb_);
xlim([0 1]);
ylim([-120 20]);
grid on

title('My Z Domain Elliptic');
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab Ellip 

[mB,mA] = ellip(n, Ap, Aa, 0.5);
H = polyval(mB, ze)./polyval(mA, ze);
Ha = abs(H); 
Hdb  = 20*log10(Ha); % Convert to dB scale
wn   = w/pi;
% Comparison of ellip to myElliptic 
figure(3)
plot(wn, Hdb); % ellip
grid on
hold on
plot(wn, Hdb_, '--r') % myElliptic
title('Z-Domain Elliptic Filter Comparison')
legend({'ellip','myElliptic'})


            %}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            End Of Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
