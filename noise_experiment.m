% Basic CSAS point-scatterer target simulator with variable noise level
%
% This is to be run as a script.
% Many figures and PNG files are produced as output,
%  titles and other annotations are made automatically,
%   and hopefully are reasonable!
% To change what's plotted, comment out what you don't want!
%
% Copyright (c) 2020,2023 Michael Robinson, Maxwell Gualtieri
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
   
% Wave phase speed (meters/second)
c = 343; % Default speed of sound in dry air

% Ambient dimension
dim = 3;

% Standoff ranges (meters)
txrange = 3;
rxrange = 3;
tx_jitter = 0.01; % (meters)
rx_jitter = 0.01; % (meters)
lock_platform_jitter = 0; % Nonzero if rx and tx jitter components are identical; rx_jitter is ignored in that case

% Scatterer center offset
target_center = [0.01,0.005]; % (meters)

% Number of measurements to collect
nlooks = 360;

% Number of frequencies to collect
nfrequencies = 200;

% Frequencies to examine (Hertz)
frequencies=linspace(0,600,nfrequencies+1);
frequencies(1)=[];

% A jitter parameter allows points to deviate from their true positions by a specified amount
jitter=0.01; % (meters, variance)

% Receiver noise level
noise_levels=linspace(0.001,0.01,20); % Variance

%% You should not need to change code below this line if you are just experimenting with parameters!

target_name={'coke_bottle','cup_open','cup_capped','pipe_open','pipe_capped'};

for noise_idx = 1:length(noise_levels),
  noise_level = noise_levels(noise_idx);
  for target_idx = 4,
    

% Data pulled from CSV file
data = dlmread([target_name{target_idx} '_scatterers.csv'],',',1,0);
scatloc=bsxfun(@plus,data(:,1:2),target_center);
scatcross=data(:,3).*exp(sqrt(-1)*data(:,4));

% Construct transmitter and receiver locations
theta_platform=linspace(0,2*pi,nlooks).';
txloc=[txrange*cos(theta_platform),txrange*sin(theta_platform)];
rxloc=[rxrange*cos(theta_platform),rxrange*sin(theta_platform)];

% Add platform jitter
tx_jitter_sig=tx_jitter*randn(size(txloc));
if( lock_platform_jitter )
  rx_jitter_sig=tx_jitter_sig;
else
  rx_jitter_sig=rx_jitter*randn(size(rxloc));
end
txloc=txloc+tx_jitter_sig;
rxloc=rxloc+rx_jitter_sig;

echos=zeros(nlooks,length(frequencies));
for i=1:length(frequencies)
  % Propagate from transmitter locations to scatterers
  tx2scat=isotropicMatrix(c/frequencies(i),txloc,scatloc);
 
  % Propagate from scatterers to transmitter
  scat2rx=isotropicMatrix(c/frequencies(i),scatloc,rxloc);

  % Store noiseless frequency response
  echos(:,i)=diag(scat2rx*bsxfun(@times,scatcross,tx2scat));
end

% Add noise
echos=echos+randn(size(echos))*noise_level;

% Save mat file
save([target_name{target_idx} "_" num2str(noise_level) ".mat"], "echos", "rxloc", "txloc", "scatloc", "noise_level")

% Save CSV file
csvwrite([target_name{target_idx} "_echos_" num2str(noise_level)  ".csv"],echos);
end
end

