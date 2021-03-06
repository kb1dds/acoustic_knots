% Testing noise and jitter robustness for sonar targets with torus knot signatures
%
% This is to be run as a script.
% Many figures and PNG files are produced as output,
%  titles and other annotations are made automatically,
%   and hopefully are reasonable!
% To change what's plotted, comment out what you don't want!
%
% Copyright (c) 2020 Michael Robinson
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
tx_jitters = [0.01]; % (meters)

% Scatterer relative rotation angle (deg)
offset_angle=28;

% Number of measurements to collect
nlooks = 200;

% Sliding window size
window_size = 20;

% Number of frequencies to collect
nfrequencies = 61;

% Pulse bandwidth
bandwidth=600;

% Frequencies to examine (Hertz)
frequencies=linspace(0,bandwidth,nfrequencies+1);
frequencies(1)=[];

fp=fopen('roc_results.csv','wt');
fprintf(fp,'bandwidth,SNR,nlooks,txrange,rxrange,jitter,eccentricity,diameter,diameter_ratio,rcs_ratio,offset,true_m,true_n,knotify_m,knotify_n,alex_poly_computed\n');

% Scatterer configuration
% Scatterers are in 2d as collections of points
% There are three such "composite" scatterers:
%  two lie at the vertices of regular polygons (with m and n sides)
%  the third is the sum of the other two
% A jitter parameter allows points to deviate from their true positions by a specified amount
m = 2; % Winding number
radius_m=1; % (meters)
n = 3; % Winding number
radius_n=1; % (meters)
jitters=[0.01]; % (meters, variance)

% Receiver noise level
noise_levels=logspace(-5,5,20); % Variance

% Number of trials for statistical test
ntrials=10;

%% You should not need to change code below this line if you are just experimenting with parameters!

% Scatterer locations
theta_m=(0:2*pi/m:2*pi).';
theta_m=theta_m(1:m);
theta_n=offset_angle*pi/180+(0:2*pi/n:2*pi).';
theta_n=theta_n(1:n);
scatloc_m=[radius_m*cos(theta_m),radius_m*sin(theta_m)];
scatloc_n=[radius_n*cos(theta_n),radius_n*sin(theta_n)];
scatloc_clean=[scatloc_m; scatloc_n];

% Construct scatterer cross sections
scatcross=ones(size(scatloc_clean,1),1);

% Construct transmitter and receiver locations
theta_platform=linspace(0,2*pi,nlooks).';
txloc_clean=[txrange*cos(theta_platform),txrange*sin(theta_platform)];
rxloc_clean=[rxrange*cos(theta_platform),rxrange*sin(theta_platform)];

correct_classifications=zeros(numel(jitters),numel(tx_jitters),numel(noise_levels),ntrials);
			    
for trial=1:1%ntrials,
  for idx_jitter=1:1%numel(jitters),
    jitter=jitters(idx_jitter);

    % Add scatterer jitter
    scatloc=scatloc_clean+jitter*randn(size(scatloc_clean));

    for idx_tx_jitter=1:numel(tx_jitters),
      tx_jitter=tx_jitters(idx_tx_jitter);
      
      % Add platform jitter
      tx_jitter_sig=tx_jitter*randn(size(txloc_clean));
      rx_jitter_sig=tx_jitter*randn(size(rxloc_clean));
  
      txloc=txloc_clean+tx_jitter_sig;
      rxloc=rxloc_clean+rx_jitter_sig;

      echos_clean=zeros(nlooks,length(frequencies));
      for i=1:length(frequencies)
        % Propagate from transmitter locations to scatterers
        tx2scat=isotropicMatrix(c/frequencies(i),txloc,scatloc);
 
        % Propagate from scatterers to transmitter
        scat2rx=isotropicMatrix(c/frequencies(i),scatloc,rxloc);

        % Store noiseless frequency response
        echos_clean(:,i)=diag(scat2rx*bsxfun(@times,scatcross,tx2scat));
      end

      for idx_noise_level=1:1%numel(noise_levels),
	noise_level=noise_levels(idx_noise_level);
	
        % Add receiver noise
        echos=echos_clean+(randn(size(echos_clean))+sqrt(-1)*randn(size(echos_clean)))*noise_level/sqrt(2);

	% SNR
	rx_snr=10*log10(sum(abs(echos_clean(:)).^2)/noise_level);

	% Update the plots
	subplot(1,3,1);
	imagesc(abs(echos))

	subplot(1,3,2);
	imagesc(abs(fft(echos,[],1)))

	% Apply sliding window aggregation
	echos_ac=echos;
	for ws=1:window_size,
	  echos=[echos, circshift(echos_ac,[ws,0])];
	end

	% Perform torus knot transform
	knf=knotify(echos,1:2*max([m,n]),1:2*max([m,n]));

	subplot(1,3,3);
	imagesc(knf);
	%pause(0.05)

	% Identify local minimum
	[jnk,idx]=min(knf(:));
	[r,c]=ind2sub(size(knf),idx);
	if( r > c ) % Standardize winding numbers on being ordered
	  r2=r;
	  r=c;
	  c=r2;
	end
	
        % Compute Alexander polynomial for the knot via external call
	csvwrite('signal.csv',[real(echos) imag(echos)]);
	system('python3 knot_invariants.py signal.csv temp.txt');
	fp2=fopen('temp.txt','rt');
	alex_poly=fgetl(fp2);
	fclose(fp2);

	fprintf(fp,'%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%d,%d,%d,%d,%s\n',bandwidth,rx_snr,nlooks,txrange,rxrange,tx_jitter,1,max(radius_m,radius_n),max(radius_m,radius_n)/min(radius_m,radius_n),0,offset_angle,m,n,r,c,alex_poly);

	% Is the local minimum the correct knot type?
	correct=((r==m)&&(c==n))||((r==n)&&(c==m));
	if correct,
	  correct_classifications(idx_jitter,idx_tx_jitter,idx_noise_level,trial)=correct_classifications(idx_jitter,idx_tx_jitter,idx_noise_level,trial)+1;
	end
      end
    end
  end
end

fclose(fp);

% Build basic ROC plot
figure
roc=squeeze(sum(sum(sum(correct_classifications,1),2),4))/(numel(jitters)*numel(tx_jitters)*ntrials);
semilogx(noise_levels,roc);
ylim([0,1]);
xlabel('Noise level');
ylabel('Correct detection rate');
