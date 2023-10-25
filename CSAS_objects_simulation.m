% Torus knots from acoustic point scatterer configurations
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
tx_jitter = 0.01; % (meters)
rx_jitter = 0.01; % (meters)
lock_platform_jitter = 0; % Nonzero if rx and tx jitter components are identical; rx_jitter is ignored in that case

% Scatterer relative rotation angle (deg)
offset_angle=28;

% Number of measurements to collect
nlooks = 200;

% Number of frequencies to collect
nfrequencies = 200;

% Frequencies to examine (Hertz)
frequencies=linspace(0,600,nfrequencies+1);
frequencies(1)=[];

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
jitter=0.01; % (meters, variance)

% Receiver noise level
noise_level=0.0001; % Variance

%% You should not need to change code below this line if you are just experimenting with parameters!

% Scatterer locations

% Coke Bottle Geometry

c = 1;

scatloc_m = [];
for i = -1:0.01:0
    scatloc_m(c,1) = i;
    scatloc_m(c,2) = -.25+.05*(cos(8*(i+1)));
    scatloc_m(c+1,1) = i;
    scatloc_m(c+1,2) = .25-.05*(cos(8*(i+1)));
    c = c+2; 
end

e = c+1;
for i = 0:0.01:.4
    scatloc_m(e,1) = i;
    scatloc_m(e,2) = .25+.02*(sin(8*(i+1)));
    scatloc_m(e+1,1) = i;
    scatloc_m(e+1,2) = -.25-.02*(sin(8*(i+1)));
    e = e+2; 
end

g = e+1;
for i = .4:0.01:.75
    scatloc_m(g,1) = i;
    scatloc_m(g,2) = .2-.2*(sqrt(i-.4));
    scatloc_m(g+1,1) = i;
    scatloc_m(g+1,2) = -.2+.2*(sqrt(i-.4));
    g = g+2; 
end

scatloc_m=scatloc_m+randn(size(scatloc_m))*jitter;

d = 1;
scatloc_n = [];
for j = -0.2:0.01:0.2
    scatloc_n(d,1) = -1;
    scatloc_n(d,2) = j;
    d = d+1; 
end





%Cup without a lid geometry

% c = 1;
% scatloc_m = [];
% for i = -.5:0.01:.5
%     scatloc_m(c,1) = i;
%     scatloc_m(c,2) = 0.2+(i+.5)*.15;
%     scatloc_m(c+1,1) = i;
%     scatloc_m(c+1,2) = -0.2-(i+.5)*.15;
%     c = c+2; 
% end
% scatloc_m=scatloc_m+randn(size(scatloc_m))*jitter;
% % scatloc_n=[radius_n*cos(theta_n),radius_n*sin(theta_n)];
% d = 1;
% scatloc_n = [];
% for j = -0.2:0.01:0.2
%     scatloc_n(d,1) = -0.5;
%     scatloc_n(d,2) = j;
%     d = d+1; 
% end

% Cup with a lid geometry

% c = 1;
% scatloc_m = [];
% for i = -.5:0.01:.5
%     scatloc_m(c,1) = i;
%     scatloc_m(c,2) = 0.2+(i+.5)*.15;
%     scatloc_m(c+1,1) = i;
%     scatloc_m(c+1,2) = -0.2-(i+.5)*.15;
%     c = c+2; 
% end
% scatloc_m=scatloc_m+randn(size(scatloc_m))*jitter;
% % scatloc_n=[radius_n*cos(theta_n),radius_n*sin(theta_n)];
% d = 1;
% scatloc_n = [];
% for j = -0.2:0.01:0.2
%     scatloc_n(d,1) = -0.5;
%     scatloc_n(d,2) = j;
%     d = d+1; 
% end
% 
% 
% e = d+1;
% for j = -0.35:0.01:0.35
%     scatloc_n(e,1) = 0.5;
%     scatloc_n(e,2) = j;
%     e = e+1; 
% end

% Pipe with no lid geometry

% c = 1;
% scatloc_m = [];
% for i = -.5:0.01:.5
%     scatloc_m(c,1) = i;
%     scatloc_m(c,2) = 0.05;
%     scatloc_m(c+1,1) = i;
%     scatloc_m(c+1,2) = -0.05;
%     c = c+2; 
% end
% scatloc_m=scatloc_m+randn(size(scatloc_m))*jitter;
% 
% d = 1;
% scatloc_n = [];
% for j = -0.05:0.01:0.05
%     scatloc_n(d,1) = 0;
%     scatloc_n(d,2) = 0;
%     scatloc_n(d+1,1) = 0;
%     scatloc_n(d+1,2) = 0;
%     d = d+2; 
% end

% Pipe with lid geometry

% c = 1;
% scatloc_m = [];
% for i = -.5:0.01:.5
%     scatloc_m(c,1) = i;
%     scatloc_m(c,2) = 0.05;
%     scatloc_m(c+1,1) = i;
%     scatloc_m(c+1,2) = -0.05;
%     c = c+2; 
% end
% scatloc_m=scatloc_m+randn(size(scatloc_m))*jitter;
% d = 1;
% scatloc_n = [];
% for j = -0.05:0.01:0.05
%     scatloc_n(d,1) = 0.5;
%     scatloc_n(d,2) = j;
%     scatloc_n(d+1,1) = -0.5;
%     scatloc_n(d+1,2) = j;
%     d = d+2; 
% end


scatloc_n=scatloc_n+randn(size(scatloc_n))*jitter;
scatloc=[scatloc_m; scatloc_n];

% Construct scatterer cross sections
scatcross_m=ones(size(scatloc_m,1),1);
scatcross_n=ones(size(scatloc_n,1),1);
scatcross=ones(size(scatloc,1),1);

% Construct transmitter and receiver locations
theta_platform=linspace(0,2*pi,nlooks).';
txloc=[txrange*cos(theta_platform),txrange*sin(theta_platform)];
rxloc=[rxrange*cos(theta_platform),rxrange*sin(theta_platform)];
txloc_m=[txrange*cos(theta_platform/m),txrange*sin(theta_platform/m)];
rxloc_m=[rxrange*cos(theta_platform/m),rxrange*sin(theta_platform/m)];
txloc_n=[txrange*cos(theta_platform/n),txrange*sin(theta_platform/n)];
rxloc_n=[rxrange*cos(theta_platform/n),rxrange*sin(theta_platform/n)];

% Add platform jitter
tx_jitter_sig=tx_jitter*randn(size(txloc));
if( lock_platform_jitter )
  rx_jitter_sig=tx_jitter_sig;
else
  rx_jitter_sig=rx_jitter*randn(size(rxloc));
end
txloc=txloc+tx_jitter_sig;
txloc_m=txloc_m+tx_jitter_sig;
txloc_n=txloc_n+tx_jitter_sig;
rxloc=rxloc+rx_jitter_sig;
rxloc_m=rxloc_m+rx_jitter_sig;
rxloc_n=rxloc_n+rx_jitter_sig;

% Construct toroidal coordinates for each transmitter and receiver location
txtoroidal=mod([theta_platform*n, theta_platform*m],2*pi);
rxtoroidal=mod([theta_platform*n, theta_platform*m],2*pi);

echos_m=zeros(nlooks,length(frequencies));
echos_n=zeros(nlooks,length(frequencies));
echos=zeros(nlooks,length(frequencies));
echos_m_torus=zeros(length(theta_platform),length(frequencies));
echos_n_torus=zeros(length(theta_platform),length(frequencies));
for i=1:length(frequencies)
  %% First scatterer: regular m-polygon with scatterers at its corners
  
  % Propagate from transmitter locations to scatterers
  tx2scat_m=isotropicMatrix(c/frequencies(i),txloc,scatloc_m);
 
  % Propagate from scatterers to transmitter
  scat2rx_m=isotropicMatrix(c/frequencies(i),scatloc_m,rxloc);

  % Store noiseless frequency response
  echos_m(:,i)=diag(scat2rx_m*bsxfun(@times,scatcross_m,tx2scat_m));

  %% Second scatterer: regular n-polygon with scatterers at its corners
  
  % Propagate from transmitter locations to scatterers
  tx2scat_n=isotropicMatrix(c/frequencies(i),txloc,scatloc_n);
 
  % Propagate from scatterers to transmitter
  scat2rx_n=isotropicMatrix(c/frequencies(i),scatloc_n,rxloc);

  % Store noiseless frequency response
  echos_n(:,i)=diag(scat2rx_n*bsxfun(@times,scatcross_n,tx2scat_n));

  %% Third scatterer: sum of the above.
  % (Note: we could have reused the code above, but didn't for clarity)
  
  % Propagate from transmitter locations to scatterers
  tx2scat=isotropicMatrix(c/frequencies(i),txloc,scatloc);
 
  % Propagate from scatterers to transmitter
  scat2rx=isotropicMatrix(c/frequencies(i),scatloc,rxloc);

  % Store noiseless frequency response
  echos(:,i)=diag(scat2rx*bsxfun(@times,scatcross,tx2scat));
  
  %% Compute overall data (for the space of targets) from which our data above is a slice
  
  % Propagate from transmitter locations to scatterers
  tx2scat_m_torus=isotropicMatrix(c/frequencies(i),txloc_m,scatloc_m);
 
  % Propagate from scatterers to transmitter
  scat2rx_m_torus=isotropicMatrix(c/frequencies(i),scatloc_m,rxloc_m);

  % Store noiseless frequency response
  echos_m_torus(:,i)=diag(scat2rx_m_torus*bsxfun(@times,scatcross_m,tx2scat_m_torus));
  
  % Propagate from transmitter locations to scatterers
  tx2scat_n_torus=isotropicMatrix(c/frequencies(i),txloc_n,scatloc_n);
 
  % Propagate from scatterers to transmitter
  scat2rx_n_torus=isotropicMatrix(c/frequencies(i),scatloc_n,rxloc_n);

  % Store noiseless frequency response
  echos_n_torus(:,i)=diag(scat2rx_n_torus*bsxfun(@times,scatcross_n,tx2scat_n_torus));
end

% Compute the full torus response
echos_torus=bsxfun(@plus,reshape(echos_n_torus,1,[],length(frequencies)),reshape(echos_m_torus,[],1,length(frequencies)));

% Add noise
echos_m=echos_m+randn(size(echos_m))*noise_level;
echos_n=echos_n+randn(size(echos_n))*noise_level;
echos=echos+randn(size(echos))*noise_level;
echos_torus=echos_torus+randn(size(echos_torus))*noise_level;

%% Plotting code below this line!

% Display the collection geometry
figure
plot(scatloc(:,1),scatloc(:,2),'b+');
hold on
plot(txloc(:,1),txloc(:,2),'r');
plot(rxloc(:,1),rxloc(:,2),'g');
plot(txloc_m(:,1),txloc_m(:,2),'k*');
plot(txloc_n(:,1),txloc_n(:,2),'r.');
axis equal
print('-dpng',['collection_geometry_m=' num2str(m) '_n=' num2str(n) '.png']);

% Display the echo data
figure
colormap gray
imagesc(frequencies,theta_platform*180/pi,abs(echos));
set(gca,'ydir','normal');
axis normal
caxis([0,0.02])
colorbar
xlabel('Frequency (Hz)')
ylabel('Look angle (deg)')
% title(['Target with ' num2str(m) '-fold symmetry'])
print('-dpng',['toroidal_flat_m=' num2str(m) '.png']);

figure
colormap gray
imagesc(abs(ifft(echos,[],2)))
set(gca,'ydir','normal');
axis normal
caxis([0,0.02])
colorbar
xlabel('Time (samples)')
ylabel('Look angle (deg)');
% title(['Target with ' num2str(n) '-fold symmetry'])
print('-dpng',['toroidal_ranges_flat_n=' num2str(n) '.png']);

figure
rangeVector=bsxfun(@minus,permute(echos,[1 3 2]),permute(echos,[3 1 2]));
dst=sqrt(sum(abs(rangeVector).^2,3));
coords=pcoa(dst,3);
plot3(coords(:,1),coords(:,2),coords(:,3),'k');
axis equal
title(['Combined target in signature space'])
print('-dpng',['toroidal_sig_m=' num2str(m) '_n=' num2str(n) '.png']);

% Saves mat file for reverse imaging
save("coke_bottle.mat", "echos", "rxloc", "txloc", "scatloc")

