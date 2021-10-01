% Composite scatter with/without trajectory distortion
%
% This is to be run as a script.
% Many figures and PNG files are produced as output,
%  titles and other annotations are made automatically,
%   and hopefully are reasonable!
% To change what's plotted, comment out what you don't want!
%
% Copyright (c) 2021 Michael Robinson
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

% Scatterer relative rotation angle (deg)
offset_angle=28;

% Number of measurements to collect
nlooks = 100;

% Number of frequencies to collect
nfrequencies = 100;

% Frequencies to examine (Hertz)
frequencies=linspace(0,600,nfrequencies+1);
frequencies(1)=[];

% Scatterer configuration
% Scatterers are in 2d as collections of points
% There are three such "composite" scatterers:
%  two lie at the vertices of regular polygons (with m and n sides)
% The n scatterer is viewed with a circular trajectory that's unevently traversed
m = 3; % Winding number
radius_m=1; % (meters)
n = 3; % Winding number
radius_n=1; % (meters)

% Receiver noise level
noise_level=0; % Variance

%% You should not need to change code below this line if you are just experimenting with parameters!

% Scatterer locations
theta_m=(0:2*pi/m:2*pi).';
theta_m=theta_m(1:m);
theta_n=offset_angle*pi/180+(0:2*pi/n:2*pi).';
theta_n=theta_n(1:n);
scatloc_m=[radius_m*cos(theta_m),radius_m*sin(theta_m)];
scatloc_n=[radius_n*cos(theta_n),radius_n*sin(theta_n)];
scatloc=[scatloc_m; scatloc_n];

% Construct scatterer cross sections
scatcross_m=ones(size(scatloc_m,1),1);
scatcross_n=ones(size(scatloc_n,1),1);
scatcross=ones(size(scatloc,1),1);

% Construct transmitter and receiver locations
theta_platform=linspace(0,2*pi,nlooks).';

% TX/RX for "m" scatterer is a circle traversed at a constant speed
txloc_m=[txrange*cos(theta_platform),txrange*sin(theta_platform)];
rxloc_m=[rxrange*cos(theta_platform),rxrange*sin(theta_platform)];

% TX/RX for "n" scatterer is a circle traversed at a variable speed
% Note that it is intentional that the plots of the "n" scatterer use
% theta_platform as an axis label, not theta_platform_n
theta_platform_n=theta_platform-0.3*(1-cos(3*theta_platform));
txloc_n=[txrange*cos(theta_platform_n),txrange*sin(theta_platform_n)];
rxloc_n=[rxrange*cos(theta_platform_n),rxrange*sin(theta_platform_n)];

echos_m=zeros(nlooks,length(frequencies));
echos_n=zeros(nlooks,length(frequencies));
echos=zeros(nlooks,length(frequencies));
for i=1:length(frequencies)
  %% First scatterer: regular m-polygon with scatterers at its corners
  
  % Propagate from transmitter locations to scatterers
  tx2scat_m=isotropicMatrix(c/frequencies(i),txloc_m,scatloc_m);
 
  % Propagate from scatterers to transmitter
  scat2rx_m=isotropicMatrix(c/frequencies(i),scatloc_m,rxloc_m);

  % Store noiseless frequency response
  echos_m(:,i)=diag(scat2rx_m*bsxfun(@times,scatcross_m,tx2scat_m));

  %% Second scatterer: regular n-polygon with scatterers at its corners
  
  % Propagate from transmitter locations to scatterers
  tx2scat_n=isotropicMatrix(c/frequencies(i),txloc_n,scatloc_n);
 
  % Propagate from scatterers to transmitter
  scat2rx_n=isotropicMatrix(c/frequencies(i),scatloc_n,rxloc_n);

  % Store noiseless frequency response
  echos_n(:,i)=diag(scat2rx_n*bsxfun(@times,scatcross_n,tx2scat_n));
end

% Add noise
echos_m=echos_m+randn(size(echos_m))*noise_level;
echos_n=echos_n+randn(size(echos_n))*noise_level;

%% Plotting code below this line!

% Display the echo data
figure
colormap gray
imagesc(frequencies,theta_platform*180/pi,abs(echos_m));
set(gca,'ydir','normal');
axis normal
caxis([0,0.002])
colorbar
xlabel('Frequency (Hz)')
ylabel('Look angle (deg)')
title(['Target with ' num2str(m) '-fold symmetry'])
print('-dpng',['toroidal_flat_m=' num2str(m) '.png']);

figure
colormap gray
imagesc(abs(ifft(echos_m,[],2)))
set(gca,'ydir','normal');
axis normal
caxis([0,0.002])
colorbar
xlabel('Time (samples)')
ylabel('Look angle (deg)');
title(['Target with ' num2str(m) '-fold symmetry'])
print('-dpng',['toroidal_ranges_flat_m=' num2str(m) '.png']);

figure
colormap gray
imagesc(frequencies,theta_platform*180/pi,abs(echos_n));
set(gca,'ydir','normal');
axis normal
caxis([0,0.002])
colorbar
xlabel('Frequency (Hz)')
ylabel('Look angle (deg)')
title(['Target with ' num2str(n) '-fold symmetry'])
print('-dpng',['toroidal_flat_n=' num2str(n) '.png']);

figure
colormap gray
imagesc(abs(ifft(echos_n,[],2)))
set(gca,'ydir','normal');
axis normal
caxis([0,0.002])
colorbar
xlabel('Time (samples)')
ylabel('Look angle (deg)');
title(['Target with ' num2str(n) '-fold symmetry'])
print('-dpng',['toroidal_ranges_flat_n=' num2str(n) '.png']);

% PCA plots of both
figure
idq=[zeros(size(echos_m,1),1); ones(size(echos_n,1),1)];
shifted_combined_echos=[echos_m, circshift(echos_m,[3,0]);
			echos_n, circshift(echos_n,[3,0])];
rangeVector=bsxfun(@minus,permute(shifted_combined_echos,[1 3 2]),
		   permute(shifted_combined_echos,[3 1 2]));
dst=sqrt(sum(abs(rangeVector).^2,3));
coords=pcoa(dst,3);
plot3(coords(idq==0,1),coords(idq==0,2),coords(idq==0,3),'k','linewidth',3);
axis equal
title(['Combined target in signature space'])
print('-dpng',['toroidal_sig_m=' num2str(m) '.png']);
figure
plot3(coords(idq==1,1),coords(idq==1,2),coords(idq==1,3),'k','linewidth',3);
axis equal
title(['Combined target in signature space'])
print('-dpng',['toroidal_sig_n=' num2str(m) '.png']);
