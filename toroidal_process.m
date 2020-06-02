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

% Scatterer relative rotation angle (deg)
offset_angle=28;

% Number of measurements to collect
nlooks = 200;

% Number of frequencies to collect
nfrequencies = 60;

% Frequencies to examine (Hertz)
frequencies=linspace(0,600,nfrequencies+1);
frequencies(1)=[];

% Scatterer configuration
m = 2; % Winding number
theta_m=(0:2*pi/m:2*pi).';
theta_m=theta_m(1:m);
radius_m=1;
n = 3; % Winding number
theta_n=offset_angle*pi/180+(0:2*pi/n:2*pi).';
theta_n=theta_n(1:n);
radius_n=1;
scatloc_m=[radius_m*cos(theta_m),radius_m*sin(theta_m)];
scatloc_n=[radius_n*cos(theta_n),radius_n*sin(theta_n)];
scatloc=[scatloc_m; scatloc_n];

%% You should not need to change code below this line if you are just experimenting with parameters!

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

% Construct toroidal coordinates for each transmitter and receiver location
txtoroidal=mod([theta_platform*n, theta_platform*m],2*pi);
rxtoroidal=mod([theta_platform*n, theta_platform*m],2*pi);

echos_m=zeros(nlooks,length(frequencies));
echos_n=zeros(nlooks,length(frequencies));
echos=zeros(nlooks,length(frequencies));
echos_m_torus=zeros(length(theta_platform),length(frequencies));
echos_n_torus=zeros(length(theta_platform),length(frequencies));
for i=1:length(frequencies)
  % Propagate from transmitter locations to scatterers
  tx2scat_m=isotropicMatrix(c/frequencies(i),txloc,scatloc_m);
 
  % Propagate from scatterers to transmitter
  scat2rx_m=isotropicMatrix(c/frequencies(i),scatloc_m,rxloc);

  % Store noiseless frequency response
  echos_m(:,i)=diag(scat2rx_m*bsxfun(@times,scatcross_m,tx2scat_m));
  
  % Propagate from transmitter locations to scatterers
  tx2scat_n=isotropicMatrix(c/frequencies(i),txloc,scatloc_n);
 
  % Propagate from scatterers to transmitter
  scat2rx_n=isotropicMatrix(c/frequencies(i),scatloc_n,rxloc);

  % Store noiseless frequency response
  echos_n(:,i)=diag(scat2rx_n*bsxfun(@times,scatcross_n,tx2scat_n));
  
  % Propagate from transmitter locations to scatterers
  tx2scat=isotropicMatrix(c/frequencies(i),txloc,scatloc);
 
  % Propagate from scatterers to transmitter
  scat2rx=isotropicMatrix(c/frequencies(i),scatloc,rxloc);

  % Store noiseless frequency response
  echos(:,i)=diag(scat2rx*bsxfun(@times,scatcross,tx2scat));
  
  %% Compute overall data (for the space of targets) from which our data is a slice
  
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

figure
colormap gray
imagesc(frequencies,theta_platform*180/pi,abs(echos));
set(gca,'ydir','normal');
axis normal
caxis([0,0.004])
colorbar
xlabel('Frequency (Hz)')
ylabel('Look angle (deg)')
title(['Combined target'])
print('-dpng',['toroidal_flat_m=' num2str(m) '_n=' num2str(n) '.png']);

figure
rangeVector=bsxfun(@minus,permute(echos,[1 3 2]),permute(echos,[3 1 2]));
dst=sqrt(sum(abs(rangeVector).^2,3));
coords=pcoa(dst,3);
plot3(coords(:,1),coords(:,2),coords(:,3),'k');
axis equal
title(['Combined target in signature space'])
print('-dpng',['toroidal_sig_m=' num2str(m) '_n=' num2str(n) '.png']);

% Display one frequency on toroidal coordinates, from which the collection is a slice
figure
colormap gray
idx=30;%length(frequencies);
imagesc(linspace(0,360,nlooks),linspace(0,360,nlooks),abs(echos_torus(:,:,idx)));
set(gca,'ydir','normal');
caxis([0,0.004]);
colorbar
xlabel(['n=' num2str(n) ' (deg)']);
ylabel(['m=' num2str(m) ' (deg)']);
title(['(m=' num2str(m) ',n=' num2str(n) ')-torus slice, f=' num2str(frequencies(idx))  ' Hz']);
print('-dpng',['toroidal_slice_f=' num2str(frequencies(idx)) '_m=' num2str(m) '_n=' num2str(n) '.png']);
hold on
plot(txtoroidal(1:end-1,1)*180/pi,txtoroidal(1:end-1,2)*180/pi,'w*','linewidth',3);
print('-dpng',['toroidal_slice_f=' num2str(frequencies(idx)) '_m=' num2str(m) '_n=' num2str(n) '_traj.png']);

% Extract collection as a slice
coll_slice=interp2(linspace(0,2*pi,nlooks),linspace(0,2*pi,nlooks),abs(echos_torus(:,:,idx)),txtoroidal(:,1),txtoroidal(:,2));
figure
plot(theta_platform*180/pi,coll_slice,'ko',theta_platform*180/pi,abs(echos(:,idx)),'k');
legend('Sliced from torus','Original');
%plot(theta_platform,coll_slice,'bo',theta_platform,abs(echos(:,idx)),theta_platform,abs(echos_m(:,idx)),theta_platform,abs(echos_n(:,idx)));
%legend('Sliced from torus','Sum',['m=' num2str(m)],['n=' num2str(n)]);
xlim([0,360]);
xlabel('Look angle (deg)');
ylabel('Cross section');
title(['(m=' num2str(m) ',n=' num2str(n) ')-torus slice versus collection cross-check plots, f=' num2str(frequencies(idx))  ' Hz']);
print('-dpng',['toroidal_crosscheck_f=' num2str(frequencies(idx)) '_m=' num2str(m) '_n=' num2str(n) '.png']);

% Draw one frequency *on* the torus with the selected slice
figure
colormap gray
sf=0.95;
[thetas,phis]=meshgrid(linspace(0,2*pi,nlooks),linspace(0,2*pi,nlooks));
xx=(2+sf*cos(phis)).*cos(thetas);
yy=(2+sf*cos(phis)).*sin(thetas);
zz=sf*sin(phis);
xsl=(2+sf*cos(txtoroidal(1:end-1,1))).*cos(txtoroidal(1:end-1,2));
ysl=(2+sf*cos(txtoroidal(1:end-1,1))).*sin(txtoroidal(1:end-1,2));
zsl=sf*sin(txtoroidal(1:end-1,1));
surf(xx,yy,zz,abs(echos_torus(:,:,idx)));
shading flat
hold on
caxis([0,0.004])
plot3(xsl,ysl,zsl,'w.')
axis equal
print('-dpng',['toroidal_3d_slice_f=' num2str(frequencies(idx)) '_m=' num2str(m) '_n=' num2str(n) '.png']);

% Estimate winding number by repeated notch filtering
figure;
imagesc(2:10,2:10,knotify(echos(:,idx),2:10,2:10));
colormap gray
xlabel('Winding number m');
xlabel('Winding number n');
title('Residual signal after removing a torus knot of a given type');
print('-dpng',['toroidal_notching_' num2str(frequencies(idx)) '_m=' num2str(m) '_n=' num2str(n) '.png']);
