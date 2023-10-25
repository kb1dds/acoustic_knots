

% Sets Axes for the image
 xaxis = -1.5:0.01:1;
 yaxis = -1:0.01:1;
 [x,y]=meshgrid(xaxis,yaxis);

% Loatds .mat file to create image out of
 load("coke_bottle.mat");
 samplePoints = [x(:),y(:)];
 object_image = tdbp_compressed(samplePoints,rxloc,txloc,0,1000,ifft(echos,[],2).');

 reshape_object_image = reshape(object_image, size(x));

%  Plots reshaped image
 figure
 imagesc(xaxis,yaxis,abs(reshape_object_image))
 colorbar
 hold on