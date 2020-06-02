function mat=isotropicMatrix(wavelength,sourcePoints,samplePoints)
% Compute the field induced by a collection of isotropic radiators
%
% function mat=isotropicMatrix(wavelength,sourcePoints,samplePoints)
%
% Input: wavelength   = operating wavelength (meters)
%        sourcePoints = location of point sources to solve (PxD, meters)
%        samplePoints = location of point where field is desired (SxD, meters)
% Output: mat = matrix that converts sources into samples (SxP)
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
  
% Wavenumber
k = 2*pi/wavelength;

% Sample points (NxP)
rangeVector=bsxfun(@minus,permute(samplePoints,[1 3 2]),permute(sourcePoints,[3 1 2]));
range=sqrt(sum(rangeVector.^2,3));

% Phase offsets (NxP) (dimensionless)
mat=exp(-sqrt(-1)*k.*range)./(4*pi*(eps+range));
mat(abs(range)<10*eps)=1; % Avoid divide-by-zero if source point overlaps sample point
