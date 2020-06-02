function coords=pcoa(dst,m)
% Compute classical Multi-dimensional scaling (Principal Coordinates Analysis)
%
% function coords=pcoa(dst,m)
%
% Input: dst = distance matrix (NxN)
%        m   = number of dimensions to use (scalar integer)
% Output: coords = matrix of coordinates (Nxm)
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

% Squared distance matrix
dst2=abs(dst).^2;

% Double centering
J=eye(size(dst2))-1/size(dst,1)*ones(size(dst,1),1)*ones(1,size(dst,1));
B=-0.5*J*dst2*J;

% Eigenanalysis
[u,sigma]=eigs(B,m);

% Construct coordinates
coords=u(:,1:m)*sqrt(sigma(1:m,1:m));
