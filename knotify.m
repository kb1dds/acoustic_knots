function ks=knotify(signal,ps,qs)
% Construct the torus knot signature for a signal
%
% function ks=knotify(signal,ps,qs)
%
% Input: signal = data input, independent columns (NxF, complex)
%        ps,qs  = winding numbers to test (Px1, Qx1)
% Output: ks = knot decomposition residuals (PxQ)
% Interpretation: the global minimum of ks is the best decomposition of the 
% signal as a torus knot
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
					     
% Preallocation
ks=zeros(numel(ps),numel(qs));

% Put signal into frequency domain, where it's easier to notch filter
signal_=fft(signal,[],1);

% Compute original signal power
ns=sum(abs(signal_(:)).^2);

for i=1:numel(ps),
  p=ps(i);
  for j=1:numel(qs),
    q=qs(j);
    
    % Construct filter tuned for winding number p
    p_idx = [p+1 size(signal_,1)-p]; % Symmetrized
    
    % Construct filter tuned for winding number q
    q_idx = [q+1 size(signal_,1)-q]; % Symmetrized

    % Apply notch filter to remove signal with winding numbers (p,q)
    signal_d=signal_;
    signal_d(p_idx,:)=0;
    signal_d(q_idx,:)=0;
    
    % Store residual as a fraction of remaining signal power
    ks(i,j)=sum(abs(signal_d(:)).^2)./ns;
  end
end
