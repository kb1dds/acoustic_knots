% Basic CSAS point-scatterer target simulator
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
   
% A jitter parameter allows points to deviate from their true positions by a specified amount
jitter=0.01; % (meters, variance)

target_name={'coke_bottle','cup_open','cup_capped','pipe_open','pipe_capped'};

for target_idx = 1:numel(target_name),

% Scatterer locations
if(target_idx == 1),
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

    d = 1;
    scatloc_n = [];
    for j = -0.2:0.01:0.2
        scatloc_n(d,1) = -1;
        scatloc_n(d,2) = j;
        d = d+1; 
    end
end 

if(target_idx == 2),
    %Cup without a lid geometry

    c = 1;
    scatloc_m = [];
    for i = -.5:0.01:.5
        scatloc_m(c,1) = i;
        scatloc_m(c,2) = 0.2+(i+.5)*.15;
        scatloc_m(c+1,1) = i;
        scatloc_m(c+1,2) = -0.2-(i+.5)*.15;
        c = c+2; 
    end
 
    d = 1;
    scatloc_n = [];
    for j = -0.2:0.01:0.2
        scatloc_n(d,1) = -0.5;
        scatloc_n(d,2) = j;
        d = d+1; 
    end
end

if(target_idx == 3),
    % Cup with a lid geometry

    c = 1;
    scatloc_m = [];
    for i = -.5:0.01:.5
        scatloc_m(c,1) = i;
        scatloc_m(c,2) = 0.2+(i+.5)*.15;
        scatloc_m(c+1,1) = i;
        scatloc_m(c+1,2) = -0.2-(i+.5)*.15;
        c = c+2; 
    end

    d = 1;
    scatloc_n = [];
    for j = -0.2:0.01:0.2
        scatloc_n(d,1) = -0.5;
        scatloc_n(d,2) = j;
        d = d+1; 
    end
 
    e = d+1;
    for j = -0.35:0.01:0.35
        scatloc_n(e,1) = 0.5;
        scatloc_n(e,2) = j;
        e = e+1; 
    end
end

if(target_idx == 4),
    % Pipe with no lid geometry

    c = 1;
    scatloc_m = [];
    for i = -.5:0.01:.5
        scatloc_m(c,1) = i;
        scatloc_m(c,2) = 0.05;
        scatloc_m(c+1,1) = i;
        scatloc_m(c+1,2) = -0.05;
        c = c+2; 
    end
 
    d = 1;
    scatloc_n = [];
    for j = -0.05:0.01:0.05
        scatloc_n(d,1) = 0;
        scatloc_n(d,2) = 0;
        scatloc_n(d+1,1) = 0;
        scatloc_n(d+1,2) = 0;
        d = d+2; 
    end
end

if(target_idx == 5),
    % Pipe with lid geometry

    c = 1;
    scatloc_m = [];
    for i = -.5:0.01:.5
        scatloc_m(c,1) = i;
        scatloc_m(c,2) = 0.05;
        scatloc_m(c+1,1) = i;
        scatloc_m(c+1,2) = -0.05;
        c = c+2; 
    end
    d = 1;
    scatloc_n = [];
    for j = -0.05:0.01:0.05
        scatloc_n(d,1) = 0.5;
        scatloc_n(d,2) = j;
        scatloc_n(d+1,1) = -0.5;
        scatloc_n(d+1,2) = j;
        d = d+2; 
    end
end


scatloc_m=scatloc_m+randn(size(scatloc_m))*jitter;
scatloc_n=scatloc_n+randn(size(scatloc_n))*jitter;
scatloc=[scatloc_m; scatloc_n];

% Construct scatterer cross sections
scatcross=ones(size(scatloc,1),1).*exp(sqrt(-1)*2*pi*randn(size(scatloc,1),1)*90/360);

fp=fopen([target_name{target_idx} '_scatterers.csv'],'wt');
fdisp(fp,'x,y,cross_mag,cross_phase');
dlmwrite(fp,[scatloc, abs(scatcross), angle(scatcross)],',');
fclose(fp);
end

