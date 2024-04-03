function Segment1D(tFun,segSet)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This demo shows how to compute the least squares on whitney 1-forms
% computed on segments on a line. Modify the parameter do_plot to
% generate or not the plot.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
%
% tFun   : selecting the function
%          1. constant, 2. linear, 3. smooth, 4. Runge 
% segSet : selecting the type of segments
%          1. equispaced, 2. random ordered 3. totally random 
%
%--------------------------------------------------------------------------
% DATES:
%--------------------------------------------------------------------------
% Written on April 03, 2024 
%
% Authors involved:
% L. Bruni Bruno, G. Elefante
%--------------------------------------------------------------------------
% LICENSE
%--------------------------------------------------------------------------
% Copyright (C) 2024-
% L. Bruni Bruno, G. Elefante
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%
% Authors:
%
% L. Bruni Bruno
% G. Elefante
%
% Date: April 04, 2024.
%--------------------------------------------------------------------------

if nargin < 2, segSet=1; end

if nargin < 1, tFun=4; end

% ------------- Parameters -------------------------
degMax = 25;
nseg = 1000;
do_plot = 1;

% Settings
f = function_test(tFun);

npti = nseg+1;
switch segSet 
    case 1
        pti = linspace(-1,1,npti);
    case 2
        rng(1981)
        pti = sort(rand(1,npti))*2-1;
    case 3
        rng(1981)
        pti = rand(1,npti)*2-1;
end

fprintf('\n\n')
fprintf('-------------------------------------------------------------\n')
fprintf('                    function type: %2d                       \n',tFun)
fprintf('                    segments type: %2d                       \n',segSet)

for deg = 1:degMax
    % Computing the f values and the
    for i = 1:nseg
        F(i,1) = integral(f,pti(i),pti(i+1))/abs(pti(i+1)-pti(i));
        for j = 1:deg+1
            W(i,j) = (pti(i+1)^j - pti(i)^j)/j/abs(pti(i+1)-pti(i));
        end
    end

    % Computing the interpolation matrix

    pti_int = linspace(-1,1,deg+2);

    for i = 1:deg+1
        fint(i,1) = integral(f,pti_int(i),pti_int(i+1))/(pti_int(i+1) - pti_int(i));
        for j = 1:deg+1
            V(i,j) = (pti_int(i+1)^j - pti_int(i)^j)/j/(pti_int(i+1) - pti_int(i));
        end
    end


    for i = 1:deg+1
        M(:,i) = pti(:).^(i-1);
    end

    % Computing the coefficients of the interpolant and of the LS approx

    coeffinterp = V\fint;
    [Q,R] = qr(W,0);
    coeffls = R\(Q'*F);

    % Error on a dense mesh of segments
    errINT(deg) = max(abs(M*coeffinterp-f(pti)'));
    errLS(deg) = max(abs(M*coeffls-f(pti)'));

    fprintf('\n \t -> deg: %3d errLS: %1.3e  errINT: %1.3e', deg, errLS(deg), errINT(deg));
end

deg = 5;
nsegSet = 10:5:200; h = 1;
pti_int = linspace(-1,1,deg+2);
for npti = nsegSet+1
    switch segSet
        case 1
            pti = linspace(-1,1,npti);
        case 2
            rng(1981)
            pti = sort(rand(1,npti))*2-1;
        case 3
            rng(1981)
            pti = rand(1,npti)*2-1;
    end
    % riempio il termine noto e la matrice per i minimi quadrati
    W = []; F = [];
    for i = 1:npti-1
        F(i,1) = integral(f,pti(i),pti(i+1))/abs(pti(i+1)-pti(i));
        for j = 1:deg+1
            W(i,j) = (pti(i+1)^j - pti(i)^j)/j/abs(pti(i+1)-pti(i));
        end
    end

    % riempio il termine noto e la matrice per l'interpolazione

    %     pti_int = linspace(-1,1,deg+2);

    V = []; fint = [];
    for i = 1:deg+1
        fint(i,1) = integral(f,pti_int(i),pti_int(i+1))/(pti_int(i+1) - pti_int(i));
        for j = 1:deg+1
            V(i,j) = (pti_int(i+1)^j - pti_int(i)^j)/j/(pti_int(i+1) - pti_int(i));
        end
    end

    M = [];
    for i = 1:deg+1
        M(:,i) = pti'.^(i-1);
    end

    coeffinterp = V\fint;
    [Q,R] = qr(W,0);
    coeffls = R\(Q'*F);

    errINT2(h) = max(abs(M*coeffinterp-f(pti)'));
    errLS2(h) = max(abs(M*coeffls-f(pti)'));
    h = h+1;
end

fprintf('\n\n')

if do_plot
    clear_figure(1);
    fig1 = figure(1);
    semilogy(1:degMax, errINT,'o-b', 1:degMax, errLS,'o-r')
    legend('Interpolated','Least squares','location','northeast')
    title('Error: interpolated vs least squares')
    xlabel('degree')
    ylabel('error')

    clear_figure(2);
    fig2 = figure(2);
    semilogy(nsegSet, errINT2,'o-b', nsegSet, errLS2,'o-r')
    legend('Interpolated','Least squares','location','southeast')
    title('Error: interpolated vs least squares')
    axis([0 200 0.2 0.5])
    xlabel('number of segments')
    ylabel('error')
end

function [f] = function_test(tFun)

switch tFun
    case 1
        f = @(x) 1+0.*x;
    case 2
        f = @(x) x;
    case 3
        f = @(x) sin(x).*exp(x)./(x.^2+x+1);
    case 4
        f = @(x) 1./(1+25*x.^2);
end

function clear_figure(nfig)

h=figure(nfig);
f_nfig=ishandle(h)&&strcmp(get(h,'type'),'figure');
if f_nfig
    clf(nfig);
end