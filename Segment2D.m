function Segment2D(tFun,segSet)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This demo shows how to compute the least squares on whitney 1-forms
% computed on segments on a triangles. Modify the parameter do_plot to
% generate or not the plot.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
%
% tFun   : selecting the function
%          1. constant, 2. smooth, 3. Runge+, 4. Runge
% segSet : selecting the type of segments
%          1. random segments, 2. small k-simplices
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

if nargin < 2, segSet=2; end
if nargin < 1, tFun=3; end

addpath('../chebfun')

fprintf('\n\n')
fprintf('-------------------------------------------------------------\n')
fprintf('                    function type : %2d                       \n',tFun)
fprintf('                    segments type : %2d                       \n',segSet)
% --------------- Parameters --------
do_plot = 1; % 0. no plot, 1. plot
MdegMax = 16;
Mdeg = 1:MdegMax;
Mev = 10000;

switch segSet 
    case 1
        X = rand(2*(3*MdegMax),2);
        X = X.*(X(:,1)+X(:,2)-1<=0) + (1-X).*(X(:,1)+X(:,2)-1>0);
    case 2
        X = equi_seg(3*MdegMax);
end

[xq,wq] = legpts(MdegMax+1);

[Fx,Fy] = function_test(tFun);


% Evaluation segments
rng(1981)
Xeval = rand(Mev,2); 
Xeval = Xeval.*(Xeval(:,1)+Xeval(:,2)-1<=0) + (1-Xeval).*(Xeval(:,1)+Xeval(:,2)-1>0);

h=1;
for M = Mdeg

% Set of the vertices of the triangle
Vert = [0, 0; 0, 1; 1,0];

if do_plot && M==5
    clear_figure(1);
    fig1 = figure(1);
    plot(Vert(:,1),Vert(:,2),'r*');
    hold on;
    for i=1:2:size(X,1)-1
    plot([X(i,1),X(i+1,1)],[X(i,2),X(i+1,2)]);
    end
end

% Computing tangents and cubature points on the segments
N = (X(2:2:end,:)-X(1:2:end-1,:))/2;

Px = X(1:2:end-1,1)+(xq'+1)/2.*(X(2:2:end,1)-X(1:2:end-1,1));
Py = X(1:2:end-1,2)+(xq'+1)/2.*(X(2:2:end,2)-X(1:2:end-1,2));

FFx = Fx(Px,Py);
FFy = Fy(Px,Py);

% Evaluating the quadrature on the segments
FF = (FFx.*N(:,1)+FFy.*N(:,2))*wq';

% Constructing the setting for the interpolation case
Xint = equi_seg_int(M);
Nint = (Xint(2:2:end,:)-Xint(1:2:end-1,:))/2;

PxInt = Xint(1:2:end-1,1)+(xq'+1)/2.*(Xint(2:2:end,1)-Xint(1:2:end-1,1));
PyInt = Xint(1:2:end-1,2)+(xq'+1)/2.*(Xint(2:2:end,2)-Xint(1:2:end-1,2));

FFxInt = Fx(PxInt,PyInt);
FFyInt = Fy(PxInt,PyInt);
FFint = (FFxInt.*Nint(:,1)+FFyInt.*Nint(:,2))*wq';

% Construction for the evaluation segments
Neval = (Xeval(2:end,:)-Xeval(1:end-1,:))./2;

PxEval = Xeval(1:end-1,1)+(xq'+1)/2.*(Xeval(2:end,1)-Xeval(1:end-1,1));
PyEval = Xeval(1:end-1,2)+(xq'+1)/2.*(Xeval(2:end,2)-Xeval(1:end-1,2));

FFxEval = Fx(PxEval,PyEval);
FFyEval = Fy(PxEval,PyEval);
FFeval = (FFxEval.*Neval(:,1)+FFyEval.*Neval(:,2))*wq';

% Selecting the multi-indices
id1 = index2D(M-1);
j1 = (id1(:,1)==0);
id2 = id1(j1,:);

% Generating the Vandermonde-like matrices
V = []; Vint = []; Veval = [];
for it = 0:2
    idx = id1;
    if it == 2
        idx = id2;
    end
    for jj = 1:size(idx,1)
        alpha = idx(jj,:);
        lambda = ((1-Px-Py).^alpha(1)).*(Px.^alpha(2)).*(Py.^alpha(3));
        lambdaInt = ((1-PxInt-PyInt).^alpha(1)).*(PxInt.^alpha(2)).*(PyInt.^alpha(3));
        lambdaEval = ((1-PxEval-PyEval).^alpha(1)).*(PxEval.^alpha(2)).*(PyEval.^alpha(3));
        switch it
            case 0
                U1 = 1-Py;
                U2 = Px;
                U1int = 1-PyInt;
                U2int = PxInt;
                U1eval = 1-PyEval;
                U2eval = PxEval;
            case 1
                U1 = Py;
                U2 = 1-Px;
                U1int = PyInt;
                U2int = 1-PxInt;
                U1eval = PyEval;
                U2eval = 1-PxEval;
            case 2
                U1 = -Py;
                U2 = Px;
                U1int = -PyInt;
                U2int = PxInt;
                U1eval = -PyEval;
                U2eval = PxEval;
        end
        V(:,end+1) = (lambda.*(U1.*N(:,1) + U2.*N(:,2)))*wq';
        Vint(:,end+1) = (lambdaInt.*(U1int.*Nint(:,1) + U2int.*Nint(:,2)))*wq';
        Veval(:,end+1) = (lambdaEval.*(U1eval.*Neval(:,1) + U2eval.*Neval(:,2)))*wq';
    end
end

% LS approximation
[Q,R] = qr(V,0);
cLS = R\(Q'*FF);
% Interpolation
cInt = Vint\FFint;

% Error in C_0 norm
SegN = vecnorm(2*Neval')';
errLS(h) = norm((FFeval-Veval*cLS)./SegN,'inf');
errINT(h) = norm((FFeval-Veval*cInt)./SegN,'inf');
h = h+1;
% STATISTICS.
fprintf('\n \t -> deg: %3d errLS: %1.3e  errINT: %1.3e', M, errLS(h-1), errINT(h-1)); 
end

if do_plot
    clear_figure(2);
    fig2 = figure(2);
    semilogy(Mdeg,errLS,'o-b')
    hold on;
    semilogy(Mdeg,errINT,'o-r')
    legend('Least squares','Interpolated','location','northeast')
    title('Error: interpolated vs least squares')
    xlabel('degree')
    ylabel('error')
end
fprintf('\n\n')

function [Fx,Fy] = function_test(tFun)

switch tFun
    case 1
        Fx = @(x,y) 1+0.*x;
        Fy = @(x,y) 1+0.*x;
    case 2
        Fx = @(x,y) exp(x+y);
        Fy = @(x,y) sin(x.*y);
    case 3
        Fx = @(x,y) 1./(1+100.*(y-1/3).^2).*(200.*(x-1/3))./(1+100.*(x-1/3).^2).^2;
        Fy = @(x,y) 1./(1+100.*(x-1/3).^2).*(200.*(y-1/3))./(1+100.*(y-1/3).^2).^2;
    case 4
        Fx = @(x,y) 1./(1+25*(x-1/3).^2);
        Fy = @(x,y) 1./(1+25*(y-1/3).^2);
end

function clear_figure(nfig)

h=figure(nfig);
f_nfig=ishandle(h)&&strcmp(get(h,'type'),'figure');
if f_nfig
    clf(nfig);
end