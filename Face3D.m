function Face3D(tFun,facSet)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This demo shows how to compute the least squares on whitney 2-forms
% computed on faces on a tetrahedron. Modify the parameter do_plot to
% generate or not the plot.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
%
% tFun   : selecting the function
%          1. constant, 2. linear, 3. smooth, 4. Runge 5. Runge v2
% facSet : selecting the type of segments
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

if nargin < 2, facSet=2; end
if nargin < 1, tFun=5; end

do_plot = 1;

% Degrees of the interpolation and ls problem
nDegMax = 8;
nDeg = 1:nDegMax;

% Generating the quadrature nodes and weights for triangles
[xyw,~]=set_xiao_gimbutas_standard(nDegMax+1);
XYZW = [xyw(:,1:2),zeros(size(xyw,1),1),xyw(:,3)];

% Generating evaluation faces
[XEval,~] = triangles3D(nDegMax+15);

% LS faces
switch facSet
    case 1
        rng(1981)
        M = 4*nchoosek(nDegMax+7,nDegMax+6);
        Xls = rand(3*M,3); XC = max(abs(Xls),[],2); XT = sum(abs(Xls),2);
        Xls = XC./XT.*Xls;
    case 2
        [Xls,~] = triangles3D(nDegMax+10); 
end

% Generating test Function
[Fx,Fy,Fz] = function_test(tFun);

fprintf('\n\n')
fprintf('-------------------------------------------------------------\n')
fprintf('                    function type : %2d                      \n',tFun)
fprintf('                    faces type    : %2d                      \n',facSet)

h = 1;
for n = nDeg

% Set of the vertices of the triangle
Vert = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];

[XX,Xint] = triangles3D(n); 

if do_plot && (n == 5)
    clear_figure(1);
    fig1 = figure(1);
    plot3(Vert(:,1),Vert(:,2),Vert(:,3),'r*');
    hold on;
    for i=1:3:size(XX,1)-2
    plot3([XX(i,1),XX(i+1,1),XX(i+2,1),XX(i,1)],[XX(i,2),XX(i+1,2),XX(i+2,2),XX(i,2)],[XX(i,3),XX(i+1,3),XX(i+2,3),XX(i,3)]);
    end
end

W1ls = Xls(2:3:end,:)-Xls(1:3:end,:);
W2ls = Xls(3:3:end,:)-Xls(1:3:end,:);

Nls = cross(W1ls,W2ls,2);


Pls = []; it = 1;
for i=1:3:size(Xls,1)
    A = [W2ls(it,:);W1ls(it,:);Nls(it,:)];
    b = Xls(i,:);
    Pls = [Pls;(XYZW(:,1:3)*A+b)'];
    it = it+1;
end

PxLS = Pls(1:3:end,:); PyLS = Pls(2:3:end,:); PzLS = Pls(3:3:end,:);

FFxLS = Fx(PxLS,PyLS,PzLS);
FFyLS = Fy(PxLS,PyLS,PzLS);
FFzLS = Fz(PxLS,PyLS,PzLS);

FFls = (FFxLS.*Nls(:,1)+FFyLS.*Nls(:,2)+FFzLS.*Nls(:,3))*XYZW(:,4);

% Quadrature points in the interpolation faces
W1int = Xint(2:3:end,:)-Xint(1:3:end,:);
W2int = Xint(3:3:end,:)-Xint(1:3:end,:);

Nint = cross(W1int,W2int,2);

Pint = []; it = 1;
for i=1:3:size(Xint,1)
    A = [W2int(it,:);W1int(it,:);Nint(it,:)];
    b = Xint(i,:);
    Pint = [Pint;(XYZW(:,1:3)*A+b)'];
    it = it+1;
end

PxINT = Pint(1:3:end,:); PyINT = Pint(2:3:end,:); PzINT = Pint(3:3:end,:);

% Evaluating the test function on the interpolation faces
FFxINT = Fx(PxINT,PyINT,PzINT);
FFyINT = Fy(PxINT,PyINT,PzINT);
FFzINT = Fz(PxINT,PyINT,PzINT);

FFint = (FFxINT.*Nint(:,1)+FFyINT.*Nint(:,2)+FFzINT.*Nint(:,3))*XYZW(:,4);

% Quadrature nodes on the evaluation faces
W1Eval = XEval(2:3:end,:)-XEval(1:3:end,:);
W2Eval = XEval(3:3:end,:)-XEval(1:3:end,:);

NEval = cross(W1Eval,W2Eval,2);

PEval = []; it = 1;
for i=1:3:size(XEval,1)
    A = [W2Eval(it,:);W1Eval(it,:);NEval(it,:)];
    b = XEval(it,:);
    PEval = [PEval;(XYZW(:,1:3)*A+b)'];
    it = it+1;
end

PxEval = PEval(1:3:end,:); PyEval = PEval(2:3:end,:); PzEval = PEval(3:3:end,:);

% Selecting the multi-indices
id1 = index3D(n-1);
j1 = (id1(:,1)==0);
id2 = id1(j1,:);


% Forming the Vandermonde-like matrix
Vls = []; VEval = []; Vint = [];
for it = 0:3
    idx = id1;
    if it == 3
        idx = id2;
    end
    for jj = 1:size(idx,1)
        alpha = idx(jj,:);
        lambdaLS = ((1-PxLS-PyLS-PzLS).^alpha(1)).*(PxLS.^alpha(2)).*(PyLS.^alpha(3)).*(PzLS.^alpha(4));
        lambdaINT = ((1-PxINT-PyINT-PzINT).^alpha(1)).*(PxINT.^alpha(2)).*(PyINT.^alpha(3)).*(PzINT.^alpha(4));
        lambdaEval = ((1-PxEval-PyEval-PzEval).^alpha(1)).*(PxEval.^alpha(2)).*(PyEval.^alpha(3)).*(PzEval.^alpha(4));
        switch it
            case 0
                U1ls = -PxLS;
                U2ls = PyLS;
                U3ls = 1-PzLS;
                U1int = -PxINT;
                U2int = PyINT;
                U3int = 1-PzINT;
                U1eval = -PxEval;
                U2eval = PyEval;
                U3eval = 1-PzEval;
            case 1
                U1ls = PxLS;
                U2ls = 1-PyLS;
                U3ls = PzLS;
                U1int = PxINT;
                U2int = 1-PyINT;
                U3int = PzINT;
                U1eval = PxEval;
                U2eval = 1-PyEval;
                U3eval = PzEval;
            case 2
                U1ls = 1-PxLS;
                U2ls = PyLS;
                U3ls = -PzLS;
                U1int = 1-PxINT;
                U2int = PyINT;
                U3int = -PzINT;
                U1eval = 1-PxEval;
                U2eval = PyEval;
                U3eval = -PzEval;
            case 3
                U1ls = PxLS;
                U2ls = -PyLS;
                U3ls = PzLS;
                U1int = PxINT;
                U2int = -PyINT;
                U3int = PzINT;
                U1eval = PxEval;
                U2eval = -PyEval;
                U3eval = PzEval;
        end
        Vls(:,end+1) = (lambdaLS.*(U1ls.*Nls(:,1) + U2ls.*Nls(:,2) + U3ls.*Nls(:,3)))*XYZW(:,4);
        Vint(:,end+1) = (lambdaINT.*(U1int.*Nint(:,1) + U2int.*Nint(:,2) + U3int.*Nint(:,3)))*XYZW(:,4);
        VEval(:,end+1) = (lambdaEval.*(U1eval.*NEval(:,1) + U2eval.*NEval(:,2) + U3eval.*NEval(:,3)))*XYZW(:,4);
    end
end

% Computing the solution of the ls problem and the interpolation prob
[Q,R] = qr(Vls,0);
cLS = R\(Q'*FFls);

cINT = Vint\FFint;

% Evaluating the function on the evaluation faces
FFxEval = Fx(PxEval,PyEval,PzEval);
FFyEval = Fy(PxEval,PyEval,PzEval);
FFzEval = Fz(PxEval,PyEval,PzEval);

FFeval = (FFxEval.*NEval(:,1)+FFyEval.*NEval(:,2)+FFzEval.*NEval(:,3))*XYZW(:,4);

% Computing the error of the interpolant and the solution of the ls problem
AreaN = vecnorm(NEval')';
errLS(h) = norm((FFeval-VEval*cLS)./AreaN,'inf');
errINT(h) = norm((FFeval-VEval*cINT)./AreaN,'inf');
h = h+1;
fprintf('\n \t -> deg: %3d errLS: %1.3e  errINT: %1.3e', n, errLS(h-1), errINT(h-1)); 
end

fprintf('\n\n')

if do_plot
    clear_figure(2);
    fig2 = figure(2);
    semilogy(nDeg,errLS,'o-b')
    hold on
    semilogy(nDeg,errINT,'o-r')
    title('Error: interpolated vs least squares')
    legend('Least squares','Interpolated','location','southwest')
    xlabel('degree')
    ylabel('error')
    xlim([0 nDegMax+1])
end


function [Fx,Fy, Fz] = function_test(tFun)

switch tFun 
    case 1
        Fx = @(x,y,z) 1+0.*x;
        Fy = @(x,y,z) 1+0.*x;
        Fz = @(x,y,z) 1+0.*x;
    case 2
        Fx = @(x,y,z) x;
        Fy = @(x,y,z) y;
        Fz = @(x,y,z) z;
    case 3
        Fx = @(x,y,z) cos(x+y);
        Fy = @(x,y,z) sin(x.*y).*sqrt(z+1);
        Fz = @(x,y,z) exp(z+y.*x);
    case 4
        Fx = @(x,y,z) 1./(1+25*(x-1/4).^2);
        Fy = @(x,y,z) 1./(1+25*(y-1/4).^2);
        Fz = @(x,y,z) 1./(1+25*(z-1/4).^2);
    case 5
        Fx = @(x,y,z) 1./(1+100.*(y-1/4).^2).*1./(1+100.*(z-1/4).^2).*(200.*(x-1/4))./(1+100.*(x-1/4).^2).^2; 
        Fy = @(x,y,z) 1./(1+100.*(x-1/4).^2).*1./(1+100.*(z-1/4).^2).*(200.*(y-1/4))./(1+100.*(y-1/4).^2).^2; 
        Fz = @(x,y,z) 1./(1+100.*(y-1/4).^2).*1./(1+100.*(x-1/4).^2).*(200.*(z-1/4))./(1+100.*(z-1/4).^2).^2; 
end

function clear_figure(nfig)

h=figure(nfig);
f_nfig=ishandle(h)&&strcmp(get(h,'type'),'figure');
if f_nfig
    clf(nfig);
end