clearvars

for gap = 1.5
    p.sim=0;
    p.gap = gap;
    Rc    = 40;      % Length of cell (um)
    p.dx  = 0.40;    % um
    p.tf  = 340;     % Length of simulation (s)
    p.dt  = 0.01;    % (s)
    p.Np  = round(Rc/p.dx);
    p.Np  = (p.Np)^2;
    %% -----------------------------------------
    % STEN parameters (using FBR model)
    % dFdt = -(a1+a2*R)*F+a3/((a4*B)^2+1)+a5
    p.a1     = 1;
    p.a2     = 37;
    p.a3     = 12;
    p.a4     = 444;
    p.a5     = 5.4e-3;
    % dBdt =
    p.b1     = 8.20;
    p.b2     = 0.06;
    p.b3     = 8880;
    % dRdt =
    p.c1     = 0.15;       %changed
    p.c2     = 0.75;       %changed
    % Diffusion coefficients
    p.DF     = 0.15;
    p.DB     = 0.15;
    p.DR     = 0.25;
    p.DZ     = 0.0125;
    p.DW     = 0.0125;
    %
    p.z1 = 0.05;
    p.z2 = 0.05;
    p.z  = 2;       %2
    p.w1 = 0.05;
    p.w2 = 0.05;
    p.w  = 1;       %5
    %
    p.noise      = 24;
    p.numsim     = uint16(p.Np);
    p.sdetype    = 'Ito';
    p.numdepvars = 5;
    % specify timeRange for the solver
    p=nullclines(p);
    F  = zeros(p.tf+1,p.Np); F(1)=p.F0;
    B  = zeros(p.tf+1,p.Np); B(1)=p.B0;
    R  = zeros(p.tf+1,p.Np); R(1)=p.R0;
    Z  = zeros(p.tf+1,p.Np); Z(1)=p.Z0;
    W  = zeros(p.tf+1,p.Np); W(1)=p.W0;
    Ts = (1/p.dt);
    %
    for k=4:5
        p.sim = p.sim+1;
        p.t0i = 0;
        p.tfi = p.tf/10;
        tic
        for ell=1:10
            p.t = p.t0i:p.dt:p.tfi;
            if ell==1
                output = SDE_euler_pai2D(p);
            else
                output = SDE_euler_pai2D(p,output(end,:));
            end
            F(p.t0i+2:p.tfi+1,:) = output(Ts+1:Ts:end,1:p.numdepvars:end);
            B(p.t0i+2:p.tfi+1,:) = output(Ts+1:Ts:end,2:p.numdepvars:end);
            R(p.t0i+2:p.tfi+1,:) = output(Ts+1:Ts:end,3:p.numdepvars:end);
            Z(p.t0i+2:p.tfi+1,:) = output(Ts+1:Ts:end,4:p.numdepvars:end);
            W(p.t0i+2:p.tfi+1,:) = output(Ts+1:Ts:end,5:p.numdepvars:end);
            p.t0i=p.tfi;
            p.tfi=p.tfi+p.tf/10;
        end
        %%
        for tt=1:p.tf+1
            F(tt,:) = noflux(F(tt,:),p.Np);
            B(tt,:) = noflux(B(tt,:),p.Np);
            R(tt,:) = noflux(R(tt,:),p.Np);
            Z(tt,:) = noflux(Z(tt,:),p.Np);
            W(tt,:) = noflux(W(tt,:),p.Np);
        end
        visualize(F,B,R,p)
        out.F = F; out.B = B; out.R = R;
        out.Z = Z; out.W = W;
        out.p = p;
        out.t = 0:p.tf;
        out.TotalF = sum(F(40:end,:));
        out.TotalB = sum(B(40:end,:));
        out.TotalR = sum(R(40:end,:))
        fname = strcat('STEN_POL_2D_',num2str(gap),'_',num2str(p.sim),'.mat');
        save(fname,'out');
        toc
    end
end


%% postprocessing

%% Visualizations
function visualize(F,B,R,p)
figure(1)
hold on,
    plot(F(:,10:20),R(:,10:20))
 hold off

[nframes,npts]=size(F);
npts = sqrt(npts);
FBR = zeros(npts,npts,3,'uint16');
figure(2);
fname ='STEN_RasGAP';

v = VideoWriter(fname,'MPEG-4');
open(v)
q=2^16-1;
for t=1:603 % nframes
    FBR(:,:,1)=q*(reshape(F(t,:),npts,npts)-min(F(:)))/(max(F(:))-min(F(:)));
    FBR(:,:,2)=q*(reshape(B(t,:),npts,npts)-min(B(:)))/(max(B(:))-min(B(:)));
    FBR(:,:,3)=q*(reshape(R(t,:),npts,npts)-min(R(:)))/(max(R(:))-min(R(:)));
    imshow(imresize(FBR,5))
    text(10,20,'Ras' ,      'Color','red',  'FontSize',16,'FontName','Arial')
    text(10,40,'PIP2',      'Color','green','FontSize',16,'FontName','Arial')
    text(10,60,'PKBA',      'Color','cyan','FontSize',16,'FontName','Arial')
    if t<=202
        text(10,80,'-20% RasGAP',      'Color','white','FontSize',16,'FontName','Arial')
        text(460,20,num2str(t),'Color','white','FontSize',16,'FontName','Arial')
    elseif t<403
        text(10,80,'WT RasGAP',      'Color','white','FontSize',16,'FontName','Arial')
         text(460,20,num2str(t-201),'Color','white','FontSize',16,'FontName','Arial')
    else
        text(10,80,'+10% RasGAP',      'Color','white','FontSize',16,'FontName','Arial')
        text(460,20,num2str(t-402),'Color','white','FontSize',16,'FontName','Arial')
    end
   drawnow
   frame = getframe(gcf);
   writeVideo(v,frame);
end
close(v);
end %function visualize(F,B,R,p)
%% Nullclines
function newp=nullclines(p)
syms F B R Z W;
newp = p;
B  = p.b1/(p.b2 + p.b3*F);
R  = p.c2*F/p.c1;
Z  = p.z2*R/p.z1;
W  = (p.w2*(p.b1/(p.b2 + p.b3*F))/p.w1);

J1 = (p.a1 + p.a2*R  + p.w*W)*F/(1+p.z*Z);
J2 = p.a3./((p.a4^2*B.^2 + 1)) + p.a5;  
eqn = -J1+J2==0;
S0 = vpasolve(eqn,F);
S0 = double(S0);
S0 = max(S0(S0==real(S0)));
newp.F0 = S0(S0>0);
newp.B0 = p.b1/(p.b2 + p.b3*newp.F0);
newp.R0 = p.c2*newp.F0/p.c1;
newp.Z0 = p.z2*newp.R0/p.z1;
newp.W0 = p.w2*newp.B0/p.w1;

figure(1)
    clf
    x = [0.0001:0.0001:0.02 0.021:0.1];
    r = p.c2*x/p.c1;
    b = p.b1./(p.b2 + p.b3*x);
    z = p.z2*r/p.z1;
    w = (p.w2.*p.b1./(p.b2 + p.b3*x))/p.w1;
    f = (p.a3./((p.a4^2*b.^2 + 1)) + p.a5).*(1+p.z*z)./x;
    f = (f-(p.a1 + p.w*w))/p.a2;
    plot(x,f,x,r,newp.F0,newp.R0,'xr')
    axis([0 0.1 0 0.1])
    drawnow
hold off
end % funnction p=nullclines(p)
%%
function xhat = SDE_euler_pai2D(p,state)

% revised 1/22/10 by Yuan:
%   1)  change SEED from 0 to sum(100*clock) as suggested, to get different
%   sequences of pseudopod-random numbers at each run
%   2)  add "param" and "input"
% revised 1/22/10 by Debojyoti:
%   1) Negative values of entities are neglected
% revised 07/27/2023 by PAI
%   1) only one input: p - should have all parameters

% Fixed stepsize Euler-Maruyama scheme for numerical solution of Ito SDE systems.
%
% usage: xhat = SDE_euler(p)
%
% IN:     p; complete vector of structural model parameters
%         p.PROBLEM; the user defined name of the current problem/experiment/example etc. (e.g. 'mySDE')
%         p.t; vector containing the equispaced simulation times sorted in ascending order.
%              It has starting simulation-time in first and ending simulation-time in last position.
%              Thus OWNTIME(i) - OWNTIME(i-1) = h, where h is the fixed stepsize
%              for the numerical intregration (i=2,3,...)
%         p.NUMDEPVARS; the number of dependent variables, i.e. the SDE dimension
%         p.NUMSIM; the number of desired simulations for the SDE numerical integration
%         p.SDETYPE; the SDE definition: must be 'Ito'
% OUTPUT: xhat; the array of the SDE approximated solution at times OWNTIME
%
% REFERENCE: [1] Kloeden, Platen and Schurz "Numerical solution of SDE through computer experiments", Springer-Verlag 2nd edition 1997
% Copyright (C) 2007, Umberto Picchini
% umberto.picchini@biomatematica.it
% http://www.biomatematica.it/Pages/Picchini.html

if(~isinteger(p.numsim) || isempty(p.numsim) || p.numsim <= 0)
    error('The number of trajectories NUMSIM must be a positive integer');
end
if(~isinteger(p.numsim) || isempty(p.numdepvars) || p.numdepvars <= 0)
    error('The number of variables NUMDEPVARS must be a positive integer');
end

if ~strcmpi(p.sdetype,'ITO')
    error('The Euler-Maruyama scheme is defined only for Ito SDE');
end

N = length(p.t);
if p.t(1)==0
    [t,xstart] = STEN_sde(p.t(1),[],'init',p);
    XVARS      = zeros(N,p.numsim*p.numdepvars);  % the predictions matrix
    xstart     = xstart';
    XVARS(1,:) = xstart((1:size(xstart,1))'*ones(1,p.numsim), :)';
else
    t          = p.t(1);
    XVARS      = zeros(N,p.numsim*p.numdepvars);  % the predictions matrix
    XVARS(1,:) = state;
end

% ugly but faster than XVARS(1,:) = repmat(xstart,1,NUMSIM);

% % Control the generation of the pseudo-random standard gaussian draws to
% % get repeatable results
% rng(0)            % PAI 07/29/2023
rng('shuffle');     % PAI 07/29/2023

for j=2:N
    % t is inherited as the starting time for this interval
    x = XVARS(j-1, :);  % the value(s) of XVARS at the start of the interval
    h = p.t(j)- t;      % the delta time (end - start) -> fixed size of the step .

    % Wiener increments generation with 'antithetic variates' variance reduction method
    Winc = zeros(1,p.numsim*p.numdepvars);
    if(mod(p.numsim*p.numdepvars,2)==0)
        Winc(1:p.numdepvars*p.numsim/2) = sqrt(h)*randn(1,p.numdepvars*p.numsim/2); % the Wiener increment(s) dWj;
        Winc(p.numdepvars  *p.numsim/2+1: end) = -Winc( 1:p.numdepvars* p.numsim/2);% retrieve the other half of the increments using the 'antithetic variates' variance reduction method;
    else
        % adjustment when (number of simulations * number of variables) is odd
        Winc(1 : round(p.numsim*p.numdepvars/2)) = sqrt(h)*randn(1,round(p.numsim*p.numdepvars/2));  % the Wiener increment(s) dWj;
        Winc(round(p.numsim*p.numdepvars/2)+1:end) = -Winc(1 : round(p.numsim*p.numdepvars/2)-1);  % retrieve the other half of the increments using the 'antithetic variates' variance reduction method;
    end
    [f,g] = STEN_sde(t,x,[],p);

    XVARS(j , :) = x + f * h + g .* Winc ;  % the Euler-Maruyama scheme for Ito SDEs with DIAGONAL noise
    ind_neglect = XVARS(j , :)<0;
    XVARS(j,ind_neglect) = XVARS(j-1,ind_neglect);
    t = p.t(j);    % now both t and j refer to the end-of-interval
end
xhat=XVARS;
end %xhat = SDE_euler_pai1D(p)
%%
function [out1,out2,out3] = STEN_sde(t,x,flag,p)

%% System equation for the excitable network: SDE: use central differencing for diffusion
% revised 3/1/2010 from RDS_sdefile.m
% SDE model definition: drift, diffusion, derivatives and initial conditions.
%
% [out1,out2,out3] = M9_sdefile(t,x,flag,bigtheta,SDETYPE,NUMDEPVARS,NUMSIM)
%
% IN:   t; working value of independent variable (time)
%       x; working value of dependent variable
%       flag; a switch, with values 'init' or otherwise
%       bigtheta; complete structural parameter vector
%       SDETYPE; the SDE definition: can be 'Ito' or 'Strat' (Stratonovich)
%       NUMDEPVARS; the number of dependent variables, i.e. the SDE dimension
%       NUMSIM; the number of desired simulations for the SDE numerical integration
% OUT:  out1; in case of flag='init' is just the initial time, otherwise it is the (vector of) SDE drift(s)
%       out2; in case of flag='init' is the initial value of the dependent variables. Otherwise it is the SDE diffusion(s)
%       out3; in case of flag='init' it is nothing. Otherwise it is the SDEs partial derivative(s) of the diffusion term

% Copyright (C) 2007, Umberto Picchini
% umberto.picchini@biomatematica.it
% http://www.biomatematica.it/Pages/Picchini.html
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%

if isempty(flag)  % NOT 'INITIAL CONDITION'
    xsplitted=cell(p.numdepvars,1);
    for i=1:p.numdepvars
        xsplitted{i} = x(i:p.numdepvars:end);
    end
    %----------------------------------------------------------------
    %  DEFINE HERE THE SDE
    %----------------------------------------------------------------

    F  = xsplitted{1};
    B  = xsplitted{2};
    R  = xsplitted{3};
    Z  = xsplitted{4};
    W  = xsplitted{5};
    %
    F = noflux(F,p.Np);
    B = noflux(B,p.Np);
    R = noflux(R,p.Np);
    Z = noflux(Z,p.Np);
    W = noflux(W,p.Np);
    %%
    switch upper(p.sdetype)
        case 'ITO'
            %------------------------------------------------------------------
            % Excitable system
            %------------------------------------------------------------------
            pol = .5;
            gap = p.gap;
            
            J1 = (gap*p.a1 + gap*p.a2*R + pol*p.w*W).*F;
            J2 = (p.a3./((p.a4^2*B.^2 + 1)) + p.a5).*(1 + pol*p.z*Z);       
            %
            J4 = p.b1;
            J3 = (p.b2 + p.b3*F).*B;
            %
            J5 = p.c1*R;
            J6 = p.c2*F;
            %
            driftF  = -J1 + J2 + p.DF./(p.dx^2)*lap2(F,p.Np);
            driftB  = -J3 + J4 + p.DB./(p.dx^2)*lap2(B,p.Np);
            driftR  = -J5 + J6 + p.DR./(p.dx^2)*lap2(R,p.Np);            
            %
            J7      = p.z1*Z;
            J8      = p.z2*R;
            driftZ  = -J7 + J8 + p.DZ./(p.dx^2)*lap2(Z,p.Np);
            %
            J9      = p.w1*W;
            J10     = p.w2*B;
            driftW  = -J9 + J10 + p.DW./(p.dx^2)*lap2(W,p.Np);

            noiseF = sqrt(J1 + J2)/p.noise;
            noiseB = sqrt(J3 + J4)/p.noise;
            noiseR = sqrt(J5 + J6)/p.noise;
            noiseZ = sqrt(J7 + J8)/p.noise;
            noiseW = sqrt(J9 + J10)/p.noise;
            %----------------------------------------------------------------
    end
    out1 = zeros(1,p.numsim*p.numdepvars);
    out1(1:p.numdepvars:end) = driftF;
    out1(2:p.numdepvars:end) = driftB;
    out1(3:p.numdepvars:end) = driftR;
    out1(4:p.numdepvars:end) = driftZ;
    out1(5:p.numdepvars:end) = driftW;

    out2 =  zeros(1,p.numsim*p.numdepvars);
    out2(1:p.numdepvars:end) = noiseF;
    out2(2:p.numdepvars:end) = noiseB;
    out2(3:p.numdepvars:end) = noiseR;
    out2(4:p.numdepvars:end) = noiseZ;
    out2(5:p.numdepvars:end) = noiseW;

    out3 = zeros(1,p.numsim*p.numdepvars);
else
    switch(flag)
        case 'init'
            out1 = t;
            %----------------------------------------------------------------
            %  DEFINE INITIAL CONDITIONS
            %----------------------------------------------------------------
            out2 = [p.F0 p.B0 p.R0 p.Z0 p.W0]; %   p.E0 p.I0 p.RR0];   
            % write here the SDE initial condition(s)
            out3 = [];
        otherwise
            error(['Unknown flag ''' flag '''.']);
    end
end
end
%%

function lapX = lap2(X,Np)
    i=2:sqrt(Np)-1;
    j=i;
    X    = reshape(X,sqrt(Np),sqrt(Np));
    lapX = zeros(size(X));
    lapX(i,j) = X(i+1,j)+X(i-1,j)+X(i,j-1)+X(i,j+1)-4*X(i,j);
    lapX = reshape(lapX,1,Np);
end

function nfX = noflux(X,Np)
    Np          = sqrt(Np);
    nfX         = reshape(X,Np,Np); 
    nfX(1,:)    = nfX(2,:); 
    nfX(Np,:)   = nfX(Np-1,:);
    nfX(:,1)    = nfX(:,2); 
    nfX(:,Np)   = nfX(:,Np-1);
    nfX         = reshape(nfX,1,Np*Np);
end