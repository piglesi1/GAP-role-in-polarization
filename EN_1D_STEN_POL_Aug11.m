clearvars
p.sim=3001
figure(1),clf
figure(2),clf

for pol = 1.25;
for gap = 1.5
    for s=1:10;
    clc
    p.sim = p.sim+1;
    p.gap = gap;
    p.pol = pol
    Rc   = 5.00;   % Radius of cell (um)
    p.dx = 0.25;   % um
    p.tf = 240;    % Length of simulation (s)
    p.dt = 0.001; % (s)
    p.Np = uint16(round(2*pi*Rc/p.dx));
    %% -----------------------------------------
    % STEN parameters (using FBR model)
    % dFdt = -(a1+a2*R)*F+a3/((a4*B)^2+1)+a5
    p.a1 = 1;
    p.a2 = 37;
    p.a3 = 12;
    p.a4 = 444;
    p.a5 = 5.4e-3;
    % dBdt =
    p.b1 = 8.20;
    p.b2 = 0.06;
    p.b3 = 8880;
    % dRdt =
    p.c1 = 0.15;  % changed
    p.c2 = 1.5;   % changed
    % Diffusion coefficients
    p.DF = 0.15;
    p.DB = 0.15;
    p.DR = 1.5;
    p.DZ = 0.0125;
    p.DW = 0.0125;
    %
    p.z1 = 0.05;
    p.z2 = 0.05;
    p.z  = 2;       %2
    p.w1 = 0.05;
    p.w2 = 0.05;
    p.w  = 1;       %5
    %
    p.noise      = 12;
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
    visualize(F,B,R,p)
    out.F = F; out.B = B; out.R = R;
    out.Z = Z; out.W = W;
    out.p = p;
    out.t = 0:p.tf;
    save(strcat('STEN_POL_1D_',num2str(p.sim)),'out');
    toc
        %%
    F=F';
    B=B';
    R=R';
    %%
    p.FBRavg=[mean(F(:)) mean(B(:)) mean(R(:)) mean(Z(:)) mean(W(:))];
    p.FBRmax=[max(F(:))  max(B(:))  max(R(:))  max(Z(:))  max(W(:))];
    p.FBRstd=[std(F(:))  std(B(:))  std(R(:))  std(Z(:))  std(W(:))];
    figure(3)
    histogram2(F(:),R(:),'DisplayStyle','tile','BinWidth',0.01)
        axis equal
        axis([0 0.25 0 0.25])
    figure(4)  
        F1=F(:,41:(p.tf-40)/2);
        F2=F(:,(p.tf-40)/2+1:p.tf);
        bw=imbinarize(F);
    bw1=imbinarize(F1);
     bw2=imbinarize(F2);
    bw1=bwfill(bw1,'holes');
    bw2=bwfill(bw2,'holes');
    stats1=regionprops(bw1,'Area','Circularity','Orientation');
    stats2=regionprops(bw2,'Area','Circularity','Orientation');
    newarea1=[stats1(:).Area];
    newarea2=[stats2(:).Area];
    idx1=find(newarea1>5);
    idx2=find(newarea2>5);
    [n1,n2]=size(F1);
    p.orient=[mean(abs([stats1(idx1).Orientation])) mean(abs([stats2(idx2).Orientation]))];
    p.circul=[mean([stats1(idx1).Circularity]) mean([stats2(idx2).Circularity])];
    p.binarea=[sum(newarea1) sum(newarea2)]; 
    p.numelements=[length(newarea1) length(newarea2)];
    p.meanarea = p.binarea./p.numelements;
    imshow(imoverlay(imadjust(F),edge(bw)))
    drawnow 
end
end
end


%% postprocessing
%% Visualizations
function visualize(F,B,R,p)
figure(1),hold
    plot(F(:,10:100),R(:,10:100))
hold off
F=out.F(41:end,:);
B=out.B(41:end,:);
R=out.R(41:end,:);
[nframes,npts]=size(F);
FB = zeros(nframes,npts,3);
FB(:,:,1)=(F-min(F(:)))/(max(F(:))-min(F(:)));
FB(:,:,2)=(B-min(B(:)))/(max(B(:))-min(B(:)));
FR = zeros(nframes,npts,3);
FR(:,:,1)=(F-min(F(:)))/(max(F(:))-min(F(:)));
FR(:,:,3)=(R-min(R(:)))/(max(R(:))-min(R(:)));
BR = zeros(nframes,npts,3);
BR(:,:,2)=(B-min(B(:)))/(max(B(:))-min(B(:)));
BR(:,:,3)=(R-min(R(:)))/(max(R(:))-min(R(:)));
figure(2);
FBR = zeros(nframes,npts,3);
FBR(:,:,1)=(F-min(F(:)))/(max(F(:))-min(F(:)));
FBR(:,:,2)=(B-min(B(:)))/(max(B(:))-min(B(:)));
FBR(:,:,3)=(R-min(R(:)))/(max(R(:))-min(R(:)));
imshow(pagetranspose(FBR))
text(10,10,'Ras' ,'Color','red','FontSize',16,'FontName','Arial')
text(10,20,'PIP2','Color','green','FontSize',16,'FontName','Arial')
text(10,30,'PKB','Color','cyan','FontSize',16,'FontName','Arial')
%%
figure(2)
    tl=tiledlayout(3,1,'TileSpacing','tight');       
    nexttile(1)
    imshow(pagetranspose(FB))
    text(10,10,'Ras' ,'Color','red','FontSize',16,'FontName','Arial')
    text(10,20,'PIP2','Color','green','FontSize',16,'FontName','Arial')
    tt='pol: %4.1f: gap: %4.1f: p.z: %4.1f';
    title(sprintf(tt,p.pol,p.gap,p.z),'FontName','Arial','FontSize',16);
    nexttile(2)
    imshow(pagetranspose(FR))
    text(10,10,'Ras','Color','red','FontSize',16,'FontName','Arial')
    text(10,20,'PKBA','Color','cyan','FontSize',16,'FontName','Arial')
    nexttile(3)
    imshow(pagetranspose(BR))
    text(10,10,'PIP2','Color','green','FontSize',16,'FontName','Arial')
    text(10,20,'PKBA','Color','cyan','FontSize',16,'FontName','Arial')
    tl.TileSpacing = 'none';
    tl.Padding = 'tight';
    drawnow
    print('-dpng',strcat('STEN_POL_1D_',num2str(p.sim)))

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
    x = [0.0001:0.0001:0.02 0.021:0.001:0.1];
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
function [out1,out2,out3] = STEN_sde(t,x,flag,p)

% System equation for the excitable network: SDE: use central differencing for diffusion
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
    %%
    switch upper(p.sdetype)
        case 'ITO'
            %------------------------------------------------------------------
            % Excitable system
            %------------------------------------------------------------------
            pol = p.pol;
            if t<(p.tf-40)/2+40
                gap = 1;
            else
                gap = p.gap;
            end
            %
            J1 = (gap*p.a1 + gap*p.a2*R+ pol*p.w*W).*F;
            J2 = (p.a3./((p.a4^2*B.^2 + 1)) + p.a5).*(1 + pol*p.z*Z);  
            %
            J4 = p.b1;
            J3 = (p.b2 + p.b3*F).*B;
            %
            J5 = p.c1*R;
            J6 = p.c2*F;
            %
            driftF  = -J1 + J2 + p.DF./(p.dx^2)*lap(F);
            driftB  = -J3 + J4 + p.DB./(p.dx^2)*lap(B);
            driftR  = -J5 + J6 + p.DR./(p.dx^2)*lap(R);            
            %
            J7      = p.z1*Z;
            J8      = p.z2*R;
            driftZ  = -J7 + J8 + p.DZ./(p.dx^2)*lap(Z);
            %
            J9      = p.w1*W;
            J10     = p.w2*B;
            driftW  = -J9 + J10 + p.DW./(p.dx^2)*lap(W);

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
function lapX = lap(X)
    lapX = circshift(X,1)+circshift(X,-1)-2*X;
end

