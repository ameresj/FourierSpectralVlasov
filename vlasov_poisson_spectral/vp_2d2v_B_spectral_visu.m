%% Spectral 2D Magnetized Vlasov-Poisson Solver
% Author: *Jakob Ameres* <http://jakobameres.com jakobameres.com>
%% Features
%
% * Fourier spectral description of density f
% * *Symplectic* (for Vlasov-Poisson) Runge Kutta *time integrator*
%
clear all; close all;
normalize=@(x) x./norm(x);
maxNumCompThreads(4);




testcase='kelvinhelmholtz';

splitting_order='strang'; %fourth10lie
splitting='scovel'; % shear std imr scovel scovel_imr
B3=16; %32

fprintf('Running: %s with B3=%d \n',splitting,B3);

% External electric field
E1_ext=[]; %E1_ext=@(x1,x2,t) 0.*x1 +0.*x2+0.*t;
E2_ext=[]; %E2_ext=@(x1,x2,t) 0.*x1 +0.*x2 +0.*t;



nviz=50;
subrot=1;

% Default parameters
% B3=10;

vmax1=4.5; vmin1=-vmax1; vmax2=vmax1; vmin2=vmin1; q=-1; m=1;
switch(testcase)
    case 'landau'
        %% Diagonal Landau Damping parameters
        eps=0.05; % Amplitude of perturbation, 0.05 for linear, 0.5 for nonlinear
        k0=normalize([1,2]).*0.5;    % Wave vector
        L=2*pi./k0; % length of domain
        q=-1;    % negative charge to mass = electrons, $\frac{q}{m}=-1$
        m=1;
        % Initial condition
        f0=@(x1,x2,v1,v2) (1+eps*cos(k0(1).*x1+k0(2).*x2)).*...
            exp(-0.5.*(v1.^2+v2.^2))/(2*pi);
        %         B3=0;
        vmax1=4.5; vmin1=-vmax1; vmax2=vmax1; vmin2=vmin1;
    case 'lambdipole'
        L=[2,2]; lambda=3.83170597020751; % J_1(lambda)=0
        U=1;     R=0.2;
        
        rho0=@(x1,x2) 2*lambda*U*cos(atan2(x2,x1)).*...
            besselj(1,lambda*hypot(x1,x2) )./...
            besselj(0,lambda*R).*( (x1.^2+x2.^2)<=R.^2);
        
        f0=@(x1,x2,v1,v2) rho0(x1-1,x2-1).*...
            exp(-0.5.*(v1.^2+v2.^2))/(2*pi);
        
        
        vmax1=4.5; vmin1=-vmax1; vmax2=vmax1; vmin2=vmin1;
        B3=10;
    case 'gaussvortex'
        L=[2,2];
        sigmax=0.1; sigmay=0.2;
        f0=@(x1,x2,v1,v2) normpdf(x1,1,sigmax).*normpdf(x2,1,sigmay).*...
            exp(-0.5.*(v1.^2+v2.^2))/(2*pi);
        
        
        vmax1=4.5; vmin1=-vmax1; vmax2=vmax1; vmin2=vmin1;
        B3=20;
        
    case 'diocotron'
        L=[2*pi,10];
        q=-1;
        m=1;
        eps=0.05;
        l=5;
        % Initial condition
        f0=@(x1,x2,v1,v2) (1+ eps.*cos(l*x1)).*(x2>=4 & x2<=5).* ...
            exp(-0.5.*(v1.^2+v2.^2))/(2*pi);
        
        
    case 'kelvinhelmholtz'
        k0=0.4;
        L=[2*pi/k0,2*pi]; % length of domain
        q=-1;
        m=1;
        nu=0.015;
        % Initial condition
        f0=@(x1,x2,v1,v2) (1+ sin(x2)+nu.*cos(k0*x1) ).*...
            exp(-0.5.*(v1.^2+v2.^2))/(2*pi);
        tmax=40; %linear till 20
        %dt=0.02;
        %         B3=10; % avoid resonances
    case 'kelvinhelmholtz_slab'
        %         epsb=0.1; % scaling
        %         B3=1/epsb;
        
        k0=0.3;
        L=[2*pi/k0,2*pi]; % length of domain
        q=-1;
        m=1;
        % Initial condition
        f0=@(x1,x2,v1,v2) (1+ sin(x2)+0.015.*sin(x2/2).*cos(k0*x1) ).*...
            exp(-0.5.*(v1.^2+v2.^2))/(2*pi);
        %         f0=@(x1,x2,v1,v2) (1+ sin(x2)+0.015.*sin(x2/2).*cos(k0*x1) ).*...
        %                            exp(-0.5.*(v1.^2+v2.^2)/epsb.^2 )/(2*pi*epsb.^2);
        %                         vmax1=vmax1*epsb;
        %                         vmin1=vmin1*epsb;
        %                         vmax2=vmax2*epsb;
        %                         vmin2=vmin2*epsb;
        
    case 'kelvinhelmholtz2'
        k0=2*pi/40;
        L=[2*pi/k0,10]; % length of domain
        q=-1;
        m=1;
        f0=@(x1,x2,v1,v2) (1+ 1.5*sech( (x2-5)/0.9).*(1+0.08.*sin(2*k0*x1))   ).*...
            exp(-0.5.*(v1.^2+v2.^2))/(2*pi);
        
end


%%%%%%%%%
% nviz=20;
% tmax=25;
%%%%%%%%%

% B3=1; % 1 10 20 40

dt=0.01*B3^2;
tmax=100;

if (B3~=0)
    dt=dt/B3;
    % tmax=tmax*B3;
    % tmax=tmax;
    tmax=tmax*B3;
end




Nx1=32; % Number of cells in spatial direction
Nv1=32;% Number of cells in velocity direciton

Nx2=Nx1;  Nv2=Nv1;

%% Define phase space mesh
v1=linspace(vmin1,vmax1,Nv1+1).';v1=v1(1:end-1);
v2=linspace(vmin2,vmax2,Nv2+1);v2=v2(1:end-1);
x1=linspace(0,L(1),Nx1+1).';x1=x1(1:end-1);
x2=linspace(0,L(2),Nx2+1);x2=x2(1:end-1);

kx1=fftshift( -Nx1/2:Nx1/2-1)'*2*pi/L(1);
kx2=fftshift( -Nx2/2:Nx2/2-1)*2*pi/L(2);
kv1=fftshift( -Nv1/2:Nv1/2-1)*2*pi/(vmax1-vmin1);
kv2=fftshift( -Nv2/2:Nv2/2-1)*2*pi/(vmax2-vmin2);
kv1=reshape(kv1,1,1,Nv1,1);
kv2=reshape(kv2,1,1,1,Nv2);

dv1=(vmax1-vmin1)/Nv1; dv2=(vmax2-vmin2)/Nv2;
dx1=L(1)/Nx1; dx2=L(2)/Nx2;
dx=dx1*dx2; dv=dv1*dv2;


%% Splitting
switch(splitting_order)
    case 'strang'
        dtA=1/2;
        dtB=1/2;
    case '2nd4lie'
        % % second order (4 Lie)
        % % Second order splitting (smaller error constant than Strang)
        y2=(2*sqrt(326)-36)^(1/3);
        alphat=(y2^2+6*y2-2)/(12*y2);
        
        dtA=[alphat, 1/2-alphat];
        dtB=flip(dtA);
    case 'fourth10lie'
        dtB=[(146+5*sqrt(19))/540, (-2+10*sqrt(19))/135, 1/5, ...
            (-23-20*sqrt(19))/270, (14-sqrt(19))/108];
        dtA=flip(dtB);
    case 'order6'
        dtw=[0.1867, 0.55549702371247839916, 0.12946694891347535806,-0.84326562338773460855];
        dtw=[dtw,1-2*sum(dtw),flip(dtw)];
        dtA=dtw/2;dtB=flip(dtw)/2;
end
% scale time step
dtA=dtA*dt;dtB=dtB*dt;


% %Build Matrix for spatial advection and all stages

switch(splitting)
case {'scovel','scovel_imr'}
    XSHIFTA1=zeros(Nx1,1,Nv1,Nv2,length(dtA));
    XSHIFTA2=zeros(1,Nx2,Nv1,Nv2,length(dtA));
    XSHIFTB1=zeros(Nx1,1,Nv1,Nv2,length(dtA));
    XSHIFTB2=zeros(1,Nx2,Nv1,Nv2,length(dtA));
    v1A1=zeros(1,1,Nv1,Nv2,length(dtA));
    v2A1=zeros(1,1,Nv1,Nv2,length(dtA));
    for jdx=1:length(dtA)
        delta=-q/m*B3;
        theta=delta*(dtA(jdx));
        XSHIFTA1(:,:,:,:,jdx)=exp(-1j*( kx1.*reshape(...
            sin(theta)*v1+ ( cos(theta)-1).*v2,1,1,Nv1,Nv2))/delta);
        XSHIFTA2(:,:,:,:,jdx)=exp(-1j*( kx2.*reshape(...
            (1-cos(theta))*v1+ sin(theta).*v2,1,1,Nv1,Nv2))/delta);
        theta=delta*(dtB(jdx));
        XSHIFTB1(:,:,:,:,jdx)=exp(-1j*( kx1.*reshape(...
            sin(theta)*v1+ (1- cos(theta)).*v2,1,1,Nv1,Nv2))/delta);
        XSHIFTB2(:,:,:,:,jdx)=exp(-1j*( kx2.*reshape(...
            (cos(theta)-1)*v1+ sin(theta).*v2,1,1,Nv1,Nv2))/delta);
        
        
        v1A1(:,:,:,:,jdx)=reshape...
            (sin(theta)*v1+ ( cos(theta)-1).*v2,1,1,Nv1,Nv2)/delta;
        v2A1(:,:,:,:,jdx)=reshape(...
            (1-cos(theta))*v1+ sin(theta).*v2,1,1,Nv1,Nv2)/delta;
        
    end
    otherwise
        XSHIFT1=zeros(Nx1,Nv1,length(dtA));
        XSHIFT2=zeros(Nx2,Nv2,length(dtA));
        for jdx=1:length(dtA)
            XSHIFT1(:,:,jdx)=exp(-1j*dtA(jdx)*(kx1(:)*v1(:)'));
            XSHIFT2(:,:,jdx)=exp(-1j*dtA(jdx)*(kx2(:)*v2(:)'));
        end
end
%Preallocate matrix for velocity advection
% VSHIFT=zeros(Nx,Nv);

%diagnostic
Nt=ceil(tmax/dt);
Epot=zeros(Nt,2);
Ekin=zeros(Nt,2);
Momentum=zeros(Nt,2);
Energy=zeros(Nt,1);
PhiL2=zeros(Nt,1);

% Set initial condition
% [XX1,XX2,VV1,VV2]=ndgrid(x1,x2,v1,v2);
% f=f0(XX1,XX2,VV1,VV2);
[XX1,XX2]=ndgrid(x1,x2); % only create 2d grid in order to save memory
F=zeros(Nx1,Nx2,Nv1,Nv2);
for vdx2=1:Nv2
    for vdx1=1:Nv1
        F(:,:,vdx1,vdx2)=f0(XX1,XX2,v1(vdx1),v2(vdx2));
    end
end




switch(splitting)
    case 'std'
        BSHIFT1A=zeros(1,1,Nv1,Nv2,length(dtA)); BSHIFT1B=zeros(1,1,Nv1,Nv2,length(dtA));
        BSHIFT2=zeros(1,1,Nv1,Nv2,length(dtA));
        for sdx=1:length(dtA)
            BSHIFT1A(:,:,:,:,sdx)=exp(1j*dtA(sdx)*q*(-reshape(v2,1,1,1,Nv2).*B3).*kv1);
            BSHIFT1B(:,:,:,:,sdx)=exp(1j*dtB(sdx)*q*(-reshape(v2,1,1,1,Nv2).*B3).*kv1);
            BSHIFT2(:,:,:,:,sdx)=exp(1j*(dtA(sdx)+dtB(sdx))*q*(reshape(v1,1,1,Nv1,1).*B3).*kv2);
            
        end
    case {'scovel','shear'}
        SHEAR1=zeros(1,1,Nv1,Nv2,length(dtA));
        SHEAR2=zeros(1,1,Nv1,Nv2,length(dtA));
        theta_rot=zeros(length(dtA),1);
        for sdx=1:length(dtA)
            theta=q/m*B3*(dtA(sdx)+dtB(sdx));
            theta_rot(sdx)=-floor((theta+pi/4)/(pi/2));
            theta=(mod(theta+pi/4,pi/2)-pi/4)/subrot;
            SHEAR1(:,:,:,:,sdx)=exp(1j*(-tan(theta/2))*reshape(v2,1,1,1,Nv2).*kv1);
            SHEAR2(:,:,:,:,sdx)=exp(1j*(sin(theta))*reshape(v1,1,1,Nv1,1).*kv2);
        end
        
    case {'scovel_imr'}
        theta=q/m*B3*(dtA+dtB);
        theta_rot=-floor((theta+pi/4)/(pi/2));
        
end



%% Fourier transform in x
tic;
switch(splitting)
    case {'scovel','scovel_imr'}
        %% Initial poisson solve
        rho=q*sum(sum(F,4),3)*dv1*dv2; %rho(x) -> integrate over v
        rho=fft(fft(rho,[],1),[],2);
        rho(1,1)=0; %remove (ion) background
        Phi=-rho./(  kx1.^2 +kx2.^2); Phi(1,1)=0;
        E1=-Phi.*(1j.*kx1);E2=-Phi.*(1j.*kx2);
        Phi=ifft(ifft(Phi,[],2),[],1,'symmetric');
        E1=ifft(ifft(E1,[],2),[],1,'symmetric');
        E2=ifft(ifft(E2,[],2),[],1,'symmetric');
        
        % Initial velocity shift
        F=fft(fft(F,[],3),[],4);
        VSHIFT=exp(1j*(dtA(1))*q*(E1.*kv1+E2.*kv2));
        F=F.*VSHIFT;
    otherwise
        F=fft(fft(F,[],1),[],2);
end


%% Visualization/Output
figrho=figure('Name', 'charge density', 'Numbertitle','off');
figv=figure('Name', 'momentum density', 'Numbertitle','off');
figx1v1=figure('Name', 'phase space (x_1,v_1)', 'Numbertitle','off');
figE=figure('Name', 'electric field', 'Numbertitle','off');

mkdir('./plots2');
% prefix=sprintf('./plots/VP2d2v_%s_%s_%s_B%d_%d_%d_', ...
%     testcase,splitting_order,splitting,B3,Nx1,Nv1 );
prefix=sprintf('./plots2/VP2d2v_%s_%s_%s_B%d_%d_%d_', ...
    testcase,splitting_order,splitting,B3,Nx1,Nv1 );
mkdir( sprintf('%splots', prefix) ); % Plots for distributions function


for tdx=1:Nt
    t=(tdx-1)*dt;
    
    if (tdx== floor(Nt*0.04) )
        eta=(Nt/tdx-1)*toc();
        fprintf('ETA: %dmin %ds -> %s \n',floor(eta/60), floor(mod(eta,60)),...
            datestr(now()+eta/(24*3600)));
    end
    
    
    %      Ekin(tdx,1)=sum(reshape(sum(fdx,4),Nv1,1).*v1(:).^2/2)*dx*dv*m;
    %             Ekin(tdx,2)=sum(reshape(sum(fdx,3),Nv2,1).*v2(:).^2/2)*dx*dv*m;
    %             Momentum(tdx,1)=sum(reshape(sum(fdx,4),Nv1,1).*v1(:))*dx*dv*m;
    %             Momentum(tdx,2)=sum(reshape(sum(fdx,3),Nv2,1).*v2(:))*dx*dv*m;
    
    
    switch(splitting)
        case {'scovel','scovel_imr'}
            for sdx=1:length(dtA)
                
                F=ifft(ifft(F,[],4),[],3,'symmetric'); % fully backtransformed
                
                %% Symmetric composed scovel with its adjoint
                % rotation in x
                F=fft(fft(F,[],1),[],2);
                verify_ampere=false;
                if (verify_ampere)
                
                    %solve Poisson
                    rho=q*sum(sum(F,4),3)*dv1*dv2; 
                    rho(1,1)=0; %remove (ion) background
                    Phi=-rho./(  kx1.^2 +kx2.^2); Phi(1,1)=0;      
                    E1=-Phi.*(1j.*kx1);E2_=-Phi.*(1j.*kx2);
        
        
%                     dE1=sum(sum( v1A1(:,:,:,:,sdx)./...
%                         (v1A1(:,:,:,:,sdx).*kx1+v2A1(:,:,:,:,sdx).*kx2)...
%                         .*F.*(XSHIFTA1(:,:,:,:,sdx) -1),4),3)*dv1*dv2;
%                        dE1=sum(sum(F.*(XSHIFTA1(:,:,:,:,sdx) -1)./kx1,4),3)*dv1*dv2;
                   
tmpv=reshape(v1,1,1,Nv1,1)./(-1j*(v1A1(:,:,:,:,sdx).*kx1+v2A1(:,:,:,:,sdx).*kx2));
tmpv(:,:,Nv1/2+1,Nv2/2+1)=1;


% squeeze(tmpv)
                       dE1=sum(sum(...
                           tmpv.*...
                           F.*(XSHIFTA1(:,:,:,:,sdx) -1),4),3)*dv1*dv2/delta^2*4;
                       
                       
%                     dE1=sum(sum( reshape(./...
%                         (v1A1(:,:,:,:,sdx).*kx1+v2A1(:,:,:,:,sdx).*kx2)...
%                         .*F.*(XSHIFTA1(:,:,:,:,sdx) -1),4),3)*dv1*dv2;
                    dE1(isnan(dE1))=0;
                    
                     dE1(:,1)=0;
                    
%                     stop
%                     0.6224*0.5193/delta
                    
%                     dv1*dv2
                
                    
                    F_=F.*XSHIFTA1(:,:,:,:,sdx);    
                     F_=F_.*XSHIFTA2(:,:,:,:,sdx);
        
                    %solve Poisson
                    rho=q*sum(sum(F_,4),3)*dv1*dv2; 
                    rho(1,1)=0; %remove (ion) background
                    Phi=-rho./(  kx1.^2 +kx2.^2); Phi(1,1)=0;      
        
                    E1_=-Phi.*(1j.*kx1);E2_=-Phi.*(1j.*kx2);
                    dE1_=(E1-E1_);
                    
                    E1=E1+dE1;
                    
                    
pcolor(abs((E1-E1_)))
                    
                    dE1_(2,2)
                    dE1(2,2)

                    
                    
                    close all;
                    
                    figure()
                    pcolor(abs(dE1_)); colorbar;
                    
                    figure()
                    pcolor(abs(dE1)); colorbar;
                    
                    
                    figure()
                     pcolor(log(abs(dE1-dE1_)));colorbar
%                     stop
                else
                F=F.*XSHIFTA1(:,:,:,:,sdx);
                F=F.*XSHIFTA2(:,:,:,:,sdx);
                    
                    
                end
                
                 
                
                
                
        
                
                
                % v x B (in the middle)
                F=ifft(ifft(F,[],2),[],1,'symmetric'); % full backtransform
                
                switch( mod(theta_rot(sdx),4) )
                    case 1
                        F=flip( permute(F,[1,2,4,3]),1+2);
                    case 2 % 180 degree
                        F=flip(flip(F,2+2),1+2);
                    case 3 % 270 degree
                        F=flip( permute(F,[1,2,4,3]),2+2);
                end
                
                
                switch(splitting)
                    case 'scovel'
                        
                        F=fft(F,[],3);
                        for rotdx=1:subrot
                            F=F.*SHEAR1(:,:,:,:,sdx);
                            F=fft(ifft(F,[],3,'symmetric'),[],4);
                            F=F.*SHEAR2(:,:,:,:,sdx);
                            F=fft(ifft(F,[],4,'symmetric'),[],3);
                            F=F.*SHEAR1(:,:,:,:,sdx);
                        end
                        F=ifft(F,[],3,'symmetric');
                        
                    case 'scovel_imr'
                        
                        F=permute(F,[3,4,1,2]);
                        
                        theta=q/m*B3*(dtA(sdx)+dtB(sdx));
                        theta=(mod(theta+pi/4,pi/2)-pi/4);
                        
                        F=imrotate(F,theta*180/pi,'bicubic','crop');
                        F=permute(F,[3,4,1,2]);
                end
                
                % adjoint rotation in x
                F=fft(fft(F,[],1),[],2);
                F=F.*XSHIFTB1(:,:,:,:,sdx);
                F=F.*XSHIFTB2(:,:,:,:,sdx);
                % -- end of scovel ---
                
                
                %% Poisson equation
                rho=q*sum(sum(F,4),3)*dv1*dv2; %rho(x) -> integrate over v
                rho(1,1)=0; %remove (ion) background
                Phi=-rho./(  kx1.^2 +kx2.^2); Phi(1,1)=0;
                E1=-Phi.*(1j.*kx1);E2=-Phi.*(1j.*kx2);
                if (sdx==1)
                    Epot(tdx,1) =(E1(:)'*E1(:))/2*dx/(Nx1*Nx2);
                    Epot(tdx,2) =(E2(:)'*E2(:))/2*dx/(Nx1*Nx2);%
                    PhiL2(tdx) =sqrt(( real(Phi(:)'*Phi(:))*dx/(Nx1*Nx2)));
                end
                E1=ifft(ifft(E1,[],2),[],1,'symmetric');
                E2=ifft(ifft(E2,[],2),[],1,'symmetric');
                F=ifft(ifft(F,[],2),[],1,'symmetric');
                
                if (sdx==1)
                    fdx=real(sum(sum(F,1),2));
                    Ekin(tdx,1)=sum(reshape(sum(fdx,4),Nv1,1).*v1(:).^2/2)*dx*dv*m;
                    Ekin(tdx,2)=sum(reshape(sum(fdx,3),Nv2,1).*v2(:).^2/2)*dx*dv*m;
                    
                    Momentum(tdx,1)=sum(reshape(sum(fdx,4),Nv1,1).*v1(:))*dx*dv*m+...
                        + sum( E2(:).*B3(:) )*dx;
                    Momentum(tdx,2)=sum(reshape(sum(fdx,3),Nv2,1).*v2(:))*dx*dv*m+...
                        + sum( -E1(:).*B3(:) )*dx;
                    Energy(tdx)=sum(Ekin(tdx,:))+sum(Epot(tdx,:));
                end
                
                % velocity shift
                F=fft(fft(F,[],3),[],4);
                if (sdx==length(dtA) && tdx==Nt)
                    tau=dtA(end); % final time step
                else
                    % wrap around
                    tau=dtB(sdx)    + dtA( mod(sdx,length(dtA))+1 );
                end
                
                F=F.*exp(1j*(tau)*q*(E1.*kv1+E2.*kv2));
            end
            
        otherwise
            for sdx=1:length(dtA)
                
                %%  H_p
                %f=fft(fft(f,[],1),[],2);
                F=bsxfun(@times,F,reshape(XSHIFT1(:,:,sdx),Nx1,1,Nv1,1));
                F=bsxfun(@times,F,reshape(XSHIFT2(:,:,sdx),1,Nx2,1,Nv2));
                
                
                
                
                
%                 stop
                
                
                
                %% Poisson equation
                rho=q*sum(sum(F,4),3)*dv1*dv2; %rho(x) -> integrate over v
                rho(1,1)=0; %remove (ion) background
                Phi=-rho./(  kx1.^2 +kx2.^2); Phi(1,1)=0;
                E1=-Phi.*(1j.*kx1);E2=-Phi.*(1j.*kx2);
                if (sdx==1)
                    Epot(tdx,1) =(E1(:)'*E1(:))/2*dx/(Nx1*Nx2);
                    Epot(tdx,2) =(E2(:)'*E2(:))/2*dx/(Nx1*Nx2);%
                    PhiL2(tdx) =sqrt(( real(Phi(:)'*Phi(:))*dx/(Nx1*Nx2)));
                end
                E1=ifft(ifft(E1,[],2),[],1,'symmetric');
                E2=ifft(ifft(E2,[],2),[],1,'symmetric');
                
                %         % add external fields
                %         if (~isempty(E1_ext)), E1=E1+E1_ext(XX1,XX2,t); end
                %         if (~isempty(E2_ext)), E2=E2+E2_ext(XX1,XX2,t); end
                %         t=t+dt*rksc(sdx);
                F=ifft(ifft(F,[],2),[],1);
                
                if (sdx==1)
                    fdx=real(sum(sum(F,1),2))*dx;
                    Ekin(tdx,1)=sum(reshape(sum(fdx,4),Nv1,1).*v1(:).^2/2)*dv*m;
                    Ekin(tdx,2)=sum(reshape(sum(fdx,3),Nv2,1).*v2(:).^2/2)*dv*m;
                    
                    Momentum(tdx,1)=sum(reshape(sum(fdx,4),Nv1,1).*v1(:))*dv*m+...
                        + sum( E2(:).*B3(:) )*dx;
                    Momentum(tdx,2)=sum(reshape(sum(fdx,3),Nv2,1).*v2(:))*dv*m+...
                        + sum( -E1(:).*B3(:) )*dx;
                    Energy(tdx)=sum(Ekin(tdx,:))+sum(Epot(tdx,:));
                end
                
                switch(splitting)
                    case {'shear','rotation','imr'}
                        
                        
                        %% H_E
                        % dV/dt = E(x)
                        F=fft(fft(F,[],3),[],4);
                        VSHIFT=exp(1j*dtA(sdx)*q*(E1.*kv1+E2.*kv2));
                        F=F.*VSHIFT;
                        %F=ifft(ifft(F,[],4),[],3);
                        
                        % Rotate in Fourier or physical space
                        % dV/dt = V x B(x)
                        
                        switch(splitting)
                            case 'imr'
                                theta=q/m*B3*(dtA(sdx)+dtB(sdx));
                                
                                F=permute(F,[3,4,1,2]);
                                F=ifft(ifft(F,[],2),[],1);
                                F=imrotate(F,theta*180/pi,'bicubic','crop');
                                F=fft(fft(F,[],1),[],2);
                                F=permute(F,[3,4,1,2]);
                            case 'shear'
                                %F=ifft(ifft(F,[],4),[],3);
                                %discrete rotation in real or Fourier space
                                switch( mod(theta_rot(sdx),4) )
                                    case 1 % 90 degree clockwise
                                        F=flip( permute(F,[1,2,2+2,1+2]),1+2);
                                    case 2 % 180 degree clockwise
                                        F=flip(flip(F,2+2),1+2);
                                    case 3 % 270 degree clockwise
                                        F=flip( permute(F,[1,2,2+2,1+2]),2+2);
                                end
                                %F=fft(F,[],3);
                                F=ifft(F,[],4);
                                for rotdx=1:subrot
                                    F=F.*SHEAR1(:,:,:,:,sdx);
                                    F=ifft(fft(F,[],4),[],3);
                                    F=F.*SHEAR2(:,:,:,:,sdx);
                                    F=ifft(fft(F,[],3),[],4);
                                    F=F.*SHEAR1(:,:,:,:,sdx);
                                    % % shear
                                    % theta=q/m*B3*(dtA(sdx)+dtB(sdx))/subrot;
                                    % F=fft(F,[],1);
                                    % F=F.*exp(1j*(-tan(theta/2))*x2.*kx1);
                                    % F=ifft(F,[],1);
                                    %
                                    % % Shear
                                    % F=fft(F,[],2);
                                    % F=F.*exp(1j*(sin(theta))*x1.*kx2);
                                    % F=ifft(F,[],2);
                                    %
                                    % % shear
                                    % F=fft(F,[],1);
                                    % F=F.*exp(1j*(-tan(theta/2))*x2.*kx1);
                                    % F=ifft(F,[],1);
                                    %  F=permute(F,[3,4,1,2]);
                                    %
                                    
                                end
                                F=fft(F,[],4);
                        end
                        %
                        %
                        
                        
                        % dV/dt = E(x)
                        %F=fft(fft(F,[],3),[],4);
                        if (dtB(sdx)==dtA(sdx))
                            F=F.*VSHIFT;
                        else
                            F=F.*exp(1j*dtB(sdx)*q*(E1.*kv1+E2.*kv2));
                        end
                        F=ifft(ifft(F,[],4),[],3);
                        
                    case 'std'
                        F=fft(F,[],3);
                        %VSHIFT=exp(1j*dtA(sdx)*q*(E1- reshape(v2,1,1,1,Nv2).*B3).*kv1);
                        VSHIFT=exp(1j*dtA(sdx)*q*(E1.*kv1));
                        F=F.*VSHIFT;
                        F=F.*BSHIFT1A(:,:,:,:,sdx);
                        F=ifft(F,[],3);
                        
                        F=fft(F,[],4);
                        %F=F.*exp(1j*(dtA(sdx)+dtB(sdx))*q*(E2+ reshape(v1,1,1,Nv1,1).*B3).*kv2);
                        F=F.*exp(1j*(dtA(sdx)+dtB(sdx))*q*(E2).*kv2);
                        F=F.*BSHIFT2;
                        F=ifft(F,[],4);
                        
                        F=fft(F,[],3);
                        if (dtB(sdx)==dtA(sdx))
                            F=F.*VSHIFT;
                        else
                            %F=F.*exp(1j*dtB(sdx)*q*(E1- reshape(v2,1,1,1,Nv2).*B3).*kv1);
                            F=F.*exp(1j*dtB(sdx)*q*(E1).*kv1);
                        end
                        F=F.*BSHIFT1B(:,:,:,:,sdx);
                        
                        F=ifft(F,[],3);
                end
                
                
                
                F=fft(fft(F,[],1),[],2);
                F=bsxfun(@times,F,reshape(XSHIFT1(:,:,end-(sdx-1)),Nx1,1,Nv1,1));
                F=bsxfun(@times,F,reshape(XSHIFT2(:,:,end-(sdx-1)),1,Nx2,1,Nv2));
                %f=ifft(ifft(f,[],2),[],1);
            end
    end
    
    
    %         %% Advection in v
    %         f=ifft(ifft(f,[],2),[],1); %back transform
    %
    %         % f is fully backtransformed, calculate diagnostics
    %         if (sdx==1)
    %             fdx=real(sum(sum(f,1),2));
    %             Ekin(tdx,1)=sum(reshape(sum(fdx,4),Nv1,1).*v1(:).^2/2)*dx*dv*m;
    %             Ekin(tdx,2)=sum(reshape(sum(fdx,3),Nv2,1).*v2(:).^2/2)*dx*dv*m;
    %             Momentum(tdx,1)=sum(reshape(sum(fdx,4),Nv1,1).*v1(:))*dx*dv*m;
    %             Momentum(tdx,2)=sum(reshape(sum(fdx,3),Nv2,1).*v2(:))*dx*dv*m;
    %             Energy(tdx)=sum(Ekin(tdx,:))+sum(Epot(tdx,:));
    %         end
    %
    %
    %
    %         %%  Advection in x
    %         f=fft(fft(f,[],1),[],2);
    %         %f=f.*exp(-1j*dt*rksd(sdx)*(kx1.*reshape(v1,1,1,Nv1,1)))...
    %         %    .*exp(-1j*dt*rksd(sdx)*(kx2.*reshape(v2,1,1,1,Nv2)));
    %
    %         f=bsxfun(@times,f,reshape(XSHIFT1(:,:,sdx),Nx1,1,Nv1,1));
    %         f=bsxfun(@times,f,reshape(XSHIFT2(:,:,sdx),1,Nx2,1,Nv2));
    %     end
    
    
    %% Visualize phase space and electric field
    if (mod(tdx-1,floor(Nt/nviz))==0)
               
        tdxviz=(tdx-1)/floor(Nt/nviz) +1;
        
        %fdx=ifft(ifft(sum(sum(f,4),3)*dv2*dv1,[],2),[],1,'symmetric');
        %pcolor(x1,x2,fdx.');
        
        % By zero padding, the Fourier representation is used for
        % interpolation below the computational grid, which
        % produces higher quality output
        padNx=max(1,floor(512/Nx1));
        x1_pad=0:dx1/padNx:L(1)-dx1/padNx;
        x2_pad=0:dx2/padNx:L(2)-dx2/padNx;
        fdx=zeros(Nx1*padNx,Nx2*padNx);
        
        padNv=max(1,floor(512/Nv1));
        v1_pad=vmin1:dv1/padNv:vmax1-dv1/padNv;
        v2_pad=vmin2:dv2/padNv:vmax2-dv2/padNv;
        fdv=zeros(Nv1*padNv,Nv2*padNv);
        
        fx1v1=zeros(Nx1*padNx,Nv1*padNv );
        
        switch(splitting)
            % Fourier transformed in v
            case {'scovel','scovel_imr'}
                fdx( (-Nx1/2:Nx1/2-1)+Nx1*padNx/2+1, ....
                    (-Nx2/2:Nx2/2-1)+Nx2*padNx/2+1)=...
                    fftshift( fft2(F(:,:,1,1))*dv2*dv1);
                %fftshift( fft2(sum(sum(ifft(ifft(F,[],4),[],3),4),3))*dv2*dv1);
                
                
                fdv( (-Nv1/2:Nv1/2-1)+Nv1*padNv/2+1, ....
                    (-Nv2/2:Nv2/2-1)+Nv2*padNv/2+1)=...
                    fftshift(squeeze(sum(sum(F,1),2)))*dx1*dx2;
                
                %             fftshift( fft2(sum(sum(ifft(ifft(F,[],4),[],3),4),3))*dv2*dv1);
                %             fdv=squeeze(sum(sum(ifft(ifft(F,[],4),[],3),1),2));
                
                fx1v1( (-Nx1/2:Nx1/2-1)+Nx1*padNx/2+1, ....
                    (-Nv1/2:Nv1/2-1)+Nv1*padNv/2+1)=...
                    fftshift(squeeze(  fft(sum(F(:,:,:,1),2)*dx2*dv2,[],1)) );
                
                %             fx1v1=squeeze( ifft( sum(F(:,:,:,1),2)*dx2*dv2,[],3,'symmetric'  ));
                
            otherwise
                fdx( (-Nx1/2:Nx1/2-1)+Nx1*padNx/2+1, ....
                    (-Nx2/2:Nx2/2-1)+Nx2*padNx/2+1)=...
                    fftshift( sum(sum(F,4),3)*dv2*dv1);
                fdv( (-Nv1/2:Nv1/2-1)+Nv1*padNv/2+1, ....
                    (-Nv2/2:Nv2/2-1)+Nv2*padNv/2+1)=...
                    fft2(ifft(ifft(squeeze(F(1,1,:,:)*dx1*dx2),[],2),1));
                
                fx1v1( (-Nx1/2:Nx1/2-1)+Nx1*padNx/2+1, ....
                    (-Nv1/2:Nv1/2-1)+Nv1*padNv/2+1)=...
                    fftshift(squeeze(  fft(sum(F(:,1,:,:),4)*dx2*dv2,[],2)) );
                
        end
        fdx=ifft(ifft(ifftshift(fdx),[],2)*padNx^2,[],1,'symmetric');
        fdv=ifft(ifft(ifftshift(fdv),[],2)*padNv^2,[],1,'symmetric');
        fx1v1=ifft(ifft(ifftshift(fx1v1),[],2)*padNv*padNx,[],1,'symmetric');
        
        figure(figrho);
        pcolor(x1_pad,x2_pad,fdx.')
        shading interp;         colormap jet; colorbar; hold on;
        xlabel('x_1');ylabel('x_2');
        title(sprintf('$$\\int f_{\\mathrm{i}}(x,v,t=%g) \\mathrm{d}v$$',...
            (tdx-1)*dt),'interpreter','latex');
        %drawnow;
        print('-dpng', sprintf('%splots/fdv_%05d.png', prefix,tdxviz));
        
        figure(figv)
        pcolor(v1_pad,v2_pad,fdv.'); shading interp;
        colormap jet; colorbar; hold on;
        xlabel('v_1');ylabel('v_2');
        title(sprintf('$$\\int f_{\\mathrm{i}}(x,v,t=%g) \\mathrm{d}x$$',...
            (tdx-1)*dt),'interpreter','latex');
        %drawnow;
        print('-dpng', sprintf('%splots/fdx_%05d.png', prefix,tdxviz));
        
        
        figure(figx1v1)
        pcolor(x1_pad,v1_pad,fx1v1.');shading interp;
        colormap jet; colorbar; hold on;
        xlabel('x_1');ylabel('v_1');
        title(sprintf('$$\\int f_{\\mathrm{i}}(x,v,t=%g) \\mathrm{d}(x_2,v_2)$$',...
            (tdx-1)*dt),'interpreter','latex');
        print('-dpng', sprintf('%splots/fx1v1_%05d.png', prefix,tdxviz));
        drawnow;
        
        

        
        
        
        % Contour plot of electrostatic potential
        figure(figE)
                
        [C_,h_]=contour(x1,x2,ifft2(Phi.','symmetric'),'Linewidth',2); colormap; hold on;
        
        grid on;
        
        unpad=round(1./min(1,(8/length(x1))));
        quiver(x1(1:unpad:end),x2(1:unpad:end),...
            E1(1:unpad:end,1:unpad:end),E2(1:unpad:end,1:unpad:end),...
            dx*unpad*0.8,'Linewidth',2);
        clabel(C_,h_); %contour labels;
        axis([x1(1),x1(end),x2(1),x2(end)])
        
        %         pcolor(x1,x2,Phi.'); shading interp; colormap jet;
        xlabel('x_1');ylabel('x_2');
        title(sprintf('$$\\vec{E}(x,t)=-\\nabla \\Phi(x,t=%g)$$',...
            (tdx-1)*dt.*[1]),'interpreter','latex');
        print('-dpng', sprintf('%splots/phi_%05d.png', prefix,tdxviz));
        drawnow;
    end
end
F=ifft(ifft(F,[],2),[],1,'symmetric'); %back transform from x
runtime=toc;



%% Discussion
set(0,'DefaultLegendFontSize',16,...
    'DefaultLegendFontSizeMode','manual',...
    'DefaultAxesFontSize', 14)
set(0,'DefaultLineMarkerSize', 10);
set(0,'DefaultLineLinewidth', 2);
time=0:dt:tmax-dt;
energy_error=abs((Energy(1)-Energy)/Energy(1));
momentum_error=sqrt(sum((Momentum-Momentum(1,:)).^2,2));


%Save a file for later
clear VSHIFT* XSHIFT* BSHIFT* F fdx fdv fx1fv1 x1_pad x2_pad v1_pad v2_pad 
save([prefix,'result.mat']);

%% Field Energys
figure('Name','Electrostatic Energy','Numbertitle','off');
semilogy(time,Epot(:,2),'-','Linewidth',2,'DisplayName','$\frac{1}{2}||E_2||^2$');
hold on;
semilogy(time,Epot(:,1),'-','Linewidth',2,'DisplayName','$\frac{1}{2}||E_1||^2$');
xlabel('time'); grid on; axis tight;hold on;
axis tight;
ylabel('electrostatic energy');
l=legend('-Dynamiclegend','Location','Best');
set(l,'Interpreter','LaTeX','FontSize',16)
print('-dpng',[prefix, 'plots/fieldenergy.png'])
print('-depsc',[prefix, 'plots/fieldenergy.eps'])


% switch testcase
% %
%
%
%  tmin=find(time>4.5*B3,1,'first');
% P=polyfit(time(tmin:end)',log(Epot(tmin:end,1)),1);
% omega=P(1);
% % P(1)=
%
% plot(time,exp(P(1)*(time)+P(2)),...
%          '--','Color',[0,0,0,0.4],'DisplayName','linear analysis');


% % 0.522337*k0
% %
%   P(1)*B3
% %
% 2*(1-k0)
%
% % 2*sqrt(1-k0^2)*k0
%  (sqrt(1-k0^2) -k0)*k0



% figure()
% semilogy(time,PhiL2,'-','Linewidth',2');
% P=polyfit(time(tmin:end)',log(PhiL2(tmin:end,1)),1); hold on;
% plot(time,exp(P(1)*(time)+P(2)),...
%          '--','Color',[0,0,0,0.4],'DisplayName','linear analysis');
% l=legend('-Dynamiclegend','Location','Best');
% set(l,'Interpreter','LaTeX','FontSize',16)
%
% P(1)*B3
% % Include decay rate for linear landau damping WITHOUT COLLISIONS
% if (eps<0.1 && norm(k0)==0.5)
%     hold on;
%     % Obtain zero of dispersion relation with dispersion_landau.m
%     omega=1.415661888604536 - 0.153359466909605i;
%     plot(time,Epot(1)*abs(exp(-1j*omega*(time-0.4))).^2,...
%         '--','Color',[0,0,0,0.4],'DisplayName','linear analysis');
%     % linear analysis with frequency
%     %plot(time,0.5*fieldenergy(1)*real(exp(-1j*omega*(time-0.4))).^2,);
%     %legend('numerical', 'linear analysis', 'linear analysis');
%
%     hold off;
% end
%
%% Kinetic energy
figure('Name','Kinetic Energy','Numbertitle','off');
plot(time,Ekin,'-','Linewidth',2); hold on;
grid on;
xlabel('time'); grid on; axis tight;
% axis tight;
ylabel('kinetic energy');
l=legend('$ \mathcal{H}_{p_1}$','$ \mathcal{H}_{p_2}$',...
    'Location','Best');
set(l,'Interpreter','LaTeX','FontSize',16)
print('-dpng',[prefix, 'plots/kineticenergy.png'])
print('-depsc',[prefix, 'plots/kineticenergy.eps'])

figure('Name','Relative Energy Error','Numbertitle','off');
semilogy(time, abs((Energy(1)-Energy)/Energy(1)));
xlabel('time'); grid on;
set(gca,'Ytick',10.^(-16:1:4))
ylabel('relative energy error');
print('-dpng',[prefix, 'plots/energy_error.png'])
print('-depsc',[prefix, 'plots/energy_error.eps'])
%
%% Momentum Error
figure('Name','Absolute Momentum Error','Numbertitle','off');
semilogy(time, abs(Momentum-Momentum(1,:)))
legend('P_1','P_2','Location','SouthEast');
xlabel('time');
title('absolute momentum error');
axis tight; grid on;
set(gca,'Ytick',10.^(-16:1:4))
print('-dpng',[prefix, 'plots/momentum_error.png'])
print('-depsc',[prefix, 'plots/momentum_error.eps'])
% close all;
