%% Spectral 2D Vlasov-Poisson Solver
% Author: *Jakob Ameres* <http://jakobameres.com jakobameres.com>
%% Features
%
% * Fourier spectral description of density f
% * *Symplectic* (for Vlasov-Poisson) Runge Kutta *time integrator*
%
clear all; close all;
normalize=@(x) x./norm(x);

testcase='bumpontail';
% External electric field
E1_ext=[]; %E1_ext=@(x1,x2,t) 0.*x1 +0.*x2+0.*t;
E2_ext=[]; %E2_ext=@(x1,x2,t) 0.*x1 +0.*x2 +0.*t;

dt=0.1; % Time step
tmax=15;
rungekutta_order=3; %Order of the runge Kutta time integrator


switch(testcase)
    case 'landau'
        %% Diagonal Landau Damping parameters
        eps=0.5; % Amplitude of perturbation, 0.05 for linear, 0.5 for nonlinear
        k0=normalize([1,2]).*0.5;    % Wave vector
        L=2*pi./k0; % length of domain
        q=-1;    % negative charge to mass = electrons, $\frac{q}{m}=-1$
        m=1;
        % Initial condition
        f0=@(x1,x2,v1,v2) (1+eps*cos(k0(1).*x1+k0(2).*x2)).*...
            exp(-0.5.*(v1.^2+v2.^2))/(2*pi);
        vmax1=4.5;
        vmin1=-4.5;
        vmax2=4.5;
        vmin2=-4.5;
        tmax=25; dt=0.1;
        
    case 'bumpontail'
        %% Diagonal Bump-on-Tail instability
        q=-1;    % negative charge to mass = electrons, $\frac{q}{m}=-1$
        m=1;
        %Nakamura, Yabe
        eps=0.001;  % We use a small amplitude because it is hard
        % to excite the correct eigenvalue, such that there
        % are spurious oscillations on top of the linear
        % growth. By choosing a small amplitude we give the
        % system more time in the linear phase to snap into
        % the right mode.
        nb=0.1; %bump size
        v0=4.5;
        bumpsigma=0.5;
        p=normalize([1,2]); % direction
        k0=p.*0.3;    % Wave vector
        
        L=2*pi./k0; % length of domain
        f0=@(x1,x2,v1,v2) (1+eps*cos(k0(1).*x1+k0(2).*x2)).*...
            ( (1-nb).*exp(-0.5.*(v1.^2+v2.^2)) + nb/bumpsigma*...
            exp(-0.5.*( (v1 - p(1)*v0).^2 + (v2-v0*p(2)).^2 )/bumpsigma.^2)  )/(2*pi);
        vmax=7.5; vmin=-4.5;
        vmax1=vmax;vmin1=vmin; vmax2=vmax;  vmin2=vmin;
        tmax=100;
        dt=0.2;
        
end



Nx1=16; % Number of cells in spatial direction
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


%% Time integration scheme for symplectic runge kutta
switch rungekutta_order
    case 1
        rksd=[1, 0]; %Symplectic Euler
        rksc=[0, 1];
    case 2
        rksd=[0.5, 0.5 ];
        rksc=[0, 1];
    case 3
        rksd=[2/3, -2/3, 1  ];
        rksc=[ 7/24, 3/4, -1/24];
    case 4
        rk4sx=real((2^(1/3) +2^(-1/3)-1)/6);
        rksd=[ 2*rk4sx +1 , -4*rk4sx-1, 2*rk4sx+1, 0];
        rksc=[ rk4sx + 0.5 , -rk4sx, -rk4sx, rk4sx +0.5];
end

% %Build Matrix for spatial advection and all stages
XSHIFT1=zeros(Nx1,Nv1,length(rksd));
XSHIFT2=zeros(Nx2,Nv2,length(rksd));
for jdx=1:length(rksd)
    XSHIFT1(:,:,jdx)=exp(-1j*dt*rksd(jdx)*(kx1(:)*v1(:)'));
    XSHIFT2(:,:,jdx)=exp(-1j*dt*rksd(jdx)*(kx2(:)*v2(:)'));
end
%Preallocate matrix for velocity advection
% VSHIFT=zeros(Nx,Nv);

%diagnostic
Nt=ceil(tmax/dt);
Epot=zeros(Nt,2);
Ekin=zeros(Nt,2);
Momentum=zeros(Nt,2);
Energy=zeros(Nt,1);

% Set initial condition
% [XX1,XX2,VV1,VV2]=ndgrid(x1,x2,v1,v2);
% f=f0(XX1,XX2,VV1,VV2);
[XX1,XX2]=ndgrid(x1,x2); % only create 2d grid in order to save memory
f=zeros(Nx1,Nx2,Nv1,Nv2);
for vdx2=1:Nv2
    for vdx1=1:Nv1
        f(:,:,vdx1,vdx2)=f0(XX1,XX2,v1(vdx1),v2(vdx2));
    end
end

%% Fourier transform in x
tic;
f=fft(fft(f,[],1),[],2);
figfx=figure('Name', 'Spatial Density', 'Numbertitle','off');
figfv=figure('Name', 'Veloctiy Density', 'Numbertitle','off');
for tdx=1:Nt
    t=(tdx-1)*dt;
    %% Loop over all stages of the symplectic Runge Kutta
    for sdx=1:length(rksd)
        
        %% Poisson equation
        rho=q*sum(sum(f,4),3)*dv1*dv2; %rho(x) -> integrate over v
        
        rho(1,1)=0; %remove (ion) background
        
        Phi=-rho./(  kx1.^2 +kx2.^2); Phi(1,1)=0;
        E1=-Phi.*(1j.*kx1);
        E2=-Phi.*(1j.*kx2);
        
        if (sdx==1)
            Epot(tdx,1) =(E1(:)'*E1(:))/2*dx/(Nx1*Nx2);
            Epot(tdx,2) =(E2(:)'*E2(:))/2*dx/(Nx1*Nx2);%
        end
        
        E1=ifft(ifft(E1,[],2),[],1,'symmetric');
        E2=ifft(ifft(E2,[],2),[],1,'symmetric');
        % add external fields
        if (~isempty(E1_ext)), E1=E1+E1_ext(XX1,XX2,t); end
        if (~isempty(E2_ext)), E2=E2+E2_ext(XX1,XX2,t); end
        t=t+dt*rksc(sdx);
        
        
        %% Advection in v
        f=ifft(ifft(f,[],2),[],1); %back transform
        
        % f is fully backtransformed, calculate diagnostics
        if (sdx==1)
            fdx=real(sum(sum(f,1),2));
            Ekin(tdx,1)=sum(reshape(sum(fdx,4),Nv1,1).*v1(:).^2/2)*dx*dv*m;
            Ekin(tdx,2)=sum(reshape(sum(fdx,3),Nv2,1).*v2(:).^2/2)*dx*dv*m;
            Momentum(tdx,1)=sum(reshape(sum(fdx,4),Nv1,1).*v1(:))*dx*dv*m;
            Momentum(tdx,2)=sum(reshape(sum(fdx,3),Nv2,1).*v2(:))*dx*dv*m;
            Energy(tdx)=sum(Ekin(tdx,:))+sum(Epot(tdx,:));
        end
        
        % Fourier transform in v
        f=fft(fft(f,[],3),[],4);
        %f=f.*exp(1j*dt*rksc(sdx)*q*E1.*kv1).*...
        %   exp(1j*dt*rksc(sdx)*q*E2.*kv2);
        f=f.*exp(1j*dt*rksc(sdx)*q*(E1.*kv1+E2.*kv2));
        
        f=ifft(ifft(f,[],4),[],3);
        
        
        %%  Advection in x
        f=fft(fft(f,[],1),[],2);
        %f=f.*exp(-1j*dt*rksd(sdx)*(kx1.*reshape(v1,1,1,Nv1,1)))...
        %    .*exp(-1j*dt*rksd(sdx)*(kx2.*reshape(v2,1,1,1,Nv2)));
        
        f=bsxfun(@times,f,reshape(XSHIFT1(:,:,sdx),Nx1,1,Nv1,1));
        f=bsxfun(@times,f,reshape(XSHIFT2(:,:,sdx),1,Nx2,1,Nv2));
    end
    
    
    %% Visualize phase space
    if (mod(tdx,floor(Nt/20))==0)
        
        figure(figfx);
        %fdx=ifft(ifft(sum(sum(f,4),3)*dv2*dv1,[],2),[],1,'symmetric');
        %pcolor(x1,x2,fdx.');
        
        
        % By zero padding, the Fourier representation is used for
        % interpolation below the computational grid, which
        % produces higher quality output
        padN=max(1,floor(512/Nx1));
        x1_pad=0:dx1/padN:L(1)-dx1/padN;
        x2_pad=0:dx2/padN:L(2)-dx2/padN;
        fdx=zeros(Nx1*padN,Nx2*padN);
        fdx( (-Nx1/2:Nx1/2-1)+Nx1*padN/2+1, ....
            (-Nx2/2:Nx2/2-1)+Nx2*padN/2+1)=fftshift(sum(sum(f,4),3)*dv2*dv1);
        fdx=ifft(ifft(ifftshift(fdx),[],2)*padN^2,[],1,'symmetric');
        pcolor(x1_pad,x2_pad,fdx.')
        shading flat;
        colormap jet; colorbar; hold on;
        title(sprintf('\\int f dv t=%g', (tdx-1)*dt));
        xlabel('x_1');ylabel('x_2');
        
        
        figure(figfv);
        padN=max(1,floor(512/Nv1));
        v1_pad=vmin1:dv1/padN:vmax1-dv1/padN;
        v2_pad=vmin2:dv2/padN:vmax1-dv2/padN;
        fdv=zeros(Nv1*padN,Nv2*padN);
        fdv( (-Nv1/2:Nv1/2-1)+Nv1*padN/2+1, ....
            (-Nv2/2:Nv2/2-1)+Nv2*padN/2+1)=...
            fftshift(fft2(squeeze( real( f(1,1,:,:)*dx1*dx2/Nx1/Nx2 ) )));
        fdv=ifft(ifft(ifftshift(fdv),[],2)*padN^2,[],1,'symmetric');
        pcolor(v1_pad,v2_pad,fdv.')
        
        shading flat;
        colormap jet; colorbar; hold on;
        title(sprintf('\\int f dx t=%g', (tdx-1)*dt));
        xlabel('v_1');ylabel('v_2');
        drawnow;
        
    end
end
f=ifft(ifft(f,[],2),[],1,'symmetric'); %back transform from x
runtime=toc;



%% Discussion
set(0,'DefaultLegendFontSize',16,...
    'DefaultLegendFontSizeMode','manual',...
    'DefaultAxesFontSize', 14)
set(0,'DefaultLineMarkerSize', 10);
set(0,'DefaultLineLinewidth', 2);
time=0:dt:tmax-dt;

%% Field Energys
figure;
semilogy(time,Epot(:,2),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_2||^2$');
hold on;
semilogy(time,Epot(:,1),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_1||^2$');
xlabel('time'); grid on; axis tight;
axis tight;
ylabel('electrostatic energy');
l=legend('-Dynamiclegend','Location','SouthEast');
set(l,'Interpreter','LaTeX','FontSize',16)

switch( testcase)
    case 'landau'
        % Include decay rate for linear landau damping WITHOUT COLLISIONS
        if (eps<0.1 && norm(k0)==0.5)
            hold on;
            % Obtain zero of dispersion relation with dispersion_landau.m
            omega=1.415661888604536 - 0.153359466909605i;
            plot(time,Epot(1)*abs(exp(-1j*omega*(time-0.4))).^2,...
                '--','Color',[0,0,0,0.4],'DisplayName','linear analysis');
            % linear analysis with frequency
            %plot(time,0.5*fieldenergy(1)*real(exp(-1j*omega*(time-0.4))).^2,);
            %legend('numerical', 'linear analysis', 'linear analysis');
            
            hold off;
        end
        
    case 'bumpontail'
%          omega= (1.0012 + 0.1981i);
%          time2=time(time<35);
%             plot(time2,Epot(1)*abs(exp(-1j*omega*(time2-7))).^2,...
%                 '--','Color',[0,0,0,0.4],'DisplayName','linear analysis');
end
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


figure('Name','Relative Energy Error','Numbertitle','off');
semilogy(time, abs((Energy(1)-Energy)/Energy(1)));
xlabel('time'); grid on;
set(gca,'Ytick',10.^(-16:1:4))
ylabel('relative energy error');


%% Momentum Error
figure('Name','Absolute Momentum Error','Numbertitle','off');
semilogy(time, abs(Momentum-Momentum(1,:)))
legend('P_1','P_2','Location','SouthEast');
xlabel('time');
title('absolute momentum error');
axis tight; grid on;
set(gca,'Ytick',10.^(-16:1:4))
