%% Spectral 2D Vlasov-Poisson Solver
% Author: *Jakob Ameres* <http://jakobameres.com jakobameres.com>
%% Features
%
% * Spectral description of density
% * *Symplectic* (for Vlasov-Poisson) Runge Kutta *time integrator*
%


%% Landau Damping parameters
eps=0.5; % Amplitude of perturbation, 0.05 for linear, 0.5 for nonlinear
k0=0.5;    % Wave vector
L=2*pi/k0; % length of domain

% Multiple species
q=[1,-1];
m=[100,1];

vth_i=0.1;

% Initial condition for each species
f0={@(x,v) (1+0*cos(k0*x))./sqrt(2*pi).*exp(-0.5.*v.^2/vth_i.^2)/vth_i,...
    @(x,v) (1+eps*cos(k0*x))./sqrt(2*pi).*exp(-0.5.*v.^2)};

Ns=length(q); % Number of species

% Background to be subtracted from the phase space plot
% background=@(x,v) 1./sqrt(2*pi).*exp(-0.5.*v.^2);
% External electric field
E_ext=@(x,t) 0.*x+0.*t;

dt=0.01; % Time step
tmax=100;
rungekutta_order=3; %Order of the runge Kutta time integrator


Nx=64; % Number of cells in spatial direction
Nv=64; % Number of cells in velocity direciton

vmax=[vth_i.^2,1]*4.5;
vmin=-vmax;


%% Define phase space mesh
v=cell(Ns,1); kv=cell(Ns,1); 
for qdx=1:Ns
v{qdx}=linspace(vmin(qdx),vmax(qdx),Nv+1);v{qdx}=v{qdx}(1:Nv);
kv{qdx}=fftshift( -Nv/2:Nv/2-1)*2*pi/(vmax(qdx)-vmin(qdx));
end

x=linspace(0,L,Nx+1).';x=x(1:end-1);
kx=fftshift( -Nx/2:Nx/2-1)'*2*pi/L;

dv=(vmax-vmin)/Nv; dx=L/Nx;



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

%Build Matrix for spatial advection and all stages
XSHIFT=zeros(Nx,Nv,Ns,length(rksd));
for jdx=1:length(rksd)
    for qdx=1:Ns
        XSHIFT(:,:,qdx,jdx)=exp(-1j*dt*rksd(jdx)*(kx*v{qdx}));
    end
end

%diagnostic
Nt=ceil(tmax/dt);
fieldenergy=zeros(Nt,Ns);
kineticenergy=zeros(Nt,Ns);
momentum=zeros(Nt,Ns);

% Set initial conditionfor each species
% [XX,VV]=ndgrid(x,v);
f=zeros(Nx,Nv,length(f0));
for qdx=1:length(f0)
    f(:,:,qdx)=f0{qdx}(XX,VV);
end

%% Fourier transform in x
f=fft(f,[],1);

figf=cell(Ns,1);
for qdx=1:Ns
    figf{qdx}=figure('Numbertitle','off','Name',...
        sprintf('Phase Space Density, q=%g, m=%g',q(qdx),m(qdx)));
end
for tdx=1:Nt
    
    %% Loop over all stages of the symplectic Runge Kutta
    for sdx=1:length(rksd)
        
        %% Poisson equation
        rho=(q(:).*dv(:))'.*squeeze(sum(f,2)); %rho(x) -> integrate over v
        
        E=rho./(-1j*kx);
        % remove constant fourier mode
        E(1,:)=0;
        if (sdx==length(rksd))
            for qdx=1:Ns
                %L2 norm of Electric field
                fieldenergy(tdx,qdx)=sum(E(:,qdx)'*E(:,qdx))/2/Nx*dx;
            end
            %EE(:,tdx)=E;
        end
        E=sum(E,2);
        E=ifft(E,'symmetric');
        E=E+E_ext(x,tdx*dt);
        %plot(x,ifft(E,'symmetric'));
        
        
        %% Advection in v
        f=ifft(f,[],1,'symmetric'); %back transform
        % f is fully backtransformed, calculate diagnostics
        if (sdx==length(rksd))
            for qdx=1:Ns
                kineticenergy(tdx,qdx)=m(qdx)*sum(sum(f(:,:,qdx),1).*v{qdx}.^2/2)*dv(qdx)*dx;
                momentum(tdx,qdx)=m(qdx)*sum(sum(f(:,:,qdx),1).*v{qdx})*dv(qdx)*dx;
            end
        end
        
        % Fourier transform in v
        f=fft(f,[],2);
        for qdx=1:length(q)
            f(:,:,qdx)=f(:,:,qdx).*exp(1j*dt*rksc(sdx)*q(qdx)/m(qdx)*E*kv{qdx});
        end
        f=ifft(f,[],2,'symmetric');
        
        
        %%  Advection in x
        f=fft(f,[],1);
        f=f.*XSHIFT(:,:,:,sdx);
        
    end
    
    
    %% Visualize phase space
    if (mod(tdx,floor(Nt/20))==0)
        for qdx=1:Ns
            figure(figf{qdx});
            pcolor(x,v{qdx},ifft(f(:,:,qdx),[],1,'symmetric')');
            shading interp;
            colormap jet; colorbar; hold on;
            title(sprintf('f_s(x,v) @ t=%g', (tdx-1)*dt));
            xlabel('x');ylabel('v');
            drawnow;
            
        end
    end
end

f=ifft(f,[],1,'symmetric');

%% Discussion
set(0,'DefaultLineLinewidth', 2);
time=(0:Nt-1)*dt;
figure('Name','Electrostatic Energy','Numbertitle','off');
semilogy(time, fieldenergy);
xlabel('time'); grid on;
ylabel('electrostatic energy');

% Include decay rate for linear landau damping WITHOUT COLLISIONS
if (eps<0.1 && k0==0.5)
    hold on;
    % Obtain zero of dispersion relation with dispersion_landau.m
    omega=1.415661888604536 - 0.153359466909605i;
    plot(time,0.5*fieldenergy(1,2)*abs(exp(-1j*omega*(time-0.4))).^2,'--');
    % linear analysis with frequency
    plot(time,0.5*fieldenergy(1)*real(exp(-1j*omega*(time-0.4))).^2,'Color',[0,0,0,0.4]);
    legend('numerical', 'linear analysis', 'linear analysis');
    
    hold off;
end

figure('Name','Kinetic Energy','Numbertitle','off');
semilogy(time, kineticenergy);
xlabel('time'); grid on;
ylabel('kinetic energy');



energy=sum(kineticenergy+fieldenergy,2);
figure('Name','Relative Energy Error','Numbertitle','off');
semilogy(time, abs((energy(1)-energy)/energy(1)));
xlabel('time'); grid on;
ylabel('relative energy error');

