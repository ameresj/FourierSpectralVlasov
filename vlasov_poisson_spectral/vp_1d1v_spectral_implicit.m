%% Spectral 2D Vlasov-Poisson Solver
% Author: *Jakob Ameres* <http://jakobameres.com jakobameres.com>
%% Features
%
% * Spectral description of density
% * *Symplectic* (for Vlasov-Poisson) Runge Kutta *time integrator*
% 
close all; clear all;

%% Landau Damping parameters
eps=0.01; % Amplitude of perturbation, 0.05 for linear, 0.5 for nonlinear
k0=0.5;    % Wave vector
L=2*pi/k0; % length of domain
qm=-1;    % negative charge to mass = electrons, $\frac{q}{m}=-1$

% Initial condition
f0=@(x,v) (1+eps*cos(k0*x))./sqrt(2*pi).*exp(-0.5.*v.^2);
% Background to be subtracted from the phase space plot
background=@(x,v) 1./sqrt(2*pi).*exp(-0.5.*v.^2);
% External electric field
E_ext=@(x,t) 0.*x+0.*t;

dt=1; % Time step
tmax=30;
rungekutta_order=3; %Order of the runge Kutta time integrator


Nx=128; % Number of cells in spatial direction
Nv=128; % Number of cells in velocity direciton
vmax=4.5;
vmin=-4.5;

%% Define phase space mesh
v=linspace(vmin,vmax,Nv+1);v=v(1:end-1);
x=linspace(0,L,Nx+1).';x=x(1:end-1);
kx=fftshift( -Nx/2:Nx/2-1)'*2*pi/L;
kv=fftshift( -Nv/2:Nv/2-1)*2*pi/(vmax-vmin);
dv=(vmax-vmin)/Nv; dx=L/Nx;

[XX,VV]=ndgrid(x,v);

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
XSHIFT=zeros(Nx,Nv,length(rksd));
for jdx=1:length(rksd)
   XSHIFT(:,:,jdx)=exp(-1j*dt*rksd(jdx)*(kx*v));
end
%Preallocate matrix for velocity advection
VSHIFT=zeros(Nx,Nv);

%diagnostic
Nt=ceil(tmax/dt);
fieldenergy=zeros(Nt,1);
kineticenergy=zeros(Nt,1);
momentum=zeros(Nt,1); 

% Poisson solver
    poissonk=1./(-1j*kx);poissonk(1)=0;
     Efun=@(f) ifft((-f(:,1)*dv).*poissonk,'symmetric');

% Set initial condition
f=f0(XX,VV);
%% Fourier transform in x
f=fft(f,[],1);f=fft(f,[],2);
figure('Name', 'Phase Space Density', 'Numbertitle','off');
for tdx=1:Nt
    
       
                           E=fft(Efun(f));
             fieldenergy(tdx)=sum(E'*E)/2/Nx*dx; %L2 norm of Electric field
             kineticenergy(tdx)=sum(ifft( (f(1,:)),'symmetric').*v.^2/2)*dv*dx;
             momentum(tdx)=sum(ifft((f(1,:)),'symmetric').*v)*dv*dx;
    
%         rho=-fft(sum(f,2)*dv);
%         E=rho./(-1j*kx); E(1)=0;
%         E=ifft(E,'symmetric');
%         Df0=v.*ifft(1j.*kx.*fft(f,[],1),[],1,'symmetric') +...
%             E.*ifft(1j.*kv.*fft(f,[],2),[],2,'symmetric');
%     f_=f;
%     
%     err=Inf;iter=0;
%     while( err>1e-20)
%            rho=-fft(sum(f,2)*dv).*0;
%         E=rho./(-1j*kx); E(1)=0;
%         E=ifft(E,'symmetric');
%         Df=(v.*ifft(1j.*kx.*fft(f_,[],1),[],1,'symmetric') +...
%             E.*ifft(1j.*kv.*fft(f_,[],2),[],2,'symmetric')+Df0)/2*(-dt) ;
%         
%         
%     err=norm(f_(:)-(f(:)+Df(:))  );
%     f_=f+Df; 
%     iter=iter+1;
%     end
%     iter
%     
%     if any(isnan(f_(:)))
%         stop
%     end
%     f=f_;

    
%             rho=(f(:,1))*dv;
%         E=rho./(-1j*kx); E(1)=0;
%         E=ifft(E,'symmetric');
%         Df0=-dt.*(fft(v.*ifft(1j*kx.*f,[],2),[],2)+...
%             fft(-E.*ifft(1j*kv.*f,[],1),[],1));





       
     Gfun=@(f,E)(fft(v.*ifft(1j*kx.*f,[],2),[],2)+...
            fft(E.*ifft(1j*kv.*f,[],1),[],1));
        
     DG=@(df) (fft(v.*ifft(1j*kx.*df,[],2),[],2)+...
                  fft(Efun(df).*ifft(1j*kv.*f,[],1),[],1)+...
                  fft(Efun(f).*ifft(1j*kv.*df,[],1),[],1));
       
    E0=Efun(f);
    f0=f;
    G0=Gfun(f,E0);        
        
    
    % Newton iterations
    err=Inf;iter=0;f_=0;
    while( err>1e-15 && iter<10)
        E1=Efun(f);
        E=(E0+E1)/2;
        %G=-dt*(  Gfun(f0,E)+Gfun(f,E))/2;
        G=-dt*(  Gfun(f0+f,E))/2;
        
        
        DG=@(df) -dt*( fft(0.5*Efun(df).*ifft(1j*kv.*f0,[],1),[],1)+...
            fft(v.*ifft(1j*kx.*df,[],2),[],2)+...
            fft(0.5*Efun(df).*ifft(1j*kv.*f,[],1),[],1)+...
            fft(E.*ifft(1j*kv.*df,[],1),[],1))/2 ;
        
        df = gmres(@(x) reshape((DG((reshape(x,Nx,Nv)))),Nx*Nv,1)-x,...
            f0(:)+G(:)-f(:),25,1e-8,10);
% rhs=        ifft2(f0+G-f);
% df = (gmres(@(x) reshape(ifft2(DG(fft2(reshape(x,Nx,Nv)))),Nx*Nv,1)-x,...
%            rhs(:),5,1e-5,10));
%        df=reshape(df,Nx,Nv); df=fft2(df);
%         
%         df = -pcg(@(x) -(reshape((DG((reshape(x,Nx,Nv)))),Nx*Nv,1)-x),...
%             f0(:)+G(:)-f(:),1e-8,1000);
%         
        df=reshape(df,Nx,Nv);
        f=f-df;
        
        err=norm(df(:))/Nx/Nv;
%         f_=f;
        %        DG=@(df) -dt*(  Gfun(f0,Efun(df)/2)+  Gfun(df,E)+ Gfun(f,Efun(df)/2))/2;
        
        %        DG=@(df) -dt*( Gfun(df,E))/2;
        % %         DG=@(df) -dt*( Gfun(f0,Efun(df)/2)+  Gfun(df,E)+ Gfun(f,Efun(df)/2) )/2;
        
        %  a=rand(Nx,Nv);
        %          b=rand(Nx,Nv);
        %         norm(Gfun(a,E)+Gfun(b,E)-Gfun(a+b,E))
        %         E1=rand(Nx,1);
        %         E2=rand(Nx,1);
        %         norm(Gfun(f,E1)+Gfun(f,E2)-Gfun(f,E1+E2))
        %
        %        norm(DG(a)+DG(b) -DG(a+b))/norm(a)
        
        %            GG = gmres(@(x) reshape((DG((reshape(x,Nx,Nv)))),Nx*Nv,1),rand(Nx*Nv,1),[],1e-10,100);
        %           GG = qmr(@(x) reshape((DG((reshape(x,Nx,Nv)))),Nx*Nv,1),rand(Nx*Nv,1),1e-10,100);
        
        %                    G=ifft2(G);  G=G(:);
        %            GG = gmres(@(x) reshape(ifft2(DG(fft2(reshape(x,Nx,Nv)))),Nx*Nv,1),G,[],1e-5,100);
        
        %            dim=2;
        %            G=ifft(G,[],dim);  G=G(:);
        %            GG = gmres(@(x) reshape(ifft(DG(fft(reshape(x,Nx,Nv),[],dim)),[],dim),Nx*Nv,1),G,[],1e-5,10);
        %     GG=fft(G,[],dim);
        % GG=G;
        
        %     err=norm(GG(:)-GG_(:)  )/Nx/Nv;
        %     f=f+reshape(GG,Nx,Nv);
        %     GG_=GG;
        iter=iter+1;
    end
% iter





%      poissonk=1./(-1j*kx);poissonk(1)=0;
%      Efun=@(f) ifft((-f(:,1)*dv).*poissonk,'symmetric');
%        
%      Gfun=@(f)(fft(v.*ifft(1j*kx.*f,[],2),[],2)+...
%             fft(Efun(f).*ifft(1j*kv.*f,[],1),[],1));
%         
%      DG=@(df) (fft(v.*ifft(1j*kx.*df,[],2),[],2)+...
%                   fft(Efun(df).*ifft(1j*kv.*f,[],1),[],1)+...
%                   fft(Efun(f).*ifft(1j*kv.*df,[],1),[],1));
%            
%     G0=-dt*Gfun(f);        
%         
%     
%     % Picard iterations
%     f_=f;    
%     err=Inf;iter=0;G_=0;
%     while( err>1e-16)
%         G=(-dt*Gfun(f_)+G0)/2;   
%         
%     err=norm(G(:)-G_(:)  )/Nx/Nv;
%     f_=f+G; 
%     G_=G;
%     iter=iter+1;
%     end
%     f=f_;
%     
    
     %% Poisson equation
%        
%         f=f+dt.*(fft(v.*ifft(1j*kx.*f,[],2),[],2) +... 
%             0.*fft(E.*ifft(-1j*kv.*f,[],1,'symmetric'),[],1));
%       
    
    
    
    
    
    
    % Crank Nicolson
    
  
    

%        %% Poisson equation
%         rho=f(:,1)*dv;
%         E=rho./(-1j*kx); E(1)=0;
%         E=ifft(E,'symmetric');
% 
%         f=f+dt.*(fft(v.*ifft(1j*kx.*f,[],2),[],2) +... 
%             0.*fft(E.*ifft(-1j*kv.*f,[],1,'symmetric'),[],1));
        
        %         E=E+E_ext(x,tdx*dt);
%         % remove constant fourier mode
%         E(1)=0; 
%         if (sdx==1)
              %E=fft(Efun(f));
%               E=fft(E);
%              fieldenergy(tdx)=sum(E'*E)/2/Nx*dx; %L2 norm of Electric field
%              kineticenergy(tdx)=sum(ifft( (f(1,:)+f0(1,:))/2,'symmetric').*v.^2/2)*dv*dx;
%              momentum(tdx)=sum(ifft((f(1,:)+f0(1,:))/2,'symmetric').*v)*dv*dx;
%             
             
          

             
%             EE(:,tdx)=E;
%         end
        %plot(x,ifft(E,'symmetric'));
    
    
    
    
    
    
%     
%     
%     
%     %% Loop over all stages of the symplectic Runge Kutta
%     for sdx=1:length(rksd)
%         
%         
%         
%         %% Advection in v
%         f=ifft(f,[],1,'symmetric'); %back transform
%         % f is fully backtransformed, calculate diagnostics
%         if (sdx==1)
%             kineticenergy(tdx)=sum(sum(f,1).*v.^2/2)*dv*dx;
%             sh_entropy(tdx)=sum(sum(f.*log(abs(f))))*dx*dv;
%             kl_entropy(tdx)=sum(sum(f.*log(abs(f./f0(XX,VV)))))*dx*dv;
%         end
%         
%         % Fourier transform in v
%         f=fft(f,[],2);
%         % Build Matrix for spatial advection
%         VSHIFT=exp(1j*dt*rksc(sdx)*E*kv);
%                     
%         f=f.*VSHIFT;
%         f=ifft(f,[],2,'symmetric');
%         
%         
%         %%  Advection in x
%         f=fft(f,[],1);
%         f=f.*XSHIFT(:,:,sdx);
%         
%     end
    
    
    %% Visualize phase space
       if (mod(tdx,floor(Nt/20))==0)
         pcolor(x,v,( real(ifft2(f))-background(XX,VV)).'); 
%         pcolor(x,v,( f-background(XX,VV)).'); 
        
%         pcolor(x,v,f')
        shading interp;
        colormap jet; colorbar; hold on;
        title(sprintf('\\deltaf @ t=%g', (tdx-1)*dt));
        xlabel('x');ylabel('v');
        drawnow;
       end
% 
end
 f=ifft(ifft(f,[],2),[],1,'symmetric');


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
plot(time,0.5*fieldenergy(1)*abs(exp(-1j*omega*(time-0.4))).^2,'--');
% linear analysis with frequency
plot(time,0.5*fieldenergy(1)*real(exp(-1j*omega*(time-0.4))).^2,'Color',[0,0,0,0.4]);
legend('numerical', 'linear analysis', 'linear analysis');

hold off;
end

figure('Name','Kinetic Energy','Numbertitle','off');
semilogy(time, kineticenergy);
xlabel('time'); grid on;
ylabel('kinetic energy');



% 
% 
energy=kineticenergy+fieldenergy;
figure('Name','Relative Energy Error','Numbertitle','off');
semilogy(time, abs((energy(1)-energy)/energy(1)));
xlabel('time'); grid on;
ylabel('relative energy error');



%% Momentum Error
figure('Name','Absolute Momentum Error','Numbertitle','off');
semilogy(time, abs(momentum-momentum(1,:)))
legend('P_1','P_2','Location','SouthEast');
xlabel('time');
title('absolute momentum error');
axis tight; grid on;
set(gca,'Ytick',10.^(-16:1:4))
% 
% 
% %Kullback-Leibler-Entropy relative to initial condition
% figure('Name','Kullback-Leibler-Entropy','Numbertitle','off');
% plot(time,kl_entropy)
% xlabel('time'); grid on;
% ylabel('entropy');
% title('Kullback-Leibler-Entropy relative to f(t=0)');
% 
% 
% %Shannon-Entropy
% figure('Name','Shannon-Entropy of f','Numbertitle','off');
% plot(time,sh_entropy)
% xlabel('time'); grid on;
% ylabel('entropy of f');
% title('Shannon-Entropy of f');
