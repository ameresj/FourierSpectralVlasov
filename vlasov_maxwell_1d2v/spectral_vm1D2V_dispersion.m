clear all; close all;


% Set number of threads to 4
maxNumCompThreads(4)
dt=0.05;      % 0.05 0.025 0.01 0.01 for (32 64 128 256)
dt2=0.05; % Diagnostic frequency

q=-1;
m=1;
N=256;
Nx=N;
Nv1=64;
Nv2=64;
Nplot=0; % Number of plots over simulation time

% splitting='strang';    %strang 2nd4lie fourth10lie
splitting='2nd4lie';    %strang 2nd4lie fourth10lie

filtration=false; % Filamentation filtration by 2/3 rule
% The Hotspots in this Code are the complex exponentials

% Noise level in order to excite modes in spatial direction of F(X),
% or the fields E2, B. E1 cannot be set since it is obtained by the 
% Poisson equation
noise=struct('X',1e-6,'E2',0,'B',1e-6);

%%
% Available test cases are:
%
% * landau - strong Landau damping
% * weibel - the Weibel instability
% * weibels - the streaming Weibel instability
%
testcase='landau'; % 'weibel','weibels','landau'


switch (testcase)
    
    case('weibel')
        % Weibel instability
        eps=1e-3;
        betar=sign(q)*1e-3; betai=0;
        k=1.25;    % Wave vector
        sigma1=0.02/sqrt(2);
        sigma2=sqrt(12)*sigma1;
        v01=0;
        v02=0;
        delta=0;
        tmax=500;
        vmax1=4.5*max(sigma1); vmin1=-vmax1;
        vmax2=4.5*max(sigma2); vmin2=-vmax2;
        
    case {'weibels','weibels2'}
        % Asymmetric Streaming Weibel instability
        sigma1=0.1/sqrt(2);
        sigma2=sigma1;
        k=0.2;
        betai=sign(q)*1e-3; betar=0;
        v01=0.5;
        v02=-0.1;
        delta=1/6.0;
        eps=0;
        tmax=400;
        vmax1=0.9;    vmin1=-0.9;
        vmax2=1.1; vmin2=-0.7;
    case ('weibels1')
        % Symmetric Streaming Weibel instability
        sigma1=0.1/sqrt(2);
        sigma2=sigma1;
        k=0.2;
        betai=sign(q)*1e-3; betar=0;
        v01=0.3;
        v02=-0.3;
        delta=0.5;
        eps=0;
        tmax=400;
        vmax1=0.9; vmin1=-0.9;
        vmax2=0.9; vmin2=-0.9;
    case('landau')
        % Strong Landau damping
        eps=0.5;
        k=0.5;
        sigma1=1;sigma2=1;
        betar=0;betai=0;v01=0;v02=0;delta=0;
        tmax=100;
        vmin1=-4.5;
        %         eps=0.5;
        %         tmax=30;
        vmin2=vmin1;
        vmax1=-vmin1; vmax2=-vmin2;
end

%% Parameters for dispersion relations
if (noise.X~=0 || noise.E2 || noise.B~=0)
    eps=0; betar=0; betai=0; 
    tmax=15;
    k=0.02;
end

L=2*pi/k;
f0=@(x,v1,v2) (1+ eps*cos(k*x)).*exp(-v1.^2./sigma1.^2/2 ).*...
    (delta*exp(-(v2-v01).^2/2/sigma2.^2) + ...
    (1-delta).*exp( -(v2-v02).^2/2/sigma2.^2))/sigma1/sigma2/(2*pi);


%%
% Splitting constants for composition methods. Set to scalar 1
% if only Strang splitting is desired
% comp_gamma=[1/(2-2^(1/3)),  -2^(1/3)/(2-2^(1/3)),  1/(2-2^(1/3))];
% comp_gamma=[pi, -2*pi, 1+pi];
switch(splitting)
    case 'strang'
        dtA=1/2;
        dtB=1/2;
    case '2nd4lie'
        % % second order (4 Lie)
        % % Second order splitting (smaller error constant than Strang)
        y2=(2*sqrt(326)-36)^(1/3);
        alphat=(y2^2+6*y2-2)/(12*y2);
        
        dtA=[alphat, 1/2-alphat];
        dtB=[1/2-alphat, alphat];
        
    case 'fourth10lie'
        dtB=[(146+5*sqrt(19))/540, (-2+10*sqrt(19))/135, 1/5, ...
            (-23-20*sqrt(19))/270, (14-sqrt(19))/108];
        dtA=flip(dtB);
    end
% scale time step
dtA=dtA*dt;dtB=dtB*dt;




%% Set up computational grids
dx=L/Nx;
dv1=(vmax1-vmin1)/Nv1;
dv2=(vmax2-vmin2)/Nv2;

x=0:dx:L-dx;
v1=vmin1:dv1:vmax1-dv1;
v2=vmin2:dv2:vmax2-dv2;

kx=fftshift(-Nx/2:1:(Nx/2-1))*2*pi/L;
kv1=fftshift(-Nv1/2:1:(Nv1/2-1))*2*pi/(vmax1-vmin1);
kv2=fftshift(-Nv2/2:1:(Nv2/2-1))*2*pi/(vmax2-vmin2);

[X, V1,V2]=ndgrid(x,v1,v2); % Spatial grid
[KX, KV1, KV2]=ndgrid(kx,kv1,kv2); % Spectral grid


%% Padding with zeros for visulization
% Ensures that a minimum resolution is used for plots
padN=ceil(128/min([Nx,Nv1,Nv2]));
x_pad=0:dx/padN:L-dx/padN;
v1_pad=vmin1:dv1/padN:vmax1-dv1/padN;
v2_pad=vmin2:dv2/padN:vmax2-dv2/padN;



%% Initial condition
F=f0(X,V1,V2);
if (noise.X~=0)
    F=bsxfun(@plus,F, noise.X*(rand(Nx,1,1)-0.5) );
end

%% Visualization
%  figure;
%  pcolor(v1,v2,squeeze(F(1,:,:)).')
%  xlabel('v_1'); ylabel('v_2');
% colorbar
% shading interp;
% drawnow;
figF=figure();
xslice=x(end);
yslice=0;
zslice=sort([v01,v02]);
clim=[0,max(F(:))*1.1];

figFv12=figure(); %slice plot at x=L/2

set(0,'DefaultLegendFontSize',16,...
    'DefaultLegendFontSizeMode','manual',...
    'DefaultAxesFontSize', 14)
set(0,'DefaultLineMarkerSize', 10);
set(0,'DefaultLineLinewidth', 1);
%prefix=['./plots/',testcase,'_'];
%prefix=['/tmp/',testcase,'_'];
mkdir(sprintf('./plots/%d',N))
prefix=sprintf('./plots/%d/%s_',N,testcase);


% Diagnostic
Nt=ceil(tmax/dt);
if (mod(dt2,dt)~=0)
   fprintf('Warning: dt2 has to be a multiple of dt\n')
   dt2=dt2-mod(dt2,dt);
end
Nt2=ceil(tmax/dt2);
Epot=zeros(Nt,2);Ekin=Epot;Bpot=zeros(Nt,1);Momentum=zeros(Nt,2);
Energy=zeros(Nt,1);

fields=struct('E1',zeros(Nx,Nt2),'E2',zeros(Nx,Nt2),'B',zeros(Nx,Nt2));

%% Initialize Fields
E1=zeros(Nx,1); E2=zeros(Nx,1); B=zeros(Nx,1);
%Initial Poisson solve
% Note that the FFT has to be normalized by Nx
RHO=q*sum(sum(fft(F,[],1),3),2)*dv1*dv2;
E1=RHO./(1j*KX(:,1,1));
E1(KX(:,1,1)==0)=0;


B(kx==k)=(betar + 1j*betai)*Nx; %Normalization of FFT by Nx
B(kx==-k)=(betar - 1j*betai)*Nx;
if (noise.B~=0)
   B=noise.B*fft((rand(Nx,1)-0.5));
   B(1)=0;
end
if (noise.E2~=0)
   E2=noise.E2*fft((rand(Nx,1)-0.5));
   E2(1)=0;
end


%% Precalculate advection operator in X
% XSHIFT=exp(-dt/2*V1.*(1j.*KX));
% XSHIFT=exp(-dt/2*V1(:,:,1).*(1j.*KX(:,:,1)));
XSHIFTA=zeros(Nx,Nv1,length(dtA));
XSHIFTB=zeros(Nx,Nv1,length(dtB));
for sdx=1:length(dtA)
    XSHIFTA(:,:,sdx)=exp(-dtA(sdx)*V1(:,:,1).*(1j.*KX(:,:,1)));
    XSHIFTB(:,:,sdx)=exp(-dtB(sdx)*V1(:,:,1).*(1j.*KX(:,:,1)));
end



% % Plot slice
% hslice = surf(linspace(0,L,100),...
%     linspace(vmin1,vmax1,100),...
%     zeros(100)+vmin2);
% % rotate(hslice,[-1,0,0],45)
% rotate(hslice,[0,1,0],atand(abs(vmin2)/L) )
% sxd = get(hslice,'XData');
% sv1d = get(hslice,'YData');
% sv2d = get(hslice,'ZData');
% shading interp
% close all;

% Reccurrence time
rtime=2*pi/(k*max([dv1,dv2,dx]));
if (rtime<=tmax)
    fprintf('WARNING: Simulation time exceeds recurrence time. \n');
    if(eps<0.2)
        
    end
end
fprintf('Recurrence occurs at t=%g\n', rtime);



tic;
% Initial transform
F=fft(F,[],2);
F=fft(F,[],3);

tdx2=0;
for tdx=1:Nt
    
    % Composition steps
    for sdx=1:length(dtA)
        
        %F=fft(F,[],2);
        %F=fft(F,[],3);
        if (sdx==1)
            U1=ifft(E1,'symmetric'); U1=reshape(U1,Nx,1,1);
            U2=ifft(E2,'symmetric'); U2=reshape(U2,Nx,1,1);
            F=F.*exp(-dtA(sdx)*q/m*1j*(U1.*KV1+U2.*KV2));
            B=B - dtA(sdx)*E2.*(1j*KX(:,1,1));
        end
        F=ifft(F,[],3);
        
        
        %% H_B
        E2 = E2 -dtA(sdx)*1j*(KX(:,1,1)).*B;
        
        %%H_p
        U3=ifft(B,'symmetric'); U3=reshape(real(U3),Nx,1,1);
        %F=fft(F,[],2);
        
        %% H_p2
        F=F.*exp(-dtA(sdx)*q/m*U3.*1j.*KV1.*V2);
        E2=E2 -dtA(sdx)*q*fft(sum(V2(:,1,:).*F(:,1,:),3))*dv2*dv1;
        F=ifft(F,[],2,'symmetric');
        
        %% H_p1
        %Magnetic field stays constant
        %U3=ifft(B,'symmetric'); U3=reshape(real(U3),Nx,1,1);
        F=fft(F,[],3);
        F=fft(F,[],1); %x-transform
        
        %H_p,A
        E1=E1+ q*sum(F(:,:,1).*(XSHIFTA(:,:,sdx)-1),2)*dv1*dv2./(1j.*KX(:,1,1)) ;
        F=F.*XSHIFTA(:,:,sdx); % advection in X
        F=ifft(F,[],1); %backtransform
        
        %H_p,L
        F=F.*exp((dtA(sdx)+dtB(sdx))*q/m*U3.*V1.*(1j.*KV2));
        
        %H_p,A
        F=fft(F,[],1); %x-transform
        E1=E1+ q*sum(F(:,:,1).*(XSHIFTB(:,:,sdx)-1),2)*dv1*dv2./(1j.*KX(:,1,1)) ;
        F=F.*XSHIFTB(:,:,sdx); % advection in X
        
        F=ifft(F,[],1); %backtransform
        F=ifft(F,[],3,'symmetric');
        
        
        %% H_p2
        F=fft(F,[],2);
        F=F.*exp(-dtB(sdx)*q/m*U3.*1j.*KV1.*V2);
        E2=E2 -dtB(sdx)*q*fft(sum(V2(:,1,:).*F(:,1,:),3))*dv2*dv1;
        
        %F=ifft(F,[],2,'symmetric');
        E2(kx==0)=0;
        E1(kx==0)=0;
        %         E2((Nx/2+1)+(-7:7) )=0;
        %B((Nx/2+1)+(-7:7))=0;
        %         E1((Nx/2+1)+(-7:7))=0;
        
        %% H_B
        E2 = E2 -dtB(sdx)*1j*(KX(:,1,1)).*B;
        
        %
        %% H_E
        %advection in V
        %F=fft(F,[],2);
        F=fft(F,[],3);
        
        U1=ifft(E1,'symmetric'); U1=reshape(U1,Nx,1,1);
        U2=ifft(E2,'symmetric'); U2=reshape(U2,Nx,1,1);
        if (sdx==length(dtA))
            F=F.*exp(-dtB(sdx)*q/m*1j*(U1.*KV1+U2.*KV2));
            B=B - dtB(sdx)*E2.*(1j*KX(:,1,1));
        else
            F=F.*exp(-(dtA(sdx+1)+dtB(sdx))*q/m*1j*(U1.*KV1+U2.*KV2));
            B=B - (dtA(sdx+1)+dtB(sdx))*E2.*(1j*KX(:,1,1));
        end
        
        %         F=ifft(F,[],3);
        %         F=ifft(F,[],2);
        
    end
    
    if filtration==true
        % Anti aliasing / Filamentation Filtration
        flx=floor((Nx/2)/3);
        flv1=floor((Nv1/2)/3);
        flv2=floor((Nv2/2)/3);
        F=fft(F,[],1);
        F(Nx/2+1-flx:Nx/2+1+flx,Nv1/2+1-flv1:Nv1/2+1+flv1,Nv2/2+1-flv2:Nv2/2+1+flv2)=0;
        F=ifft(F,[],1);
        B(Nx/2+1-flx:Nx/2+1+flx)=0;
        E1(Nx/2+1-flx:Nx/2+1+flx)=0;
        E2(Nx/2+1-flx:Nx/2+1+flx)=0;
    end
    
    %     %diffusion
    %     D=1e-5;
    %     F=F.*exp(-D*dt*1j*(KV1.^2+KV2.^2));
    
    
    % Diagnostic
    Epot(tdx,1) =(E1'*E1)*dx/Nx/2;
    Epot(tdx,2) =(E2'*E2)*dx/Nx/2;
    Bpot(tdx)=(B'*B)*dx/Nx/2;
    
    % F is Fourier transformed in v1 and v2
    Ekin(tdx,1) = squeeze(ifft(sum(F(:,:,1)*dv2,1)*dx,'symmetric'))*v1'.^2*dv1/2*m;
    Ekin(tdx,2) = squeeze(ifft(sum(F(:,1,:)*dv1,1)*dx,'symmetric')).'*v2'.^2*dv2/2*m;
    Momentum(tdx,1)=squeeze(ifft(sum(F(:,:,1)*dv2,1)*dx,'symmetric'))*v1'*dv1/2*m...
                                    + real( E2'*B )*dx/Nx;
    Momentum(tdx,2)=squeeze(ifft(sum(F(:,1,:)*dv1,1)*dx,'symmetric')).'*v2'*dv2/2*m...
                                     -real( E1'*B)*dx/Nx;
    Energy(tdx)=sum(Ekin(tdx,:))+sum(Epot(tdx,:))+sum(Bpot(tdx));
    
    if (mod((tdx-1)*dt,dt2)==0)
        tdx2=tdx2+1;
        fields.E1(:,tdx2)=E1/Nx; %ifft(E1,'symmetric');
        fields.E2(:,tdx2)=E2/Nx; %ifft(E2,'symmetric');
        fields.B(:,tdx2)=B/Nx; %ifft(B,'symmetric');
    end
    
    fprintf( 'E_pot=(%g,%g), %05.2f%%\n',Epot(tdx,:),tdx/Nt*100 )
    if  mod(tdx-1,floor(Nt/Nplot))==0
        
        % Pad for backtransform for higher resolution
        F2=fft(F,[],1);
        F2=fftshift(F2);
        
        Fpad=zeros(Nx*padN,Nv1*padN,Nv2*padN);
        Fpad( (-Nx/2:Nx/2-1)+Nx*padN/2+1, (-Nv1/2:Nv1/2-1)+Nv1*padN/2+1,...
            (-Nv2/2:Nv2/2-1)+Nv2*padN/2 +1 )=F2;
        Fpad=ifftshift(Fpad)*padN^3; 
        %if false
        F2=ifft(Fpad,[],1); 
        
        F2=ifft(F2,[],2);
        F2=ifft(F2,[],3,'symmetric');
        F2=abs(F2);
        
        
        figure(figF);
        h=slice(x_pad,v1_pad,v2_pad,...
            permute(abs(F2),[2,1,3]),xslice,yslice,zslice);       
        caxis(clim)
        colormap jet;
        shading interp;
        lighting gouraud;axis tight; caxis(clim)
        xlabel('x'); ylabel('v_1'); zlabel('v_2');
        title(sprintf('t=%07.2f', (tdx-1)*dt))
        view(-45,30);
        set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
        alpha('color');
        alphamap('rampup')
        alphamap('increase',.03)
        % The following is a
        c = colorbar();
        % Manually flush the event queue and force MATLAB to render the colorbar
        % necessary on some versions
        drawnow;
        % Get the color data of the object that correponds to the colorbar
        cdata = c.Face.Texture.CData;
        % Change the 4th channel (alpha channel) to the given alphamap from (0-255)
        cdata(end,:) = uint8(floor(alphamap()*255));
        % Ensure that the display respects the alpha channel
        c.Face.Texture.ColorType = 'truecoloralpha';
        % Update the color data with the new transparency information
        c.Face.Texture.CData = cdata;
        drawnow;
        print(figF,'-dpng',sprintf('%sdensity3d_%07d.png',prefix,tdx))         
        
        figure(figFv12) %slice plot at x=0
        pcolor(v1_pad,v2_pad,squeeze(F2(1,:,:)).')
        shading interp;  colormap jet; colorbar;
        caxis(clim)
        xlabel('v_1'); ylabel('v_2')
        drawnow;
        title(sprintf('x=0, t=%07.2f', (tdx-1)*dt))
        
        print(figFv12,'-dpng',sprintf('%sdensity2d_%07d.png',prefix,tdx))
        
    end
    
end
runtime=toc;
time=(0:dt:tmax-dt);


%% Test Gauss Law
F=ifft(F,[],3);
F=ifft(F,[],2,'symmetric');


RHO=q*fft(sum(sum(F,3),2),[],1)*dv1*dv2;
E1_=RHO./(1j*kx');
E1_(1)=0; % constant background
E1_(Nx/2+1)=E1(Nx/2+1); %get rid off the highest mode and its anoying
% assymmetry

if filtration==true
    E1_(Nx/2+1-flx:Nx/2+1+flx)=0;
end

gauss_err=max(abs(E1(2:end)-E1_(2:end)));

%% Discussion
% Save all results to file for later
save([prefix,'result.mat']);


% clear all; close all;
% load('plots/128/landau_result.mat'); close all;
% plots/128/weibels1_result.mat
% plots/128/weibel_result.mat  plots/128/weibels2_result.mat


%% Gauss Error
fprintf('Relative Gauss error: %g\n', gauss_err/min(abs(E1(2:end))));
fprintf('Gauss error: %02.3g\n', gauss_err);

%% Field Energys
time=0:dt:tmax-dt;
figure;
semilogy(time,Epot(:,2),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_2||^2$');
hold on;
semilogy(time,Epot(:,1),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_1||^2$');
semilogy(time,Bpot,'Linewidth',1, 'DisplayName','$\frac{1}{2}||B_3||^2$')
xlabel('time'); grid on; axis tight;
axis tight;
title('electromagnetic energy');
axis([time(1),time(end),1e-12,1.1*max([Epot(:);Bpot(:)])])
switch(testcase)
    case 'weibel'
        omega=0.02784;
        timelin=time(time<250);
        plot(timelin,0.5*Bpot(1)*exp(2*omega*(timelin-30)),'--','Linewidth',2,...
            'Color',[0,0,0,0.7],'DisplayName','lin. analysis');
        
    case {'weibels','weibels2'}
        timelin=time(time<100);
        omega=0.03;
        plot(timelin,1e-6*exp(2*omega*(timelin)),'--','Linewidth',2,...
            'Color',[0,0,0,0.7],'DisplayName','lin. analysis');
        axis([time(1),time(end),1e-10,max([Epot(:);Bpot(:)])]);
        
    case 'landau'
        if (eps<0.1)
            omega=1.415661888604536 - 0.153359466909605i;
            plot(time,0.5*Epot(1,1)*abs(exp(-1j*omega*(time-0.4))).^2);
        end
        
        
end
l=legend('-Dynamiclegend','Location','SouthEast');
set(l,'Interpreter','LaTeX','FontSize',16)
print('-dpng',[prefix, 'fieldenergy.png'])


%% Kinetic energy
time=0:dt:tmax-dt;
figure;
plot(time,Ekin,'-','Linewidth',2); hold on;
grid on;
xlabel('time'); grid on; axis tight;
% axis tight;
title('kinetic energy');
l=legend('$ \mathcal{H}_{p_1}$','$ \mathcal{H}_{p_2}$',...
    'Location','SouthEast');
set(l,'Interpreter','LaTeX','FontSize',16)
print('-dpng',[prefix, 'kineticenergy.png'])



%% Energy Error
figure;
Energy=sum(Epot,2)+Bpot+sum(Ekin,2);

% %
% plot(time,sum(Ekin,2)); hold on;
% plot(time,sum(Ekin,2)+Bpot+Epot(:,1))

semilogy(time, abs((Energy-Energy(1)))/Energy(1),'Linewidth',1) ;hold on;
grid on; axis tight;
title('relative energy error');
xlabel('time');
ylabel('$|\mathcal{H}(t)-\mathcal{H}(0)|/\mathcal{H}(0)$',...
    'interpreter','latex')
set(gca,'Ytick',10.^(-16:1))
print('-dpng',[prefix, 'energy_error.png'])


%% Momentum Error
figure;
semilogy(time, abs(Momentum-Momentum(2,:)))
legend('P_1','P_2','Location','SouthEast');
xlabel('time');
title('absolute momentum error');
axis tight; grid on;
set(gca,'Ytick',10.^(-16:2:1))
print('-dpng',[prefix, 'momentum_error.png'])



%% Dispersion Relation by post-processing
time2=0:dt2:tmax-dt2;
omega2=(-Nt2/2:Nt2/2-1)*2*pi/tmax;
kx2=fftshift(kx);

%% Spatial field harmonics
figure('Name', 'Spatial harmonics of E_1');
semilogy(time2, abs(fields.E1)); grid on;
xlabel('time'); ylabel('$|\hat{E}_1(t)|$','interpreter','latex')
figure('Name', 'Spatial harmonics of E_2');
semilogy(time2, abs(fields.E2)); grid on;
xlabel('time'); ylabel('$|\hat{E}_2(t)|$','interpreter','latex')
figure('Name', 'Spatial harmonics of B');
semilogy(time2, abs(fields.B)); grid on;
xlabel('time'); ylabel('$|\hat{B}(t)|$','interpreter','latex')



%% Matrix pencil
% The time interval can be changed in order to obtain a dispersion
% relation for the nonlinear state
tlimits=[0.1*tmax,tmax]; % Time interval    
Nr=6;    %Number of supsected branches
roots.E1=zeros(Nr,Nx/2-2);roots.E2=zeros(Nr,Nx/2-2);roots.B=zeros(Nr,Nx/2-2);
for kdx=2:Nx/2-1
    tdx=(time2>=tlimits(1) & time2<=tlimits(2));
    [om,D]=fMatPen(fields.E1(kdx,tdx).', Nr);
    roots.E1(:,kdx-1)=(-1j*D+om)/dt2;
    
    [om,D]=fMatPen(fields.E2(kdx,tdx).', Nr);
    roots.E2(:,kdx-1)=(-1j*D+om)/dt2;
    
    [om,D]=fMatPen(fields.B(kdx,tdx).', Nr);
    roots.B(:,kdx-1)=(-1j*D+om)/dt2;
end
roots.KX=repmat(kx(2:Nx/2-1),Nr,1);



figure();
scatter( roots.KX(:), real(roots.E1(:)),50,...
           imag(roots.E1(:)),'filled'); 
cb=colorbar; ylabel(cb,'growth rate')
caxis([-max(abs(imag(roots.E1(:)))),max(abs(imag(roots.E1(:))))]);
colormap jet;
xlabel('k'); ylabel('frequency \omega');
title('Dispersion roots for E_1');
grid on;

figure();
scatter( roots.KX(:), real(roots.E2(:)),50,...
           imag(roots.E2(:)),'filled'); 
cb=colorbar; ylabel(cb,'growth rate')
caxis([-max(abs(imag(roots.E2(:)))),max(abs(imag(roots.E2(:))))]);
colormap jet;
xlabel('k'); ylabel('frequency \omega');
title('Dispersion roots for E_2');
grid on;

figure();
scatter( roots.KX(:), real(roots.B(:)),50,...
           imag(roots.B(:)),'filled'); 
cb=colorbar; ylabel(cb,'growth rate')
caxis([-max(abs(imag(roots.B(:)))),max(abs(imag(roots.B(:))))]);
colormap jet;
xlabel('k'); ylabel('frequency \omega');
title('Dispersion roots for B');
grid on;



%% Time and spatial field harmonics
fields2.E1=fftshift(fft(fftshift(fields.E1,1),[],2),2)/Nt2*dt2; 
fields2.E2=fftshift(fft(fftshift(fields.E2,1),[],2),2)/Nt2*dt2; 
fields2.B=fftshift(fft(fftshift(fields.B,1),[],2),2)/Nt2*dt2; 

figure('Name','Dispersion E_1')
pcolor(kx2,omega2,log(abs( fields2.E1)).' );
shading interp;
colormap hot; colorbar;
xlabel('k'); ylabel('\omega');
axis([kx2(Nx/2+2),kx2(end),-3,3])

figure('Name','Dispersion E_2')
pcolor(kx2,omega2,log(abs( fields2.E2)).' );
shading interp;
colormap hot; colorbar;
xlabel('k'); ylabel('\omega');
axis([kx2(Nx/2+2),kx2(end),-3,3])

figure('Name','Dispersion B')
pcolor(kx2,omega2,log(abs( fields2.B)).' );
shading interp;
colormap hot; colorbar;
xlabel('k'); ylabel('\omega');
axis([kx2(Nx/2+2),kx2(end),-3,3])












