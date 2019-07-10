clear all; close all;
dt=0.05;
% %dt=0.05;

tmax=150;

q=-1;
m=1;
N=32;
Nx=N;
Nv1=N;
Nv2=N;
Nplot=4; % Number of plots over simulation time



filtration=false; % Filamentation filtration by 2/3 rule
% The Hotspots in this Code are the complex exponentials

%%
% Available test cases are:
%
% * landau - strong Landau damping
% * weibel - the Weibel instability
% * weibels - the streaming Weibel instability
%
testcase='weibel'; % 'weibel','landau','weibels'
Hp1='exact'; %split, exact


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
        %tmax=500;
        tmax=150;
        vmax1=2.5*max(sigma1); vmin1=-vmax1;
        vmax2=2.5*max(sigma2); vmin2=-vmax2;
        
    case ('weibels')
        % Streaming Weibel instability
        sigma1=0.1/sqrt(2);
        sigma2=sigma1;
        k=0.2;
        betai=sign(q)*1e-3; betar=0;
        v01=0.5;
        v02=-0.1;
        delta=1/6.0;
        eps=0;
        tmax=400;
        vmin1=-0.3;
        vmax1=0.3;
        vmin2=-0.5;
        vmax2=1.1;
        
        
    case('landau')
        % Strong Landau damping
        eps=0.5;
        k=0.5;
        sigma1=1;sigma2=1;
        betar=0;betai=0;v01=0;v02=0;delta=0;
        tmax=100;
        vmin1=-4.5; vmin2=vmin1;
        vmax1=-vmin1; vmax2=-vmin2;
end

L=2*pi/k;
f0=@(x,v1,v2) (1+ eps*cos(k*x)).*exp(-v1.^2./sigma1.^2/2 ).*...
    (delta*exp(-(v2-v01).^2/2/sigma2.^2) + ...
    (1-delta).*exp( -(v2-v02).^2/2/sigma2.^2))/sigma1/sigma2/(2*pi);


%%
% Splitting constants for composition methods. Set to scalar 1
% if only Strang splitting is desired
% comp_gamma=[1/(2-2^(1/3)),  -2^(1/3)/(2-2^(1/3)),  1/(2-2^(1/3))];
comp_gamma=1;
% comp_gamma=[pi, -2*pi, 1+pi];

%% Set up computational grids
dx=L/Nx;
dv1=(vmax1-vmin1)/Nv1;
dv2=(vmax2-vmin2)/Nv2;

x=0:dx:L-dx;
v1=vmin1:dv1:vmax1-dv1;
v2=vmin2:dv2:vmax2-dv2;

kx=fftshift(-Nx/2:1:Nx/2-1)*2*pi/L;
kv1=fftshift(-Nv1/2:1:Nv1/2-1)*2*pi/(vmax1-vmin1);
kv2=fftshift(-Nv2/2:1:Nv2/2-1)*2*pi/(vmax2-vmin2);

[X, V1,V2]=ndgrid(x,v1,v2); % Spatial grid
[KX, KV1, KV2]=ndgrid(kx,kv1,kv2); % Spectral grid

%% Initial condition
F=f0(X,V1,V2);
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
Epot=zeros(Nt,2);Ekin=Epot;Bpot=zeros(Nt,1);Momentum=zeros(Nt,2);
Energy=zeros(Nt,1);

%% Initialize Fields
E1=zeros(Nx,1); E2=zeros(Nx,1); B=zeros(Nx,1);
%Initial Poisson solve
% Note that the FFT has to be normalized by Nx
RHO=q*sum(sum(fft(F,[],1),3),2)*dv1*dv2;
E1=RHO./(1j*KX(:,1,1));
E1(KX(:,1,1)==0)=0;


B(kx==k)=(betar + 1j*betai)*Nx; %Normalization of FFT by Nx
B(kx==-k)=(betar - 1j*betai)*Nx;

%% Precalculate advection operator in X
XSHIFT=zeros(Nx,Nv1,length(comp_gamma));
for jdx=1:length(comp_gamma)
    switch(Hp1)
        case 'split'
            deltat=dt/2*comp_gamma(jdx);
        case 'exact'
            deltat=dt*comp_gamma(jdx);
    end
    XSHIFT(:,:,jdx)=exp(-deltat*V1(:,:,1).*(1j.*KX(:,:,1)));
end

%% Matrix DFT
if (Hp1=='exact')
    %FFT in x as matrix
    DFTx=fliplr(vander( exp(-1j*2*pi/Nx*(0:Nx-1).') ));
    IDFTx=fliplr(vander( exp(1j*2*pi/Nx*(0:Nx-1).') ))/Nx;
    %IDFTx=ifft( eye(Nx),[],1,'symmetric');
    DFx=IDFTx*diag(1j*kx)*DFTx;
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
for tdx=1:Nt
    
    
    % Composition steps of the second order Strang splitting
    for splitdx=1:length(comp_gamma)
        deltat=comp_gamma(splitdx)*dt;
        %% H_E
        %advection in V
        %      F=fft(F,[],2);
        %      F=fft(F,[],3);
        
        U1=ifft(E1,'symmetric'); U1=reshape(U1,Nx,1,1);
        U2=ifft(E2,'symmetric'); U2=reshape(U2,Nx,1,1);
        F=F.*exp(-deltat/2*q/m*1j*(U1.*KV1+U2.*KV2));
        
        B=B - deltat/2*E2.*(1j*KX(:,1,1));
        
        F=ifft(F,[],3);
        F=ifft(F,[],2,'symmetric');
        
        
        %% H_B
        E2 = E2 -deltat/2*1j*(KX(:,1,1)).*B;
        
        %% H_p2
        F=fft(F,[],2);
        U3=ifft(B,'symmetric'); U3=reshape(real(U3),Nx,1,1);
        F=F.*exp(-deltat/2*q/m*U3.*1j.*KV1.*V2);
        F=ifft(F,[],2,'symmetric');
        
        E2=E2 -deltat/2*q*fft(sum(sum(V2.*F,3),2))*dv2*dv1;
        E2(KX(:,1,1)==0)=0;
        
        
         %% H_p1
        switch(Hp1)
            case 'split' % Symmetric strang splitting
                %H_p,A
                F=fft(F,[],1); %x-transform
                E1=E1+ q*sum(sum(F,3).*(XSHIFT(:,:,splitdx)-1),2)*dv1*dv2./(1j.*KX(:,1,1)) ;
                E1(KX(:,1,1)==0)=0;
                F=F.*reshape(XSHIFT(:,:,splitdx),Nx,Nv1,1); % advection in X
                F=ifft(F,[],1,'symmetric'); %backtransform
                
                
                %H_p,L
                U3=ifft(B,'symmetric'); U3=reshape(real(U3),Nx,1,1);
                F=fft(F,[],3);
                F=F.*exp(deltat*q/m*U3.*V1.*(1j.*KV2));
                F=ifft(F,[],3,'symmetric');
                
                %H_p,A
                F=fft(F,[],1); %x-transform
                E1=E1+ q*sum(sum(F,3).*(XSHIFT(:,:,splitdx)-1),2)*dv1*dv2./(1j.*KX(:,1,1)) ;
                E1(KX(:,1,1)==0)=0;
                F=F.*reshape(XSHIFT(:,:,splitdx),Nx,Nv1,1); % advection in X
                F=ifft(F,[],1,'symmetric'); %backtransform
                
            case 'exact'               
                E1=E1+ q*sum(fft(sum(F,3),[],1).*...
                      (XSHIFT(:,:,splitdx)-1),2)*dv1*dv2./(1j.*KX(:,1,1)) ;
                E1(KX(:,1,1)==0)=0;
                
                %Advection
                U3=ifft(B,'symmetric');
                F=fft(F,[],3);
                for v2dx=1:Nv2
                    
                    
                    expA_=expm(-deltat*dv1*(DFx-q/m*diag(U3).*(1j.*kv2(v2dx)) ));
                    expB_=expm(+deltat*dv1*(DFx-q/m*diag(U3).*(1j.*kv2(v2dx)) ));
                    expA=expA_;expB=expB_;
                    for v1dx=1:Nv1/2-1
                    
                        
                    F(:,Nv1/2+1+v1dx,v2dx)=expA*F(:,Nv1/2+1+v1dx,v2dx);
                    expA=expA*expA_;
                    F(:,Nv1/2+1-v1dx,v2dx)=expB*F(:,Nv1/2+1-v1dx,v2dx);
                    expB=expB*expB_;
                    
                    
                    end
                    F(:,1,v2dx)=expB*F(:,1,v2dx);
                end
                F=ifft(F,[],3,'symmetric');
                %fprintf('t=%g\n', tdx*dt);
        end
        
        
   
        
        
        
        %% H_p2
        F=fft(F,[],2);
        U3=ifft(B,'symmetric'); U3=reshape(real(U3),Nx,1,1);
        F=F.*exp(-deltat/2*q/m*U3.*1j.*KV1.*V2);
        F=ifft(F,[],2,'symmetric');
        
        E2=E2 -deltat/2*q*fft(sum(sum(V2.*F,3),2))*dv2*dv1;
        E2(KX(:,1,1)==0)=0;
        
        %% H_B
        E2 = E2 -deltat/2*1j*(KX(:,1,1)).*B;
        
        %     RHO=q*sum(sum(fft(F,[],1),3),2)*dv1*dv2;
        %     gauss_err=max(abs(E1(2:end)-RHO(2:end)./(1j*KX(2:end,1,1))))/max(abs(E1(2:end)));
        %     fprintf('Gauss error: %g\n', gauss_err);
        %
        %% H_E
        %advection in V
        F=fft(F,[],2);
        F=fft(F,[],3);
        
        U1=ifft(E1,'symmetric'); U1=reshape(U1,Nx,1,1);
        U2=ifft(E2,'symmetric'); U2=reshape(U2,Nx,1,1);
        F=F.*exp(-deltat/2*q/m*1j*(U1.*KV1+U2.*KV2));
        B=B - deltat/2*E2.*(1j*KX(:,1,1));
        
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
    % Diagnostic
    Epot(tdx,1) =(E1'*E1)*dx/Nx/2;
    Epot(tdx,2) =(E2'*E2)*dx/Nx/2;
    Bpot(tdx)=(B'*B)*dx/Nx/2;
    
    
    F2=ifft(F,[],2);
    F2=ifft(F2,[],3,'symmetric');
    Ekin(tdx,1) = 0.5*m*sum(F2(:).*V1(:).^2)*dx*dv1*dv2;
    Ekin(tdx,2) = 0.5*m*sum(F2(:).*V2(:).^2)*dx*dv1*dv2;
    Momentum(tdx,1)=m*sum(F2(:).*V1(:))*dx*dv1*dv2 + real( E2'*B )*dx/Nx;
    Momentum(tdx,2)=m*sum(F2(:).*V2(:))*dx*dv1*dv2 -real( E1'*B)*dx/Nx;
   
    
    Energy(tdx)=sum(Ekin(tdx,:))+sum(Epot(tdx,:))+sum(Bpot(tdx));
    
    % F is Fourier transformed in v1 and v2
    %     Ekin(1,tdx) = squeeze(ifft(sum(F(:,:,1)*dv2,1)*dx,'symmetric'))*v1'.^2*dv1/2;
    %     Ekin(2,tdx) = squeeze(ifft(sum(F(:,1,:)*dv1,1)*dx,'symmetric')).'*v2'.^2*dv2/2
    
    %     size(sum(F(:,1,:)*dv1,1)*dx)
    
    %     Epot(1,tdx) =(E1'*E1)*(2*pi/L)*pi/2/Nx^2;
    %     Epot(2,tdx) =(E2'*E2)*(2*pi/L)*pi/2/Nx^2;
    %     Bpot(tdx)=(B'*B)*(2*pi/L)*pi/2/Nx^2;
    
    
    
    
%         if (mod(tdx,10)==0)
%     
%             surf(squeeze(X(1,:,:)), squeeze(V1(1,:,:)), ...
%                 squeeze(V2(1,:,:)), squeeze(F(1,:,:)))
%                 view(140,40)
%     
%            hold on;
%             surf(squeeze(X(:,1,:)), squeeze(V1(:,1,:)), ...
%                 squeeze(V2(:,1,:)), squeeze(F(:,1,:)))
%               hold on;
%             surf(squeeze(X(:,:,Nv2/2)), squeeze(V1(:,:,Nv2/2-1)), ...
%                 squeeze(V2(:,:,Nv2/2)), squeeze(F(:,:,Nv2/2-1)))
%     %
%     %         surf(squeeze(X(:,1,:)), squeeze(V1(:,1,:)), ...
%     %             squeeze(V2(:,1,:)), squeeze(F(:,1,:)))
%             xlabel('x'); ylabel('v_1'); zlabel('v_2');
%     
%     %         shading interp;
%     
%     
%     %         hx = slice(permute(X,[2 1 3]),...
%     %             permute(V1,[2 1 3]),...
%     %             permute(V2,[2 1 3]),permute(F,[2 1 3]), sxd,sv1d,sv2d);
%     %         hx.FaceColor = 'interp';
%     %         hx.EdgeColor = 'none';
%     %         axis tight;
%               shading interp;   axis tight
%              hold off;             colorbar;
%     
%                  drawnow;
%     
%         end
    
    
    % hv1 = slice(x,y,z,v,[],max(y(:)),[]);
    % hv1.FaceColor = 'interp';
    % hv1.EdgeColor = 'none';
    %
    % hv2 = slice(x,y,z,v,[],[],zmin);
    % hv2.FaceColor = 'interp';
    % hv2.EdgeColor = 'none';
    
    %     surf(X(:,:,1),V1(:,:,1)
    
    %     pcolor( X(:,:,1),V1(:,:,1),sum(F,3)/Nv2*(vmax2-vmin2) ); shading interp;
    %
    %
    % F2=ifft(F,[],2,'symmetric');
    % F2=real(ifft(F2,[],3,'symmetric'));
    %
    % pcolor(v1,v2,squeeze(real(F(1,:,:))).')
    % xlabel('v_1'); ylabel('v_2');
    % colorbar
    % shading interp;
    % drawnow;
    
    
    if mod(tdx-1,floor(Nt/Nplot))==0
        %F2=ifft(F,[],2,'symmetric');
        %F2=abs(ifft(F2,[],3,'symmetric'));
        
        
        
        figure(figF);
        h=slice(permute(X,[2,1,3]),permute(V1,[2,1,3]),permute(V2,[2,1,3]),...
            permute(abs(F2),[2,1,3]),xslice,yslice,zslice);
        %         caxis(clim)
        colormap jet;
        shading interp;
        lighting gouraud;axis tight; caxis(clim)
        xlabel('x'); ylabel('v_1'); zlabel('v_2');
        title(sprintf('t=%f', (tdx-1)*dt))
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
        print('-dpng',sprintf('%sdensity3d_%07d.png',prefix,tdx))
        
        
        
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
fprintf('Relative Gauss error: %g\n', gauss_err/min(abs(E1(2:end))));
fprintf('Gauss error: %02.3g\n', gauss_err);

%% Discussion

%Save a file for later
% save([prefix,'result.mat']);


%% Field Energys
time=0:dt:tmax-dt;
figure;
semilogy(time,Epot(:,2),'-.','Linewidth',1,'DisplayName','$\frac{1}{2}||E_2||^2$');
hold on;
semilogy(time,Epot(:,1),'-.','Linewidth',1,'DisplayName','$\frac{1}{2}||E_1||^2$');
semilogy(time,Bpot,'Linewidth',1, 'DisplayName','$\frac{1}{2}||B_3||^2$')
xlabel('time'); grid on; axis tight;
axis tight;
title('electromagnetic energy');
% axis([time(1),time(end),1e-10,max([Epot(:);Bpot(:)])])
switch(testcase)
    case 'weibel'
        omega=0.02784;
        timelin=time(time<250);
        plot(timelin,0.5*Bpot(1)*exp(2*omega*(timelin-30)),'--','Linewidth',2,...
            'Color',[0,0,0,0.7],'DisplayName','lin. analysis');
        
    case 'weibels'
        timelin=time(time<100);
        omega=0.03;
        plot(timelin,1e-6*exp(2*omega*(timelin)),'--','Linewidth',2,...
            'Color',[0,0,0,0.7],'DisplayName','lin. analysis');
       axis([time(1),time(end),1e-10,max([Epot(:);Bpot(:)])]); 

    case 'landau'
        if (eps<0.1)
            % omega=1.415661888604536 - 0.153359466909605i;
            % plot(time,0.5*fieldenergy(1)*abs(exp(-1j*omega*(time-0.4))).^2);
        end
        
        
end
l=legend('-Dynamiclegend','Location','SouthEast');
set(l,'Interpreter','LaTeX','FontSize',16)
print('-dpng',[prefix, 'fieldenergy.png'])


%% Kinetic energy
time=0:dt:tmax-dt;
figure;
plot(time,Ekin,'-.','Linewidth',2); hold on;
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
semilogy(time, abs(Momentum-Momentum(1,:)))
legend('P_1','P_2','Location','SouthEast');
xlabel('time');
title('absolute momentum error');
axis tight; grid on;
set(gca,'Ytick',10.^(-16:2:1))
 print('-dpng',[prefix, 'momentum_error.png'])


%
%
% figure();
% xslice=x(end);
% yslice=0;
% zslice=[v01,v02];
% slice(permute(X,[2,1,3]),permute(V1,[2,1,3]),permute(V2,[2,1,3]),...
%          permute(real(F),[2,1,3]),xslice,yslice,zslice);
%  alphamap('rampup');
%
% shading interp;
% lighting gouraud;
% axis tight;
% colormap jet;
% xlabel('x'); ylabel('v_1'); zlabel('v_2');


%
% figure();
% xslice=x(end);
% yslice=0;
% zslice=[v01,0,v02];
% h=slice(permute(X,[2,1,3]),permute(V1,[2,1,3]),permute(V2,[2,1,3]),...
%          permute(abs(F),[2,1,3]),xslice,yslice,zslice);
%
% colormap jet;
%
% shading interp;
% lighting gouraud;
% axis tight;
% xlabel('x'); ylabel('v_1'); zlabel('v_2');
%
% view(-50,40);
% set(h,'EdgeColor','none',...
% 'FaceColor','interp',...
% 'FaceAlpha','interp')
% alpha('color');
% alphamap('rampup')
% alphamap('increase',.05)
% colormap jet;
% % The following is a
% c = colorbar();
% % Manually flush the event queue and force MATLAB to render the colorbar
% % necessary on some versions
% drawnow;
% % Get the color data of the object that correponds to the colorbar
% cdata = c.Face.Texture.CData;
% % Change the 4th channel (alpha channel) to the given alphamap from (0-255)
% cdata(end,:) = uint8(floor(alphamap()*255));
% % Ensure that the display respects the alpha channel
% c.Face.Texture.ColorType = 'truecoloralpha';
% % Update the color data with the new transparency information
% c.Face.Texture.CData = cdata;


%
% cb=colorbar
%
% annotation('textbox',...
%     cb.Position,...
%     'FitBoxToText','off',...
%     'FaceAlpha',0.5,...
%     'EdgeColor',[1 1 1],...
%     'BackgroundColor',[1 1 1]);
% colormap parula;
% h=colorbar;
% facealpha(h)
% h=colorbar;
% cb = findobj(gca,'Type','axes','Tag','Colorbar');
% alpha(h,0.1)
% alpha(cbIm,0.1)

% alphamap('decrease',.1)


%
%
%
% figure();
% scatter3(X(:),V1(:),V2(:),abs(F(:)))
%
% alphamap('increase',abs(min(F(:)))*4)
% figure;
% alphamap('vdown')
% am = get(gcf,'Alphamap');
% plot(am)

% pcolor(squeeze(real(F(end/2,:,:))))
% contourslice(X,V1,V2,real(F),1:4,[],0)
%
% contourslice(real(F(1:8:end,:,:)),x,v1,v2,8)
% xlabel('x')
% figure()
% vol3d('cdata',real(F),'texture','3D');
% view(45,15);
% % daspect(1./voxel_size2);
% axis tight;axis off;
%
% camlight; camlight(-90,-10); camlight(180,-10);lighting phong;
% alphamap('rampup');
% alphamap(0.1 .* alphamap);
%  shading interp;

%
%  xslice = L/2;
% yslice = (vmax1-vmin1)/2+vmin1;
% zslice = (vmax2-vmin2)/2+vmin2;

% xslice=linspace(0,x(end),4);
% yslice=[];
% zslice=[];
% slice(permute(X,[2,1,3]),permute(V1,[2,1,3]),permute(V2,[2,1,3]),...
%         permute(real(F),[2,1,3]),xslice,yslice,zslice);
% slice(permute(X(:,:,1:end/2),[2,1,3]),permute(V1(:,:,1:end/2),[2,1,3]),...
%        permute(V2(:,:,1:end/2),[2,1,3]),...
%         permute(real(F(:,:,1:end/2)),[2,1,3]),xslice,yslice,zslice); hold on;
% slice(permute(X(:, end/2:end ,end/2:end) ,[2,1,3]),permute(V1(:,end/2:end ,end/2:end),[2,1,3]),...
%           permute(V2(:, end/2:end ,end/2:end),[2,1,3]),...
%         permute(real(F(:, end/2:end ,end/2:end)),[2,1,3]),xslice,yslice,zslice);
% slice(permute(X,[2,1,3]),permute(V1,[2,1,3]),permute(V2,[2,1,3]),...
%         permute(real(F),[2,1,3]),x(end),[],[]);
%
%
% shading interp;
% lighting phong;


%
%
% DATA=real(F);




%
%
% [x,y,z,v] = flow;
%
% figure;
% xslice = 5;
% yslice = 0;
% zslice = 0;
% slice(x,y,z,v,xslice,yslice,zslice);
%
%
% view(3);
% axis on;
% grid on;
% light;
% lighting phong;
% camlight('left');
% shading interp;
%
%


%
%
%  figure;
%  data=real(F)
% patch(isocaps(X,V1,V2,data,.5),...
%    'FaceColor','interp','EdgeColor','none');
% p1 = patch(isosurface(X,V1,V2,data,.7),...
%    'FaceColor','blue','EdgeColor','none');
% % isonormals(data,p1);
% view(3);
% axis vis3d tight
% camlight left
% colormap('parula');
% lighting gouraud
% camlight;
%  camlight(-90,-10);
% % camlight(180,-10);

%
%
% semilogy(time,Epot.'); hold on;
% semilogy(time,Bpot.')

%
%
% %  semilogy(time,Epot(1,:).'); hold on;
%
%
%
% semilogy(time,Epot(1,:).'); hold on;
% semilogy(time,Bpot.')
%
% fieldenergy=Epot(1,:);
%
% % Include decay rate for linear landau damping WITHOUT COLLISIONS
% % if (eps<0.1 && kx==0.5)
% hold on;
% % Obtain zero of dispersion relation with dispersion_landau.m
% omega=1.415661888604536 - 0.153359466909605i;
% plot(time,0.5*fieldenergy(1)*abs(exp(-1j*omega*(time-0.4))).^2);
% % linear analysis with frequency
% plot(time,0.5*fieldenergy(1)*real(exp(-1j*omega*(time-0.4))).^2);
% legend('numerical', 'linear analysis', 'linear analysis');
%
% hold off;
% % end
