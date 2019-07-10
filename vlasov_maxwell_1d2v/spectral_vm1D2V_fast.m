clear all;close all;


% Set number of threads to 4
maxNumCompThreads(4)
dt=0.05;      % 0.05 0.025 0.01 0.01 for (32 64 128 256)
tmax=150;

q=-1;
m=1;
N=64;
Nx=N;
Nv1=N;
Nv2=N;
Nplot=800; % Number of plots over simulation time 40, Inf for all
cmap=parula; %jet 


splitting='strang';    %strang 2nd4lie fourth10lie
% splitting='2nd4lie';    %strang 2nd4lie fourth10lie

filtration=false; % Filamentation filtration by 2/3 rule
% The Hotspots in this Code are the complex exponentials

%%
% Available test cases are:
%
% * landau - strong Landau damping
% * weibel - the Weibel instability
% * weibels - the streaming Weibel instability
%
testcase='weibels'; % 'weibel','weibels','landau'
Hp='split'; %split, exact, exact_split, split2
Hp_steps=1; %substeps


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
        %vmax1=0.8;    vmin1=-0.8;
        %vmax2=1.05; vmin2=-0.55;
        vmax1=0.3;    vmin1=-0.3;
        vmax2=1.05; vmin2=-0.55;
        %         vmax1=0.6;    vmin1=-vmax1;
        %         vmax2=0.8; vmin2=-0.6;
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

x=reshape(0:dx:L-dx,Nx,1,1);
v1=reshape(vmin1:dv1:vmax1-dv1,1,Nv1,1); 
v2=reshape(vmin2:dv2:vmax2-dv2,1,1,Nv2);

kx=reshape(fftshift(-Nx/2:1:(Nx/2-1))*2*pi/L,Nx,1,1);
kv1=reshape(fftshift(-Nv1/2:1:(Nv1/2-1))*2*pi/(vmax1-vmin1),1,Nv1,1);
kv2=reshape(fftshift(-Nv2/2:1:(Nv2/2-1))*2*pi/(vmax2-vmin2),1,1,Nv2);

[X, V1,V2]=ndgrid(x,v1,v2); % Spatial grid
% [KX, KV1, KV2]=ndgrid(kx,kv1,kv2); % Spectral grid



% Zero padded grid for plotting
padN=max(1,floor(256/Nx));
x_pad=0:dx/padN:L-dx/padN;
v1_pad=vmin1:dv1/padN:vmax1-dv1/padN;
v2_pad=vmin2:dv2/padN:vmax2-dv2/padN;



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

figFv12=figure(); %slice plot at x=L/2

set(0,'DefaultLegendFontSize',16,...
    'DefaultLegendFontSizeMode','manual',...
    'DefaultAxesFontSize', 14)
set(0,'DefaultLineMarkerSize', 10);
set(0,'DefaultLineLinewidth', 1);
%prefix=['./plots/',testcase,'_'];
%prefix=['/tmp/',testcase,'_'];
mkdir(sprintf('./plots/%d',N))
prefix=sprintf('./plots/%d/%s_%s',N,testcase,Hp);


% Diagnostic
Nt=ceil(tmax/dt); if isinf(Nplot), Nplot=Nt; end
Epot=zeros(Nt,2);Ekin=Epot;Bpot=zeros(Nt,1);Momentum=zeros(Nt,2);
Energy=zeros(Nt,1);

%% Initialize Fields
E1=zeros(Nx,1); E2=zeros(Nx,1); B=zeros(Nx,1);
%Initial Poisson solve
% Note that the FFT has to be normalized by Nx
RHO=q*sum(sum(fft(F,[],1),3),2)*dv1*dv2;
E1=RHO./(1j*kx);
E1(kx==0)=0;


B(kx==k)=(betar + 1j*betai)*Nx; %Normalization of FFT by Nx
B(kx==-k)=(betar - 1j*betai)*Nx;


%% Precalculate advection operator in X
switch(Hp)
    case {'split','split2'}
        XSHIFTA=zeros(Nx,Nv1,length(dtA));
        XSHIFTB=zeros(Nx,Nv1,length(dtB));
        for sdx=1:length(dtA)
            XSHIFTA(:,:,sdx)=exp(-dtA(sdx)*v1.*(1j.*kx));
            XSHIFTB(:,:,sdx)=exp(-dtB(sdx)*v1.*(1j.*kx));
        end
    case 'exact_split'
        XSHIFT=zeros(Nx,Nv1,length(dtA));
        XSHIFT05=zeros(Nx,Nv1,length(dtA));
        XSHIFT1=zeros(Nx,Nv1,length(dtA));
        for sdx=1:length(dtA)
            XSHIFT(:,:,sdx)=exp(-(dtA(sdx)+dtB(sdx))*v1.*(1j.*kx));
            XSHIFT05(:,:,sdx)=exp(-(dtA(sdx)+dtB(sdx))/Hp1_steps/2*v1.*(1j.*kx));
            XSHIFT1(:,:,sdx)=exp(-(dtA(sdx)+dtB(sdx))/Hp1_steps*v1.*(1j.*kx));
        end
    case {'exact','order8','order6'}
        XSHIFT=zeros(Nx,Nv1,length(dtA));
        for sdx=1:length(dtA)
            XSHIFT(:,:,sdx)=exp(-(dtA(sdx)+dtB(sdx))*v1.*(1j.*kx));
        end
end



%% Matrix DFT
if (strcmp(Hp,'exact'))
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
    
    for sdx=1:length(dtA)
        
        %F=fft(F,[],2);
        %F=fft(F,[],3);
        if (sdx==1)
            U1=ifft(E1,'symmetric'); U1=reshape(U1,Nx,1,1);
            U2=ifft(E2,'symmetric'); U2=reshape(U2,Nx,1,1);
            F=F.*exp(-dtA(sdx)*q/m*1j*(U1.*kv1+U2.*kv2));
            B=B - dtA(sdx)*E2.*(1j*kx);
        end
        F=ifft(F,[],3);
        
        
        %% H_B
        E2 = E2 -dtA(sdx)*1j*kx.*B;
        
        %%H_p
        U3=ifft(B,'symmetric'); U3=reshape(real(U3),Nx,1,1);
        %F=fft(F,[],2);
        
        if ~strcmp(Hp,'split2')
            %% H_p2
            F=F.*exp(-dtA(sdx)*q/m*U3.*1j.*kv1.*V2);
            E2=E2 -dtA(sdx)*q*fft(sum(V2(:,1,:).*F(:,1,:),3))*dv2*dv1;
        end
        F=ifft(F,[],2,'symmetric');
        
        %% H_p1
        switch(Hp)
            case 'split' % Symmetric strang splitting
                %Magnetic field stays constant
                U3=ifft(B,'symmetric'); U3=reshape(real(U3),Nx,1,1);
                F=fft(F,[],3);
                F=fft(F,[],1); %x-transform
                
                %H_p,A
                E1=E1+ q*sum(F(:,:,1).*(XSHIFTA(:,:,sdx)-1),2)*dv1*dv2./(1j.*kx);
                F=F.*XSHIFTA(:,:,sdx); % advection in X
                F=ifft(F,[],1); %backtransform
                
                %H_p,L
                F=F.*exp((dtA(sdx)+dtB(sdx))*q/m*U3.*v1.*(1j.*kv2));
                
                %H_p,A
                F=fft(F,[],1); %x-transform
                E1=E1+ q*sum(F(:,:,1).*(XSHIFTB(:,:,sdx)-1),2)*dv1*dv2./(1j.*kx);
                F=F.*XSHIFTB(:,:,sdx); % advection in X
                
                F=ifft(F,[],1); %backtransform
                F=ifft(F,[],3,'symmetric');
                
            case 'split2'
                %Rotation in Hp
                F=fft(F,[],1); %x-transform
                E1=E1+ q*sum(sum(F,3).*(XSHIFTA(:,:,sdx)-1),2)*dv1*dv2./(1j.*kx) ;
                Hp_OP=(1./(v1.*1j.*kx));
                Hp_OP(kx==0,:,:)=1;Hp_OP(:,v1==0,:)=1;
                Hp_OP=Hp_OP.*(XSHIFTA(:,:,sdx)-1);
                 E2=E2+ q*sum(sum(F(:,:,:).*v2,3).*Hp_OP,2)*dv1*dv2;                   
                F=F.*XSHIFTA(:,:,sdx); % advection in X
                F=ifft(F,[],1,'symmetric'); %backtransform
                
                % Rotation in B
                U3=ifft(B,'symmetric'); U3=reshape(real(U3),Nx,1,1);
                theta=q/m*U3*(dtA(sdx)+dtB(sdx));
                SHEAR1=exp(1j*(-tan(theta/2)).*reshape(v2,1,1,Nv2).*...
                    reshape(kv1,1,Nv1,1));
%                 SHEAR2=exp(1j*(sin(theta)).*reshape(v1,1,Nv1,1).*...
%                     reshape(kv2,1,1,Nv2));
%                 F=fft(F,[],2);
%                 F=F.*SHEAR1;
%                 F=ifft(fft(F,[],3),[],2);
%                 F=F.*SHEAR2;
%                 F=ifft(fft(F,[],2),[],3);
%                 F=F.*SHEAR1;
%                 F=ifft(F,[],2,'symmetric');
                
                F=fft(F,[],2);
                F=F.*SHEAR1;
                F=fft(ifft(F,[],2,'symmetric'),[],3);
                F=F.*exp(1j*(sin(theta)).*v1.*kv2);
                F=fft(ifft(F,[],3,'symmetric'),[],2);
                F=F.*SHEAR1;
                F=ifft(F,[],2,'symmetric');
                
                
                
                F=fft(F,[],1); %x-transform
                E1=E1+ q*sum(sum(F,3).*(XSHIFTA(:,:,sdx)-1),2)*dv1*dv2./(1j.*kx) ;
                Hp_OP=(1./(v1.*1j.*kx));
                Hp_OP(kx==0,:,:)=1;Hp_OP(:,v1==0,:)=1;
                Hp_OP=Hp_OP.*(XSHIFTB(:,:,sdx)-1);
                E2=E2+ q*sum(sum(F.*v2,3).*Hp_OP,2)*dv1*dv2;                   
                F=F.*XSHIFTB(:,:,sdx); % advection in X
                F=ifft(F,[],1,'symmetric'); %backtransform
                
                
            case 'exact_split'
                U3=ifft(B,'symmetric'); U3=reshape(real(U3),Nx,1,1);
                E1=E1+ q*sum(fft(sum(F,3),[],1).*...
                    (XSHIFT(:,:,sdx)-1),2)*dv1*dv2./(1j.*kx) ;
                E1(kx==0)=0;
                VSHIFT=exp((dtA(sdx)+dtB(sdx))/Hp1_steps*q/m*U3.*v1.*(1j.*kv2));
                F=fft(F,[],3);
                F=fft(F,[],1); %x-transform
                F=F.*XSHIFT05(:,:,sdx); % advection in X
                
                F=ifft(F,[],1); %backtransform
                for expdx=1:Hp1_steps
                    F=F.*VSHIFT;
                    if (expdx<Hp1_steps)
                        F=fft(F,[],1); %x-transform
                        F=F.*XSHIFT1(:,:,sdx); % advection in X
                        F=ifft(F,[],1); %backtransform
                    end
                end
                F=fft(F,[],1); %x-transform
                F=F.*XSHIFT05(:,:,sdx); % advection in X
                F=ifft(F,[],1); %backtransform
                F=ifft(F,[],3,'symmetric');
            case 'exact'
                U3=ifft(B,'symmetric');
                
                U3=reshape(real(U3),Nx,1,1);
                E1=E1+ q*sum(fft(sum(F,3),[],1).*...
                    (XSHIFT(:,:,sdx)-1),2)*dv1*dv2./(1j.*kx) ;
                E1(kx==0)=0;
                
                %Advection
                U3=ifft(B,'symmetric');
                F=fft(F,[],3);
                for v2dx=1:Nv2
                    
                    deltat=dtA(sdx)+dtB(sdx);
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
                
            case 'order8'
                U3=ifft(B,'symmetric');
                
                U3=reshape(real(U3),Nx,1,1);
                E1=E1+ q*sum(fft(sum(F,3),[],1).*...
                    (XSHIFT(:,:,sdx)-1),2)*dv1*dv2./(1j.*kx) ;
                E1(kx==0)=0;
                
                dtw=[25/194, 0.58151408710525096243,-0.41017537146985013753, ....
                    0.18514693571658773265,-0.40955234342085141934,0.14440594108001204106,...
                    0.27833550039367965131, 0.31495668391629485789];
                dtw=[dtw,1-2*sum(dtw),flip(dtw)]*(dtA(sdx)+dtB(sdx));
                
                F=fft(F,[],3);
                F=fft(F,[],1); %x-transform
                
                
                %H_p,A
                
                F=F.*exp(-dtw(1)/2*v1.*(1j.*kx)); % advection in X
                
                for sdx2=1:length(dtw)
                    
                    %H_p,L
                    F=ifft(F,[],1); %backtransform
                    F=F.*exp( dtw(sdx2)*q/m*U3.*v1.*(1j.*kv2));
                    F=fft(F,[],1); %x-transform
                    %H_p,A
                    if (sdx2<length(dtw))
                        F=F.*exp(-(dtw(sdx2)+dtw(sdx2+1) )/2*...
                            v1.*(1j.*kx));  % advection in X
                    end
                    
                end
                
                %H_p,A
                F=F.*exp(-dtw(end)/2*v1.*(1j.*kx)); % advection in X
                
                F=ifft(F,[],1); %backtransform
                F=ifft(F,[],3,'symmetric');
                
                
                
        end
        
        F=fft(F,[],2);
        if ~strcmp(Hp,'split2')
            %% H_p2
            F=F.*exp(-dtB(sdx)*q/m*U3.*1j.*kv1.*v2);
            E2=E2 -dtB(sdx)*q*fft(sum(V2(:,1,:).*F(:,1,:),3))*dv2*dv1;
        end
        
        %F=ifft(F,[],2,'symmetric');
        E2(kx==0)=0;
        E1(kx==0)=0;
        %         E2((Nx/2+1)+(-7:7) )=0;
        %B((Nx/2+1)+(-7:7))=0;
        %         E1((Nx/2+1)+(-7:7))=0;
        
        %% H_B
        E2 = E2 -dtB(sdx)*1j*(kx).*B;
        
        %
        %% H_E
        %advection in V
        %F=fft(F,[],2);
        F=fft(F,[],3);
        
        U1=ifft(E1,'symmetric'); U1=reshape(U1,Nx,1,1);
        U2=ifft(E2,'symmetric'); U2=reshape(U2,Nx,1,1);
        if (sdx==length(dtA))
            F=F.*exp(-dtB(sdx)*q/m*1j*(U1.*kv1+U2.*kv2));
            B=B - dtB(sdx)*E2.*(1j*kx);
        else
            F=F.*exp(-(dtA(sdx+1)+dtB(sdx))*q/m*1j*(U1.*kv1+U2.*kv2));
            B=B - (dtA(sdx+1)+dtB(sdx))*E2.*(1j*kx);
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
    %     F=F.*exp(-D*dt*1j*(Kv1.^2+kv2.^2));
    
    
    % Diagnostic
    Epot(tdx,1) =(E1'*E1)*dx/Nx/2;
    Epot(tdx,2) =(E2'*E2)*dx/Nx/2;
    Bpot(tdx)=(B'*B)*dx/Nx/2;
    
    % F is Fourier transformed in v1 and v2
    Ekin(tdx,1) = sum(ifft(sum(F(:,:,1)*dv2,1)*dx,'symmetric').*v1.^2)*dv1/2*m;
    Ekin(tdx,2) = sum(ifft(sum(F(:,1,:)*dv1,1)*dx,'symmetric').*v2.^2)*dv2/2*m;
    Momentum(tdx,1)=sum(ifft(sum(F(:,:,1)*dv2,1)*dx,'symmetric').*v1)*dv1/2*m...
        + real( E2'*B )*dx/Nx;
    Momentum(tdx,2)=sum(ifft(sum(F(:,1,:)*dv1,1)*dx,'symmetric').*v2)*dv2/2*m...
        -real( E1'*B)*dx/Nx;
    Energy(tdx)=sum(Ekin(tdx,:))+sum(Epot(tdx,:))+sum(Bpot(tdx));
    
    
    fprintf( 'E_pot=(%g,%g), %05.2f%%\n',Epot(tdx,:),tdx/Nt*100 )
    if  mod(tdx-1,floor(Nt/Nplot))==0
        
        %if false
        F2=ifft(F,[],2);
        F2=ifft(F2,[],3,'symmetric');
        
        %         F2=fft(F2,Nv2*padN,3);F2=fft(F2,Nv1*padN,2);F2=fft(F2,Nx*padN,1);
        %         F2=fft(F2,[],3);F2=fft(F2,[],2);F2=fft(F2,[],1);
        %         F2=ifft(F2,Nx*padN,1);F2=ifft(F2,Nv1*padN,2);
        %         F2=ifft(F2,Nv2*padN,3,'symmetric');
        
        
        F2=abs(F2);
        
        
        
        figure(figF);
        h=slice(permute(x,[2,1,3]),permute(v1,[2,1,3]),permute(v2,[2,1,3]),...
            permute(abs(F2),[2,1,3]),xslice,yslice,zslice);
        %         caxis(clim)
        colormap(cmap);
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
        pcolor(v1(:),v2(:),squeeze(F2(1,:,:)).')
        shading interp;  colormap(cmap); colorbar;
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
E1_=RHO./(1j*kx);
E1_(1)=0; % constant background
E1_(Nx/2+1)=E1(Nx/2+1); %get rid off the highest mode and its anoying
% assymmetry

if filtration==true
    E1_(Nx/2+1-flx:Nx/2+1+flx)=0;
end

gauss_err=max(abs(E1(2:end)-E1_(2:end)));

%% Discussion
% save(prefix
%Save a file for later
save([prefix,'result.mat']);


% clear all; close all;
% load('plots/128/landau_result.mat'); close all;
%  load('plots/32/weibel_result.mat'); close all;
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






stop
%% Splitting error
%F0=F;


% Exact reference value
U3=ifft(B,'symmetric');
DFTx=fliplr(vander( exp(-1j*2*pi/Nx*(0:Nx-1).') ));
IDFTx=fliplr(vander( exp(1j*2*pi/Nx*(0:Nx-1).') ))/Nx;
%IDFTx=ifft( eye(Nx),[],1,'symmetric');
DFx=IDFTx*diag(1j*kx)*DFTx;
tic;
for idx=1:5
    
    F_ref=fft(F0,[],3);
    for v2dx=1:Nv2
        for v1dx=1:Nv1
            F_ref(:,v1dx,v2dx)=expm(-dt*v1(v1dx)*(DFx-q/m*diag(U3).*(1j.*kv2(v2dx)) ))...
                *F_ref(:,v1dx,v2dx);
        end
    end
    F_ref=ifft(F_ref,[],3,'symmetric');
end
expm_time=toc/5;

tic;
for idx=1:20
    
    F=fft(F0,[],3);
    for v2dx=1:Nv2
        deltat=dt;
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
end
expm_time2=toc/20;

% norm(F(:)-F_ref(:))/norm(F_ref(:))

idx=1; clear spliterror_strang spliterror_lie spliterror_fourth10lie ...
    lie_time strang_time fourth10lie_time;
sub_steps_test=2.^(0:5);
for idx=1:length(sub_steps_test)
    sub_steps=sub_steps_test(idx);
    
    
    % Second order Strang splitting
    tic;
    VSHIFT=exp(dt/sub_steps*q/m*U3.*v1.*(1j.*kv2));
    XSHIFT1=exp(-(dt/sub_steps)*v1.*(1j.*kx));
    lie_time(idx)=toc;
    XSHIFT05=exp(-(dt/sub_steps/2)/2*v1.*(1j.*kx));
    
    F=fft(F0,[],3);
    F=fft(F,[],1); %x-transform
    F=F.*XSHIFT05; % advection in X
    F=ifft(F,[],1); %backtransform
    for expdx=1:sub_steps
        F=F.*VSHIFT;
        if (expdx<sub_steps)
            F=fft(F,[],1); %x-transform
            F=F.*XSHIFT1; % advection in X
            F=ifft(F,[],1); %backtransform
        end
    end
    F=fft(F,[],1); %x-transform
    F=F.*XSHIFT05; % advection in X
    F=ifft(F,[],1); %backtransform
    F=ifft(F,[],3,'symmetric');
    strang_time(idx)=toc;
    spliterror_strang(idx)=norm(F(:)-F_ref(:))/norm(F_ref(:));
    
    tic;
    % Lie splitting
    F=fft(F0,[],3);
    for expdx=1:sub_steps
        F=fft(F,[],1); %x-transform
        F=F.*XSHIFT1; % advection in X
        F=ifft(F,[],1); %backtransform
        F=F.*VSHIFT;
    end
    F=ifft(F,[],3,'symmetric');
    lie_time(idx)=lie_time(idx)+toc;
    
    spliterror_lie(idx)=norm(F(:)-F_ref(:))/norm(F_ref(:));
    
    
end


for idx=1:length(sub_steps_test)
    sub_steps=sub_steps_test(idx);
    tic;
    dtA=[(146+5*sqrt(19))/540, (-2+10*sqrt(19))/135, 1/5, ...
        (-23-20*sqrt(19))/270, (14-sqrt(19))/108]*dt/sub_steps;
    
    dtB=flip(dtA);
    
    dtB=flip(dtB); dtA=flip(dtA);
    XSHIFTA=zeros(Nx,Nv1,length(dtA));
    %XSHIFTB=zeros(Nx,Nv1,length(dtB));
    VSHIFT=zeros(Nx,Nv1,Nv2,length(dtB));
    for sdx=1:length(dtA)
        XSHIFTA(:,:,sdx)=exp(-dtA(sdx)*v1.*(1j.*kx));
        %XSHIFTB(:,:,sdx)=exp(-dtB(sdx)*v1.*(1j.*kx));
        VSHIFT(:,:,:,sdx)=exp((dtA(sdx)+dtB(sdx))*q/m*U3.*v1.*(1j.*kv2));
    end
    F=fft(F0,[],3);
    F=fft(F,[],1); %x-transform
    for expdx=1:sub_steps
        
        
        for sdx=1:length(dtA)
            %H_p,A
            F=F.*XSHIFTA(:,:,sdx); % advection in X
            F=ifft(F,[],1); %backtransform
            
            %H_p,L
            F=F.*VSHIFT(:,:,:,sdx);
            
            %H_p,A
            F=fft(F,[],1); %x-transform
            %F=F.*XSHIFTB(:,:,sdx); % advection in X
            F=F.*XSHIFTA(:,:,length(dtA)-sdx+1); % advection in X
            
        end
    end
    F=ifft(F,[],1); %backtransform
    F=ifft(F,[],3,'symmetric');
    fourth10lie_time(idx)=toc;
    
    spliterror_fourth10lie(idx)=norm(F(:)-F_ref(:))/norm(F_ref(:));
    tic;
end
% F=F0;

dtw=[25/194, 0.58151408710525096243,-0.41017537146985013753, ....
    0.18514693571658773265,-0.40955234342085141934,0.14440594108001204106,...
    0.27833550039367965131, 0.31495668391629485789];
dtw=[dtw,1-2*sum(dtw)];

clear spliterror_order8;
for idx=1:length(sub_steps_test)
    sub_steps=sub_steps_test(idx);
    tic;
    dtw=[25/194, 0.58151408710525096243,-0.41017537146985013753, ....
        0.18514693571658773265,-0.40955234342085141934,0.14440594108001204106,...
        0.27833550039367965131, 0.31495668391629485789];
    dtw=[dtw,1-2*sum(dtw)]*dt/sub_steps;
    %dtw=dt/sub_steps;
    XSHIFT=zeros(Nx,Nv1,length(dtw));
    VSHIFT=zeros(Nx,Nv1,Nv2,length(dtw));
    for sdx=1:length(dtw)
        XSHIFT(:,:,sdx)=exp(-dtw(sdx)/2*v1.*(1j.*kx));
        VSHIFT(:,:,:,sdx)=exp( dtw(sdx)*q/m*U3.*v1.*(1j.*kv2));
    end
    F=fft(F0,[],3);
    F=fft(F,[],1); %x-transform
    
    for expdx=1:sub_steps
        
        for sdx=[ 1:length(dtw), length(dtw)-1:-1:1]
            
            %H_p,A
            F=F.*XSHIFT(:,:,sdx); % advection in X
            F=ifft(F,[],1); %backtransform
            
            %H_p,L
            F=F.*VSHIFT(:,:,:,sdx);
            
            %H_p,A
            F=fft(F,[],1); %x-transform
            F=F.*XSHIFT(:,:,sdx); % advection in X
            
            
        end
        
    end
    F=ifft(F,[],1); %backtransform
    F=ifft(F,[],3,'symmetric');
    order8_time(idx)=toc;
    
    spliterror_order8(idx)=norm(F(:)-F_ref(:))/norm(F_ref(:));
end

clear spliterror_order6 order6_time;
for idx=1:length(sub_steps_test)
    sub_steps=sub_steps_test(idx);
    tic;
    dtw=[0.1867, 0.55549702371247839916, 0.12946694891347535806,-0.84326562338773460855];
    
    dtw=[dtw,1-2*sum(dtw)]*dt/sub_steps;
    %dtw=dt/sub_steps;
    XSHIFT=zeros(Nx,Nv1,length(dtw));
    VSHIFT=zeros(Nx,Nv1,Nv2,length(dtw));
    for sdx=1:length(dtw)
        XSHIFT(:,:,sdx)=exp(-dtw(sdx)/2*v1.*(1j.*kx));
        VSHIFT(:,:,:,sdx)=exp( dtw(sdx)*q/m*U3.*v1.*(1j.*kv2));
    end
    F=fft(F0,[],3);
    F=fft(F,[],1); %x-transform
    
    for expdx=1:sub_steps
        
        for sdx=[ 1:length(dtw), length(dtw)-1:-1:1]
            
            %H_p,A
            F=F.*XSHIFT(:,:,sdx); % advection in X
            F=ifft(F,[],1); %backtransform
            
            %H_p,L
            F=F.*VSHIFT(:,:,:,sdx);
            
            %H_p,A
            F=fft(F,[],1); %x-transform
            F=F.*XSHIFT(:,:,sdx); % advection in X
            
            
        end
        
    end
    F=ifft(F,[],1); %backtransform
    F=ifft(F,[],3,'symmetric');
    order6_time(idx)=toc;
    
    spliterror_order6(idx)=norm(F(:)-F_ref(:))/norm(F_ref(:));
end




set(0,'DefaultLineLinewidth', 2);
figure()
loglog(sub_steps_test,spliterror_lie); hold on;  grid on;
loglog(sub_steps_test,spliterror_strang,'-x')
loglog(sub_steps_test,spliterror_fourth10lie,'-o');
loglog(sub_steps_test,spliterror_order6,'-o');
loglog(sub_steps_test,spliterror_order8,'-o');

legend('lie','strang','fourth10lie','order6','order8','Location','SouthWest')
xlabel('sub steps');
ylabel('relative error in f');
axis tight;
print('-dpng',[prefix, 'Hp1_approx_error.png'])


figure()
loglog(sub_steps_test,lie_time); hold on; grid on;
plot(sub_steps_test,strang_time,'-x')
plot(sub_steps_test,fourth10lie_time,'-o')
loglog(sub_steps_test,order6_time,'-o');
loglog(sub_steps_test,order8_time,'-o');
plot(sub_steps_test,expm_time2*ones(length(sub_steps_test),1),'k--');
axis tight;
legend('lie','strang','fourth10lie','order6','order8','expm','Location','NorthWest')
xlabel('sub steps'); ylabel('wall time');
print('-dpng',[prefix, 'Hp1_approx_time.png'])



% Acurracy over time
figure()
semilogy(sub_steps_test,lie_time.*spliterror_lie); hold on;
plot(sub_steps_test,strang_time.*spliterror_strang); hold on;
plot(sub_steps_test,fourth10lie_time.*spliterror_fourth10lie); hold on;
plot(sub_steps_test,order6_time.*spliterror_order6); hold on;
plot(sub_steps_test,order8_time.*spliterror_order8); hold on;

plot(sub_steps_test,expm_time2*ones(length(sub_steps_test),1)*1e-14,'k--'); hold on;
set(gca,'Ytick',10.^(-16:2:1))
% print('-dpng',[prefix, 'momentum_error.png'])
% %





%
% ref=load('plots/128/weibels_result.mat');
% semilogy(ref.time,ref.Epot(:,2),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_2||^2$');
% hold on;
% semilogy(ref.time,ref.Epot(:,1),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_1||^2$');
% semilogy(ref.time,ref.Bpot,'Linewidth',1, 'DisplayName','$\frac{1}{2}||B_3||^2$')
%
%
% figure()
% semilogy(ref.time,abs(ref.Epot(:,2)-Epot(:,2))/norm(ref.Epot(:,2)),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_2||^2$'); hold on;
% semilogy(ref.time,abs(ref.Epot(:,2)-Epot(:,2))/norm(ref.Epot(:,2)),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_1||^2$');
% semilogy(ref.time,abs(ref.Bpot-Bpot)/norm(ref.Bpot),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||B_3||^2$');
%
% hold on;
% semilogy(ref.time,ref.Epot(:,1),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_1||^2$');
% semilogy(ref.time,ref.Bpot,'Linewidth',1, 'DisplayName','$\frac{1}{2}||B_3||^2$')
%
%
% figure()
% semilogy(time, abs((Energy-Energy(1)))/Energy(1),'Linewidth',1) ;hold on;
% semilogy(time, abs((ref.Energy-ref.Energy(1)))/ref.Energy(1),'Linewidth',1) ;hold on;

