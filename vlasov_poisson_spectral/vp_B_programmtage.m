
close all; clear all;
set(0,'DefaultLegendFontSize',16,...
    'DefaultLegendFontSizeMode','manual',...
    'DefaultAxesFontSize', 14)
set(0,'DefaultLineMarkerSize', 10);
set(0,'DefaultLineLinewidth', 2);

% prefix=@(method,B)sprintf('plots/VP2d2v_kelvinhelmholtz_strang_%s_B%d_16_32_result.mat',method,B);
prefix=@(method,B)sprintf('plots/VP2d2v_kelvinhelmholtz_strang_%s_B%d_32_32_result.mat',method,B);



mkdir('./programmtage');
pltdir='./programmtage/';

fig1=figure();
fig2=figure();

for B=[1,2,4,8,16,32]
    

scovel=load(prefix('scovel',B));
shear=load(prefix('shear',B));
std=load(prefix('std',B));
% imr=load(prefix('imr',B));

calc_momentum_error=@(data) sqrt(sum((data.Momentum-data.Momentum(1,:)).^2,2));

scovel.momentum_error=calc_momentum_error(scovel);
shear.momentum_error=calc_momentum_error(shear);
std.momentum_error=calc_momentum_error(std);

% scovel_imr=load(prefix('scovel_imr',B));

% scovel_imr=load(sprintf('plots/VP2d2v_kelvinhelmholtz_strang_scovel_imr_B%d_16_32_result.mat',B));

% Dispersion relation
k0=scovel.k0;
%gamma=sqrt(1-k0^2)*k0;
gamma=2*(1-k0)*k0;


% % semilogy(scovel.time/B,scovel.Epot(:,1),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_1||^2$'); hold on;
%  figure()
%  semilogy(std.time/B,std.energy_error,'--','DisplayName','splitting'); hold on;
%  semilogy(shear.time/B,shear.energy_error,'-.','DisplayName','exp.Boris'); hold on;
%  semilogy(scovel.time/B,scovel.energy_error,'-','DisplayName','Scovel'); hold on;
%  %semilogy(imr.time/B,imr.energy_error,'-','DisplayName','imr'); hold on;
%  xlabel('fluid time');  grid on;  %axis tight;
%  ylabel('relative energy error');
% legend('Location','SouthEast');
% print('-depsc', sprintf('%senergy_error_B%d.eps',pltdir,B))
% print('-dpng', sprintf('%senergy_error_B%d.png',pltdir,B))
%  title(sprintf('B_3=%d',B));
% 
% 
% %  figure()
% %  semilogy(std.time/B,std.kine,'-','DisplayName','splitting'); hold on;
% %  semilogy(shear.time/B,shear.momentum_error,'-','DisplayName','exp.Boris'); hold on;
% %  semilogy(scovel.time/B,scovel.momentum_error,'-','DisplayName','Scovel'); hold on;
% %  %semilogy(imr.time/B,imr.energy_error,'-','DisplayName','imr'); hold on;
% %  title(sprintf('B_3=%d',B));
% %  xlabel('fluid time');  grid on;  %axis tight;
% %  ylabel('relative energy error');
% 
%  figure()
%  semilogy(std.time/B,std.momentum_error,'--','DisplayName','splitting'); hold on;
%  semilogy(shear.time/B,shear.momentum_error,'-.','DisplayName','exp.Boris'); hold on;
%  semilogy(scovel.time/B,scovel.momentum_error,'-','DisplayName','Scovel'); hold on;
%  %semilogy(imr.time/B,imr.energy_error,'-','DisplayName','imr'); hold on;
%  xlabel('fluid time');  grid on;  %axis tight;
%  ylabel('momentum error');
% legend('Location','SouthEast');
% print('-depsc', sprintf('%smomentum_error_B%d.eps',pltdir,B))
% print('-dpng', sprintf('%smomentum_error_B%d.png',pltdir,B))
%  title(sprintf('B_3=%d',B));


figure()
% % semilogy(std.time/B,std.Epot(:,1),'--','DisplayName','splitting'); hold on;
% semilogy(shear.time/B,shear.Epot(:,1),'-.','DisplayName','exp.Boris'); hold on;
semilogy(scovel.time/B,scovel.Epot(:,1),'-','DisplayName','Scovel'); hold on;
 xlabel('fluid time');  grid on;  %axis tight;
 ylabel('$\frac{1}{2}||E_1||^2$','Interpreter','latex');
time=std.time; 
   time=time(time<20);
la=semilogy(time,      0.02.*exp(gamma*time),'-','Linewidth',3,'Color',0.8*[1,1,1],'DisplayName','lin.approx.');
uistack(la,'bottom');
legend('Location','SouthEast'); 




 
figure()
% semilogy(std.time/B,std.Epot(:,1),'--','DisplayName','splitting'); hold on;
semilogy(shear.time/B,shear.Epot(:,1),'-.','DisplayName','exp.Boris'); hold on;
semilogy(scovel.time/B,scovel.Epot(:,1),'-','DisplayName','Scovel'); hold on;
%    semilogy(imr.time/B,imr.Epot(:,1),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_1||^2$'); hold on;
% semilogy(scovel_imr.time/B,scovel_imr.Epot(:,1),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_1||^2$'); hold on;
% semilogy(std.time/B,std.Epot(:,2),'-','DisplayName','splitting'); hold on;
%   semilogy(shear.time/B,shear.Epot(:,2),'-.','DisplayName','exp.Boris'); hold on;
%   semilogy(scovel.time/B,scovel.Epot(:,2),'--','DisplayName','Scovel'); hold on;


 xlabel('fluid time');  grid on;  %axis tight;
 ylabel('$\frac{1}{2}||E_1||^2$','Interpreter','latex');
time=std.time; 
   time=time(time<20);
la=semilogy(time,      0.02.*exp(gamma*time),'-','Linewidth',3,'Color',0.8*[1,1,1],'DisplayName','lin.approx.');
uistack(la,'bottom');
legend('Location','SouthEast'); 

print('-depsc', sprintf('%sEpot1_B%d.eps',pltdir,B))
print('-dpng', sprintf('%sEpot1_B%d.png',pltdir,B))
title(sprintf('B_3=%d',B));

%  figure(fig1);
%   semilogy(scovel.time/B,scovel.Epot(:,1),'-','DisplayName',sprintf('B_3=%d',B)); hold on;
%  xlabel('fluid time');  grid on;  %axis tight;
%  ylabel('$\frac{1}{2}||E_1||^2$','Interpreter','latex');
%    legend('Location','SouthEast');
% 
% 
%  figure(fig2);
%   semilogy(scovel.time/B,scovel.energy_error,'-','DisplayName',sprintf('B_3=%d',B)); hold on;
%  xlabel('fluid time');  grid on;  %axis tight;
%  ylabel('relative energy error');
%    legend('Location','SouthEast');



% semilogy(shear.time/B,shear.Epot(:,1),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_1||^2$'); hold on;

% end
% 
% 
% figure(fig1);
% la=semilogy(time, 0.015.*exp(gamma*time),'-','Linewidth',5,...
%     'Color',0.7*[1,1,1],'DisplayName','lin.approx.');
% uistack(la,'bottom');
% lg=legend(); set(lg,'Fontsize',12);
% print('-depsc', sprintf('%sscovel_Epot1_B%d.eps',pltdir,B))
% print('-dpng', sprintf('%sscovel_Epot1_B%d.png',pltdir,B))
%  figure(fig2);
% print('-depsc', sprintf('%sscovel_energy_error_B%d.eps',pltdir,B))
% print('-dpng', sprintf('%sscovel_energy_error_B%d.png',pltdir,B))
% 
% 
% 
% 
% % Scovel 
% close all;
% 
% for B=[1,32]
%     figure()
% for method={'strang', '2nd4lie', 'fourth10lie'}
% 
%    
% prefix=@(method,B)sprintf('plots/VP2d2v_kelvinhelmholtz_%s_scovel_B%d_16_16_result.mat',method,B);
% 
% data=load(prefix(method{1},B));
%  semilogy(data.time(1:10:end)/B,data.energy_error(1:10:end),'-','DisplayName',method{1},'Linewidth',1); hold on;
% 
% 
% end
% grid on;
% end
% legend()
% 
% 
% 
% % 
% % 
% % 
% % semilogy(shear.time,shear.Epot(:,1),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_1||^2$'); hold on;
% % semilogy(scovel.time,scovel.Epot(:,1),'-','Linewidth',1,'DisplayName','$\frac{1}{2}||E_1||^2$');