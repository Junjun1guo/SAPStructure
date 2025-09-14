clear;clc;
waves =1:40;
delta = load(['dt.txt']);
T1= load(['T1.txt']);
Damp = 0.05;
for i = 1:1:36 %case
    T = T1(i,1);
    for j = 1:1:length(waves)
        wave = waves(j);
        wave2str = num2str(wave);
        acc = load([wave2str,'.acc']);
%         地震动输入的时间间隔
        dt = delta(wave);
        [sa,sv,sd] = spectrasa(dt,acc,T,Damp);
        SAw(i,j) = sa;
        SVw(i,j) = sv;
        SDw(i,j) = sd;
    end
end
SF=SAw/0.1;
for i = 1:36
    for j=1:40
        SFFinal((i-1)*40+j,1)=SF(i,j);
    end
end

% SA=SAw';
% SaTarget=SAw/0.1;
% [jc,w]= max(SAw,[],2);
% Tg=0.01*w;
% figure();
% for j = 1:1:length(waves)
%     plot(T,SAw(j,:)./9.8);
%     hold on;
% end

% for j = 1:1:length(waves)
%     plot(T,SVw(j,:));
%     hold on;
% end
% legend('25','11','35','34','12','17');
% % set(gca,'Xscale','log');
% set(gca, 'fontsize',12,'FontName','Monospaced','YAxisLocation','left');
% set(gcf,'unit','centimeters','position',[10 5 7.5 6]);
% saveas(gcf,['spec','-log','-all','.pdf']);
% 
% scale = zeros(40,1);
% cnt = 0;
% 
% for j = 1:1:length(waves)
%     wave = waves(j);
%     cnt = cnt+1;
%     scale(wave) = SAw(j,1)./9.8;
% end
% adjust = 0.1./scale;

