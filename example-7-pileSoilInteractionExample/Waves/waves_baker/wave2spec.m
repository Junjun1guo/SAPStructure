clear;
waves =1:40;
delta = load(['dt.txt']);
Damp = 0.05;
% 反应谱的间隔
T = 0.7;
[SAw, SVw, SDw] = deal(zeros(length(waves),length(T)));
for j = 1:1:length(waves)
    for i = 1:1:length(T)
        wave = waves(j);
        wave2str = num2str(wave);
%         if length(wave2str) == 1
%             wave2str = ['0',wave2str];
%         end
        acc = load([wave2str,'.acc']);
%         地震动输入的时间间隔
        dt = delta(wave);
        t = T(i);
%         单位是m/s^2
%         G = 9.8;
        [sa,sv,sd] = spectrasa(dt,acc,t,Damp);
        SAw(j,i) = sa;
        SVw(j,i) = sv;
        SDw(j,i) = sd;
    end
end
SA=SAw';
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

