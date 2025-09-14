clc
clear
% m controls the number of gms
% tic;
g=9.81;
dt=textread('dt.txt');
[m,n]=size(dt);
mkdir ('dispTH');
% SF=textread('SF_for_Sa10.txt');
for i=1:31

    % load acc history
    % acc=0;
    % acc=load([num2str(i),'.acc']);
    acc=load([num2str(i),'.acc']);    
    [accnum,l]=size(acc);
    t=dt(i)*accnum;
    gmtime=0:dt(i):t;
    
    % PGA(g)
    PGA(i,1)=max(abs(acc));
    disp(['PGA of ',num2str(i),' Wave is Obtain'])
    
    % PGV(cm/s)
    velo=0;
    for j=1:(accnum-1),
        velo(j+1)=(acc(j)+acc(j+1))/2*dt(i)*100*g+velo(j);
    end
    PGV(i,1)=max(abs(velo));
    disp(['PGV of ',num2str(i),' Wave is Obtain'])
%     
%     % PGD(cm)
%     displ=0;
%     for j=1:(accnum-1),
%         displ(j+1)=(velo(j)+velo(j+1))/2*dt(i)+displ(j);
%     end
%     PGD(i,1)=max(abs(displ));
%     disp(['PGD of ',num2str(i),' Wave is Obtain'])
%     dispTH=0;
%     dispTH=displ./100.; % unit m
%     fid = fopen(['dispTH\',num2str(i),'.disp'],'wt'); 
%     fprintf(fid,'%g\n',dispTH);   
%     fclose(fid);
%     disp(['dispTH of ',num2str(i),' Wave is Recorded'])
%  
%     % PGV/PGA (s)
%     VAratio(i,1)=PGV(i,1)/100./(PGA(i,1)*g);
%     disp(['PGV/PGA of ',num2str(i),' Wave is Obtain'])
%     
    %Arias Intensity (Ia) (m/sec)
    npyt=1;
    nwxw=1;
    npp=1;
    Ia(i,1)=sum((dt(i).*((acc*9.81).^2)))*pi/2/9.81;
    disp(['Ia of ',num2str(i),' Wave is Obtain'])
% 
%     % percentages of Ia corresponding time
%     IaD5=0;
%     IaD75=0;
%     IaD95=0;
%     
%     for j=1:accnum,
%         IaIa(j)=sum((dt(i).*((acc(1:j)*9.81).^2)))*pi/2/9.81;
%         % D5
%         if IaIa(j)>=0.05*Ia(i,1) & npyt==1,
%             npyt=npyt+1;
%             IaD5=j;
%         end
%         % D75
%         if IaIa(j)>=0.75*Ia(i,1) & nwxw==1,
%             nwxw=nwxw+1;
%             IaD75=j;
%         end
%         % D95
%         if IaIa(j)>=0.95*Ia(i,1) & npp==1,
%             npp=npp+1;
%             IaD95=j;
%         end
%     end
%     
%     %A95 parameter (g)
%     %accpx=sort(acc);
%     %for j=1:accnum,
%     %    Iapx(j)=sum((dt(i).*((accpx(1:j)*9.81).^2)))*pi/2/9.81;
%     %    if Iapx(j)>=0.95*Ia(i,1)
%     %        A95(i,1)=accpx(j);
%     %    end
%     %end
%     %disp(['A95 of ',num2str(i),' Wave is Obtain'])
%     
%     % significant duration (sec)
%     SD(i,1)=(IaD95-IaD5)*dt(i);
%     disp(['Significant Duration of ',num2str(i),' Wave is Obtain'])
%     
%     % RMS of acceleration (g)
%     SigniDurAcc=0;
%     SigniDurAcc=acc(IaD5:IaD95)';
%     Arms(i,1)=sqrt(sum(dt(i).*(SigniDurAcc.^2))/SD(i,1));
%     disp(['Arms of ',num2str(i),' Wave is Obtain'])
%     
%     % RMS of velocity  (cm/s)
%     SigniDurVel=0;
%     SigniDurVel=velo(IaD5:IaD95);    
%     Vrms(i,1)=sqrt(sum(dt(i).*(SigniDurVel.^2))/SD(i,1));
%     disp(['Vrms of ',num2str(i),' Wave is Obtain'])
%     
%     % RMS of displacement (cm)
%     SigniDurDisp=0;
%     SigniDurDisp=displ(IaD5:IaD95); 
%     Drms(i,1)=sqrt(sum(dt(i).*(SigniDurDisp.^2))/SD(i,1));
%     disp(['Drms of ',num2str(i),' Wave is Obtain'])
%     
%     % Characteristic Intensity (Ic)[g^1.5/S^0.5]
%     Ic(i,1)=sqrt(Arms(i).^3)*sqrt(SD(i,1));
%     disp(['Ic of ',num2str(i),' Wave is Obtain'])
% 
%     % Shaking intensity rate (m/s2)
%     IaD5D75(i,1)=(IaD75-IaD5)*dt(i);
%     SIR(i,1)=0.7*Ia(i,1)/IaD5D75(i,1);
%     disp(['SIR of ',num2str(i),' Wave is Obtain'])
% 
%     % Specific Energy Density (SED) (cm2/sec)
%     SED(i,1)=sum(dt(i)*(velo.^2));
%     disp(['SED of ',num2str(i),' Wave is Obtain'])
%     
    % Cumulative Absolute Velocity (CAV) (cm/sec)
    CAV(i,1)=sum(dt(i).*abs(acc(1:end).*g))*100;
    disp(['CAV of ',num2str(i),' Wave is Obtain'])
% 
%     % Cumulative Absolute Velocity with Threshold of 5 cm/s2 (CAV5) (cm/sec)
%     ABSacc=0;
%     ABSacc=abs(acc.*g.*100); % cm/s2
%     acc5=0;
%     for j=1:accnum    
%         if ABSacc(j,1) < 5
%            acc5(j,1) = 0;
%         else
%            acc5(j,1)=ABSacc(j,1);
%         end
%     end
%     CAV5(i,1)=sum(dt(i).*acc5(1:end)); % cm/s
%     disp(['CAV5 of ',num2str(i),' Wave is Obtain'])
% 
%     % Cumulative absolute displacement (CAD) (cm)
%     CAD(i,1)=sum(dt(i).*abs(velo(1:end)));
%     disp(['CAD of ',num2str(i),' Wave is Obtain'])
% 
%     % Effective Design Acceleration (EDA) (g)
%     N=4;
%     fc=9;
%     wc=fc/50;
%     [Bb,Aa]=butter(N,wc,'low');
%     dtlb1=filter(Bb,Aa,acc);
%     EDA(i,1)=max(abs(dtlb1));
%     disp(['EDA of ',num2str(i),' Wave is Obtain'])
% 
%     % Fajfar intensity (FI = PGV*(SD^0.25)) (cm/s^0.75)
%     FI(i,1)=PGV(i,1)*(SD(i,1)^0.25);
%     
%     % Sustained maximum acceleration (SMA) (g)
%     accabspx=0;
%     accabspx=sort(abs(acc));
%     SMA(i,1)=accabspx(end-2);
%     disp(['SMA of ',num2str(i),' Wave is Obtain'])
%     
%     % Sustained maximum velocity (SMV) (cm/sec)
%     velopx=0;
%     velopx=sort(abs(velo));
%     SMV(i,1)=velopx(end-2);
%     disp(['SMV of ',num2str(i),' Wave is Obtain'])
%     
    % Sa02(5%) (g)  
    T1=0.2;
    [sa,sv,sd]=spectrasa(dt(i),acc,T1,0.05);
    Sa02(i,1)=sa;
    disp(['Sa02 of ',num2str(i),' Wave is Obtain'])

    % Sa10(5%) (g)
    Tn=1.0;
    [sa,sv,sd]=spectrasa(dt(i),acc,Tn,0.05);
    Sa10(i,1)=sa;
    disp(['Sa10 of ',num2str(i),' Wave is Obtain'])
%    
%     % Sd02 (cm)
%     [sa,sv,sd]=spectrasa(dt(i),acc,T1,0.05);
%     Sd02(i,1)=sd*g*100;
%     disp(['Sd02 of ',num2str(i),' Wave is Obtain'])
%     
    % Sa20(5%) (g)  
    T2=2;
    [sa,sv,sd]=spectrasa(dt(i),acc,T2,0.05);
    Sa20(i,1)=sa;
    disp(['Sa20 of ',num2str(i),' Wave is Obtain'])
%     
%     % Sv20 (cm/sec)
%     [sa,sv,sd]=spectrasa(dt(i),acc,T2,0.05);
%     Sv20(i,1)=sv*g*100;
%     disp(['Sv20 of ',num2str(i),' Wave is Obtain'])
%     
%     % Sd20 (cm)
%     [sa,sv,sd]=spectrasa(dt(i),acc,T2,0.05);
%     Sd20(i,1)=sd*g*100;
%     disp(['Sd20 of ',num2str(i),' Wave is Obtain'])
% 
%     % Effective peak acceleration (EPA) (g)
%     ty=1;
%     ega=0;
%     for T=0.1:dt(i):0.5,
%         [sa,sv,sd]=spectrasa(dt(i),acc,T,0.05);
%         ega(ty)=sa;
%         ty=ty+1;
%     end
%     EPA(i,1)=mean(ega)/2.5;
% 
%     % Effective peak velocity (EPV) (cm/s)
%     egv=0;
%     Tv=1.0;
%     [sa,sv,sd]=spectrasa(dt(i),acc,Tv,0.05);
%     EPV(i,1)=sv*g*100/2.5;
%     
    %Housner intensity (cm)
    tu=1;
    for T=0.1:dt(i):2.5,
        [sa,sv,sd]=spectrasa(dt(i),acc,T,0.05);
        ihsum(tu)=sv*dt(i);
        tu=tu+1;
    end
    HI(i,1)=sum(ihsum)*100*g;
    ihsum=[];
    disp(['HI of ',num2str(i),' Wave is Obtain'])
% 
%     T3=0.8;
%     [sa,sv,sd]=spectrasa(dt(i),acc,T3,0.05);
%     Sa08(i,1)=sa;
%     disp(['Sa08 of ',num2str(i),' Wave is Obtain'])
%     disp(['====== ',num2str(i),' Wave Properties is Obtain ======'])
end
 IMWaves=[PGA PGV Ia CAV HI Sa02 Sa10 Sa20];
% save('F:\ValidationForML\FullModel\post\IMWave.mat','IMWaves');

% IM=[PGA SMA CAV CAV5 Arms Ia Ic SIR EDA PGV SMV CAD Vrms FI SED PGD Drms VAratio SD HI EPA EPV Sa02 Sa20 Sa10 Sv20 Sd02 Sd20  Sa08];
% % save IM.txt -ascii IM;
% toc;

