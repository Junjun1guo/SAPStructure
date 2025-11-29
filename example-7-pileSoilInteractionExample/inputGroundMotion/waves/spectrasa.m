function [sa,sv,sd]=spectrasa(dt,acc,T,Damp)
    %Sa(5%)
         j=1;
         [accnum,l]=size(acc);
         Dt=dt;
         Accelerate=acc;
         count=accnum;
         Displace=zeros(1,count); 
         Velocity=zeros(1,count); 
         AbsAcce=zeros(1,count); 

         Frcy=2*pi/T ; 
         DamFrcy=Frcy*sqrt(1-Damp*Damp); 
         e_t=exp(-Damp*Frcy*Dt);
         s=sin(DamFrcy*Dt);
         c=cos(DamFrcy*Dt);
         A=zeros(2,2);
         A(1,1)=e_t*(s*Damp/sqrt(1-Damp*Damp)+c);
         A(1,2)=e_t*s/DamFrcy;
         A(2,1)=-Frcy*e_t*s/sqrt(1-Damp*Damp);
         A(2,2)=e_t*(-s*Damp/sqrt(1-Damp*Damp)+c);
         d_f=(2*Damp^2-1)/(Frcy^2*Dt); 
         d_3t=Damp/(Frcy^3*Dt);
         B=zeros(2,2);
         B(1,1)=e_t*((d_f+Damp/Frcy)*s/DamFrcy+(2*d_3t+1/Frcy^2)*c)-2*d_3t;
         B(1,2)=-e_t*(d_f*s/DamFrcy+2*d_3t*c)-1/Frcy^2+2*d_3t;
         B(2,1)=e_t*((d_f+Damp/Frcy)*(c-Damp/sqrt(1-Damp^2)*s)-(2*d_3t+1/Frcy^2)*(DamFrcy*s+Damp*Frcy*c))+1/(Frcy^2*Dt);
         B(2,2)=e_t*(1/(Frcy^2*Dt)*c+s*Damp/(Frcy*DamFrcy*Dt))-1/(Frcy^2*Dt);
         for ii=1:(count-1) 
              Displace(ii+1)=A(1,1)*Displace(ii)+A(1,2)*Velocity(ii)+B(1,1)*Accelerate(ii)+B(1,2)*Accelerate(ii+1);
              Velocity(ii+1)=A(2,1)*Displace(ii)+A(2,2)*Velocity(ii)+B(2,1)*Accelerate(ii)+B(2,2)*Accelerate(ii+1);
              AbsAcce(ii+1)=-2*Damp*Frcy*Velocity(ii+1)-Frcy^2*Displace(ii+1);
         end
         MDis(j)=max(abs(Displace));
         MVel(j)=max(abs(Velocity));
         if T==0.0
            MAcc(j)=max(abs(Accelerate));
         else
            MAcc(j)=max(abs(AbsAcce));
         end
         Displace=zeros(1,count);
         Velocity=zeros(1,count);
         AbsAcce=zeros(1,count);

         sa=MAcc;
         sv=MVel;
         sd=MDis;