%{
    This code simulates cavity detuning and accelerating gradient for a 650 MHz cavity. 
    The parameters used are apporoximation of the real values. In this configuration the 
    maximum dt possible is 5e-8, anything larger will result in numerical instabilities. 

%}

close all
clear all

tic; 

tmax = 3;
dt = 0.05e-6;
t = (0:round(tmax/dt)-1).'*dt;
w = 2*pi*[0:numel(t)/2-1 -numel(t)/2:-1].'/(numel(t)*dt);

Qext=3*10^6;
Q0=2.7*10^10;
beta =Q0/Qext;
QL=Q0/(1+beta); 
%650 MHz cavity Low beta 
w0=2*pi*650*10^6;  %Units in Hz
RL=(0.5*341)*QL; %units in Ohm 
whalf= w0/(2*QL);
wT = whalf/(1+1/beta);
Leff=1.038; %Units in Meters
tau=2*QL/w0;  %Units in seconds 
Amp=(25*10^6 *Leff)/(RL); 
forward=0.5*25*Leff*ones(size(t,1),1); % now in units of Voltage
%forward is the current of the generator
% pre-allocate variable size 
V=zeros(size(t)); %Cavity voltage
DW_m=zeros(size(t,1),10);   %Angular frequency per mode, 10 in this case
DW=zeros(size(t));          %Total angular frequency
wdot=zeros(size(t,1),10);   %derivate of mode angular frequency
Eacc=zeros(size(t));        %Accelerating gradient

% Terms for the LFD
%Dw''+(2/tau_m)*Dw'+O^2_m*Dw=driving term 
O_m=2*pi*[157 182 189 215 292 331 380 412 462 471]; % Modes of the cavity 
O_m=O_m.^2; 
k_p=2*pi*2.9636*[.02 .03 .055 .55 .085 .029 .052 .19 .075 .095]; % The sum of he array is 1.181, coupling of LFD
K_LFD=sum(k_p)/2/pi; %Static LFD 
G=k_p.*O_m; %
tau_m=2*[56.8 113.2 70.38 25 300 304.2 409.46 305.11 202.19 205.86]*10^-3; %time constant of the mechanical modes

% Terms for the piezos
%Dw''+(2/tau_m)*Dw'+O^2_m*Dw=driving term 
k_piezo=2*pi*16.9348*[.02 .03 .055 .55 .085 .029 .052 .19 .075 .095]; %coupling of the piezos
G_piezo=k_piezo.*O_m;

% Terms for the piezos
%Dw''+(2/tau_m)*Dw'+O^2_m*Dw=driving term 
k_micro=2*pi*64.0982*[.02 .03 .055 .55 .085 .029 .052 .19 .075 .095]; 
G_micro=k_micro.*O_m;

%Terms for the ANC 
mu = 1e-9; 
eta = 1e-9;
gamma = 0.01; 
I_m = zeros(4,1);
Q_m = zeros(4,1);  
phi_m  = zeros(4,1);
omega_m = 2*pi*[5; 10; 20; 50];
piezo_all = zeros(size(t,1),1);
Im_k = I_m; % hold the previous value
Qm_k = Q_m; 

% Driving terms for the second order detuning: d^2/dt^2(delta Omega_m) 
% Mainly used for debugging, 
m_term = zeros(size(t,1),10);
p_term = zeros(size(t,1),10);
V_term = zeros(size(t,1),10);

microphonics =  0.01*sin(2*pi*5*t) + 0.01*sin(2*pi*10*t)+ 0.01*sin(2*pi*20*t)+ 0.01*sin(2*pi*50*t); % uses same frequency as the omega_m
V(1)=2*forward(1)+ 0j;

for k=1:1:size(forward,1)-1 
    % First Order cavity equation
    f1=((-whalf+DW(k)*1i)*V(k)+2*whalf*(forward(k))); 
    V(k+1)=V(k)+ dt*f1;
    Eacc(k)=abs((V(k)-V(1)))/Leff; 
    
    %2nd order Detuning Equation (Mechanical Mode)
    m_term(k,:) = G_micro*microphonics(k);
    p_term(k,:) = G_piezo*piezo_all(k);
    V_term(k,:) = -G*(Eacc(k))^2;
    f4=(V_term(k,:) + p_term(k,:) + m_term(k,:) -(2./tau_m).*wdot(k,:)-O_m.*DW_m(k,:));
    wdot(k+1,:)=wdot(k,:)+dt*f4;
    DW_m(k+1,:)=DW_m(k,:)+dt*wdot(k,:);
    DW(k+1,1)=sum(DW_m(k+1,:));
    
    %Feedback Code using Least Mean Square algorithm (LMS), note that there is no delay
    %---------------
    I_m = Im_k -mu*DW(k)*cos(omega_m*t(k)-phi_m); %k+1
    Q_m = Qm_k -mu*DW(k)*sin(omega_m*t(k)-phi_m); %k+1
    phi_m  = phi_m-eta*DW(k)*(Im_k.*sin(omega_m*t(k)-phi_m)-Qm_k.*cos(omega_m*t(k)-phi_m)); %k+1

    piezo_m = I_m.*cos(omega_m*t(k+1)) + Q_m.*sin(omega_m*t(k+1)); %k+1
    piezo_all(k+1)=sum(piezo_m); % at k
    %detuning(k)=detuning(k)+piezo_all(k);
    Im_k = I_m;
    Qm_k = Q_m;

end

detuning  = DW/2/pi;

figure(10)
yyaxis left
plot(1000*t,abs(V/(Leff)),'LineWidth',2)
xlabel('Time [ms]')
ylabel('Eacc [MV/m]')

yyaxis right
plot(1000*t,detuning,'LineWidth',2)
xlabel('Time [ms]')
ylabel('\Deltaf [Hz]')


legend('Eacc','Detuning')
legend('boxoff')

figure(2)
plot(t,sum(m_term,2)/2/pi)
hold on 
plot(t,sum(p_term,2)/2/pi)
plot(t,sum(V_term,2)/2/pi)

legend('microphonics term','piezo term')
legend('boxoff')


figure(3)
plot(t,piezo_all)
hold on 
plot(t,microphonics)

elapsedTime = toc;
disp(['Execution time: ', num2str(elapsedTime), ' seconds']);