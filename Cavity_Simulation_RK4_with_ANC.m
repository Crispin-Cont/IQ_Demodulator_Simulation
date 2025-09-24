close all
clear all

tic; 

tmax = 0.5;
dt =0.07e-6;
t = (0:round(tmax/dt)-1).'*dt;
w = 2*pi*[0:numel(t)/2-1 -numel(t)/2:-1].'/(numel(t)*dt);

Qext=3*10^6;
Q0=2.7*10^10;
beta =Q0/Qext;
QL=Q0/(1+beta); 
%1.3 GHz Cavity 
%{
w0=2*pi*1300*10^6;  %Units in Hz
RL=(0.5*1080)*QL; %units in Ohm 
 %}

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
V=zeros(size(t));
DW_m=zeros(size(t,1),10);
DW=zeros(size(t));
wdot=zeros(size(t,1),10);
Eacc=zeros(size(t)); 

% Terms for the LFD
%Dw''+(2/tau_m)*Dw'+O^2_m*Dw=driving term 
O_m=2*pi*[157 182 189 215 292 331 380 412 462 471];
O_m=O_m.^2; 
k_p=2*pi*2.9636*[.02 .03 .055 .55 .085 .029 .052 .19 .075 .095]; % The sum of he array is 1.181
K_LFD=sum(k_p)/2/pi; %Static LFD 
G=k_p.*O_m; %
tau_m=2*[56.8 113.2 70.38 25 300 304.2 409.46 305.11 202.19 205.86]*10^-3; 

% Terms for the piezos
%Dw''+(2/tau_m)*Dw'+O^2_m*Dw=driving term 
k_piezo=2*pi*16.9348*[.02 .03 .055 .55 .085 .029 .052 .19 .075 .095]; 
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


microphonics =  0.01*sin(2*pi*5*t) + 0.01*sin(2*pi*10*t)+ 0.01*sin(2*pi*20*t)+ 0.01*sin(2*pi*50*t);
Im_k = I_m; % hold the previous value
Qm_k = Q_m; 
V(1)=2*forward(1)+ 0j;
y=zeros(2*10+2,1); % hold term for DW, dotDW, and cavity voltage
y(21)=V(1);
for k=1:1:size(forward,1)-1 
    
    k1 = dt*cavity_dynamics(y,V(1),Leff,forward(k),whalf,G_micro,G_piezo,G,tau_m,O_m,microphonics(k),piezo_all(k));
    k2 = dt*cavity_dynamics(y+0.5*dt*k1,V(1),Leff,forward(k),whalf,G_micro,G_piezo,G,tau_m,O_m,microphonics(k),piezo_all(k));
    k3 = dt*cavity_dynamics(y+0.5*dt*k2,V(1),Leff,forward(k),whalf,G_micro,G_piezo,G,tau_m,O_m,microphonics(k),piezo_all(k));
    k4 = dt*cavity_dynamics(y+dt*k3,V(1),Leff,forward(k),whalf,G_micro,G_piezo,G,tau_m,O_m,microphonics(k),piezo_all(k));

    y = y + (1/6)*(k1+2*k2+2*k3+k4);
    DW(k+1) = sum(y(1:10));
    V(k+1) = y(21)+ y(22)*1j;
    
    %Feedback Code 
    I_m = Im_k -mu*DW(k)*cos(omega_m*t(k)-phi_m); %k+1
    Q_m = Qm_k -mu*DW(k)*sin(omega_m*t(k)-phi_m); %k+1
    phi_m  = phi_m-eta*DW(k)*(Im_k.*sin(omega_m*t(k)-phi_m)-Qm_k.*cos(omega_m*t(k)-phi_m)); %k+1

    piezo_m = I_m.*cos(omega_m*t(k+1)) + Q_m.*sin(omega_m*t(k+1)); %k
    piezo_all(k+1)=sum(piezo_m); % at k
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


figure(3)
plot(t,piezo_all)
hold on 
plot(t,microphonics)

elapsedTime = toc;
disp(['Execution time: ', num2str(elapsedTime), ' seconds']);


function f=cavity_dynamics(y,Vconst,Leff,VF,whalf,G_micro,G_piezo,G,tau_m,O_m,microphonics,piezo_all)
    DW_modes = y(1:10);
    dotDW_modes = y(11:20);
    vcavity_real = y(21);
    vcavity_imag = y(22);
    VC = vcavity_real + vcavity_imag*1j;
    DW = sum(DW_modes);
    Eacc = (VC-Vconst)/Leff;
    dvdt=((-whalf+DW*1i)*VC+2*whalf*VF);
    dvdt_real = real(dvdt);
    dvdt_imag = imag(dvdt);

    m_term = G_micro*microphonics;
    p_term = G_piezo*piezo_all;
    V_term = -G*abs(Eacc)^2;
    d_dotDW_modes_dt=(V_term' + p_term' + m_term' -(2./tau_m').*dotDW_modes-O_m'.*DW_modes);
    

    f = [dotDW_modes;d_dotDW_modes_dt;dvdt_real;dvdt_imag]; 
end 