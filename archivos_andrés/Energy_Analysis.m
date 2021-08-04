%
% Se realiza en análisis del consumo de energía del sistema.
% Para esto se tiene en cuneta la energía que consume en un ciclo
% y con ésto se conoce la enegía que consume en 1 seg.
% Se tomó la funcioin (k/2)*x^2 para la enegía que adquiere la barra
% y se analizó para un ciclo completo en cada una de las frecuencias.
%
% Sería recomendable analizar la energía que se pierde en cada una
% de las frecuencias.
%
 
%% Convert Flex Sensor to meters  time_topic_3, flexion2
flex = 4.1199*flexion2;
flex = flex-mean(flex);
flex = flex./100;

time_flex = time_topic_3;

% Filter of the signal
windowSize = 2; 
b = (1/windowSize)*ones(1,windowSize);
a = 2;
flex_filtered = filter(b,a,flex);

figure(1)
plot(time_flex,flex)
hold on
plot(time_flex,flex_filtered)
xlabel('Time (s)');
ylabel('Amplitud (m)')

[freq_flex,P1] = Frequency_Spectrum(time_flex,flex)

figure(2)
plot(freq_flex*2*pi,P1)
xlabel('Frequency (rad/s)')
ylabel('Amplitud');

%% Angular Frequencies
% angle_pos = angule_position;
% time_pos = time_topic_4;


%% Flexion Energy  time_flex, flex
% kstiff = 18.11;
beam_properties
k_formulation = (3*E*Inertia)/(L^3); 
Flex_Energy_Test1 = (k_formulation/2)*(flex.^2);

%flexion
clear F
% To have a window of 16 seg (min_t+4) and (max_t-1)
min_t = 4.4;
max_t = 25.4;       %4.4+21
add_t = 21;

i = 1;
j = 1;

k = 1;
for h=1:length(Flex_Energy_Test1)
    if time_flex(h)>=min_t+4 &&  time_flex(h)<=max_t-1
        F(i) = Flex_Energy_Test1(h);
        flex_rad(i,j) = flex(h);
        time_rad(i,j) = time_flex(h);
        i = i + 1;
        flag = true;
    elseif time_flex(h)>max_t && flag==true
        flag = false;
        Frms(k) = rms(F);
        clear F
        min_t = min_t+add_t;       
        max_t = max_t+add_t;
        k = k + 1;
        i = 1;
        j = j + 1;
    end
end

omega = [12.9548,14.0495,15.0661,16.0045,16.995,18.015,19.0281,20.0707]';
%omega = [13,14,15,16,17,18,19,20]';
T = 2*pi./omega;
scale_to_seg = 1./T;

Elastic_Energy_1 = Frms(1:8)'.*scale_to_seg;   % Energy used in 1 seg
figure(3)
plot(omega,Elastic_Energy_1)
title('Elastic Energy Stored')
xlabel('Frequency (rad/s)');
ylabel('Energy (J)');


%% Energy Consumed by Voltage and Current 
% time_topic_1, voltage    time_topic_3, current
% the firt measurement it's at 6 seg
volt = voltage;
volt_time = time_topic_1;
cur = current;
cur_time = time_topic_3;
% figure(4)
% plot(volt_time,volt)

%voltage
clear V
% To have a window of 16 seg (min_t+4) and (max_t-1)
min_t = 4.4;
max_t = 25.4;       %4.4+21
add_t = 21;
i = 1;
k = 1;
flag = false;
for h=1:length(volt)
%     volt_time(i)
    if volt_time(h)>=min_t+4 &&  volt_time(h)<=max_t-1
        V(i) = volt(h);
        i = i + 1;
        flag = true;
    elseif volt_time(h)>max_t && flag==true
        i = 1;
        flag = false;
        Vrms(k) = rms(V);
        k = k + 1;
        clear V
        min_t = min_t+add_t;
        max_t = max_t+add_t;
    end
end

% current
clear I
% To have a window of 16 seg (min_t+4) and (max_t-1)
min_t = 4.4;
max_t = 25.4;       %4.4+21
add_t = 21;
i = 1;
k = 1;
for h=1:length(cur)
    if cur_time(h)>=min_t+4 &&  cur_time(h)<=max_t-1
        I(i) = cur(h);
        i = i + 1;
        flag = true;
    elseif cur_time(h)>max_t && flag==true
        i = 1;
        flag = false;
        Irms(k) = rms(I);
        k = k + 1;
        clear I
        min_t = min_t+add_t;
        max_t = max_t+add_t;
    end
end

power = Vrms(1:8).*(Irms(1:8)/1000);                  % Power consumtion in one cycle

omega = [12.9548,14.0495,15.0661,16.0045,16.995,18.015,19.0281,20.0707]';
%omega = [13,14,15,16,17,18,19,20]';
T = 2*pi./omega;
scale_to_seg = 1./T;

Electrical_Energy = power(1:8)'.*scale_to_seg;   % Energy used in 1 seg
figure(4)
plot(omega,Electrical_Energy)
title('Electrica Energy Consumed')
xlabel('Frequency (rad/s)');
ylabel('Energy (J)');

%% E_out/E_in test1
% 
% E_in_out = Elastic_Energy_1./Electrical_Energy;
% 
% figure(5)
% plot(omega,E_in_out)
% title('Relation between the accumulated energy and the given energy')
% ylabel('Energy Output / Energy Input')
% xlabel('Frequency (rad/s)')

%% Dissipated Energy  time_rad, flex_rad
c =  0.0085;

% 13 rad/s
flex_rad(1:end-1,1)=flex_rad(1:end-1,1)-mean(flex_rad(1:end-1,1));
[D,flex_pos1]=findpeaks(flex_rad(1:end-1,1),'MinPeakHeight',0.044,'MinPeakDistance',10);
% plot(time_rad(1:end-1,1),flex_rad(1:end-1,1),time_rad(flex_pos1,1),D,'or') 
Pmean(1) = mean(flex_rad(flex_pos1,1));
Tmean(1) = time_rad(end-1,1)-time_rad(1,1);

% Pmean1_Vector = ones(length(time_rad(1:end-1,1)),1)*Pmean1;
% hold on
% plot(time_rad(1:end-1,1),flex_rad(1:end-1,1))
% hold on
% plot(time_rad(1:end-1,1),Pmean1_Vector)

deltaW_(1) = omega(1)*pi*c*(Pmean(1)^2);  % This is interval time 16 aprox 
% deltaW(1)  = deltaW_1_1/Tmean(1);

%  14 rad/s
flex_rad(:,2)=flex_rad(:,2)-mean(flex_rad(:,2));
[D,flex_pos2]=findpeaks(flex_rad(:,2),'MinPeakHeight',0.044,'MinPeakDistance',10);
% plot(time_rad(1:end,2),flex_rad(1:end,2),time_rad(flex_pos2,2),D,'or') 
Pmean(2) = mean(flex_rad(flex_pos2,2));
Tmean(2) = time_rad(end,2)-time_rad(1,2);

deltaW_(2) = omega(2)*pi*c*(Pmean(2)^2);  % This is interval time 16 aprox 
% deltaW(2) = deltaW_2_1/Tmean(2);
 
%  15 rad/s
flex_rad(:,3)=flex_rad(:,3)-mean(flex_rad(:,3));
[D,flex_pos3]=findpeaks(flex_rad(:,3),'MinPeakHeight',0.044,'MinPeakDistance',10);
% plot(time_rad(1:end-1,3),flex_rad(1:end-1,3),time_rad(flex_pos3,3),D,'or') 
Pmean(3) = mean(flex_rad(flex_pos3,3));
Tmean(3) = time_rad(end-1,3)-time_rad(1,3);

deltaW_(3) = omega(3)*pi*c*(Pmean(3)^2);  % This is interval time 16 aprox 
% deltaW(3) = deltaW_3_1/Tmean(3);

% %  16 rad/s
flex_rad(:,4)=flex_rad(:,4)-mean(flex_rad(:,4));
[D,flex_pos4]=findpeaks(flex_rad(:,4),'MinPeakHeight',0.044,'MinPeakDistance',10);
% plot(time_rad(1:end,4),flex_rad(1:end,4),time_rad(flex_pos4,4),D,'or') 
Pmean(4) = mean(flex_rad(flex_pos4,4));
Tmean(4) = time_rad(end,4)-time_rad(1,4);

deltaW_(4) = omega(4)*pi*c*(Pmean(4)^2);  % This is interval time 16 aprox 
% deltaW(4) = deltaW_4_1/Tmean(4);

%  17 rad/s
flex_rad(:,5) = flex_rad(:,5)-mean(flex_rad(:,5));
[D,flex_pos5]=findpeaks(flex_rad(:,5),'MinPeakHeight',0.044,'MinPeakDistance',10);
% plot(time_rad(1:end,5),flex_rad(1:end,5),time_rad(flex_pos5,5),D,'or') 
Pmean(5) = mean(flex_rad(flex_pos5,5));
Tmean(5) = time_rad(end,5)-time_rad(1,5);
% Pmean5_Vector = ones(length(time_rad(1:end,5)),1)*Pmean5;
% hold on
% plot(time_rad(1:end,5),Pmean5_Vector)

deltaW_(5) = omega(5)*pi*c*(Pmean(5)^2);  % This is interval time 16 aprox 
% deltaW(5) = deltaW_5_1/Tmean(5);

%  18 rad/s
flex_rad(:,6) = flex_rad(:,6)-mean(flex_rad(:,6));
[D,flex_pos6]=findpeaks(flex_rad(:,6),'MinPeakHeight',0.044,'MinPeakDistance',10);
% plot(time_rad(1:end,6),flex_rad(1:end,6),time_rad(flex_pos6,6),D,'or') 
Pmean(6) = mean(flex_rad(flex_pos6,6));
Tmean(6) = time_rad(end,6)-time_rad(1,6);

deltaW_(6) = omega(6)*pi*c*(Pmean(6)^2);  % This is interval time 16 aprox 
% deltaW(6) = deltaW_6_1/Tmean(6);

% 19 rad/s
flex_rad(:,7) = flex_rad(:,7)-mean(flex_rad(:,7));
[D,flex_pos7]=findpeaks(flex_rad(:,7),'MinPeakHeight',0.020,'MinPeakDistance',10);
% plot(time_rad(1:end,7),flex_rad(1:end,7),time_rad(flex_pos7,7),D,'or') 
Pmean(7) = mean(flex_rad(flex_pos7,7));
Tmean(7) = time_rad(end,7)-time_rad(1,7);

deltaW_(7) = omega(7)*pi*c*(Pmean(7)^2);  % This is interval time 16 aprox 
% deltaW(7) = deltaW_7_1/Tmean(7);

% 20 rad/s
flex_rad(:,8) = flex_rad(:,8)-mean(flex_rad(:,8));
[D,flex_pos8]=findpeaks(flex_rad(:,8),'MinPeakHeight',0.015,'MinPeakDistance',8);
% plot(time_rad(1:end,8),flex_rad(1:end,8),time_rad(flex_pos8,8),D,'or') 
Pmean(8) = mean(flex_rad(flex_pos8,8));
Tmean(8) = time_rad(end,8)-time_rad(1,8);

deltaW_(8) = omega(8)*pi*c*(Pmean(8)^2);  % This is interval time 16 aprox 
% deltaW(8) = deltaW_8_1/Tmean(8);

beam_properties
k_formulation = (3*E*Inertia)/(L^3); 
Flex_Energy_Test2 = (k_formulation/2)*(Pmean.^2); 

omega = [12.9548,14.0495,15.0661,16.0045,16.995,18.015,19.0281,20.0707]';
%omega = [13,14,15,16,17,18,19,20]';
T = (2*pi./omega);
scale_to_seg = 1./T;

% By 2* because it's two parts of the beam
Elastic_Energy_2 = 2*Flex_Energy_Test2'.*scale_to_seg;   % Energy used in 1 seg
deltaW = 2*deltaW_'.*scale_to_seg; 

figure(7)
plot(omega,Elastic_Energy_2)
title('Elastic Energy Stored')
xlabel('Frequency (rad/s)');
ylabel('Energy (J)');

%% E_out/E_in test2

E_in_out = Elastic_Energy_2./Electrical_Energy;
E_in_out_with_damping = (Elastic_Energy_2-deltaW)./Electrical_Energy;

figure(10)
plot(omega,E_in_out)
title('Relation between the accumulated energy and the given energy')
ylabel('Energy Output / Energy Input')
xlabel('Frequency (rad/s)')

figure(11)
grid on
hold on
plot(omega,E_in_out_with_damping,'--*')
axis([12 21 0 2.5])
for p=1:8
 text(omega(p),E_in_out_with_damping(p),num2str(round(E_in_out_with_damping(p),3)))
end
% title('Relation between the accumulated energy minus energy Dissipated with respect the given energy')
% ylabel('Elastic Energy Stored - Dissipated Energy/Electrica Energy Consumed (J)')
ylabel('Energy Relation')
xlabel('Frequency (rad/s)')

figure(12)
grid on
hold on
plot(omega,deltaW)
ylabel('Energy Dissipated')
xlabel('Frequency (rad/s)')



figure(13)
grid on
hold on
plot(omega,Electrical_Energy,'--*')
for p=1:8
 text(omega(p),Electrical_Energy(p),num2str(round(Electrical_Energy(p),3)))
end
% title('Electrica Energy Consumed')
xlabel('Frequency (rad/s)');
ylabel('Electrical Energy(J)');

figure(14)
grid on
hold on
plot(omega,Elastic_Energy_2,'--*')
for p=1:8
 text(omega(p),Elastic_Energy_2(p),num2str(round(Elastic_Energy_2(p),3)))
end
% title('Elastic Energy Stored')
xlabel('Frequency (rad/s)');
ylabel('Elastic Energy (J)');

figure(15)
grid on
hold on
plot(omega,deltaW,'--*')
for p=1:8
 text(omega(p),deltaW(p),num2str(round(deltaW(p),3)))
end
% title('Dissipated Energy')
ylabel('Dissipated Energy (J)');
xlabel('Frequency (rad/s)');

% figure(13)
% subplot(1,3,1)
% grid on
% hold on
% plot(omega,Electrical_Energy,'--*')
% for p=1:8
%  text(omega(p),Electrical_Energy(p),num2str(round(Electrical_Energy(p),2)))
% end
% % title('Electrica Energy Consumed')
% xlabel('Frequency (rad/s)');
% ylabel('Electrica Energy Consumed (J)');
% subplot(1,3,2)
% grid on
% hold on
% plot(omega,Elastic_Energy_2,'--*')
% for p=1:8
%  text(omega(p),Elastic_Energy_2(p),num2str(round(Elastic_Energy_2(p),2)))
% end
% % title('Elastic Energy Stored')
% xlabel('Frequency (rad/s)');
% ylabel('Elastic Energy Stored (J)');
% subplot(1,3,3)
% grid on
% hold on
% plot(omega,deltaW,'--*')
% for p=1:8
%  text(omega(p),deltaW(p),num2str(round(deltaW(p),2)))
% end
% % title('Dissipated Energy')
% ylabel('Dissipated Energy (J)');
% xlabel('Frequency (rad/s)');
% 
% %% deltaW divied by Energy stored
% 
% figure(7)
% WW = deltaW'./Elastic_Energy;
% plot(omega,WW)
% title('Relation between the dissipated energy and the elastic energy is transformed in work')
% xlabel('Frequency (rad/s)')
% ylabel('Dissipated Energy / Energy Used')
