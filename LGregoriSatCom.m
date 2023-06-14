%% SAT COMM ASSIGNEMENT 
% Author: Ludovico Gregori URN 6778145
% Date: 15/12/2022

% Formatting
clearvars; close all; clc; format longG;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesFontSize', 16);

% DATA & FORMULAS to calculate faster every value
d2r = pi/180;
r2d = 180/pi;

%% COMPLETE LINK BUDGET ANALISYS (TABLE1)

%% elevation and true azimuth angles, range

R_earth = 6378.137*10^3; %km to meters, earth mean radius
h = 35786*10^3; %km to meters, geo altitude
sigma = R_earth/(R_earth+h)

psi_e = 1.8; %degrees East, target longitude
psi_s = 9.0; %degrees East, satellite longitude
psi_es = (psi_e - psi_s)*d2r

teta_e = 46.07 %degree N, target latitude
teta_e_rad = teta_e*d2r;
beta = acos(cos(teta_e_rad)*cos(psi_es));
beta_deg = beta*r2d

elevation = atan2(cos(beta)-sigma,sin(beta));
elevation_deg = elevation * r2d

azimuth = atan2(tan(psi_es),sin(teta_e_rad));
azimuth_deg = azimuth*r2d + 180

range = h*sqrt(1+0.4199*(1-cos(beta))) %distance target-satellite

%% free space path loss value
ligthspeed = 3*10^8; %m/s 
frequency = 20*10^9; %GHz
lambda = ligthspeed/frequency;

%free space path loss 
FSPL = 20*log10(4*pi*range/lambda) %dB
FSPL2 = 32.45 + 20*log10(20000) + 20*log10(range/1000) %dB

%% spreading loss
spreading_loss = 10*log10(4*pi*(range)^2) %dB

%% antenna gain 
antenna_diameter = 5.1324; %meters, diameter
sat_area = pi*(antenna_diameter/2)^2;
eta_area = 0.7;
Ae = eta_area * sat_area;
Gsat = (4*pi*Ae/(lambda^2));

GsatdB = 10*log10(Gsat);
GsatdB = 20*log10(antenna_diameter*pi*20/0.3) + 10*log10(70/100) %dBi

%% PFD power flux density, clear sky
EIRPsat = 58.600; %dBW, table line 17
EIRP2 = EIRPsat + 3; %dB 
atm_loss = 0.49; %dB, table line 31
feeder_loss = 0.3; %dB,
PFD_clearsky = EIRPsat - spreading_loss - atm_loss - feeder_loss %dBWm^-2
PFD_clearsky2 = EIRP2 - spreading_loss - atm_loss - feeder_loss %dBWm^-2

%% noise temperature and noise figure 

NF = 2; %dB noise figure from LNB, paper
noise_temp = 290*(10^(NF/10)-1) %Kelvin

%% system noise temperature
alpha = 10^(-feeder_loss/10) ;
antenna_temp = 43; %Kelvin, pag.8
receiver_temp = noise_temp;
feeder_temp = 290; %Kelvin
system_temp = alpha*antenna_temp + (1-alpha)*feeder_temp + receiver_temp %check elevation site temperature in case add here

%% figure of merit
%G/T referred to LNA 
G_T_sat = GsatdB - 10*log10(system_temp)  %dB

%% symbol rate, information rate, BW
m = 3;
R_coding = 0.75; %FEC rate, code rate
symbol_rate = 50*10^6; %MBaud
roll_fact = 0.2;

coding_rate = m * symbol_rate %baud per sec
information_rate = R_coding * coding_rate %bit per sec, bit rate
occupied_BW = (1+roll_fact)*symbol_rate %Hz

%area of the receiver antenna
semimajor_axis = 0.77/2; %m
semiminor_axis = 0.72/2; %m
receiver_area = pi*semiminor_axis*semimajor_axis; %m2
efficiency_receiver = 0.70; %eta for effective area, user's antenna efficiency line 6 TABLE 1
eff_receiver_area = efficiency_receiver * receiver_area;
Greceiver = (4*pi*eff_receiver_area)/(lambda^2);
Greceiver_dB = 10*log10(Greceiver)

averagediameter = (0.72+0.77)/2;
averagearea = pi*averagediameter^2/4;
Greceiver2dB = 20*log10((averagediameter)*pi*20/0.3)+10*log10(0.7)

%% C_N  required
eta_required = 2.228;
Eb_N0required = 4.41; %dB
C_N_table_required = 7.12; %dB 
Boltz_cost = -228.6 ; %dBW/KelvinHz
k_boltz = 10^-(Boltz_cost/10);
pointing_loss = 0.55; %dB
losses = FSPL2 + atm_loss + pointing_loss; %dB
T_rx = 290; %Kelvin, average value
Margin_table = 1; %dB
r=0.75;

G_T_rx = Greceiver_dB - 10*log10(system_temp)
C_N_required = Eb_N0required - 10*log10(1+roll_fact) + 10*log10(m) + 10*log10(r) %dB
C_Ntot = EIRPsat - losses + Greceiver_dB - Boltz_cost - 10*log10(system_temp*occupied_BW); %dB
C_Ntot2 = EIRP2 - losses + Greceiver_dB - Boltz_cost - 10*log10(system_temp*occupied_BW); %dB

%% C/No thermal
%feeder loss included into G_T term.
C_No = EIRPsat - losses + G_T_rx - Boltz_cost %dB/Hz
C_No2 = EIRP2 - losses + G_T_rx - Boltz_cost %dB/Hz

%% C/Io 

test_long = 1.8; %deg, Limousin user's longitude
sat_long = 9.0; %deg, satellite longitude
deltat = 46.07*d2r; %rad, user's latitude 
deltaL = abs((test_long - sat_long)*d2r); %rad, longitude difference 

lambda_test = acos(cos(deltat)*cos(deltaL)); %rad, 
Az = acos(sin(deltat)/sin(lambda_test)); %rad, Azimuth
eta = atan2(0.151294*sin(lambda_test),(1-0.151294*cos(lambda_test))); %deg, nadir angle
etadeg = eta*r2d;
test_site_X_deg = -eta*sin(Az)*180/pi
test_site_Y_deg = eta*cos(Az)*180/pi

beam19_X = -0.83752;
beam19_Y = 6.934528;

%off axis angle
%alpha_angle = abs(r_x-r_beam19) %deg, length b

%beta angle between test site and interference beams
beam6_X = -0.93968;
beam6_Y = 6.668402;
beam9_X = -1.11908;
beam9_Y = 6.889933;
beam12_X = -0.65812;
beam12_Y = 6.712996;
beam25_X = -1.01692;
beam25_Y = 7.156057;
beam28_X = -0.55596;
beam28_Y = 6.979122;
beam32_X = -0.73536;
beam32_Y = 7.200652;

alpha_angle  = sqrt((beam19_X-test_site_X_deg)^2 + (beam19_Y-test_site_Y_deg)^2)  %deg, length 
beta_angle6  = sqrt((beam6_X-test_site_X_deg)^2  + (beam6_Y-test_site_Y_deg)^2)   %deg, length 
beta_angle9  = sqrt((beam9_X-test_site_X_deg)^2  + (beam9_Y-test_site_Y_deg)^2)   %deg, length 
beta_angle12 = sqrt((beam12_X-test_site_X_deg)^2 + (beam12_Y-test_site_Y_deg)^2)  %deg, length 
beta_angle25 = sqrt((beam25_X-test_site_X_deg)^2 + (beam25_Y-test_site_Y_deg)^2)  %deg, length 
beta_angle28 = sqrt((beam28_X-test_site_X_deg)^2 + (beam28_Y-test_site_Y_deg)^2)  %deg, length 
beta_angle32 = sqrt((beam32_X-test_site_X_deg)^2 + (beam32_Y-test_site_Y_deg)^2)  %deg, length 

% alpha_angle = 4.6673e-04;
% beta_angle6  = 0.2848;
% beta_angle9  = 0.2853;
% beta_angle12 = 0.2393;
% beta_angle25 = 0.2435;
% beta_angle28 = 0.2848;
% beta_angle32 = 0.2853;

gain_alpha  = ((0)+((alpha_angle)-(0))*((-0.03148)-(0))/((0.01)-(0)));%dB
gain_beam6  = ((-23.2723)+((beta_angle6)-(0.28))*((-22.5581)-(-23.2723))/((0.29)-(0.28))) %dB
gain_beam9  = ((-23.2723)+((beta_angle9)-(0.28))*((-22.5581)-(-23.2723))/((0.29)-(0.28))) %dB
gain_beam12 = ((-23.2723)+((beta_angle12)-(0.28))*((-22.5581)-(-23.2723))/((0.29)-(0.28))) %dB
gain_beam25 = ((-23.2723)+((beta_angle25)-(0.28))*((-22.5581)-(-23.2723))/((0.29)-(0.28))) %dB
gain_beam28 = ((-23.2723)+((beta_angle28)-(0.28))*((-22.5581)-(-23.2723))/((0.29)-(0.28))) %dB
gain_beam32 = ((-23.2723)+((beta_angle32)-(0.28))*((-22.5581)-(-23.2723))/((0.29)-(0.28))) %dB

psi_e = 1.8; %degrees East, target longitude
teta_e_rad = 46.07*d2r; %rad, target latitude

%beam6
psi_beam6 = 1.304963; %degrees East, beam 6 longitude
psi_es_beam6 = (psi_e - psi_beam6)*d2r %longitude difference, deltaL used in the figures
beta_beam6 = acos(cos(teta_e_rad)*cos(psi_es_beam6)); %discover this angle!!!!!
range_beam6 = h*sqrt(1+0.4199*(1-cos(beta_beam6))) %distance target-satellite
%free space path loss beam 6
FSPL_beam6 = 20*log10((4*pi*range_beam6)/(lambda)) %dB

%beam9
psi_beam9 = -0.57647; %degrees East, beam 9 longitude
psi_es_beam9 = (psi_e - psi_beam9)*d2r %longitude difference, deltaL used in the figures
beta_beam9 = acos(cos(teta_e_rad)*cos(psi_es_beam9)); %discover this angle!!!!!
range_beam9 = h*sqrt(1+0.4199*(1-cos(beta_beam9))) %distance target-satellite
%free space path loss beam 6
FSPL_beam9 = 20*log10(4*pi*range_beam9/lambda) %dB

%beam12
psi_beam12 = 3.583562; %degrees East, beam 12 longitude
psi_es_beam12 = (psi_e - psi_beam12)*d2r %longitude difference, deltaL used in the figures
beta_beam12 = acos(cos(teta_e_rad)*cos(psi_es_beam12)); %discover this angle!!!!!
range_beam12 = h*sqrt(1+0.4199*(1-cos(beta_beam12))) %distance target-satellite
%free space path loss beam 6
FSPL_beam12 = 20*log10(4*pi*range_beam12/lambda) %dB

%beam25
psi_beam25 = -0.20506; %degrees East, beam 25 longitude
psi_es_beam25 = (psi_e - psi_beam25)*d2r %longitude difference, deltaL used in the figures
beta_beam25 = acos(cos(teta_e_rad)*cos(psi_es_beam25)); %discover this angle!!!!!
range_beam25 = h*sqrt(1+0.4199*(1-cos(beta_beam25))) %distance target-satellite
%free space path loss beam 6
FSPL_beam25 = 20*log10(4*pi*range_beam25/lambda) %dB

%beam28
psi_beam28 = 4.187885; %degrees East, beam 28 longitude
psi_es_beam28 = (psi_e - psi_beam28)*d2r %longitude difference, deltaL used in the figures
beta_beam28 = acos(cos(teta_e_rad)*cos(psi_es_beam28)); %discover this angle!!!!!
range_beam28 = h*sqrt(1+0.4199*(1-cos(beta_beam28))) %distance target-satellite
%free space path loss beam 6
FSPL_beam28 = 20*log10(4*pi*range_beam28/lambda) %dB

%beam32
psi_beam32 = 2.299807; %degrees East, beam 32 longitude
psi_es_beam32 = (psi_e - psi_beam32)*d2r %longitude difference, deltaL used in the figures
beta_beam32 = acos(cos(teta_e_rad)*cos(psi_es_beam32)); %discover this angle!!!!!
range_beam32 = h*sqrt(1+0.4199*(1-cos(beta_beam32))) %distance target-satellite
%free space path loss beam 6
FSPL_beam32 = 20*log10(4*pi*range_beam32/lambda) %dB

%C/I contributes, not dBs
Galpha_FSPL  = 10^((gain_alpha-FSPL)/10)
Gbeta_FSPL_6 = 10^((gain_beam6-FSPL_beam6)/10)
Gbeta_FSPL_9 = 10^((gain_beam9-FSPL_beam9)/10)
Gbeta_FSPL_12= 10^((gain_beam12-FSPL_beam12)/10)
Gbeta_FSPL_25 = 10^((gain_beam25-FSPL_beam25)/10)
Gbeta_FSPL_28 = 10^((gain_beam28-FSPL_beam28)/10)
Gbeta_FSPL_32 = 10^((gain_beam32-FSPL_beam32)/10)

%C/I tot clear sky condition
C_I_tot = Galpha_FSPL / (Gbeta_FSPL_6 + Gbeta_FSPL_9 + Gbeta_FSPL_12 + Gbeta_FSPL_25 + Gbeta_FSPL_28 + Gbeta_FSPL_32) %linear

C_ItotnoFSPL = 10^(gain_alpha/10) / (10^((gain_beam6 + gain_beam9 + gain_beam12 + gain_beam25 + gain_beam28 + gain_beam32)/10)) %dB

C_I_tot_dB = 10*log10(C_I_tot) %dB
C_I0_dB = C_I_tot_dB + 10*log10(occupied_BW) %dB Hz

%C/(N0+I0)tot
C_N0I0_tot = 10*log10(1/((10^-(C_No/10)+10^-(C_I0_dB/10)))) %dB/Hz
C_N0I02_tot = 10*log10(1/((10^-(C_No2/10)+10^-(C_I0_dB/10)))) %dB/Hz

%C/(N+I) tot
C_NI_tot = C_N0I0_tot - 10*log10(occupied_BW) %dB 
C_NI2_tot = C_N0I02_tot - 10*log10(occupied_BW) %dB 

%C/I adjacent satellite interference
C_Iadj0 = 91.28; %dBHz
C_Iadj  = C_Iadj0 - 10*log10(occupied_BW)

%Eb/NoIo tot
Eb_N0I0_tot  = C_NI_tot + 10*log10(1+roll_fact) - 10*log10(m) - 10*log10(r) %dB, equation 20a page 27
Eb_N0I02_tot = C_NI2_tot + 10*log10(1+roll_fact) - 10*log10(m) - 10*log10(r) %dB, equation 20a page 27

%Eb/NI tot 
Eb_NI_tot  = Eb_N0I0_tot - 10*log10(occupied_BW);
Eb_NI2_tot = Eb_N0I02_tot - 10*log10(occupied_BW);

%Margin
Eb_N0required = 4.41; %dB
Margin  = Eb_N0I0_tot - Eb_N0required  %dB
Margin2 = Eb_N0I02_tot - Eb_N0required %dB

%channel capacity noise and interference
C = occupied_BW*log2(1 + 10^(C_NI_tot/10)) %bits per input symbol, baud
C2 = occupied_BW*log2(1 + 10^(C_NI2_tot/10)); %bits per input symbol, baud

%% rain attenuation

%effective rain height
height_rain = 5.0 - 0.075*(teta_e - 23) %km

%slant-path length
height_site = 0.087; %km
slant_length = (height_rain - height_site)/sin(elevation) %km

%horizontal projection
horizontal_proj = slant_length*cos(elevation) %km

%rain intensity R1percent
R001 = 40.5; %mm/h
L0 = 35*exp(-0.015*R001)

%reduction factor
r001 = 1/(1+horizontal_proj/L0) %adimensional

%% specific attenuation
%need to calculate k and alpha_rain coefficients before obtaining specific
%attenuation value
k_h = 0.0751;
k_v  = 0.0691;
alpha_rain_h = 1.099;
alpha_rain_v = 1.065;
pol_tilt_angle = 45; %deg , teta line 15, TABLE 1
pol_tilt_rad = pol_tilt_angle * d2r; %rad
k_total = (k_h+k_v + (k_h - k_v)*(cos(elevation)^2)*cos(2*pol_tilt_rad))/2 ;
alpha_rain = (k_h*alpha_rain_h + k_v*alpha_rain_v + (k_h*alpha_rain_h - k_v*alpha_rain_v)*(cos(elevation)^2)*cos(2*pol_tilt_rad))/(2*k_total);
gammaR = k_total*(R001)^alpha_rain %dB/km

%% rain attenuation percentage 
A001 = gammaR * slant_length * r001 %dB
p=0.3; %availability percentage 
beta_rain = 0; %latitude station 46.07 deg > 36 deg. 
A03 = A001*(p/0.01)^-(0.655+0.033*log(p)-0.045*log(A001)-beta_rain*(1-p)*sin(elevation)) %dB

%% increase noise temp and total system noise temp
loss_precip = A03;
increase_temp = alpha*275*(1-10^(-loss_precip/10)) %Kelvin
total_temp_rain = system_temp + increase_temp 
        
%new values of all the parameter functions of the noise temperature system
atm_rain_loss = 1.24; %dB 
rain_losses = FSPL2 + atm_rain_loss + pointing_loss + A03; 
PFD_rain = EIRPsat - spreading_loss - atm_rain_loss  - feeder_loss
PFD_rain2 = EIRP2 - spreading_loss - atm_rain_loss  - feeder_loss

G_T_rx_rain = Greceiver_dB - 10*log10(total_temp_rain)
C_Ntot_rain = EIRPsat - rain_losses + Greceiver_dB - Boltz_cost - 10*log10(total_temp_rain*occupied_BW) %dB
C_Ntot2_rain = EIRP2 - rain_losses + Greceiver_dB - Boltz_cost - 10*log10(total_temp_rain*occupied_BW) %dB

C_No_rain = EIRPsat - rain_losses + G_T_rx - Boltz_cost %dB/Hz
C_No2_rain = EIRP2 - rain_losses + G_T_rx - Boltz_cost %dB/Hz

C_I0_rain = C_I0_dB %dB Hz

%C/I adjacent satellite interference
C_Iadj0 = 91.28; %dBHz
C_Iadj_rain = C_Iadj0 - 10*log10(occupied_BW)

CoverN0I0_rain = 10*log10(1/((10^-(C_No_rain/10)+10^-(C_I0_rain/10)))) %dBHz
CoverN0I02_rain = 10*log10(1/((10^-(C_No2_rain/10)+10^-(C_I0_rain/10)))) %dBHz

CoverNI_rain = CoverN0I0_rain - 10*log10(occupied_BW) %dB
CoverNI2_rain = CoverN0I02_rain - 10*log10(occupied_BW) %dB

Eb_N0I0_rain = CoverNI_rain + 10*log10(1+roll_fact) - 10*log10(m) - 10*log10(r) %dB
Eb_N0I02_rain = CoverNI2_rain + 10*log10(1+roll_fact) - 10*log10(m) - 10*log10(r) %dB

%Margin
Eb_N0required = 4.41; %dB
Margin_rain = Eb_N0I0_rain - Eb_N0required  %dB
Margin2_rain = Eb_N0I02_rain - Eb_N0required  %dB

%channel capacity in the rain
C_rain = occupied_BW*log2(1 + 10^(CoverNI_rain/10)) %bits per input symbol, baud
C2_rain = occupied_BW*log2(1 + 10^(CoverNI2_rain/10)); %bits per input symbol, baud


%% part4 data rate with the 3dB increase for the sat EIRP 
% EIRPsat = 58.6 dB
% Gsat = 8.088281539355804e+05 , linear
% information_rate =   112500000 bit per sec
% EIRP2 = EIRPsat + 3; %dB 
% EIRP (dBWatt) = Power (dbWatt) + Gain (dBi)  

% power at the satellite antenna 
power1dB = EIRPsat - GsatdB %dBW
power1_linear = 10^(power1dB/10) %Watt
power2_dB = EIRP2 - GsatdB %dBW
power2_linear = 10^(power2_dB/10) %Watt

%energy per bit
Eb1 = power1_linear/information_rate %Watt sec = Joule 
Eb2 = power2_linear/information_rate %Watt sec = Joule

% C = occupied_BW*log2(1 + C_N_tot)
C1_noise = occupied_BW*log2(1 + 10^(C_Ntot/10))
C2_noise = occupied_BW*log2(1 + 10^(C_Ntot2/10))
C1rain_noise = occupied_BW*log2(1 + 10^(CoverNI_rain/10))
C2rain_noise = occupied_BW*log2(1 + 10^(CoverNI2_rain/10))

% C = occupied_BW*log2(1 + C_NI_tot)
C1_noise_inter = C
C2_noise_inter = C2
C1rain_noise_inter = C_rain
C2rain_noise_inter = C2_rain

%bit rate Rb
% coding_rate = m * symbol_rate %baud per sec
% information_rate = R_coding * coding_rate %bit per sec, bit rate
% occupied_BW = (1+roll_fact)*symbol_rate %Hz

spectral_eff = information_rate/occupied_BW;

bit_rate1 = spectral_eff*C1_noise /(log2(1 + 10^(C_Ntot/10))) %bit per sec, equal to information_rate, GOOD.

bit_rate2 = 10^((C_NI2_tot-Eb_NI2_tot)/10) 

spectral_eff2 = bit_rate2/occupied_BW 


