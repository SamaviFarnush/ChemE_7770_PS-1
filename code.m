clc
clear all

%control term
inducer_concentration=[0.0001:0.0001:10];               %mM
K_binding=0.3;                                          %mM
n=1.5;
weight_1=0.26;
weight_2=300;
binding_function=(inducer_concentration.^n)./(K_binding^n+inducer_concentration.^n);                        %non-dimensional
control_function=(weight_1+weight_2.*binding_function)./(1+weight_1+weight_2.*binding_function);            %control term

%Specific rate of transcription
max_transcription_rate=42;                                   %elongation_constant
k_elongation=max_transcription_rate/3075;     
k_initiation=1/42;
cell_volume=0.75*10^-15;                                     %L/cell
total_ribosome_R_XT=0.25*4600*10^6/(6.023*10^23*cell_volume);%micro-M
gene_concentration_Gp=2500*10^6/(6.023*10^23*cell_volume);   %micro-M 
                                                     
saturation_constant_Kxp=(80-50)/(38-10)*k_initiation;               
time_constant_tau_Xp=k_elongation/k_initiation;                               
thermodynamic_limit_r_x=k_elongation*total_ribosome_R_XT*gene_concentration_Gp/(saturation_constant_Kxp*time_constant_tau_Xp+(time_constant_tau_Xp+1)*gene_concentration_Gp);
rate_r_xhat=thermodynamic_limit_r_x.*control_function;

%mRNA concentration
doubling_time_td=30;                                %min
specific_growth_rate=log(2)/(60*doubling_time_td);  %1/s
mRNA_halflife=4;                                    %min
degradation_rate=log(2)/(mRNA_halflife*60);         %/s 
dry_weight=281*10^-15;                              %gDW/cell
mRNA_concentration=(rate_r_xhat./(degradation_rate+specific_growth_rate))*cell_volume/dry_weight*10^3; %nmol/gDW

%Plot
semilogx(inducer_concentration,mRNA_concentration,'-');
xlabel('Inducer (mmol/L)');
ylabel('mRNA concentration (nano-mol/gDW)');
%ylim([0 10]);
title('PS-1');
grid on