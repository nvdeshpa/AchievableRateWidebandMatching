clear all;
close all;

%% Universal constants
c= 3*10^8;  %% speed of light
kb= 1.380649*10^(-23); %% Boltzmann constant
T=290; %% Temperature in Kelvin

%% Other constants
R=50; %% real part of reference impedance
Z_0=R; %% reference impedance
Dtxrx=500; %% distance between Tx and Rx
G=1.5; %% Gain of Chu's antenna
f_c= 7*10^9; %% center frequency
lambda_c= c/f_c; %% center wavelength

%% Varying parameter
BW_arr=[0.1:0.1:0.6]*f_c;  %% range of bandwidths

%% Achievable rate for different cases
achievable_rate_SISO_ideal_arr= zeros(length(BW_arr),1);  %% Ideal -> No Bode/Fano theory
achievable_rate_SISO_optimized_arr= zeros(length(BW_arr),1); %% Shannon+ Bode/Fano theory
achievable_rate_SISO_freq_flat_arr= zeros(length(BW_arr),1); %% Frequency-flat matching
achievable_rate_SISO_no_MN_arr= zeros(length(BW_arr),1); %% no matching
achievable_rate_SISO_fixed_Q_center_match_arr=zeros(length(BW_arr),1); %% conjugate matching at center frequency
achievable_rate_SISO_optimized_circuit_arr=zeros(length(BW_arr),1); %% Circuit optimized using ADS for Shannon+ Bode/Fano bound
achievable_rate_SISO_freq_flat_circuit_arr=zeros(length(BW_arr),1); %% Circuit optimized using ADS for frequency-flat matching

%% What to plot?
%% 1 indicates to plot, 0 indicated no plot

plot_SNR_optimized=1;
plot_SNR_optimized_ckt=1;
plot_SNR_no_MN=1;
plot_SNR_freq_flat=0;
plot_SNR_freq_flat_ckt=1;
plot_SNR_fixed_Q_center_match=1;

optimized_ckt_ADS_file_available=1;
ff_ckt_ADS_file_available=1;

for bwi=1:length(BW_arr)
    
    BW= BW_arr(bwi);
    f_min= f_c-BW/2; %% smallest frequency in the band
    f_max= f_c+ BW/2; %% highest frequency in the band
    f_discrete= [f_min: (f_max-f_min)/500: f_max];
    sigma_sq= kb*T*BW; %% Noise variance
    Ps=250*10^(-3); %% Signal power
    Ps_by_sigma_sq=Ps/sigma_sq;

    
    
    %% Simulation of Chu's antenna 
    a= lambda_c/10;  %% size of Chu's antenna
    %% Define scattering parameter as function of frequency over a specified bandwidth
    S_T=  @(f) 1./(2*a^2/c^2.*(1j*2*pi*f).^2  + 2*a/c.*(1j*2*pi*f)+1);
    z11= @(f)  1./((1j*2*pi.*f).*a/(c*R))+ 1./(1./((1j*2*pi.*f).*a*R/c)+ 1/R);
    Z_RT=  @(f) c./(2*pi.*f*Dtxrx).*G.*real(z11(f));
    %% Compute |1- S_T|^2 as function of frequency over a specified bandwidth
    Tx_trans_T= @(f) (abs(1- S_T(f))).^2.*(Z_RT(f)).^2/Z_0^2./(abs(1+ z11(f)/Z_0)).^2;
    %% Compute effective SNR and equivalent SNR
    PSD= @(f) Ps_by_sigma_sq.*((f<=f_max)&(f>=f_min));
    SNR_eff=  @(f) Tx_trans_T(f).*PSD(f);
    SNR_eq= @(f) SNR_eff(f)./(1-abs(S_T(f)).^2);
    
    % Compute ideal Shannon achievable rate
    Transmission_coeff_ideal= @(f)  (f<inf)&(f> 0);  
    
    [achievable_rate_SISO_ideal, SE_ideal, SNR_ideal]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_ideal, f_min, f_max);
    achievable_rate_SISO_ideal_arr(bwi)=achievable_rate_SISO_ideal;
    %% Plot ideal SE as function of frequency
    figure
    
    plot(f_discrete/10^9, 10*log10(abs(SNR_ideal(f_discrete))), '--', 'LineWidth',2 , 'DisplayName','Upper bound (Shannon)');
    hold on
    
    % Bode fano constraints
    B1=2*a/c;
    B2= 4/3*a^3/c^3;
    f1= @(f) f.^(-2)/(2*pi^2);
    f2= @(f) f.^(-4)/(8*pi^4);
    % Optimized transmission coefficient as function of frequency and 2 Lagrange parameters
    
    Tau= @(f, mu_1, mu_2) max(0, (1+log(2)./SNR_eq(f).*(mu_1.*f1(f)+ mu_2.*f2(f)))./(1- log(2)*(mu_1.*f1(f)+ mu_2.*f2(f))));
    
    % Constraint equations
    %% Determine optimal mu1 and mu2
    %% Fix mu 1 = 0 and run bisection search on mu 2
    mu_1=0;
    mu_2_max=0;
    SNR_eq_max= max(SNR_eq(f_discrete));
    mu_2_min= -SNR_eq_max/log(2)/f2(f_c);
    Max_iter_b_search=1000;
    converged_flag_mu2search=0;
    C1_max_mu2= integral(  @(f) (f1(f).*log(1./(1-Tau(f, mu_1, mu_2_max)))), 0,Inf, 'Waypoints', [f_min, f_max], 'AbsTol', 10^(-13), 'RelTol',10^(-10)) ;
    C2_max_mu2= integral(  @(f) (f2(f).*log(1./(1-Tau(f, mu_1, mu_2_max)))), 0, Inf,'Waypoints', [f_min, f_max], 'AbsTol', 10^(-13), 'RelTol',10^(-10)) ;
    C1_min_mu2= integral(  @(f) (f1(f).*log(1./(1-Tau(f, mu_1, mu_2_min)))), 0,Inf, 'Waypoints', [f_min, f_max], 'AbsTol', 10^(-13), 'RelTol',10^(-10)) ;
    C2_min_mu2= integral(  @(f) (f2(f).*log(1./(1-Tau(f, mu_1, mu_2_min)))), 0, Inf,'Waypoints', [f_min, f_max], 'AbsTol', 10^(-13), 'RelTol',10^(-10)) ;
    for iter=1:Max_iter_b_search
        mu_2_mid= (mu_2_min+mu_2_max)/2;
        C1_mid_mu2= integral(  @(f) (f1(f).*log(1./(1-Tau(f, mu_1, mu_2_mid)))), 0,Inf, 'Waypoints', [f_min, f_max], 'AbsTol', 10^(-13), 'RelTol',10^(-10)) ;
        C2_mid_mu2= integral(  @(f) (f2(f).*log(1./(1-Tau(f, mu_1, mu_2_mid)))), 0, Inf,'Waypoints', [f_min, f_max], 'AbsTol', 10^(-13), 'RelTol',10^(-10)) ;
        if (C1_mid_mu2 -B1<0 ) && (C2_mid_mu2-B2<0)
            mu_2_min=mu_2_mid;
            if (abs(C2_mid_mu2 -B2)/B2 < 10^(-12))
                mu_2_converged=mu_2_mid;
                %fprintf("mu 2 search has converged")
                converged_flag_mu2search=1;
                break;
            end
        else 
            mu_2_max= mu_2_mid;
        end
    end
    mu_1_fixed=0;
    if converged_flag_mu2search==1
        Transmission_coeff_opt_mu2opt= @(f) Tau(f, mu_1_fixed, mu_2_converged);
        [achievable_rate_SISO_optimized_mu2opt, SE_optimized]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_opt_mu2opt, f_min, f_max);
    end
    %% Fix mu 2 = 0 and run bisection search on mu 1
    mu_2=0;
    mu_1_max=0;
    SNR_eq_max= max(SNR_eq(f_discrete));
    mu_1_min= -SNR_eq_max/log(2)/f1(f_c);
    Max_iter_b_search=5000;
    converged_flag_mu1search=0;
    C1_max_mu1= integral(  @(f) (f1(f).*log(1./(1-Tau(f, mu_1_max, mu_2)))), 0,Inf, 'Waypoints', [f_min, f_max], 'AbsTol', 10^(-13), 'RelTol',10^(-10)) ;
    C2_max_mu1= integral(  @(f) (f2(f).*log(1./(1-Tau(f, mu_1_max, mu_2)))), 0, Inf,'Waypoints', [f_min, f_max], 'AbsTol', 10^(-13), 'RelTol',10^(-10)) ;
    C1_min_mu1= integral(  @(f) (f1(f).*log(1./(1-Tau(f, mu_1_min, mu_2)))), 0,Inf, 'Waypoints', [f_min, f_max], 'AbsTol', 10^(-13), 'RelTol',10^(-10)) ;
    C2_min_mu1= integral(  @(f) (f2(f).*log(1./(1-Tau(f, mu_1_min, mu_2)))), 0, Inf,'Waypoints', [f_min, f_max], 'AbsTol', 10^(-13), 'RelTol',10^(-10)) ;
    for iter=1:Max_iter_b_search
        mu_1_mid= (mu_1_min+mu_1_max)/2;
        C1_mid_mu1= integral(  @(f) (f1(f).*log(1./(1-Tau(f, mu_1_mid, mu_2)))), 0,Inf, 'Waypoints', [f_min, f_max], 'AbsTol', 10^(-13), 'RelTol',10^(-10)) ;
        C2_mid_mu1= integral(  @(f) (f2(f).*log(1./(1-Tau(f, mu_1_mid, mu_2)))), 0, Inf,'Waypoints', [f_min, f_max], 'AbsTol', 10^(-13), 'RelTol',10^(-10)) ;
        if (C1_mid_mu1 -B1<0 ) && (C2_mid_mu1-B2<0)
            mu_1_min=mu_1_mid;
            if (abs(C1_mid_mu1 -B1)/B1 < 10^(-12))
                mu_1_converged=mu_1_mid;
                %fprintf("mu 1 search has converged")
                converged_flag_mu1search=1;

                break;
            end
        else 
            mu_1_max= mu_1_mid;
        end
    end
    mu_2_fixed=0;
    if converged_flag_mu1search==1
        Transmission_coeff_opt_mu1opt= @(f) Tau(f, mu_1_converged, mu_2_fixed);
        [achievable_rate_SISO_optimized_mu1opt, SE_optimized]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_opt_mu1opt, f_min, f_max);
    end

    if converged_flag_mu1search==1 && converged_flag_mu2search==1
        if achievable_rate_SISO_optimized_mu1opt > achievable_rate_SISO_optimized_mu2opt
            mu_1_optimal= mu_1_converged;
            mu_2_optimal= mu_2_fixed;
        else
            mu_1_optimal= mu_1_fixed;
            mu_2_optimal= mu_2_converged;
        end
    elseif converged_flag_mu1search==1
        mu_1_optimal= mu_1_converged;
        mu_2_optimal= mu_2_fixed;
    elseif converged_flag_mu2search==1
        mu_1_optimal= mu_1_fixed;
        mu_2_optimal= mu_2_converged;
    else
        fprintf("No convergence")
    end


    Transmission_coeff_opt= @(f) Tau(f, mu_1_optimal, mu_2_optimal);
    quantize_tx_coeff(f_discrete, Transmission_coeff_opt);
    [achievable_rate_SISO_optimized, SE_optimized, SNR_optimized]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_opt, f_min, f_max);
    
    if plot_SNR_optimized==1
        plot(f_discrete/10^9, 10*log10(abs(SNR_optimized(f_discrete))), 'LineWidth',2 , 'DisplayName', 'Upper bound (Shannon and Bode/Fano)');
    end
    achievable_rate_SISO_optimized_arr(bwi)=achievable_rate_SISO_optimized;
    clear mu_1_optimal
    clear mu_2_optimal
    %% ADS circuit
    if f_c== 7*10^9 || f_c==5*10^9
        if f_c==7*10^9
            tdfread("SISO_7GHz_responses/fc7G_BW"+ num2str((BW/f_c*10))+ ".txt");
        elseif f_c== 5*10^9
            tdfread("SISO_5GHz_responses/fc5G_BW"+ num2str((BW/f_c*10))+ ".txt");
        end
        achievable_rate_SISO_optimized_circuit_arr(bwi)=trapz(freq, log2(1+ SNR_eq(freq).*Eq_Tx));
        SE_optimized_ckt= log2(1+ SNR_eq(freq).*Eq_Tx);
        SNR_optimized_ckt = SNR_eq(freq).*Eq_Tx;
        Eq_Tx_optimized_ckt=Eq_Tx;
    else
        fprintf("Circuit response for "+ num2str(f_c/10^9) + " GHz is not available");
        optimized_ckt_ADS_file_available=0;
    end
    

    if plot_SNR_optimized_ckt==1  && optimized_ckt_ADS_file_available==1
        plot(freq/10^9, 10*log10(abs(SNR_optimized_ckt)), '-.', 'LineWidth',2 ,'DisplayName', 'Optimal transmission (circuit in ADS)');
    end
    %% Frequency flat matching which satisfies the Bode Fano constraints
    %% Boxcar matching: reference paper: Diversity Limits of Compact Broadband Multi-Antenna Systems
    reflection_cf_1= exp(-4*a/c*pi^2./(integral( @(f) f.^(-2)  , f_min, f_max) ));
    reflection_cf_2= exp(-32*pi^4*a^3/(3*c^3)./(integral(@(f) f.^(-4)  , f_min, f_max)  ));
    reflection_cf_freq_flat= max(reflection_cf_1, reflection_cf_2);
    Transmission_coeff_freq_flat= @(f)  (1-reflection_cf_freq_flat)*((f<=f_max)&(f>=f_min)); 

     
    [achievable_rate_SISO_freq_flat, SE_freq_flat, SNR_freq_flat]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_freq_flat, f_min, f_max);
    if plot_SNR_freq_flat==1
        plot(f_discrete/10^9, 10*log10(abs(SNR_freq_flat(f_discrete))), 'DisplayName','Frequency-flat transmission coefficient');
    end
    achievable_rate_SISO_freq_flat_arr(bwi)= achievable_rate_SISO_freq_flat;
    if f_c== 7*10^9 
        
        tdfread("SISO_7GHz_responses/ff_fc7G_BW"+ num2str((BW/f_c*10))+ ".txt")
        
        achievable_rate_SISO_freq_flat_circuit_arr(bwi)=trapz(freq, log2(1+ SNR_eq(freq).*Eq_Tx));
        Eq_Tx_ff_ckt= Eq_Tx;
        SE_freq_flat_ckt= log2(1+ SNR_eq(freq).*Eq_Tx);
        SNR_freq_flat_ckt = SNR_eq(freq).*Eq_Tx;
        if plot_SNR_freq_flat_ckt==1
            plot(freq/10^9, 10*log10(abs(SNR_freq_flat_ckt)), '-.', 'LineWidth',2 ,'DisplayName', 'Frequency-flat transmission (circuit in ADS)');
        end
    else
        fprintf("Circuit response for "+ num2str(f_c/10^9) + " GHz is not available");
        ff_ckt_ADS_file_available=0;
    end
    
    

    [achievable_rate_SISO_no_MN, SE_no_MN, SNR_no_MN]=compute_achievable_rate_SISO(SNR_eff, Transmission_coeff_ideal, f_min, f_max);
    achievable_rate_SISO_no_MN_arr(bwi)=achievable_rate_SISO_no_MN;
 
    if plot_SNR_no_MN==1
        plot(f_discrete/10^9, 10*log10(abs(SNR_no_MN(f_discrete))), ':', 'LineWidth',2 ,'DisplayName','No matching');
    end
    
    [tau_fixed_Q_center_match]= fixed_Q_center_match(a, Z_0, f_c, S_T);
    [achievable_rate_SISO_fixed_Q_center_match, SE_fixed_Q_center_match,  SNR_fixed_Q_center_match]=compute_achievable_rate_SISO(SNR_eq, tau_fixed_Q_center_match, f_min, f_max);
    
    if plot_SNR_fixed_Q_center_match==1
        plot(f_discrete/10^9, 10*log10(abs(SNR_fixed_Q_center_match(f_discrete))), ':','LineWidth',2 ,'DisplayName','Conjugate matching at center frequency');
    end
    achievable_rate_SISO_fixed_Q_center_match_arr(bwi)= achievable_rate_SISO_fixed_Q_center_match;
    legend
    grid on
    xlabel('Frequency (in GHz)')
    ylabel('Signal to noise ratio (in dB)')
    xlim([f_min/10^9, f_max/10^9]);
    title("Signal to noise ratio versus frequency for a bandwidth of " + num2str(BW/10^9)+ " GHz");
    legend('Location','best')
    if optimized_ckt_ADS_file_available==1 || ff_ckt_ADS_file_available==1
        %% Plot transmission coefficient as function of frequency
        figure
        plot(freq/10^9, Transmission_coeff_opt(freq),'LineWidth',2 ,'DisplayName', 'Optimal transmission (theoretical)');
        hold on 
        if optimized_ckt_ADS_file_available==1
            plot(freq/10^9,  Eq_Tx_optimized_ckt, '--', 'LineWidth',2 ,'DisplayName', 'Optimal transmission (circuit in ADS)');
        end   
    
        plot(freq/10^9, Transmission_coeff_freq_flat(freq),'LineWidth',2 ,'DisplayName', 'Frequency-flat transmission (theoretical)');
        if ff_ckt_ADS_file_available==1
            plot(freq/10^9, Eq_Tx_ff_ckt, '--', 'LineWidth',2 ,'DisplayName', 'Frequency-flat transmission (circuit in ADS)');
        end
        plot(freq/10^9, tau_fixed_Q_center_match(freq), 'LineWidth',2 ,'DisplayName', 'Conjugate matching');
        xlim([f_min/10^9, f_max/10^9]);
        grid on
        legend('Location','best')
        xlabel('Frequency (in GHz)')
        ylabel('Transmission coefficient (linear scale)')
        title("Transmission coefficient versus frequency for a bandwidth of "+num2str(BW/10^9)+ " GHz");
    end
end

figure

plot(BW_arr/10^9,achievable_rate_SISO_ideal_arr/10^9, '-o', 'MarkerSize',12, 'LineWidth',2 , 'DisplayName', 'Upper bound (Shannon)'  )
hold on
if plot_SNR_optimized==1
    plot(BW_arr/10^9, achievable_rate_SISO_optimized_arr/10^9, '-x', 'MarkerSize',12,  'LineWidth',2,  'DisplayName','Upper bound (Shannon and Bode/Fano)');
end
if optimized_ckt_ADS_file_available==1
    plot(BW_arr/10^9, achievable_rate_SISO_optimized_circuit_arr/10^9, '--sq', 'MarkerSize',12,  'LineWidth',2, 'DisplayName','Optimal transmission (circuit in ADS)');
end
if plot_SNR_freq_flat==1
    plot(BW_arr/10^9, achievable_rate_SISO_freq_flat_arr/10^9, '-<', 'MarkerSize',12, 'LineWidth',2, 'DisplayName','Frequency-flat transmission coefficient');
end
if ff_ckt_ADS_file_available==1
    plot(BW_arr/10^9, achievable_rate_SISO_freq_flat_circuit_arr/10^9, '--*','MarkerSize',12,  'LineWidth',2, 'DisplayName','Frequency-flat transmission (circuit in ADS)');
end
if plot_SNR_no_MN==1
    plot(BW_arr/10^9,achievable_rate_SISO_no_MN_arr/10^9,'-.^' , 'MarkerSize',12, 'LineWidth',2, 'DisplayName', 'No matching' )
end
if plot_SNR_fixed_Q_center_match==1
    plot(BW_arr/10^9,achievable_rate_SISO_fixed_Q_center_match_arr/10^9, '--v', 'MarkerSize',12,  'LineWidth',2, 'DisplayName', 'Conjugate matching' )
end
legend('Location','best')
grid on
ylabel('Achievable rate (Gbits/sec)')
xlabel('Bandwidth (in GHz)')
title("Achievable rate as function of bandwidth for center frequency " + num2str(f_c/10^9)+ " GHz")


function [achievable_rate_SISO, SE, SNR_final]= compute_achievable_rate_SISO(SNR_eq, Transmission_coeff, f_min, f_max)

SE= @(f) log2(1+ SNR_eq(f).*Transmission_coeff(f));
SNR_final= @(f) SNR_eq(f).*Transmission_coeff(f);
achievable_rate_SISO= integral(SE, f_min, f_max);
end


function [tau_fixed_Q_center_match]= fixed_Q_center_match(a, R, f_c, S_eq)
c= 3*10^8;  %% speed of light
omega_c=2*pi*f_c;

R_eff= R*(1- 1/(1+omega_c^2*a^2/c^2));
C_eff= a/(c*R)*(1+omega_c^2*a^2/c^2);
Q= sqrt(R/R_eff-1);
C_m=Q/(omega_c*R);
L_res= 1/(omega_c^2*C_eff);
L1= Q*R_eff/omega_c;
L_m= L1+ L_res;
z11= @(f) 1./(1j*(2*pi.*f)*C_m);
z12= @(f) 1./(1j*(2*pi.*f)*C_m);
z21= @(f) 1./(1j*(2*pi.*f)*C_m);
z22= @(f) 1j*(2*pi.*f)*L_m+ 1./(1j*(2*pi.*f)*C_m);
delta_z= @(f) (z22(f) +R).*(z11(f)+R)- z12(f).*z21(f);
S_21= @(f) 2*z12(f)*R./delta_z(f);
S_22= @(f) ( (z11(f) +R).*(z22(f)-R) - z12(f).*z21(f))./delta_z(f);
tau_fixed_Q_center_match= @(f) abs(S_21(f)).^2./(abs(1- S_22(f).*S_eq(f))).^2.*(1-(abs(S_eq(f))).^2);

end

function quantize_tx_coeff(f_arr, tau)

sub_groups_num=7;
total_points=length(f_arr);
tau_lower=zeros(sub_groups_num,1);
tau_upper=zeros(sub_groups_num,1);
for n=1:sub_groups_num
    if n< sub_groups_num
        f_sub= f_arr((n-1)*floor(total_points/sub_groups_num)+1: (n)*floor(total_points/sub_groups_num) );
        
    else
        f_sub= f_arr((n-1)*floor(total_points/sub_groups_num)+1: end);
    end
    f_min= min(f_sub);
    f_max= max(f_sub);
    tau_lower(n)= tau(f_min);
    tau_upper(n)= tau(f_max);
    %fprintf("%f %f  %f %f ", f_min/10^9, f_max/10^9,  tau_lower(n), tau_upper(n))
end

end