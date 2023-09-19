clear all;
close all;
%% Universal constants
c=3*10^8;  %% speed of light
kb= 1.380649*10^(-23); %% Boltzmann constant
T=290;  %% Temperature in Kelvin

%% Other constants
R=50; %% reference impedance
Z_0=R; %% reference impedance
f_c= 7*10^9;  %% center frequency
lambda_c= c/f_c;  %% center wavelength


BW_arr=[0.1:0.1:0.6]*f_c; %% range of bandwidths


%% Achievable rate for different cases
achievable_rate_MISO_ideal_arr=zeros(length(BW_arr),1);   %% Ideal -> No Bode/Fano theory
achievable_rate_MISO_optimized_arr=zeros(length(BW_arr),1);  %% Shannon+ Bode/Fano theory
achievable_rate_MISO_fixed_Q_center_match_arr=zeros(length(BW_arr),1); %% conjugate matching at center frequency
achievable_rate_MISO_no_MN_arr=zeros(length(BW_arr),1); %% no matching
achievable_rate_MISO_freq_flat_arr=zeros(length(BW_arr),1); %% Frequency-flat matching
achievable_rate_MISO_optimized_circuit_arr=zeros(length(BW_arr),1); %% Circuit optimized using ADS for Shannon+ Bode/Fano bound
achievable_rate_MISO_freq_flat_circuit_arr=zeros(length(BW_arr),1); %% Circuit optimized using ADS for frequency-flat matching

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

%% Beamforming mode
mode=1;  %% 1 is odd and 0 is even
if mode==1
    theta_ff=pi/2; %% 0 indicates broadside and pi/2 indicated endfire
elseif mode==0
    theta_ff=0; %% 0 indicates broadside and pi/2 indicated endfire
end




for bwi=1:length(BW_arr)
    BW= BW_arr(bwi); 
    f_min= f_c-BW/2;
    f_max= f_c+ BW/2;
    f_discrete= [f_min: (f_max-f_min)/500: f_max];
    
    %% Simulation of electric Chu's antenna array of 2 elements
    a= lambda_c/15;
    d=0.5*lambda_c;
    
    %% Define impedance parameters in terms of the  frequency 'f'
    %% self impedances
    z11= @(f)  1./((1j*2*pi.*f).*a/(c*R))+ 1./(1./((1j*2*pi.*f).*a*R/c)+ 1/R);
    z22= @(f)  1./((1j*2*pi.*f).*a/(c*R))+ 1./(1./((1j*2*pi.*f).*a*R/c)+ 1/R);
    
    %% Mutual impedances
    beta_angle=  pi/2;
    gamma_angle= pi/2;
    z21= @(f) -3*sqrt(real(z11(f)).*real(z22(f))).*(0.5*sin(beta_angle)*sin(gamma_angle)*(1./((1j*2*pi.*f).*d/c) +  1./((1j*2*pi.*f).*d/c).^2  +  1./((1j*2*pi.*f).*d/c).^3  ) ...
       + cos(gamma_angle)*cos(beta_angle)*(  1./((1j*2*pi.*f).*d/c).^2  +  1./((1j*2*pi.*f).*d/c).^3   ) ).*exp(-(1j*2*pi.*f).*d/c);
    z12= @(f) z21(f);
    
    %% Convert impedance matrix to scattering matrix
    delta_z= @(f) (z22(f) +R).*(z11(f)+R)- z12(f).*z21(f);
    S_21= @(f) 2*z12(f)*R./delta_z(f);
    S_12= @(f) S_21(f);
    S_22= @(f) ( (z11(f) +R).*(z22(f)-R) - z12(f).*z21(f))./delta_z(f);
    S_11= @(f) ( (z11(f) -R).*(z22(f)+R) - z12(f).*z21(f))./delta_z(f);
    
    %% Fit to a rational model
    f_fit_freq= [10^7:10^7:10*10^9];
    S_T_fit_freq= create_S_matrix_at_fit_freq(S_11, S_12, S_21, S_22, f_fit_freq);
    S_T_fit_func_1= rationalfit(f_fit_freq, S_T_fit_freq);
    sobj=sparameters(S_T_fit_freq, f_fit_freq);
    S_T_fit_func = makepassive(S_T_fit_func_1,sobj);
    shape_vec= size(S_T_fit_func(1,1).A);
    Np_fit=shape_vec(1);
    S_11_rational=S_T_fit_func(1,1).D;
    S_12_rational=S_T_fit_func(1,2).D;
    S_21_rational=S_T_fit_func(2,1).D;
    S_22_rational=S_T_fit_func(2,2).D;
    C_11=S_T_fit_func(1,1).C;
    C_12=S_T_fit_func(1,2).C;
    C_21=S_T_fit_func(2,1).C;
    C_22=S_T_fit_func(2,2).C;
    A_11=S_T_fit_func(1,1).A;
    A_12=S_T_fit_func(1,2).A;
    A_21=S_T_fit_func(2,1).A;
    A_22=S_T_fit_func(2,2).A;
    s=tf('s');
    for n=1:Np_fit
        S_11_rational = S_11_rational + C_11(n)/(s-A_11(n));
        S_12_rational = S_12_rational + C_12(n)/(s-A_12(n));
        S_21_rational = S_21_rational + C_21(n)/(s-A_21(n));
        S_22_rational = S_22_rational + C_22(n)/(s-A_22(n));
    end
    S_11_rational=S_11_rational*exp(-s*S_T_fit_func(1,1).delay);
    S_12_rational=S_12_rational*exp(-s*S_T_fit_func(1,2).delay);
    S_21_rational=S_21_rational*exp(-s*S_T_fit_func(2,1).delay);
    S_22_rational=S_22_rational*exp(-s*S_T_fit_func(2,2).delay);
    %% Verify the fitting
    [num_mat_11, den_mat_11 ]= create_num_den_coeff(S_11_rational);
    S_11_f=@(f) create_S_rational_func_f(num_mat_11, den_mat_11, f);
    [num_mat_12, den_mat_12 ]= create_num_den_coeff(S_12_rational);
    S_12_f=@(f) create_S_rational_func_f(num_mat_12, den_mat_12, f);
    [num_mat_21, den_mat_21 ]= create_num_den_coeff(S_21_rational);
    S_21_f=@(f) create_S_rational_func_f(num_mat_21, den_mat_21, f);
    [num_mat_22, den_mat_22 ]= create_num_den_coeff(S_22_rational);
    S_22_f=@(f) create_S_rational_func_f(num_mat_22, den_mat_22, f);

    
    
    %% Equivalent load expression
    theta_bf=2*pi.*f_c.*d*sin(theta_ff)/c;
    phi_bf=exp(1j*theta_bf);
    S_eq_rational= (S_11_rational+phi_bf*S_12_rational + phi_bf*S_21_rational + phi_bf^2*S_22_rational)/2;
    
    
    
    %% Wireless channel for MISO system
    Dtxrx=500;
    
    G=1.5;
    Z_RT_LOS_11= @(f)  c./(2*pi.*f*Dtxrx).*G.*real(z11(f));
    Z_RT_LOS_12= @(f)  c./(2*pi.*f*Dtxrx).*G.*real(z11(f)).*exp(-1j*2*pi.*f.*d*sin(theta_ff)/c);
    Z_RT_LOS   = @(f) [Z_RT_LOS_11(f), Z_RT_LOS_12(f) ].';
    %% Scattering matrix of equivalent load using even and odd mode beamformers
    S_eq= @(f) 0.5*(S_11(f)+ phi_bf*S_12(f)+ phi_bf*S_21(f)+ phi_bf^2*S_22(f));
    Z_eq= @(f) R*(1+ S_eq(f))./(1- S_eq(f));
    sigma_sq= kb*T*BW;
    Ps=250*10^(-3);
    Ps_by_sigma_sq=Ps/sigma_sq;
    
    PSD= @(f) Ps_by_sigma_sq.*((f<=f_max)&(f>=f_min));
    
    
    Tx_trans_T = @(f)  abs((Z_RT_LOS_11(f)).*(1-S_11(f)- phi_bf*S_12(f) )/sqrt(2)   + (Z_RT_LOS_12(f)).*(-S_21(f)+ phi_bf*(1- S_22(f)) )/sqrt(2)  ).^2/Z_0^2./(abs(1+ z11(f)/Z_0)).^2;
    SNR_eff=  @(f) Tx_trans_T(f).*PSD(f);
    SNR_eq= @(f) SNR_eff(f)./(1-abs(S_eq(f)).^2); 
    %% Reduce the model order
    S_rational= S_eq_rational;
    if mode==0
        [S_rational,info] = balred(S_rational, 2);
    elseif mode==1
        [S_rational,info] = balred(S_rational, 4);
    end
    [~, z]  =pzmap(S_rational);
   
    
    [num_load,den_load] = tfdata(S_rational);
    
    num_mat=cell2mat(num_load(1,1));
    den_mat=cell2mat(den_load(1,1));
    
    num_mat2=num_mat.*((-1).^([1:1:length(num_mat)]-1));
    den_mat2=den_mat.*((-1).^([1:1:length(den_mat)]-1));
   
    syms s
    S_rational_poly1 = poly2sym(num_mat,s)/poly2sym(den_mat,s);
    S_rational_poly2 = poly2sym(num_mat2,s)/poly2sym(den_mat2,s);
    
    
    S_load=vpasolve(S_rational_poly1*S_rational_poly2 -1==0, s);
    S_eq_sim= zeros(length(f_discrete), 1);
    for i=1:length(f_discrete)
        s=1j*2*pi*f_discrete(i);
        S_eq_sim(i)=subs(S_rational_poly1);
    end
    
    S_pos= double(S_load( (real(S_load)>0) ));  
    %% Bode Fano upper bounds
    S_rational_eval_at_S_pos=zeros(length(S_pos),1);
    B_BF=zeros(length(S_pos),1);
    for ei=1:length(S_pos)
        s=S_pos(ei);
        S_rational_eval_at_S_pos(ei)= subs(S_rational_poly1);
        Ni=1;
        Di=1;
        for z_even_i=1:length(z)
            Ni= Ni*(S_pos(ei)+ z(z_even_i));
            Di= Di*(S_pos(ei)- z(z_even_i));
        end
        B_BF(ei)=-pi/2*(log(abs(S_rational_eval_at_S_pos(ei)*Ni/Di )));
    end
    
    B1=B_BF(1);
    B2=B_BF(2);
    f_load1= @(f) pi/2*real(1./(S_pos(1) - 1j*2*pi.*f)   + 1./(S_pos(1) + 1j*2*pi.*f));
    f_load2= @(f) pi/2*real(1./(S_pos(2) - 1j*2*pi.*f)   + 1./(S_pos(2) + 1j*2*pi.*f));
     
    sum_mu_i_f_i=   @(f, mu_1, mu_2) mu_1.*f_load1(f)  +  mu_2.*f_load2(f);
    Tau_LOAD= @(f, mu_1, mu_2) max(0, (1+log(2)./SNR_eq(f).*( sum_mu_i_f_i(f, mu_1, mu_2)   ))./(1- log(2)*(sum_mu_i_f_i(f, mu_1, mu_2) ) ));
    
    Transmission_coeff_ideal= @(f)  (f<inf)&(f> 0);  
        
    [achievable_rate_MISO_ideal, SE_ideal,  SNR_ideal]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_ideal, f_min, f_max);
    achievable_rate_MISO_ideal_arr(bwi)=achievable_rate_MISO_ideal;
    
    %% Plot ideal SE as function of frequency
    figure
    plot(f_discrete/10^9, 10*log10(abs(SNR_ideal(f_discrete))), '--', 'LineWidth',2 , 'DisplayName','Upper bound (Shannon)');
    hold on
    
    % Constraint equations
    %% Determine optimal mu1 and mu2
    %% Fix mu 1 = 0 and run bisection search on mu 2
    mu_1=0;
    mu_2_max=0;
    SNR_eq_max= max(SNR_eq(f_discrete));
    mu_2_min= -SNR_eq_max/log(2)/f_load2(f_c);
    Max_iter_b_search=1000;
    converged_flag_mu2search=0;
    C1_max_mu2= integral(  @(f) (f_load1(f).*log(1./(1-Tau_LOAD(f, mu_1, mu_2_max)))), 0,10^20, 'Waypoints', [f_min, f_max]) ;
    C2_max_mu2= integral(  @(f) (f_load2(f).*log(1./(1-Tau_LOAD(f, mu_1, mu_2_max)))), 0, 10^20,'Waypoints', [f_min, f_max]) ;
    C1_min_mu2= integral(  @(f) (f_load1(f).*log(1./(1-Tau_LOAD(f, mu_1, mu_2_min)))), 0,10^20, 'Waypoints', [f_min, f_max]) ;
    C2_min_mu2= integral(  @(f) (f_load2(f).*log(1./(1-Tau_LOAD(f, mu_1, mu_2_min)))), 0, 10^20,'Waypoints', [f_min, f_max]) ;
    for iter=1:Max_iter_b_search
        mu_2_mid= (mu_2_min+mu_2_max)/2;
        C1_mid_mu2= integral(  @(f) (f_load1(f).*log(1./(1-Tau_LOAD(f, mu_1, mu_2_mid)))), 0,10^20, 'Waypoints', [f_min, f_max]) ;
        C2_mid_mu2= integral(  @(f) (f_load2(f).*log(1./(1-Tau_LOAD(f, mu_1, mu_2_mid)))), 0, 10^20,'Waypoints', [f_min, f_max]) ;
        if (C1_mid_mu2 - B1<0 ) && (C2_mid_mu2- B2<0)
            mu_2_min=mu_2_mid;
            if (abs(C2_mid_mu2 - B2)/B2 < 10^(-6))
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
        Transmission_coeff_opt_mu2opt= @(f) Tau_LOAD(f, mu_1_fixed, mu_2_converged);
        [achievable_rate_SISO_optimized_mu2opt, SE_optimized]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_opt_mu2opt, f_min, f_max);
    end
    %% Fix mu 2 = 0 and run bisection search on mu 1
    mu_2=0;
    mu_1_max=0;
    SNR_eq_max= max(SNR_eq(f_discrete));
    mu_1_min= -SNR_eq_max/log(2)/f_load1(f_c);
    Max_iter_b_search=1000;
    converged_flag_mu1search=0;
    C1_max_mu1= integral(  @(f) (f_load1(f).*log(1./(1-Tau_LOAD(f, mu_1_max, mu_2)))), 0,10^20, 'Waypoints', [f_min, f_max]) ;
    C2_max_mu1= integral(  @(f) (f_load2(f).*log(1./(1-Tau_LOAD(f, mu_1_max, mu_2)))), 0, 10^20,'Waypoints', [f_min, f_max]) ;
    C1_min_mu1= integral(  @(f) (f_load1(f).*log(1./(1-Tau_LOAD(f, mu_1_min, mu_2)))), 0,10^20, 'Waypoints', [f_min, f_max]) ;
    C2_min_mu1= integral(  @(f) (f_load2(f).*log(1./(1-Tau_LOAD(f, mu_1_min, mu_2)))), 0, 10^20,'Waypoints', [f_min, f_max]) ;
    for iter=1:Max_iter_b_search
        mu_1_mid= (mu_1_min+mu_1_max)/2;
        C1_mid_mu1= integral(  @(f) (f_load1(f).*log(1./(1-Tau_LOAD(f, mu_1_mid, mu_2)))), 0,10^20, 'Waypoints', [f_min, f_max]) ;
        C2_mid_mu1= integral(  @(f) (f_load2(f).*log(1./(1-Tau_LOAD(f, mu_1_mid, mu_2)))), 0, 10^20,'Waypoints', [f_min, f_max]) ;
        if (C1_mid_mu1 - B1<0 ) && (C2_mid_mu1-B2<0)
            mu_1_min=mu_1_mid;
            if (abs(C1_mid_mu1 -B1)/B1 < 10^(-6)) 
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
        Transmission_coeff_opt_mu1opt= @(f) Tau_LOAD(f, mu_1_converged, mu_2_fixed);
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
        %fprintf("No convergence")
    end

    
    Transmission_coeff_opt= @(f) Tau_LOAD(f, mu_1_optimal, mu_2_optimal);
    quantize_tx_coeff(f_discrete, Transmission_coeff_opt);
    [achievable_rate_MISO_optimized, SE_optimized, SNR_optimized]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_opt, f_min, f_max);
    
    plot(f_discrete/10^9, 10*log10(abs(SNR_optimized(f_discrete))), 'LineWidth',2 , 'DisplayName', 'Upper bound (Shannon and Bode/Fano)');
    achievable_rate_MISO_optimized_arr(bwi)=achievable_rate_MISO_optimized;
    clear mu_1_optimal
    clear mu_2_optimal
    %quantize_tx_coeff(f_discrete, Transmission_coeff_opt);
    %% ADS circuit
    if mode==0
        tdfread("MISO_7G_even/B"+ num2str((BW/f_c*10))+ "_7Gfc_even_MISO.txt")
    elseif mode==1
        tdfread("MISO_7G_odd/B"+ num2str((BW/f_c*10))+ "_7Gfc_odd_MISO.txt")
    end
    Eq_Tx_optimized_ckt=Eq_Tx;
    achievable_rate_MISO_optimized_circuit_arr(bwi)=trapz(freq, log2(1+ SNR_eq(freq).*Eq_Tx));
    SE_optimized_ckt= log2(1+ SNR_eq(freq).*Eq_Tx);
    SNR_optimized_ckt = SNR_eq(freq).*Eq_Tx;
    
    plot(freq/10^9, 10*log10(abs(SNR_optimized_ckt)), '-.', 'LineWidth',2 ,'DisplayName', 'Optimal transmission (circuit in ADS)');
    

    %% Frequency flat matching which satisfies the Bode Fano constraints
    %% Boxcar matching: reference paper: Diversity Limits of Compact Broadband Multi-Antenna Systems
    reflection_cf_1= exp(-B1./(integral( @(f) f_load1(f)  , f_min, f_max) ));
    reflection_cf_2= exp(-B2./(integral(@(f) f_load2(f)  , f_min, f_max)  ));
    reflection_cf_freq_flat= max(reflection_cf_1, reflection_cf_2);
    Transmission_coeff_freq_flat= @(f)  (1-reflection_cf_freq_flat)*((f<=f_max)&(f>=f_min)); 
    [achievable_rate_MISO_freq_flat, SE_freq_flat, SNR_freq_flat]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_freq_flat, f_min, f_max);
    
    %plot(f_discrete/10^9, 10*log10(abs(SNR_freq_flat(f_discrete))), 'DisplayName','Frequency-flat transmission coefficient')
    achievable_rate_MISO_freq_flat_arr(bwi)= achievable_rate_MISO_freq_flat;
    if mode==0
        tdfread("MISO_7G_even/B"+ num2str((BW/f_c*10))+ "_7Gfc_even_ff_MISO.txt")
    elseif mode==1
        tdfread("MISO_7G_odd/B"+ num2str((BW/f_c*10))+ "_7Gfc_odd_ff_MISO.txt")
    end
    Eq_Tx_ff_ckt= Eq_Tx;
    achievable_rate_MISO_freq_flat_circuit_arr(bwi)=trapz(freq, log2(1+ SNR_eq(freq).*Eq_Tx));
    SE_freq_flat_ckt= log2(1+ SNR_eq(freq).*Eq_Tx);
    SNR_freq_flat_ckt = SNR_eq(freq).*Eq_Tx;
    plot(freq/10^9, 10*log10(abs(SNR_freq_flat_ckt)), '-.', 'LineWidth',2 ,'DisplayName', 'Frequency-flat transmission (circuit in ADS)');

    %% No matching
    [achievable_rate_MISO_no_MN, SE_no_MN, SNR_no_MN]=compute_achievable_rate_SISO(SNR_eff, Transmission_coeff_ideal, f_min, f_max);
    achievable_rate_MISO_no_MN_arr(bwi)=achievable_rate_MISO_no_MN;
    
    plot(f_discrete/10^9, 10*log10(abs(SNR_no_MN(f_discrete))), ':', 'LineWidth',2 ,'DisplayName','No matching');
    
    %% Conjugate matching
    [tau_fixed_Q_center_match]= fixed_Q_center_match(Z_eq, R, f_c, S_eq);
    [achievable_rate_MISO_fixed_Q_center_match, SE_fixed_Q_center_match,   SNR_fixed_Q_center_match]=compute_achievable_rate_SISO(SNR_eq, tau_fixed_Q_center_match, f_min, f_max);
    plot(f_discrete/10^9, 10*log10(abs(SNR_fixed_Q_center_match(f_discrete))), ':','LineWidth',2 ,'DisplayName','Conjugate matching at center frequency');
    achievable_rate_MISO_fixed_Q_center_match_arr(bwi)= achievable_rate_MISO_fixed_Q_center_match;
    legend
    grid on
    xlabel('Frequency (in GHz)')
    ylabel('Signal to noise ratio (in dB)')
    xlim([f_min/10^9, f_max/10^9]);
    title("Bandwidth ="+num2str(BW/10^6)+ " MHz"+ ", normalized antenna size="+num2str(a/lambda_c))
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
plot(BW_arr/10^9,achievable_rate_MISO_ideal_arr/10^9, '-o', 'MarkerSize',12, 'LineWidth',2 , 'DisplayName', 'Upper bound (Shannon)'  )
hold on
plot(BW_arr/10^9, achievable_rate_MISO_optimized_arr/10^9, '-x', 'MarkerSize',12,  'LineWidth',2,  'DisplayName','Upper bound (Shannon and Bode/Fano)');
plot(BW_arr/10^9, achievable_rate_MISO_optimized_circuit_arr/10^9, '--sq', 'MarkerSize',12,  'LineWidth',2, 'DisplayName','Optimal transmission (circuit in ADS)');
%plot(BW_arr/10^9, achievable_rate_MISO_freq_flat_arr/10^9, '-<', 'MarkerSize',12, 'LineWidth',2, 'DisplayName','Frequency-flat transmission coefficient');
plot(BW_arr/10^9, achievable_rate_MISO_freq_flat_circuit_arr/10^9, '--*','MarkerSize',12,  'LineWidth',2, 'DisplayName','Frequency-flat transmission (circuit in ADS)');
plot(BW_arr/10^9,achievable_rate_MISO_no_MN_arr/10^9,'-.^' , 'MarkerSize',12, 'LineWidth',2, 'DisplayName', 'No matching' )
plot(BW_arr/10^9,achievable_rate_MISO_fixed_Q_center_match_arr/10^9, '--v', 'MarkerSize',12,  'LineWidth',2, 'DisplayName', 'Conjugate matching' )
legend('Location','best')
grid on
ylabel('Achievable rate (Gbits/sec)')
xlabel('Bandwidth (in GHz)')
title('Achievable rate as function of bandwidth for center frequency 7 GHz')

function S_T_fit_freq= create_S_matrix_at_fit_freq(S_11, S_12, S_21, S_22, f_fit_freq)

S_T_fit_freq=zeros(2,2, length(f_fit_freq));
for fi=1:length(f_fit_freq)
    fit_freq= f_fit_freq(fi);
    S_T_fit_freq(1,1, fi)= S_11(fit_freq);
    S_T_fit_freq(1,2, fi)= S_12(fit_freq);
    S_T_fit_freq(2,1, fi)= S_21(fit_freq);
    S_T_fit_freq(2,2, fi)= S_22(fit_freq);
end

end

function S_rational_func_f= create_S_rational_func_f(num_coeff, den_coeff, f)

N_num_coeff=length(num_coeff);
N_den_coeff= length(den_coeff);

N=0;
for ni=1:N_num_coeff
    N= N+ ((1j*2*pi.*f).^(N_num_coeff- ni))*double(num_coeff(ni));
end

D=0;
for di=1:N_den_coeff
    D= D+ ((1j*2*pi.*f).^(N_den_coeff- di))*double(den_coeff(di));
end
S_rational_func_f=N./D;
end

function [num, den]= create_num_den_coeff(S)

[num_cell,den_cell] = tfdata(S);
num=cell2mat(num_cell);
den=cell2mat(den_cell);

end

function [achievable_rate_SISO, SE, SNR_final]= compute_achievable_rate_SISO(SNR_eq, Transmission_coeff, f_min, f_max)

SE= @(f) log2(1+ SNR_eq(f).*Transmission_coeff(f));
SNR_final= @(f) SNR_eq(f).*Transmission_coeff(f);
achievable_rate_SISO= integral(SE, f_min, f_max);
end


function [tau_fixed_Q_center_match]= fixed_Q_center_match(Z, R, f_c, S_eq)
omega_c=2*pi*f_c;
R_eff= real(Z(f_c));
C_eff= -1./(omega_c*imag(Z(f_c)));
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
    fprintf("%f %f  %f %f ", f_min/10^9, f_max/10^9,  tau_lower(n), tau_upper(n))
end

end