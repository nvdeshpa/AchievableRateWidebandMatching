clear all;
close all;
rng default
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
plot_SNR_no_MN=0;
plot_SNR_freq_flat=0;
plot_SNR_freq_flat_ckt=0;
plot_SNR_fixed_Q_center_match=0;

optimized_ckt_ADS_file_available=1;
ff_ckt_ADS_file_available=1;
no_rx_impairment=1;

%% Beamforming mode
% mode=0;  %% 1 is odd and 0 is even
% if mode==1
%     theta_ff=pi/2; %% 0 indicates broadside and pi/2 indicated endfire
% elseif mode==0
%     theta_ff=0; %% 0 indicates broadside and pi/2 indicated endfire
% end
% %     %% Simulation of electric Chu's antenna array of 2 elements
% %     a= lambda_c/15;
% %     d=0.5*lambda_c;

%% Wireless channel for MISO system
Dtxrx=500;
theta_ff=pi/4;
include_multipath=1;
if include_multipath==1
    L=8;
    alpha_L= 1/sqrt(2)*(randn(1, L)  +1j*randn(1,L));
    delta=pi/20;
    theta_L = unifrnd(theta_ff-delta, theta_ff+delta, 1, L);
else
    L=0;
    alpha_L=0;
    theta_L =0;
end


%% Simulation of dipole antenna array of N elements
d=0.5*lambda_c;
l=lambda_c/4;
a=lambda_c/500; %% radius of thin wire antenna
%% Fit to a rational model
f_fit_freq= [4*10^9:10^7:10*10^9];
%N=4; %% number of antennas
N=10; %% number of antennas
initial_compute=1;


%% Scattering parameters from MATLAB Antenna toolbox
width_d=0.1*l;


if initial_compute==1
    [sobj, S_T_fit_freq, Directivity_linear_scale]= scattering_parameters_linear_dipole_array(l, width_d, d, f_fit_freq, N);
    save("N10_d_pt5_l_lambda_by_4_w_l_by_10", "S_T_fit_freq", "Directivity_linear_scale");
else
    data_mat= load('N8_d_pt5_l_lambda_by_4_w_l_by_10.mat');
    S_T_fit_freq=data_mat.S_T_fit_freq;
    Directivity_linear_scale=data_mat.Directivity_linear_scale;
end
Z_T_fit_freq = s2z(S_T_fit_freq,R);
%% Equivalent load expression
theta_bf=2*pi.*f_c.*d*sin(theta_ff)/c;
phi_bf=exp(1j*theta_bf);
f_c_index= (length(f_fit_freq)+1)/2;
beamformer= create_beamformer_optimal_discrete_set(N, f_c, f_c_index, S_T_fit_freq, Z_T_fit_freq, Dtxrx, d, l, theta_ff,  Z_0, Directivity_linear_scale, include_multipath, alpha_L, theta_L );
%beamformer= [1; phi_bf;  phi_bf^2; phi_bf^3]/sqrt(4);
%beamformer= [1; phi_bf;  ]/sqrt(2);
S_eq=  create_equivalent_load_N_dipole_antennas( f_fit_freq, beamformer, S_T_fit_freq);
rfwrite(S_eq,f_fit_freq,'N8_d_pt5_l_lambda_by_4_w_l_by_10_theta_0.s2p')
Z_eq=  R*(1+ S_eq)./(1- S_eq);
S_T_fit_func_1= rationalfit(f_fit_freq, S_eq);
sobj=sparameters(S_eq, f_fit_freq);
S_T_fit_func = makepassive(S_T_fit_func_1,sobj);
S_eq_rational = create_rational_form_equivalent_load(S_T_fit_func, 1, 1);
%% Verify the fitting
[num_mat_11, den_mat_11 ]= create_num_den_coeff(S_eq_rational);
S_f_validation=@(f) create_S_rational_func_f(num_mat_11, den_mat_11, f);

%% Reduce the model order
S_rational= S_eq_rational;
% reduced_order=4;
% [S_rational,info] = balred(S_rational, reduced_order);

[p, z]  =pzmap(S_rational);
[num_load,den_load] = tfdata(S_rational);

num_mat=cell2mat(num_load(1,1));
den_mat=cell2mat(den_load(1,1));

num_mat2=num_mat.*((-1).^([1:1:length(num_mat)]-1));
den_mat2=den_mat.*((-1).^([1:1:length(den_mat)]-1));

syms s
S_rational_poly1 = poly2sym(num_mat,s)/poly2sym(den_mat,s);
S_rational_poly2 = poly2sym(num_mat2,s)/poly2sym(den_mat2,s);


S_load=vpasolve(S_rational_poly1*S_rational_poly2 -1==0, s);
S_eq_sim= zeros(length(f_fit_freq), 1);
for i=1:length(f_fit_freq)
    s=1j*2*pi*f_fit_freq(i);
    S_eq_sim(i)=subs(S_rational_poly1);
end

% s_grid_real=[-5*10^10:10^9: 10^9 ];
% s_grid_imag=[-10^11: 10^9:10^11];
% s_grid=zeros(length(s_grid_real), length(s_grid_imag));
% s_grid_decision=zeros(length(s_grid_real), length(s_grid_imag));
% for real_i=1:length(s_grid_real)
%     for imag_i = 1:length(s_grid_imag)
%         s= s_grid_real(real_i)+1j*s_grid_imag(imag_i);
%         s_grid(real_i, imag_i)=subs(S_rational_poly1);
%         if  abs(s_grid(real_i, imag_i))<=1
%             s_grid_decision(real_i, imag_i)=1;
%         end
%     end
% end

S_pos= double(S_load( (real(S_load)>0) ));  
%% Bode Fano upper bounds
S_rational_eval_at_S_pos=zeros(length(S_pos),1);
B_BF=zeros(length(S_pos),1);
for ei=1:length(S_pos)
    s=S_pos(ei);
%     if abs(real(s))/abs(imag(s))<10^(-8)  
%         %% Bound 1
%         Pinv_sum=0;
%         for p_i=1:length(p)
%             Pinv_sum=Pinv_sum+ 1/(p(p_i) - 1j*abs(imag(s)));
%         end
%         Zinv_sum=0;
%         for z_i=1:length(z)
%             Zinv_sum= Zinv_sum+ 1/(z(z_i)+ 1j*abs(imag(s)));
%         end
%         B_BF(ei)=-pi/2*(Pinv_sum+Zinv_sum );
%     else
        %% Bound 2
        S_rational_eval_at_S_pos(ei)= subs(S_rational_poly1);
        Ni=1;
        Di=1;
        for z_even_i=1:length(z)
            Ni= Ni*(S_pos(ei)+ z(z_even_i));
            Di= Di*(S_pos(ei)- z(z_even_i));
        end
        B_BF(ei)=-pi/2*(log(abs(S_rational_eval_at_S_pos(ei)*Ni/Di )));
    %end
end


num_con=length(S_pos);
f_load =  create_f_load_vector(S_pos,f_fit_freq,num_con);


%sum_mu_i_f_i =  @( mu_arr) create_sum_mu_i_f_i(mu_arr, f_fit_freq, S_pos, num_con);



for bwi=1:length(BW_arr)
    BW= BW_arr(bwi); 
    f_min= f_c-BW/2;
    f_max= f_c+ BW/2;
    f_c_index= (length(f_fit_freq)+1)/2;
    f_res= f_fit_freq(2)- f_fit_freq(1);
    f_max_index=f_c_index+ ceil(BW/2/f_res);
    f_min_index=f_c_index-ceil(BW/2/f_res);
    f_discrete= f_fit_freq(f_min_index:1:f_max_index);
    f_len=length(f_discrete);
    S_T_freq =  S_T_fit_freq(:,:, f_min_index:1:f_max_index);
    Z_T_freq = Z_T_fit_freq(:,:, f_min_index:1:f_max_index);
    sum_mu_i_f_i =  @( mu_arr) create_sum_mu_i_f_i(mu_arr, f_discrete, f_load(:, f_min_index:1:f_max_index ));
    %f_discrete= [f_min: (f_max-f_min)/500: f_max];  
    
    
    
    
    
    
    
    sigma_sq= kb*T*BW;
    
    
    NoiseFig=2;
    
    sigma_sq_rx_noise= kb*T*BW*(NoiseFig-1);
    Ps=250*10^(-3);
    Ps_by_sigma_sq=Ps/sigma_sq;
    tau_rx=0.7;
    if no_rx_impairment==1
        Noise=kb*T;
    else
        Noise=(kb*T)*(NoiseFig-1 + tau_rx )/tau_rx;
    end
    %Ps_by_sigma_sq_rx_impaired=tau_rx*Ps/(sigma_sq*tau_rx +sigma_sq_rx_noise);
    
    %PSD= @(f) Ps_by_sigma_sq.*((f<=f_max)&(f>=f_min));
    %PSD =  Ps_by_sigma_sq.*ones(1, f_len);
    
   
    
    Tx_trans_T = create_Tx_trans_T_N_dip(f_discrete, S_T_freq, Z_T_freq,  N, Dtxrx, d, l, theta_ff,  beamformer, Z_0, Directivity_linear_scale(f_min_index:1:f_max_index), include_multipath, alpha_L, theta_L);
    
    
    
    
    [PSD_water_fill, PSD_equal]= water_fill_allocation(Tx_trans_T./(1-abs(S_eq(f_min_index:1:f_max_index)).^2), Ps, Noise, f_discrete);
    equal_power_allocation=1;
    if equal_power_allocation==1
        PSD =PSD_equal;
    else
        PSD =PSD_water_fill;
    end
    SNR_eff=   Tx_trans_T.*PSD;
    %SNR_eff_water_fill=   Tx_trans_T.*PSD_water_fill;
    SNR_eq= SNR_eff./(1-abs(S_eq(f_min_index:1:f_max_index)).^2); 
    %SNR_eq_water_fill= SNR_eff_water_fill./(1-abs(S_eq(f_min_index:1:f_max_index)).^2); 
    Tau_LOAD= @(mu_arr) max(0, (1+log(2)./SNR_eq.*( sum_mu_i_f_i(mu_arr)   ))./(1- log(2)*(sum_mu_i_f_i(mu_arr) ) ));
    
    
    Transmission_coeff_ideal= ones(1, f_len);  
        
    [achievable_rate_MISO_ideal, SE_ideal,  SNR_ideal]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_ideal, f_discrete);
    %[achievable_rate_MISO_ideal_water_fill, SE_ideal_water_fill,  SNR_ideal_water_fill]=compute_achievable_rate_SISO(SNR_eq_water_fill, Transmission_coeff_ideal, f_discrete);
    achievable_rate_MISO_ideal_arr(bwi)=achievable_rate_MISO_ideal;
    
    %% Plot ideal SE as function of frequency
    figure
    plot(f_discrete/10^9, 10*log10(abs(SNR_ideal)), '--', 'LineWidth',2 , 'DisplayName','Upper bound (Shannon)');
    hold on
    
    
    %mu_arr_opt = compute_optimal_mu_arr(num_con, f_load(:, f_min_index:1:f_max_index ), B_BF, Tau_LOAD, SNR_eq, f_min, f_max, f_discrete, f_c);
    mu_arr_opt = compute_optimal_mu_arr(num_con, f_load(:, f_min_index:1:f_max_index ), B_BF, Tau_LOAD, SNR_eq, f_discrete);
    
    Transmission_coeff_opt=  Tau_LOAD(mu_arr_opt);
    quantize_tx_coeff(f_discrete, Transmission_coeff_opt, f_min_index, f_max_index );
    [achievable_rate_MISO_optimized, SE_optimized, SNR_optimized]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_opt, f_discrete);
    
    plot(f_discrete/10^9, 10*log10(abs(SNR_optimized)), 'LineWidth',2 , 'DisplayName', 'Upper bound (Shannon and Bode/Fano)');
    achievable_rate_MISO_optimized_arr(bwi)=achievable_rate_MISO_optimized;
    clear mu_1_optimal
    clear mu_2_optimal
   
    %% ADS circuit
%     if mode==0
%         tdfread("MISO_7G_even/B"+ num2str((BW/f_c*10))+ "_7Gfc_even_MISO.txt")
%     elseif mode==1
%         tdfread("MISO_7G_odd/B"+ num2str((BW/f_c*10))+ "_7Gfc_odd_MISO.txt")
%     end
    tdfread("N8_d_pt5_l_lambda_by_4_w_l_by_10_theta_pi_by_4_bw_PT" + num2str((BW/f_c*10))+"_Eq_Tx.txt")
    Eq_Tx_optimized_ckt=Eq_Tx;
    achievable_rate_MISO_optimized_circuit_arr(bwi)=trapz(freq(f_min_index:1:f_max_index), log2(1+ (SNR_eq.').*Eq_Tx(f_min_index:1:f_max_index)));
    SE_optimized_ckt= log2(1+ (SNR_eq.').*Eq_Tx(f_min_index:1:f_max_index));
    SNR_optimized_ckt = (SNR_eq.').*Eq_Tx(f_min_index:1:f_max_index);
    
    plot(freq(f_min_index:1:f_max_index)/10^9, 10*log10(abs(SNR_optimized_ckt)), '-.', 'LineWidth',2 ,'DisplayName', 'Optimal transmission (circuit in ADS)');
%     

    %% Frequency flat matching which satisfies the Bode Fano constraints
    %% Boxcar matching: reference paper: Diversity Limits of Compact Broadband Multi-Antenna Systems

    %Transmission_coeff_freq_flat= @(f) compute_tau_freq_flat(B_BF, f_load, num_con , f_min, f_max, f);
    Transmission_coeff_freq_flat=compute_tau_freq_flat(B_BF, f_load(:, f_min_index:1:f_max_index ), num_con ,  f_discrete);
    [achievable_rate_MISO_freq_flat, SE_freq_flat, SNR_freq_flat]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_freq_flat, f_discrete);
    
    plot(f_discrete/10^9, 10*log10(abs(SNR_freq_flat)), 'DisplayName','Frequency-flat transmission coefficient')
    achievable_rate_MISO_freq_flat_arr(bwi)= achievable_rate_MISO_freq_flat;
%     if mode==0
%         tdfread("MISO_7G_even/B"+ num2str((BW/f_c*10))+ "_7Gfc_even_ff_MISO.txt")
%     elseif mode==1
%         tdfread("MISO_7G_odd/B"+ num2str((BW/f_c*10))+ "_7Gfc_odd_ff_MISO.txt")
%     end
    tdfread("N8_d_pt5_l_lambda_by_4_w_l_by_10_theta_pi_by_4_bw_PT" + num2str((BW/f_c*10))+ "_Eq_Tx_ff.txt")
    Eq_Tx_ff_ckt= Eq_Tx;
    achievable_rate_MISO_freq_flat_circuit_arr(bwi)=trapz(freq(f_min_index:1:f_max_index), log2(1+ (SNR_eq.').*Eq_Tx(f_min_index:1:f_max_index)));
    SE_freq_flat_ckt= log2(1+ (SNR_eq.').*Eq_Tx(f_min_index:1:f_max_index));
    SNR_freq_flat_ckt = (SNR_eq.').*Eq_Tx(f_min_index:1:f_max_index);
    plot(freq(f_min_index:1:f_max_index)/10^9, 10*log10(abs(SNR_freq_flat_ckt)), '-.', 'LineWidth',2 ,'DisplayName', 'Frequency-flat transmission (circuit in ADS)');

    %% No matching
    [achievable_rate_MISO_no_MN, SE_no_MN, SNR_no_MN]=compute_achievable_rate_SISO(SNR_eff, Transmission_coeff_ideal, f_discrete);
    achievable_rate_MISO_no_MN_arr(bwi)=achievable_rate_MISO_no_MN;
    
    plot(f_discrete/10^9, 10*log10(abs(SNR_no_MN)), ':', 'LineWidth',2 ,'DisplayName','No matching');
    
    %% Conjugate matching
    [tau_fixed_Q_center_match]= fixed_Q_center_match(Z_eq, R, f_c, S_eq, f_c_index, f_discrete, f_min_index, f_max_index);
    [achievable_rate_MISO_fixed_Q_center_match, SE_fixed_Q_center_match,   SNR_fixed_Q_center_match]=compute_achievable_rate_SISO(SNR_eq, tau_fixed_Q_center_match, f_discrete);
    plot(f_discrete/10^9, 10*log10(abs(SNR_fixed_Q_center_match)), ':','LineWidth',2 ,'DisplayName','Conjugate matching at center frequency');
    achievable_rate_MISO_fixed_Q_center_match_arr(bwi)= achievable_rate_MISO_fixed_Q_center_match;
    legend
    grid on
    xlabel('Frequency (in GHz)')
    ylabel('Signal to noise ratio (in dB)')
    xlim([f_min/10^9, f_max/10^9]);
    title("Bandwidth ="+num2str(BW/10^6)+ " MHz"+ ", normalized antenna size="+num2str(a/lambda_c))
    legend('Location','best')
% 
%     if optimized_ckt_ADS_file_available==1 || ff_ckt_ADS_file_available==1
%         %% Plot transmission coefficient as function of frequency
%         figure
%         plot(freq/10^9, Transmission_coeff_opt(freq.'),'LineWidth',2 ,'DisplayName', 'Optimal transmission (theoretical)');
%         hold on 
%         if optimized_ckt_ADS_file_available==1
%             plot(freq/10^9,  Eq_Tx_optimized_ckt, '--', 'LineWidth',2 ,'DisplayName', 'Optimal transmission (circuit in ADS)');
%         end   
%     
%         plot(freq/10^9, Transmission_coeff_freq_flat(freq.'),'LineWidth',2 ,'DisplayName', 'Frequency-flat transmission (theoretical)');
%         if ff_ckt_ADS_file_available==1
%             plot(freq/10^9, Eq_Tx_ff_ckt, '--', 'LineWidth',2 ,'DisplayName', 'Frequency-flat transmission (circuit in ADS)');
%         end
%         plot(freq/10^9, tau_fixed_Q_center_match(freq.'), 'LineWidth',2 ,'DisplayName', 'Conjugate matching');
%         xlim([f_min/10^9, f_max/10^9]);
%         grid on
%         legend('Location','best')
%         xlabel('Frequency (in GHz)')
%         ylabel('Transmission coefficient (linear scale)')
%         title("Transmission coefficient versus frequency for a bandwidth of "+num2str(BW/10^9)+ " GHz");
%     end
end




figure
plot(BW_arr/10^9,achievable_rate_MISO_ideal_arr/10^9, '-o', 'MarkerSize',12, 'LineWidth',2 , 'DisplayName', 'Upper bound (Shannon)'  )
hold on
plot(BW_arr/10^9, achievable_rate_MISO_optimized_arr/10^9, '-x', 'MarkerSize',12,  'LineWidth',2,  'DisplayName','Upper bound (Shannon and Bode/Fano)');
plot(BW_arr/10^9, achievable_rate_MISO_optimized_circuit_arr/10^9, '--sq', 'MarkerSize',12,  'LineWidth',2, 'DisplayName','Optimal transmission (circuit in ADS)');
plot(BW_arr/10^9, achievable_rate_MISO_freq_flat_arr/10^9, '-<', 'MarkerSize',12, 'LineWidth',2, 'DisplayName','Frequency-flat transmission coefficient');
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

function [achievable_rate_SISO, SE, SNR_final]= compute_achievable_rate_SISO(SNR_eq, Transmission_coeff, f_arr)

% SE= @(f) log2(1+ SNR_eq(f).*Transmission_coeff(f));
% SNR_final= @(f) SNR_eq(f).*Transmission_coeff(f);
% achievable_rate_SISO= integral(SE, f_min, f_max);
SE= log2(1+ SNR_eq.*Transmission_coeff);
SNR_final=  SNR_eq.*Transmission_coeff;
achievable_rate_SISO= trapz(f_arr, SE);

end


function [tau_fixed_Q_center_match]= fixed_Q_center_match(Z, R, f_c, S_eq, f_c_index, f_arr, f_min_index, f_max_index)
omega_c=2*pi*f_c;
R_eff= real(Z(f_c_index));
C_eff= -1./(omega_c*imag(Z(f_c_index)));
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
tau_fixed_Q_center_match=  abs(S_21(f_arr)).^2./(abs(1- S_22(f_arr).*S_eq(f_min_index:1:f_max_index))).^2.*(1-(abs(S_eq(f_min_index:1:f_max_index))).^2);

end
function quantize_tx_coeff(f_arr, tau, f_min_index, f_max_index)

sub_groups_num=7;
total_points=length(f_arr);
tau_lower=zeros(sub_groups_num,1);
tau_upper=zeros(sub_groups_num,1);
for n=1:sub_groups_num
    if n< sub_groups_num
        f_min_i=(n-1)*floor(total_points/sub_groups_num)+1;
        f_max_i=(n)*floor(total_points/sub_groups_num);
        f_sub= f_arr(f_min_i: f_max_i );
        
    else
        f_min_i=(n-1)*floor(total_points/sub_groups_num)+1;
        f_max_i=length(f_arr);
        f_sub= f_arr((n-1)*floor(total_points/sub_groups_num)+1: end);
        
    end
    [f_min]= min(f_sub);
    [f_max]= max(f_sub);
    tau_lower(n)= tau(f_min_i);
    tau_upper(n)= tau(f_max_i);
    fprintf("%f %f  %f %f \n", f_min/10^9, f_max/10^9,  tau_lower(n), tau_upper(n))
end

end


function S = create_scattering_matrix_for_N_Chu_antennas_ULA(N, d , a, f_fit_freq)

R=50;
c=3*10^8;
Z_mat=zeros(N, N , length(f_fit_freq));
S= zeros(N,N, length(f_fit_freq));
for fi=1:length(f_fit_freq)
     fit_freq= f_fit_freq(fi);
     Z_s= 1./((1j*2*pi.*fit_freq).*a/(c*R))+ 1./(1./((1j*2*pi.*fit_freq).*a*R/c)+ 1/R); %% self impedance
     d_arr= d*[1:1:N-1];
     Z_m=  -1.5*abs(real(Z_s)).*(1./((1j*2*pi.*fit_freq).*d_arr/c) +  1./((1j*2*pi.*fit_freq).*d_arr/c).^2  +  1./((1j*2*pi.*fit_freq).*d_arr/c).^3).*exp(-(1j*2*pi.*fit_freq).*d_arr/c); 
     Z_vec= [Z_s, Z_m];
     Z_mat(:, :, fi)= toeplitz(Z_vec, Z_vec);
     S(:,:,fi)= inv(Z_mat(:, :, fi) + R*eye(N) )*(Z_mat(:, :, fi) - R*eye(N));
end


end

function Z= MI_sbs_dip(d2lambda, l2lambda)

mu_0=4*pi*10^(-7);
epsilon_0=8.854*10^(-12); 
eta=sqrt(mu_0/ epsilon_0);

u0=2*pi*d2lambda;
u1=2*pi*(sqrt(d2lambda.^2+ l2lambda^2) + l2lambda);
u2=2*pi*(sqrt(d2lambda.^2+ l2lambda^2) - l2lambda);

R= eta./(4*pi).*(2*cosint(u0)-cosint(u1)-cosint(u2));
X= -eta./(4*pi).*(2*sinint(u0)- sinint(u1)-sinint(u2));
Z=R+1j*X;
end

function Z= self_impedance_dip(k0, l, a)
C=0.5772; %% Euler constant 
epsilon=8.854*10^(-12); 
c=3*10^8; %% speed of light
R= (  C+ log(k0*l) - cosint(k0*l)   +0.5*sin(k0*l)*(sinint(2*k0*l)  - 2*sinint(k0*l)) +0.5*cos(k0*l)*(C+log(k0*l/2) ...
     + cosint(2*k0*l)  -2*cosint(k0*l))   )/(2*pi*epsilon*c);

X= ( 2*sinint(k0*l) + cos(k0*l)*(2*sinint(k0*l)- sinint(2*k0*l))- ...
    sin(k0*l)*(2*cosint(k0*l)  - cosint(2*k0*l)-   cosint(2*k0*a^2/l))  )/(4*pi*epsilon*c);

Z=R+1j*X;
end

function [S, Z_mat] = create_scattering_matrix_for_N_dipole_antennas_ULA(N, d , l, a, f_fit_freq, R)

%R=50;
c=3*10^8;
Z_mat=zeros(N, N , length(f_fit_freq));
S= zeros(N,N, length(f_fit_freq));
k_arr=2*pi*f_fit_freq/c;
for fi=1:length(f_fit_freq)
     fit_freq= f_fit_freq(fi);
     k=k_arr(fi);
     Z_s= self_impedance_dip(k, l, a);
     %Z_s= 1./((1j*2*pi.*fit_freq).*a/(c*R))+ 1./(1./((1j*2*pi.*fit_freq).*a*R/c)+ 1/R); %% self impedance
     d_arr= d*[1:1:N-1];
     Z_m= MI_sbs_dip(d_arr/c*fit_freq,l/c*fit_freq);
     %Z_m=  -1.5*abs(real(Z_s)).*(1./((1j*2*pi.*fit_freq).*d_arr/c) +  1./((1j*2*pi.*fit_freq).*d_arr/c).^2  +  1./((1j*2*pi.*fit_freq).*d_arr/c).^3).*exp(-(1j*2*pi.*fit_freq).*d_arr/c); 
     Z_vec= [Z_s, Z_m];
     Z_mat(:, :, fi)= toeplitz(Z_vec, Z_vec);
     S(:,:,fi)= inv(Z_mat(:, :, fi) + R*eye(N) )*(Z_mat(:, :, fi) - R*eye(N));
end


end

function S_eq_rational = create_rational_form_equivalent_load(S_T_fit_func, N, beamformer)
shape_vec= size(S_T_fit_func(1,1).A);
Np_fit=shape_vec(1);
C= zeros(N,N, Np_fit);
A= zeros(N,N, Np_fit);
s=tf('s');
S_rational=zeros(N)*s;


for n1=1:N
    for n2=1: N
        C(n1, n2, :)= S_T_fit_func(n1,n2).C;
        %C(n2, n1, :)= C(n1, n2, :);
        A(n1, n2, :)= S_T_fit_func(n1,n2).A;
        %A(n2, n1, :)= A(n1, n2, :);
        
        
        
        for n=1:Np_fit
             S_rational(n1, n2)= S_rational(n1,n2)+ C(n1, n2, n)/(s- A(n1,n2,n));
        end
        %S_rational(n1, n2)=S_rational(n2, n1);
    end
end


% C_11=S_T_fit_func(1,1).C;
% C_12=S_T_fit_func(1,2).C;
% C_21=S_T_fit_func(2,1).C;
% C_22=S_T_fit_func(2,2).C;
% A_11=S_T_fit_func(1,1).A;
% A_12=S_T_fit_func(1,2).A;
% A_21=S_T_fit_func(2,1).A;
% A_22=S_T_fit_func(2,2).A;
% 
% for n=1:Np_fit
%     S_rational= S_rational+ C(:, :, n)./(s- A(:,:,n));
% %     S_11_rational = S_11_rational + C_11(n)/(s-A_11(n));
% %     S_12_rational = S_12_rational + C_12(n)/(s-A_12(n));
% %     S_21_rational = S_21_rational + C_21(n)/(s-A_21(n));
% %     S_22_rational = S_22_rational + C_22(n)/(s-A_22(n));
% end

%S_eq_rational= (S_rational(1,1)+phi_bf*S_rational(1,2) + phi_bf*S_rational(2,1) + phi_bf^2*S_rational(2,2))/2;
S_rational=S_rational+S_T_fit_func.D;
S_eq_rational=  beamformer.'*S_rational*beamformer;
end

function S_eq= create_equivalent_load_N_antennas(N, d, a, f_arr, beamformer)
f_length= length(f_arr);
S_eq=zeros(1, f_length);
S_T_f=create_scattering_matrix_for_N_Chu_antennas_ULA(N, d , a, f_arr);
for fi=1:f_length
    S_eq(fi)= beamformer.'*S_T_f(:,:,fi)*beamformer;
end

end

function S_eq= create_equivalent_load_N_dipole_antennas(f_arr, beamformer, S_T_f)
f_length= length(f_arr);
S_eq=zeros(1, f_length);

for fi=1:f_length
    S_eq(fi)= beamformer.'*S_T_f(:,:,fi)*beamformer;
end

end



function Tx_trans_T = create_Tx_trans_T(f_arr,  N, Dtxrx, d, theta_ff, G, a, beamformer, Z_0, R, c)

f_length= length(f_arr);
Tx_trans_T=zeros(f_length,1);

Z_RT_LOS = create_LOS_channel(f_arr, N, Dtxrx, d, theta_ff, G , a, R);
S_T_freq = create_scattering_matrix_for_N_Chu_antennas_ULA(N, d , a, f_arr);

for fi=1:f_length
    freq_i= f_arr(fi);
    Z_s= 1./((1j*2*pi.*freq_i).*a/(c*R))+ 1./(1./((1j*2*pi.*freq_i).*a*R/c)+ 1/R);
    Tx_trans_T(fi)= abs(beamformer.'*(eye(N)-S_T_freq(:,:, fi))*Z_RT_LOS(fi,:).')^2/Z_0^2./(abs(1+ Z_s/Z_0)).^2;    
end
Tx_trans_T=Tx_trans_T.';
end

function Tx_trans_T = create_Tx_trans_T_N_dip(f_arr, S_T_freq, Z_T_freq,  N, Dtxrx, d, l, theta_ff, beamformer, Z_0, directivity_linear_scale, include_multipath, alpha_L, theta_L)

f_length= length(f_arr);
Tx_trans_T=zeros(f_length,1);

%Z_RT_LOS = create_LOS_channel(f_arr, N, Dtxrx, d, theta_ff, G , a, R);
Z_RT_LOS = create_LOS_channel_N_dip(f_arr, Z_T_freq,  N, Dtxrx, d, theta_ff,  l , directivity_linear_scale, include_multipath, alpha_L, theta_L);
%S_T_freq = create_scattering_matrix_for_N_Chu_antennas_ULA(N, d , a, f_arr);
%S_T_freq =create_scattering_matrix_for_N_dipole_antennas_ULA(N, d , l, a, f_arr);

for fi=1:f_length   
    Z_s= Z_T_freq(1,1, fi);   
    Tx_trans_T(fi)= abs(beamformer.'*(eye(N)-S_T_freq(:,:, fi))*Z_RT_LOS(fi,:).')^2/Z_0^2./(abs(1+ Z_s/Z_0)).^2;    
end
Tx_trans_T=Tx_trans_T.';
end

function Z_RT_X = create_LOS_channel(f_arr, N, Dtxrx, d, theta_ff, G, a, R , include_multipath)

c=3*10^8;
f_length= length(f_arr);
Z_RT_LOS=zeros(f_length, N);
Z_RT_NLOS=zeros(f_length, N);
K=10^0.9;
random_multipath= 1/sqrt(K)*(1/sqrt(2)*(randn(1,N) +1j*randn(1,N)));
for fi=1:f_length
    freq_i= f_arr(fi);
    Z_s= 1./((1j*2*pi.*freq_i).*a/(c*R))+ 1./(1./((1j*2*pi.*freq_i).*a*R/c)+ 1/R);
    Z_RT_LOS(fi, :)= c./(2*pi.*freq_i*Dtxrx).*G.*real(Z_s)*exp(-1j*2*pi.*freq_i.*d*sin(theta_ff)/c*[0:1:N-1]);
    %Z_RT_NLOS(fi, :)= c./(2*pi.*freq_i*Dtxrx).*G.*real(Z_s)*(exp(-1j*2*pi.*freq_i.*d*sin(theta_ff)/c*[0:1:N-1])  + 1/sqrt(K)*(1/sqrt(2)*(randn(1,N) +1j*randn(1,N))) );
    Z_RT_NLOS(fi, :)= c./(2*pi.*freq_i*Dtxrx).*G.*real(Z_s)*(exp(-1j*2*pi.*freq_i.*d*sin(theta_ff)/c*[0:1:N-1])  + random_multipath );
end


if include_multipath==1
    Z_RT_X=Z_RT_NLOS;
else
    Z_RT_X=Z_RT_LOS;
end

end


function Z_RT_X = create_LOS_channel_N_dip(f_arr, Z_T_freq, N, Dtxrx, d, theta_ff,  l , Directivity_linear_scale, include_multipath, alpha_L, theta_L)

c=3*10^8;
f_length= length(f_arr);
Z_RT_LOS=zeros(f_length, N);
Z_RT_NLOS=zeros(f_length, N);
%K=10^0.9;
%random_multipath= 1/sqrt(K)*(1/sqrt(2)*(randn(1,N) +1j*randn(1,N)));
%k_arr=2*pi*f_arr/c;
%F=@(t, k) ((cos(k*l/2.*cos(t))-cos(k*l/2))./sin(t)).^2;

% L=8;
% alpha_L= 1/sqrt(2)*(randn(1, L)  +1j*randn(1,L));
% delta=pi/20;
% theta_L = unifrnd(theta_ff-delta, theta_ff+delta, 1, L);
for fi=1:f_length
    freq_i= f_arr(fi);
    %k=k_arr(fi);
    %Z_s= self_impedance_dip(k, l, a);
    Z_s= Z_T_freq(1,1, fi);
    %Z_s= 1./((1j*2*pi.*freq_i).*a/(c*R))+ 1./(1./((1j*2*pi.*freq_i).*a*R/c)+ 1/R);
    %G= sqrt(F(pi/2, k));
    %G=1.5;
    G= (Directivity_linear_scale(fi));
    Z_RT_LOS(fi, :)= c./(2*pi.*freq_i*Dtxrx).*G.*real(Z_s)*exp(-1j*2*pi.*freq_i.*d*sin(theta_ff)/c*[0:1:N-1]);
    %Z_RT_NLOS(fi, :)= c./(2*pi.*freq_i*Dtxrx).*G.*real(Z_s)*(exp(-1j*2*pi.*freq_i.*d*sin(theta_ff)/c*[0:1:N-1])  + 1/sqrt(K)*(1/sqrt(2)*(randn(1,N) +1j*randn(1,N))) );
    random_multipath= alpha_L*exp(-1j*2*pi.*freq_i.*d*sin(theta_L.')/c*[0:1:N-1]);
    Z_RT_NLOS(fi, :)= c./(2*pi.*freq_i*Dtxrx).*G.*real(Z_s)*(exp(-1j*2*pi.*freq_i.*d*sin(theta_ff)/c*[0:1:N-1]) + random_multipath  );
end

if include_multipath==1
    Z_RT_X=Z_RT_NLOS;
else
    Z_RT_X=Z_RT_LOS;
end

end

function f_load = create_f_load_vector(S_pos, f_arr, S_pos_length)
f_length= length(f_arr);

f_load=zeros(S_pos_length, f_length);

for fi=1:f_length
    freq_i=f_arr(fi);
    for si=1:S_pos_length
%         if abs(real(S_pos(si)))/abs(imag(S_pos(si)))<10^(-8)  
%             f_load(si, fi)= pi/2*(1/(abs(imag(S_pos(si))) -   2*pi.*freq_i)^2   + 1/(abs(imag(S_pos(si))) +   2*pi.*freq_i)^2);
%         else
            f_load(si, fi)= pi/2*real(1./(S_pos(si) - 1j*2*pi.*freq_i)   + 1./(S_pos(si) + 1j*2*pi.*freq_i));
       % end
    end
        
end

end


function  sum_mu_i_f_i = create_sum_mu_i_f_i(mu_arr, f_arr, f_load)
f_length= length(f_arr);
%f_load =  create_f_load_vector(S_pos, f_arr, S_pos_length);
sum_mu_i_f_i  = zeros(1, f_length);
for fi=1:f_length
    sum_mu_i_f_i(fi)= mu_arr*f_load(:, fi);    
end


end


function mu_arr_opt = compute_optimal_mu_arr(num_con, f_load, B_BF, Tau_LOAD, SNR_eq, f_discrete)

achievable_rate_MISO_optimized_mu_i_opt=zeros(1, num_con);
achievable_rate_MISO_optimized_mu_best=0;
mu_best= zeros(1, num_con);
for nci=1:num_con
    mu_max= zeros(1, num_con);
    mu_min= zeros(1, num_con);
    SNR_eq_max= max(SNR_eq);
    %f_load_i= @(f) compute_f_load_i(f_load, nci, f);
    f_load_i= f_load(nci, :);
    mu_min(nci)= -SNR_eq_max/log(2)/f_load_i((length(f_discrete)+1)/2);
    %C_mu_max = compute_C_mu_arr(mu_max, Tau_LOAD, f_load, num_con , f_min, f_max);
    %C_mu_min = compute_C_mu_arr(mu_min, Tau_LOAD, f_load, num_con , f_min, f_max);
    
    C_mu_max = compute_C_mu_arr(mu_max, Tau_LOAD, f_load, num_con , f_discrete);
    C_mu_min = compute_C_mu_arr(mu_min, Tau_LOAD, f_load, num_con , f_discrete);
    Max_iter_b_search=1000;
    converged_flag_mu_i_search=0;
    
    for iter=1:Max_iter_b_search
        mu_mid= (mu_min+mu_max)/2;
        %C_mu_mid = compute_C_mu_arr(mu_mid, Tau_LOAD, f_load, num_con , f_min, f_max);
        C_mu_mid = compute_C_mu_arr(mu_mid, Tau_LOAD, f_load, num_con , f_discrete);
        if C_mu_mid - B_BF < 0
            mu_min=mu_mid;
            %if (abs(C_mu_mid(nci) -B_BF(nci))/B_BF(nci) < 10^(-6)) 
            %if (abs((C_mu_mid(nci) -B_BF(nci))*mu_mid(nci)) < 10^(-2)) 
            if (abs((C_mu_mid(nci) -B_BF(nci))) < 10^(-10)) 
                mu_converged=mu_mid;
                
                %fprintf("mu 1 search has converged")
                converged_flag_mu_i_search=1;
                break;
            end
        else 
            mu_max= mu_mid;
        end          
    end
    if converged_flag_mu_i_search==1
        Transmission_coeff_opt_mu_i_opt= Tau_LOAD(mu_converged);
        [achievable_rate_MISO_optimized_mu_i_opt(nci), ~]=compute_achievable_rate_SISO(SNR_eq, Transmission_coeff_opt_mu_i_opt, f_discrete);
        if achievable_rate_MISO_optimized_mu_i_opt(nci)>achievable_rate_MISO_optimized_mu_best
            achievable_rate_MISO_optimized_mu_best=achievable_rate_MISO_optimized_mu_i_opt(nci);
            mu_best= mu_converged;
        end
    end
end
mu_arr_opt=mu_best;
end


function C_mu_arr = compute_C_mu_arr(mu_arr, Tau_LOAD, f_load, num_con , f_arr)
C_mu_arr=zeros(num_con,1);


for nci=1:num_con
    %f_load_i= @(f) compute_f_load_i(f_load, nci, f);
    f_load_i = f_load(nci, :);
    %C_mu_arr(nci)=integral(  @(f) (f_load_i(f).*log(1./(1-Tau_LOAD(f, mu_arr)))), 0,10^20, 'Waypoints', [f_min, f_max]) ;
    C_mu_arr(nci)= trapz(f_arr, f_load_i.*log(1./(1-Tau_LOAD(mu_arr))) );
end

end

function f_load_i = compute_f_load_i(f_load, nci, f_arr)

F_load= f_load(f_arr);
f_load_i= F_load(nci, :);

end

function  tau_freq_flat= compute_tau_freq_flat(B_BF, f_load, num_con ,  f_arr)
reflection_cf_i=zeros(num_con,1);
for nci=1:num_con
    f_load_i= f_load(nci, :);
    %reflection_cf_i(nci)= exp(-B_BF(nci)./(integral( @(f) f_load_i(f)  , f_min, f_max) ));
    reflection_cf_i(nci)= exp(-B_BF(nci)./(trapz( f_arr,  f_load_i  ) ));
end

reflection_cf_max= max(reflection_cf_i);
%tau_freq_flat=  (1-reflection_cf_max)*((f_arr<=f_max)&(f_arr>=f_min)); 
tau_freq_flat=  (1-reflection_cf_max)*(ones(1, length(f_arr))); 
end

function [sobj, S_T, D_linear]= scattering_parameters_linear_dipole_array(length_d, width_d, dx, freq_range, N_ant)

dipole_length=length_d;
dipole_width=width_d;
d = dipole('Length',dipole_length,'Width',dipole_width);

dipole_array = linearArray;
%dipole_array=phased.ULA(N_ant, dx, 'Element', d);
dipole_array.Element = d;
dipole_array.NumElements = N_ant;
dipole_array.ElementSpacing = dx;

%D=directivity(dipole_array, freq_range, 0);
D_dbi=pattern(d, freq_range, 0, pi/2);
D_linear= 10.^(D_dbi/10);
sobj = sparameters(dipole_array,freq_range);

S_T= sobj.Parameters;

end

function [PSD_water_fill, PSD_equal]= water_fill_allocation(Tx_trans_T, Ps, Noise, f_discrete)

%PSD_water_fill=zeros(1, length(f_discrete));
SNR= Tx_trans_T/Noise;
eta_max= max(SNR);
eta_min= 0 ;
Max_iter=500;
f_len=length(f_discrete);

for iter=1:Max_iter
    eta_mid= (eta_max+eta_min)/2;
    constraint= trapz(f_discrete, ( max(0, 1/eta_mid- 1./SNR)));
    %constraint= (f_discrete(2)-f_discrete(1))*sum(( max(0, 1/eta_mid- 1./SNR)));
    if constraint > Ps
        eta_min= eta_mid;
    else
        eta_max=eta_mid;
        
    end
    if abs(constraint-Ps)< 10^(-15)
            PSD_water_fill=max(0, 1/eta_mid- 1./SNR)/Noise ;
            PSD_equal= (ones(1, f_len)*constraint/(f_discrete(f_len)- f_discrete(1)))/Noise;
            %PSD_equal= trapz(f_discrete,PSD_water_fill )/length(f_discrete)*ones(1, length(f_discrete));
            break;
    end
end
end


function beamformer= create_beamformer_optimal_discrete_set(N, f_c, f_c_index, S_T_fit_freq, Z_T_fit_freq, Dtxrx, d, l, theta_ff,  Z_0, Directivity_linear_scale, include_multipath, alpha_L, theta_L )

beamformer_discrete_set= ((2*(dec2bin(2^N-1:-1:0)-'0')-1).')/sqrt(N);
Tx_trans_T_arr=zeros(2^N,1);
for i=1:2^N
    beamformer_i=beamformer_discrete_set(:, i);
    S_T_freq= S_T_fit_freq(:, :, f_c_index);
    Z_T_freq= Z_T_fit_freq(:, :, f_c_index);
    
    Tx_trans_T_arr(i) = create_Tx_trans_T_N_dip(f_c, S_T_freq, Z_T_freq,  N, Dtxrx, d, l, theta_ff,  beamformer_i, Z_0, Directivity_linear_scale(f_c_index), include_multipath, alpha_L, theta_L);
    
    
end

[Max_Tx_trans_T_arr, Tx_trans_T_arr_i]= max(Tx_trans_T_arr);
beamformer=beamformer_discrete_set(:,Tx_trans_T_arr_i);

end