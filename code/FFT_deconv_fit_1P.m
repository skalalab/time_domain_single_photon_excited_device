function [Tm, A0n, T1, A1n, T2, A2n, IRF, decay_orig_offset_rmvd, Chi_Square_red_opt, Chi_Square_cyan_opt, f_cut_ind_opt, decay_orig_norm, decay_deconv_final_opt, decay_reconv_opt, fit_opt_XY, decay_reconv2_opt] = FFT_deconv_fit_1P(time, decay_orig, SHG, Shift, fig)

L = size(decay_orig,1);
fs = 1/(1e-9*time(2));  %time in nanoseconds converted to seconds - this gives the sampling frequency

decay_orig_fft = fft(decay_orig, L+1);  % (L+1)-point fft of the decay

offset_start_frac = 0.968;               %%% calculate offset (i.e., ambient light signal) from the tail end of the decay
offset_end_frac = 1.00;                 %%% 

offset = mean(decay_orig(round(offset_start_frac*L):round(offset_end_frac*L)));    % the background level is calculated from the end of the decay histogram

decay = decay_orig - offset;

decay_orig_offset_rmvd = decay;
decay_orig_norm = decay_orig/sum(decay_orig);  % normalize the decay by its sum
offset_frac = offset/sum(decay);               % offset as a fraction of the decay sum (to be added to the reconvolved decay later)
decay = decay/sum(decay);                      % normalize the offset-removed decay by its sum
[~, peak_ind] = max(decay);

%%

IRF = SHG/max(SHG(:));
IRFinterp = interp1((1:1:size(IRF,1))', IRF, (1:0.1:size(IRF,1))');  % 10 fold interpolation of the IRF decay histogram
IRFinterp = circshift(IRFinterp, round(10*Shift));      % apply desired shift
IRF = interp1((1:1:size(IRFinterp,1))', IRFinterp, (1:10:size(IRFinterp,1))');  % 1/10 fold - bring back to original time resolution

%%
warning('off');
IRF_fft = fft(IRF, L+1);
decay_fft = fft(decay, L+1);

decay_deconv_fft = decay_fft./IRF_fft;

%   find the cut-off frequency by comparing fft values with their expected
%   standard deviation per freq bin (freq where they meet has cutt-off SNR=1)
%   Based on the following paper:
%   Andre, J. C., et al. "Applications of fast Fourier transform to deconvolution in single photon counting." Journal of Physical Chemistry 83.17 (1979): 2285-2294.

FR = real(decay_fft);           FI = imag(decay_fft);
GR = real(IRF_fft);             GI = imag(IRF_fft);
DR = real(decay_deconv_fft);    DI = imag(decay_deconv_fft);

var_FR(1:L/2+1, 1) = 0.5*(FR(1) + FR(1:2:L+1));
var_FI(1:L/2+1, 1) = 0.5*(FR(1) - FR(1:2:L+1));

var_GR(1:L/2+1, 1) = 0.5*(GR(1) + GR(1:2:L+1));
var_GI(1:L/2+1, 1) = 0.5*(GR(1) - GR(1:2:L+1));

var_DR(1:L/2+1, 1) =    ((GR(1:L/2+1)./(GR(1:L/2+1).^2+GI(1:L/2+1).^2)).*sqrt(var_FR)).^2 + ...
                        ((GI(1:L/2+1)./(GR(1:L/2+1).^2+GI(1:L/2+1).^2)).*sqrt(var_FI)).^2 + ...
                        ((FR(1:L/2+1)./(GR(1:L/2+1).^2+GI(1:L/2+1).^2)-2*DR(1:L/2+1).*GR(1:L/2+1)./(GR(1:L/2+1).^2+GI(1:L/2+1).^2)).*sqrt(var_GR)).^2 +...
                        ((FI(1:L/2+1)./(GR(1:L/2+1).^2+GI(1:L/2+1).^2)+2*DR(1:L/2+1).*GI(1:L/2+1)./(GR(1:L/2+1).^2+GI(1:L/2+1).^2)).*sqrt(var_GI)).^2;

var_DI(1:L/2+1, 1) =    ((GI(1:L/2+1)./(GR(1:L/2+1).^2+GI(1:L/2+1).^2)).*sqrt(var_FR)).^2 + ...
                        ((GR(1:L/2+1)./(GR(1:L/2+1).^2+GI(1:L/2+1).^2)).*sqrt(var_FI)).^2 + ...
                        ((FI(1:L/2+1)./(GR(1:L/2+1).^2+GI(1:L/2+1).^2)-2*DI(1:L/2+1).*GR(1:L/2+1)./(GR(1:L/2+1).^2+GI(1:L/2+1).^2)).*sqrt(var_GR)).^2 +...
                        ((FR(1:L/2+1)./(GR(1:L/2+1).^2+GI(1:L/2+1).^2)-2*DI(1:L/2+1).*GI(1:L/2+1)./(GR(1:L/2+1).^2+GI(1:L/2+1).^2)).*sqrt(var_GI)).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_cut_ind1 = find(abs(sqrt(var_DR)./DR(1:L/2+1))>5, 1, 'first');   % default >5 
f_cut_ind2 = find(abs(sqrt(var_DI)./DI(1:L/2+1))>5, 1, 'first');   % default >5 
f_cut_ind3 = 10;                                                   % manually-set maximum allowed number of freq. components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('on');
f_cut_ind_limit = min([f_cut_ind1, f_cut_ind2, f_cut_ind3]);

Chi_Square_cyan_opt = Inf;

%   FOR loop: cycle through several frequency cut-off points and compare Chi^2 to pick the best cut-off frequency
for f_cut_ind = 2:f_cut_ind_limit
    fs = 1/(1e-9*time(2));  %time in nanoseconds
    f = fs*(0:L)'/L;
    
    A_Trail = -2*pi*f(f_cut_ind)*(real(decay_deconv_fft(f_cut_ind)))/(imag(decay_deconv_fft(f_cut_ind)));
    B_Trail = ((real(decay_deconv_fft(f_cut_ind)))/A_Trail)*(A_Trail^2+4*pi^2*f(f_cut_ind)^2);
    warning('off');
    decay_deconv_fft_Exp_Trail(1:f_cut_ind-1, 1) = decay_deconv_fft(1:f_cut_ind-1, 1);
    decay_deconv_fft_Exp_Trail(f_cut_ind:(L/2+1), 1) = B_Trail./(A_Trail+2*pi*1i*f(f_cut_ind:(L/2+1)));
    decay_deconv_fft_Exp_Trail((L/2+2):(L+1), 1) = flip(conj(decay_deconv_fft_Exp_Trail(2:(L/2+1))));
    warning('on');
    decay_deconv_Exp_Trail = ifft(decay_deconv_fft_Exp_Trail, 'symmetric');         % inverse Fourier transform
    decay_deconv_Exp_Trail = decay_deconv_Exp_Trail/sum(decay_deconv_Exp_Trail);
    
    decay_deconv_final = decay_deconv_Exp_Trail(1:size(time,1));

%%  Reconvolve the extracted decay with the IRF and calculate Chi_sqaure for goodness of agreement with measured decay
   
    decay_reconv = conv(decay_deconv_Exp_Trail,IRF); decay_reconv = decay_reconv + offset_frac*sum(decay_reconv(1:L));
    decay_reconv = decay_reconv/sum(decay_reconv(1:L));
    
    
    Chi_Square_red  = (max(decay_orig(:))/max(decay_orig_norm(:)))*...
                      sum(((decay_reconv(1:L)-decay_orig_norm(1:L)).^2)...
                          ./(decay_orig_norm(1:L)+eps))/L;
   
   
%%  Exponential fitting of deconvolved decay trace and reconvolution with the IRF
    warning('off');
    
    b0 = [0.02;-2;0.02;-0.4];    % starting parameter set                                                         %%%
    lb = [0; -10; 0; -10];                                                                                        %%%
    ub = [1000; -0.1; 1000; -0.1];                                                                                %%%
    
%%% Bi-exponential decay fit
    f2 = fit(time(2:L-peak_ind)-time(1),decay_deconv_final(2:L-peak_ind),'exp2', 'StartPoint', b0, 'Lower', lb, 'Upper', ub);   %%%
% %     f2 = fit(time(2:L-peak_ind)-time(1),decay_deconv_final(2:L-peak_ind),'exp1');   %%% single-exponential decay
    b = [offset_frac*sum(f2(time(2:L-peak_ind))); f2.a; -f2.b; f2.c; -f2.d];        
% %     b = [offset_frac*sum(f2(time(2:L-peak_ind))); f2.a; -f2.b];                     %%% single-exponential decay

    decay_reconv2 = conv(cat(1,decay_deconv_final(1),f2(time(2:end))),IRF); decay_reconv2 = decay_reconv2 + offset_frac*sum(decay_reconv2(1:L)); 
%
    decay_reconv2_norm = decay_reconv2(1:L)/sum(decay_reconv2(1:L));
    
                      
    Chi_Square_cyan = (max(decay_orig(:))/max(decay_orig_norm(:)))*...
                      sum(((decay_reconv2_norm(1:L)-decay_orig_norm(1:L)).^2)...
                          ./(decay_orig_norm(1:L)+eps))/L;

    if Chi_Square_cyan < Chi_Square_cyan_opt
        Chi_Square_cyan_opt = Chi_Square_cyan;
        Chi_Square_red_opt = Chi_Square_red;
        decay_deconv_final_opt = decay_deconv_final;
        decay_reconv_opt = decay_reconv;
        b_opt = b;
        decay_reconv2_opt = decay_reconv2_norm;
        f_cut_ind_opt = f_cut_ind;
        f_opt = f2;
    end
end
%%
%%%
fit_opt_XY = [time(2:L-peak_ind), f_opt(time(2:L-peak_ind))];

%%% Single-exponential decay fit
%{
T1 = 1/b_opt(3); A1n = b_opt(2)/(b_opt(1)+b_opt(2));
T2 = 0; A2n = 0;

A0n = b_opt(1)/(b_opt(1)+b_opt(2));
Tm = T1;
%}
%%% Bi-exponential decay fit
%
if b_opt(3) >= b_opt(5)
    T1 = 1/b_opt(3); A1n = b_opt(2)/(b_opt(1)+b_opt(2)+b_opt(4));
    T2 = 1/b_opt(5); A2n = b_opt(4)/(b_opt(1)+b_opt(2)+b_opt(4));
    
else
    T1 = 1/b_opt(5); A1n = b_opt(4)/(b_opt(1)+b_opt(2)+b_opt(4));
    T2 = 1/b_opt(3); A2n = b_opt(2)/(b_opt(1)+b_opt(2)+b_opt(4));
end
A0n = b_opt(1)/(b_opt(1)+b_opt(2)+b_opt(4));
Tm = A1n*T1+A2n*T2;
%

warning('on');


end
