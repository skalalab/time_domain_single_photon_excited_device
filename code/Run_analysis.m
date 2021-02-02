close all;
clearvars;

filedir = 'C:\Users\ksamimi\Desktop\temp\delete\2020_11_01_Tcells_repeat2\Dish1\Dish1_decays\';
% filedir = 'C:\Users\ksamimi\Desktop\temp\delete\2020_11_01_Tcells_repeat2\Dish2\Dish2_decays\';

BckgrndDecay = zeros(200,1);
% BckgrndDecay = 0.3*BckgrndDecay_Dish1;
% BckgrndDecay = 0.4*BckgrndDecay_Dish2;


all_decays = [];
bad_decays = [];


c = [255,   0,      0;
     255,   165,    0;
     0,     255,    0;
     0,     0,      255;
     255,   0,    255]./255;            % colors

DecayFiles = dir(strcat(filedir, 'd_cfd___TCell_*.csv'));
% DecayFiles = dir(strcat(filedir, 'd_cfd___RBC_*.csv'));
% DecayFiles = dir(strcat(filedir, 'd_cfd___Bckgrnd_*.csv'));




Shift = 0.0;  % in steps of 0.1 time bin


time = (0:(20/199):20)';        % time vector (0 to 20) - in nanoseconds [ns]
SHG = circshift(SHG1, 0);  % IRF measured from second harmonic generation or a fast decaying standard sample.


%%
phasor_freq = 0.050;                    % in Giga Hertz [GHz]

scrsz = get(0,'ScreenSize');
fig = figure('Position',[round(1.1*scrsz(3)) round(0.1*scrsz(4)) round(0.8*scrsz(3)) round(0.8*scrsz(4))]);
% NADH_fits={};
filecount=0;
for file = DecayFiles'
    filecount = filecount+1;
    data = importdata(strcat(filedir, file.name));
    decay1 = data(:,1) - BckgrndDecay;     % Subtract dish's background fluorescence


    all_decays = cat(2, all_decays, decay1);
    decay = circshift(decay1, 0);
    
    
    L = size(decay,1);
    if sum(decay)==0, decay=ones(size(decay)); end
    [Tm, A0n, T1, A1n, T2, A2n, IRF, decay_orig_offset_rmvd, Chi2_red, Chi2_cyan, f_cutoff_ind,...
        decay_orig_norm, decay_deconv_final_opt, decay_reconv_opt, fit_opt_XY, decay_reconv2_opt] = FFT_deconv_fit_1P(time, decay, SHG, Shift, fig);


    [M, Phi, G, S] = phasor_calc(phasor_freq, time, decay_orig_offset_rmvd, SHG);       % Calculate phasor coordinates and correct using the IRF
    


    figure(fig), plot(time, decay_orig_norm, 'LineWidth', 2, 'Color', 'b'); hold on; ax=gca;ax.FontSize=25;ax.LineWidth=2; ylim([-0.005 max(IRF(:))/30*sum(IRF)]); ylabel('Norm. Count', 'Interpreter', 'latex', 'FontSize', 25); xlabel('time [ns]', 'Interpreter', 'latex', 'FontSize', 25);
    figure(fig), plot(time, IRF*(max(decay_orig_norm)/max(IRF)), 'LineWidth', 2, 'Color', 'g');         % plot IRF
    figure(fig), plot(time, decay_deconv_final_opt, 'LineWidth', 2, 'Color', 'm');                      % plot IRF-deconvolved decay
    figure(fig), plot(time, decay_reconv_opt(1:L), 'LineWidth', 2, 'Color', 'r');                       % reconvolution with IRF should reproduce the measured decay
    figure(fig), plot (fit_opt_XY(:,1), fit_opt_XY(:,2), 'Color', 'k');                                 % plot bi-exponential fit of the deconvolved decay
    figure(fig), plot (time, decay_reconv2_opt, 'LineWidth', 1.5, 'Color', 'c');                        % reconvolution of the bi-exponential fit with IRF
        
    disp(char(strcat('Chi2_red ='," ",num2str(Chi2_red,'%05.4f'), "  ", 'Chi2_cyn ='," ",num2str(Chi2_cyan,'%05.4f'), "  ",...
        'A1% ='," ",num2str(100*A1n/(A1n+A2n),'%04.2f'), "  ", 'Tm ='," ",num2str(1000*Tm,'%04.0f'),'[ps]', "   ",...
        'T1 ='," ",num2str(1000*T1,'%04.0f'),'[ps]', "  ",  'T2 ='," ",num2str(1000*T2,'%04.0f'),'[ps]', "  ", 'f_cutoff_ind =',...
        " ",num2str(f_cutoff_ind), " ", 'Shift ='," ",num2str(Shift,'%5.1f'), " ", 'num_iter ='," ",num2str(num_iter), " ", 'count=',num2str(sum(decay)))));

    pause(0.01)
        
    NADH_fits{filecount,1} = cat(2, [1000*Tm, A0n, 1000*T1, (A1n/(A1n+A2n)), 1000*T2, (A2n/(A1n+A2n))], [M, Phi, G, S, Shift, sum(decay), Chi2_cyan]);
    hold off;
    
    if (Tm>10)||(A2n<0.05)||(abs(A0n)>0.02)||(Chi2_cyan>15)
        warning(char(strcat('Check fit for cell in', " ", file.name, "  ",'A0n=', num2str(A0n), "  ",'Tm=', num2str(Tm),'[ns]', " ", 'count=',num2str(sum(decay)))));
        bad_decays = cat(2, bad_decays, decay);
    end

end
NADH_fits(cellfun(@isempty, NADH_fits)) = {NaN(size(NADH_fits{1,1}))};

for i = 1:size(NADH_fits,1)
    NADH_fits_array(i,:) = NADH_fits{i};
end

%%%%%

NADH_fits_array_Dish1 = NADH_fits_array;
% NADH_fits_array_Dish2 = NADH_fits_array;


%%
%%%%%
NADH_Activated = NADH_fits_array_Dish1;
NADH_Naive = NADH_fits_array_Dish2;
%%%%%

%% Remove outliers based on photon count, phasor, or lifetime parameters
ind = (round(NADH_Activated(:,3),2)==100)|(round(NADH_Activated(:,3),2)>=600)|(NADH_Activated(:,12)<=50000);
mask = ones(size(NADH_Activated,1), 1);
mask(ind) = NaN;
mask = repmat(mask,[1 size(NADH_Activated,2)]);
NADH_Activated = NADH_Activated.*mask;

ind = (round(NADH_Naive(:,3),2)==100)|(round(NADH_Naive(:,3),2)>=600)|(NADH_Naive(:,12)<=30000);
mask = ones(size(NADH_Naive,1), 1);
mask(ind) = NaN;
mask = repmat(mask,[1 size(NADH_Naive,2)]);
NADH_Naive = NADH_Naive.*mask;


%%    Phasor Plot
XLim=[0 1];
YLim=[0 0.55];
u = (0:0.01:100)';
w = 2*pi*phasor_freq*[0.5; 1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
scrsz = get(0,'ScreenSize');
f1 = figure('Position',[round(1.1*scrsz(3)) round(0.1*scrsz(4)) round(0.8*scrsz(3)) round(0.8*scrsz(4))]);
hold on; ax=gca;ax.FontSize=25;ax.LineWidth=2; ax.XAxisLocation='origin'; ax.YAxisLocation='origin'; axis equal; ylabel('$S$', 'Interpreter', 'latex', 'FontSize', 35); xlabel('$G$', 'Interpreter', 'latex', 'FontSize', 35); ax.XLim=XLim; ax.YLim=YLim;

scatter(NADH_Activated(:,9),NADH_Activated(:,10),10, 'r','filled'); hold on;  % phasor at w

scatter(NADH_Naive(:,9),NADH_Naive(:,10),10, 'b','filled'); hold on;  % phasor at w

alpha(0.5);

plot(1./(1+u.^2),u./(1+u.^2),'color','k');
scatter(1./(1+w.^2),w./(1+w.^2),'k', 'filled');

txt = sprintf('0.5'); text((0.973*ax.XLim(2)),(0.3*ax.YLim(2)+0.02*mean(YLim)), txt,'HorizontalAlignment','left', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 18);
txt = sprintf('1 ns'); text((0.91*ax.XLim(2)),(0.5*ax.YLim(2)+0.09*mean(YLim)), txt,'HorizontalAlignment','left', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 18);
txt = sprintf('3 ns'); text((0.52*ax.XLim(2)),(0.5*ax.YLim(2)+0.88*mean(YLim)), txt,'HorizontalAlignment','left', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 18);
txt = sprintf('10 ns'); text((0.04*ax.XLim(2)),(0.5*ax.YLim(2)+0.09*mean(YLim)), txt,'HorizontalAlignment','left', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 18);

%% Plot Tm distributions
fig = figure('Position',[round(1.1*scrsz(3)) round(0.1*scrsz(4)) round(0.8*scrsz(3)) round(0.8*scrsz(4))]); hold on;
histfit(reshape(NADH_Naive(:,1), [], 1),20,'kernel'); histfit(reshape(NADH_Activated(:,1), [], 1),20,'kernel'); alpha(0.6); ax=gca;ax.FontSize=25;ax.LineWidth=2;
legend('Quiescent mean lifetime histogram', 'Quiescent Tm distribution fit', 'Activated mean lifetime histogram', 'Activated Tm distribution fit', 'Interpreter', 'tex', 'Location', 'northeast')
ylabel('Number of cells', 'Interpreter', 'latex', 'FontSize', 25); xlabel('Tm [ps]', 'Interpreter', 'latex', 'FontSize', 25);

%% Plot T1 distributions
fig = figure('Position',[round(1.1*scrsz(3)) round(0.1*scrsz(4)) round(0.8*scrsz(3)) round(0.8*scrsz(4))]); hold on;
histfit(reshape(NADH_Naive(:,3), [], 1),20,'kernel'); histfit(reshape(NADH_Activated(:,3), [], 1),20,'kernel'); alpha(0.6); ax=gca;ax.FontSize=25;ax.LineWidth=2;
legend('Quiescent short component histogram', 'Quiescent T1 distribution fit', 'Activated short component histogram', 'Activated T1 distribution fit', 'Interpreter', 'tex', 'Location', 'northeast')
ylabel('Number of cells', 'Interpreter', 'latex', 'FontSize', 25); xlabel('T1 [ps]', 'Interpreter', 'latex', 'FontSize', 25);

%% Plot T2 distributions
fig = figure('Position',[round(1.1*scrsz(3)) round(0.1*scrsz(4)) round(0.8*scrsz(3)) round(0.8*scrsz(4))]); hold on;
histfit(reshape(NADH_Naive(:,5), [], 1),20,'kernel'); histfit(reshape(NADH_Activated(:,5), [], 1),20,'kernel'); alpha(0.6); ax=gca;ax.FontSize=25;ax.LineWidth=2;
legend('Quiescent long component histogram', 'Quiescent T2 distribution fit', 'Activated long component histogram', 'Activated T2 distribution fit', 'Interpreter', 'tex', 'Location', 'northeast')
ylabel('Number of cells', 'Interpreter', 'latex', 'FontSize', 25); xlabel('T2 [ps]', 'Interpreter', 'latex', 'FontSize', 25);

legend('Activated', 'Quiescent', 'Interpreter', 'tex', 'Location', 'northeast')

%% Plot Alpha1% distributions
fig = figure('Position',[round(1.1*scrsz(3)) round(0.1*scrsz(4)) round(0.8*scrsz(3)) round(0.8*scrsz(4))]); hold on;
XLim=[60 100];
histfit(100*reshape(NADH_Naive(:,4), [], 1),20,'kernel'); histfit(100*reshape(NADH_Activated(:,4), [], 1),20,'kernel'); alpha(0.6);
legend('Quiescent short fraction histogram', 'Quiescent \alpha_1% distribution fit', 'Activated short fraction histogram', 'Activated \alpha_1% distribution fit', 'Interpreter', 'tex', 'Location', 'northwest')
ax=gca;ax.FontSize=25;ax.LineWidth=2; ax.XLim=XLim; ax.YLim=[0 27]; 
ylabel('Number of cells', 'Interpreter', 'latex', 'FontSize', 25); xlabel('\alpha_1%', 'Interpreter', 'tex', 'FontSize', 25);

