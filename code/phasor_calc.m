function [M, Phi, G, S] = phasor_calc(f, time, decay, IRF)
w = 2*pi*f;
G_decay = sum(decay.*cos(w*time))/sum(decay(:));
S_decay = sum(decay.*sin(w*time))/sum(decay(:));

G_IRF = sum(IRF.*cos(w*time))/sum(IRF(:));
S_IRF = sum(IRF.*sin(w*time))/sum(IRF(:));

GS = [G_IRF, S_IRF; -S_IRF, G_IRF]*[G_decay; S_decay]/(G_IRF^2 + S_IRF^2);
% GS = [0.7312, 0.6381; -0.6381, 0.7312]*[G_decay; S_decay];                  % Coumarin6 correction


G = GS(1);
S = GS(2);
Phi = pi-atan2(S,-G);
M = sqrt(G.^2 + S.^2);
end