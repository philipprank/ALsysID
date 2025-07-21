%%% Frequency Response Visualization %%%

% 1) Systemkoeffizienten
b = [0, -0.05];
a = [1, -0.95];
% b = [0,  0.1, -0.2];    % Zähler: 0·z^0 +0.1·z⁻¹ −0.2·z⁻²
% a = [1, -1.5,  0.9];    % Nenner: 1 −1.5·z⁻¹ +0.9·z⁻²

% Abtastzeit / -frequenz
Ts   = 1;
fs   = 1/Ts;

% Frequenzraster und Frequenzgang (einseitig 0…π)
nfft = 512;
[H, omega] = freqz(b, a, nfft, 'half');

% Bereich für die gemeinsame X-Achse
omega_lim = [0, pi];

%% 2) Ensemble-mittlere PSD des Eingangssignals
% Workspace muss enthalten: u_opt_loc (1400×1×20)

% Monte-Carlo-Parameter
M      = size(u_opt_loc,3);    % z.B. 20
L      = size(u_opt_loc,1);    % z.B. 1400
k      = 100;                  % Periodenlänge
window = hann(k);              % Hann-Window über eine Periode
nover  = floor(0.5 * k);       % 50 % Überlapp
Pxx    = zeros(nfft/2+1, M);   % Speicher (257×M)

% PSD per Lauf mit pwelch und nfft=512
for i = 1:M
    u_k = squeeze(u_opt_loc(:,:,i));  % 1400×1
    u_k = u_k - mean(u_k);            % DC entfernen
    [Pxx(:,i), f] = pwelch(u_k, window, nover, nfft, fs);
end
Pxx_avg = mean(Pxx, 2);
omega_psd = 2*pi * f / fs;  % = [0 … π] rad/sample

%% 3) Gemeinsamer Plot mit zwei Y-Achsen
figure; clf;

% linke Achse: Amplitudengang
yyaxis left
plot(omega, abs(H), 'b-', 'LineWidth', 1.5);
ylabel('|H(e^{j\omega})|', 'Color', 'b');
set(gca, 'YColor', 'b');
xlim(omega_lim);
grid on

% rechte Achse: PSD
yyaxis right
plot(omega_psd, 10*log10(Pxx_avg), 'r-', 'LineWidth', 1.5);
ylabel('PSD [dB/Hz]', 'Color', 'r');
set(gca, 'YColor', 'r');
xlim(omega_lim);

% Gemeinsame Beschriftungen
xlabel('\omega [rad/sample]');
legend({'|H(e^{j\omega})|','PSD von u\_opt\_loc'}, ...
       'Location','northeast');
title('System-Amplitude und Eingang-PSD gegen \omega');


%% 2) System 4: MIMO ARX(2,2)
% Zustandsraum (Ts=1):
%   x = [y1(t-1); y2(t-1); y1(t-2); y2(t-2); u1(t-1); u2(t-1); u1(t-2); u2(t-2)]
A4 = [
  0.7007  -0.7007   0.1363   0.4396   0.1100  -0.1200   0.1300  -0.1400;
  0.7007   0.7007   0.2413  -0.3996   0.1500  -0.1600   0.1700  -0.1800;
  1.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000;
  0.0000   1.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000;
  0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000;
  0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000;
  0.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000;
  0.0000   0.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000
];

B4 = [
    zeros(4,2);
    eye(2);
    zeros(2,2)
    ];

% Ausgabe-Matrix (beide Outputs aus den ersten 2 Zuständen)
C4 = [eye(2), zeros(2,6)];
D4 = zeros(2,2);

% Zustandsraum-System anlegen
sys4 = ss(A4, B4, C4, D4, 1);

% Frequenzraster
omega4 = linspace(0, pi, 512);

% Frequenzantwort (2×2×512)
H4 = freqresp(sys4, omega4);

% Für jeden ω den größten Singulärwert bestimmen
sv4 = zeros(size(omega4));
for k = 1:length(omega4)
    Hk    = squeeze(H4(:,:,k));
    sv4(k) = max(svd(Hk));
end

% Plot des maximalen Singulärwerts
figure;
plot(omega4, sv4, 'LineWidth', 1.5)
grid on
xlabel('\omega [rad/sample]')
ylabel('max \sigma\{H_{sys4}(e^{j\omega})\}')
title('Maximaler Singulärwert-Frequenzgang von System 4 (MIMO ARX(2,2))')

%% ------------------------------------------------------------------------
%% Vollständiges Skript inkl. Durchschnitts-PSD aller optimierten Inputs
%% ------------------------------------------------------------------------

A = [
  0.7007  -0.7007   0.1363   0.4396   0.1100  -0.1200   0.1300  -0.1400;
  0.7007   0.7007   0.2413  -0.3996   0.1500  -0.1600   0.1700  -0.1800;
  1.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000;
  0.0000   1.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000;
  zeros(4,8);
  0.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000;
  0.0000   0.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000
];
B = [ zeros(4,2); eye(2); zeros(2,2) ];
C = [ eye(2), zeros(2,6) ];
D = zeros(2,2);
sys4 = ss(A,B,C,D,1);   % Ts = 1 s

%% 2) Lade dein adaptiv optimiertes Signal --------------------------------
% Nach deiner Monte-Carlo-Schleife musst Du im Workspace haben:
%   u_adapt  % Dimension T×2×M

if ~exist('u_opt_loc','var')
    error('Variable u_opt_loc nicht gefunden. Bitte Simulation zuvor ausführen.');
end
[T, du, M] = size(u_opt_loc);
if du~=2
    error('Erwarte genau 2 Eingangskanäle in u_opt_loc (du=2).');
end

%% 3) Mittlere PSD aller Kanäle über alle M Runs --------------------------
Fs   = 1;       % Abtastrate [Hz]
Nf   = 1024;    % FFT-Punkte
Psum = zeros(Nf/2+1,1);

for m = 1:M
    for ch = 1:du
        u_ch = squeeze(u_opt_loc(:,ch,m));
        Pxx  = pwelch(u_ch, [], [], Nf, Fs, 'onesided');
        Psum = Psum + Pxx;
    end
end

% Mittelwert bilden und normalisieren
Pavg   = Psum / (du*M);
Pavg_n = Pavg / max(Pavg);

% Frequenzvektor für PSD
f_psd = linspace(0, Fs/2, Nf/2+1)';

%% 4) Singulärwert-Frequenzgang des MIMO-Systems ---------------------------
omega = linspace(0, pi, 512);
H4    = freqresp(sys4, omega);
sv    = zeros(size(omega));

for k = 1:numel(omega)
    Hk    = squeeze(H4(:,:,k));
    sv(k) = max(svd(Hk));
end

% Umrechnung in Hz: f = omega/pi * (Fs/2)
f_sv = omega/pi * (Fs/2);

%% 5) Overlay-Plot auf gemeinsamer f-Achse in Hz --------------------------
figure;
yyaxis left
plot(f_sv, sv, 'LineWidth',1.5)
ylabel('max\,\sigma\{H_{sys4}(e^{j\omega})\}')
ylim([0 max(sv)*1.1])

yyaxis right
plot(f_psd, Pavg_n, '--', 'LineWidth',1.5)
ylabel('Normierte mittlere PSD\{u\}')
ylim([0 1.1])

xlabel('f [Hz]')
title('MIMO-System (max. Singulärwert) & mittlere PSD aller adaptiven Inputs')
legend('max \sigma\{H\}','PSD_{avg}\{u\}','Location','best')
grid on
xlim([0 Fs/2])

