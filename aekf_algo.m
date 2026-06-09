ocv_soc = load("ocv_soc_haritasi.mat");
hppc_verisi = load("thevenin_parametreleri.mat");
drive_cycle = load("10degC\Drive Cycles\03-27-17_09.06 10degC_US06_Pan18650PF.mat");

% hppc'den elde ettiğimiz veriler
hppc_akim = double(hppc_verisi.I);
hppc_akim = hppc_akim(:);
hppc_R0 = double(hppc_verisi.R0);
hppc_R0 = hppc_R0(:);
hppc_R1 = double(hppc_verisi.R1);
hppc_R1 = hppc_R1(:);
hppc_C1 = double(hppc_verisi.C1);
hppc_C1 = hppc_C1(:);

%c20 discharge verisinden elde ettiğimiz veriler
soc_haritasi = ocv_soc.soc_haritasi; 
ocv_haritasi = ocv_soc.ocv_haritasi;

toplam_kapasite = 3.00;
hppc_ah = double(hppc_verisi.Ah);
hppc_soc = (toplam_kapasite - abs(hppc_ah)) / toplam_kapasite;
hppc_soc = hppc_soc(:);

% thevenin parametrelerini akım ve soc değerine göre düzenleme
harita_R0 = scatteredInterpolant(hppc_akim, hppc_soc, hppc_R0, 'linear', 'nearest');
harita_R1 = scatteredInterpolant(hppc_akim, hppc_soc, hppc_R1, 'linear', 'nearest');
harita_C1 = scatteredInterpolant(hppc_akim, hppc_soc, hppc_C1, 'linear', 'nearest');

% drive cycle'dan simülasyon için aldığımız veriler
gercek_akim_verisi = drive_cycle.meas.Current;
gercek_voltaj_verisi = drive_cycle.meas.Voltage;

N = length(gercek_akim_verisi); 
soc_tahminleri = zeros(N, 1);
voltaj_tahminleri = zeros(N, 1);

% aekf ilk durum için belirlediğimiz değerler
dt = 1;
x_k = [1.0; 0.0];
P_k = [0.1, 0; 0, 0.1];
Q = [0.00001, 0; 0, 0.001];
R = 0.05;
alfa = 0.98;

for k = 1:N
    I_anlik = gercek_akim_verisi(k);
    soc_eski = x_k(1);
    
    % haritadan anlık veriyi çekiyoruz
    R0_anlik = harita_R0(I_anlik, soc_eski);
    R1_anlik = harita_R1(I_anlik, soc_eski);
    C1_anlik = harita_C1(I_anlik, soc_eski);
    tau_anlik = R1_anlik * C1_anlik;
    
    A = [1, 0; 0, exp(-dt / tau_anlik)];
    B = [dt / 10800; R1_anlik * (1 - exp(-dt / tau_anlik))];
    
    %durum uzay ifadesi
    x_k_tahmin = A * x_k + B * I_anlik;  
    P_k_tahmin = A * P_k * A' + Q;       

    soc_tahmin = max( 0.002, min( 0.998, x_k_tahmin(1) ) );
    V1_tahmin = x_k_tahmin(2);
    
    ocv_tahmin = interp1(soc_haritasi, ocv_haritasi, soc_tahmin, 'pchip', 'extrap');
    
    V_tahmin = ocv_tahmin - V1_tahmin + I_anlik * R0_anlik;
    
    delta_soc = 0.001;
    ocv_arti = interp1(soc_haritasi, ocv_haritasi, soc_tahmin + delta_soc, 'pchip', 'extrap');
    ocv_eksi = interp1(soc_haritasi, ocv_haritasi, soc_tahmin - delta_soc, 'pchip', 'extrap');
    turev_ocv = (ocv_arti - ocv_eksi) / (2 * delta_soc);
    
    H = [turev_ocv, -1];
    
    V_gercek = gercek_voltaj_verisi(k);
    y_k = V_gercek - V_tahmin;
   
    R_gecici = alfa * R + (1 - alfa) * (y_k^2 - H * P_k_tahmin * H');
    R = max(0.015, min(0.2, R_gecici)); 
    
    S = H * P_k_tahmin * H' + R; 
    K = P_k_tahmin * H' / S; % kalman kazanci hesabi
    
    x_k = x_k_tahmin + K * y_k;
    
    x_k(1) = max(0.0, min(1.05, x_k(1))); 
    
    P_k = (eye(2) - K * H) * P_k_tahmin;
    P_k = (P_k + P_k') / 2; 
    
    soc_tahminleri(k) = x_k(1);
    voltaj_tahminleri(k) = V_tahmin;
end

zaman_ekseni = (0:N-1) * dt;

figure('Name', 'SOC Tahmini');
plot(zaman_ekseni, soc_tahminleri * 100, 'LineWidth', 2);
title('AEKF ile Gercek Zamanli Doluluk Orani Tahmini');
xlabel('Zaman saniye');
ylabel('Doluluk Orani %');
grid on;

figure('Name', 'Terminal Voltaj Karsilastirmasi');
plot(zaman_ekseni, gercek_voltaj_verisi, 'b', 'LineWidth', 1.5);
hold on;
plot(zaman_ekseni, voltaj_tahminleri, 'r--', 'LineWidth', 1.5);
title('Gercek ve Tahmini Terminal Voltaji');
xlabel('Zaman saniye');
ylabel('Voltaj V');
legend('Gercek Olcum', 'Filtre Tahmini');
grid on;

hata_vektoru = gercek_voltaj_verisi - voltaj_tahminleri;
rmse_degeri = sqrt(mean(hata_vektoru.^2));
disp('Voltaj Tahmin Hatasi RMSE:');
disp(rmse_degeri);


