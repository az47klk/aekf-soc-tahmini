clear; close; close all;

hppc = load("./25degC/5 pulse disch/03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat");
c20_ocv = load("./25degC/C20 OCV and 1C discharge tests_start_of_tests/05-08-17_13.26 C20 OCV Test_C20_25dC.mat");
drive_cycle = load("./25degC/Drive cycles/03-20-17_01.43 25degC_US06_Pan18650PF.mat");

dI = diff (hppc.meas.Current);
kesilme_anlari = find(dI > 0.5);
N_darbe = length(kesilme_anlari);

% Degerleri tutacagimiz bos dizileri donguden once olusturuyoruz
R0 = zeros(N_darbe, 1);
R1 = zeros(N_darbe, 1);
C1 = zeros(N_darbe, 1);
tau = zeros(N_darbe, 1);
I = zeros(N_darbe, 1);
V_p = zeros(N_darbe, 1);

fprintf("\n\n\n\n%80.0f ===== THEVENIN MODEL PARAMETRELERI ===== 0\n", 0.0);
for i=1:length(kesilme_anlari)
    lokal_kesilme_ani = kesilme_anlari(i);

    if i==length(kesilme_anlari) ikinci_kesilme_ani = length(hppc.meas.Voltage);
    else ikinci_kesilme_ani = kesilme_anlari(i+1); end
    
    V_eski = hppc.meas.Voltage(lokal_kesilme_ani);
    V_yeni = hppc.meas.Voltage(lokal_kesilme_ani + 2);
    I_akim = hppc.meas.Current(lokal_kesilme_ani);
    
    delta_V = V_yeni - V_eski;   
    
    R_0 = delta_V / abs(I_akim); % burada thevenin modelinin R0 direncini (yani pilin ic direnci de diyebiliriz) bulduk.

    V_tepe = max(hppc.meas.Voltage(lokal_kesilme_ani:ikinci_kesilme_ani));
    Exp_aralik = V_tepe - V_yeni;
    V_hedef = V_eski + (0.632 * Exp_aralik); % 0.632 degeri (1-e^(-1)) den gelmektedir.( V(t) = V_max * (1-e^(-t/tau)) )
    
    V_tarama = hppc.meas.Voltage(lokal_kesilme_ani:ikinci_kesilme_ani);
    t_tarama = hppc.meas.Time(lokal_kesilme_ani:ikinci_kesilme_ani);

    hedef_index = find(V_tarama >= V_hedef, 1, "first");

    tau_suresi = t_tarama(hedef_index) - t_tarama(1);

    R_1 = Exp_aralik / abs(I_akim);

    C_1 = tau_suresi/R_1;
    
    Ah(i) = hppc.meas.Ah(lokal_kesilme_ani + hedef_index);
    R0(i) = R_0;
    R1(i) = R_1;
    C1(i) = C_1;
    tau(i) = tau_suresi;
    I(i) = round(abs(I_akim), 1);
    V_p(i) = V_tepe;
    fprintf("%2.i. Tn: %9.3fs,   Tn+1: %9.3fs,   V_p: %6.4fV,   I_a: %10.6fA,   R_0: %7.4fmOhm,   R_1: %8.4fmOhm,   C_1: %9.4fF,   tau: %8.6f,   h_i: %6.0f,   Ah: %f \n",i, hppc.meas.Time(lokal_kesilme_ani),  hppc.meas.Time(ikinci_kesilme_ani),V_tepe, I_akim, R_0*1000, R_1*1000, C_1, tau_suresi, hedef_index+lokal_kesilme_ani, hppc.meas.Ah(lokal_kesilme_ani + hedef_index));
end

save("thevenin_parametreleri.mat", "Ah", "R0", "R1", "C1", "tau", "I","V_p");

% interpolasyon ile ocv-soc haritasını çıkarma
% bunu yapıyoruz çünkü c20 verilerimiz pürüzlü ve biraz gürültülü.
zaman_c20 = c20_ocv.meas.Time(1:1247);
akim_c20 = c20_ocv.meas.Current(1:1247);
voltaj_c20 = c20_ocv.meas.Voltage(1:1247);

toplam_kapasite_As = abs(trapz(zaman_c20, akim_c20));
toplam_kapasite_Ah = toplam_kapasite_As / 3600;

cekilen_yuk_dizisi = cumtrapz(zaman_c20, abs(akim_c20));
soc_c20_dizisi = 1 - cekilen_yuk_dizisi / toplam_kapasite_As;

[soc_tekil, benzersiz_indeksler] = unique(soc_c20_dizisi, 'stable');
voltaj_tekil = voltaj_c20(benzersiz_indeksler);

soc_haritasi = linspace(0, 1, 1000);
ocv_haritasi = interp1(soc_tekil, voltaj_tekil, soc_haritasi, 'pchip');
fprintf("\n\n\n\n%20.0f ===== OCV-SOC DEGERLERI ===== 0\n", 0.0);
for i=1:1000
    fprintf("%4.i. SOC: %6.5f, OCV: %6.5f\n", i, soc_haritasi(i), ocv_haritasi(i));
end

save("ocv_soc_haritasi.mat", "soc_haritasi", "ocv_haritasi");

figure('Name', 'OCV-SOC Interpolasyon Haritasi', 'Position', [100 100 800 600]);
plot(soc_tekil, voltaj_tekil, 'b', 'LineWidth', 3);
hold on;
plot(soc_haritasi, ocv_haritasi, 'r--', 'LineWidth', 2);
xlabel('SOC Doluluk Orani');
ylabel('OCV Voltaj');
legend('Gercek C20 Verisi', 'pchip Interpolasyon Haritasi', 'Location', 'northwest');
grid on;



