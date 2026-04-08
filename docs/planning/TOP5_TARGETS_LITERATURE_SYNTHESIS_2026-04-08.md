# Rekomendasi Top-5 Calibration Targets (Unified VSD 0D)

## Tujuan Dokumen
Dokumen ini merapikan daftar target kalibrasi dan rujukan yang Anda kumpulkan untuk konteks VSD pre-surgery, dengan fokus pada metrik yang:
1. klinis bermakna,
2. terukur secara praktis,
3. sensitif terhadap parameter model 0D.

## Ringkasan Rekomendasi Top-5 (Pre-Surgery Priority)
1. Qp/Qs
2. PAP_mean (mean pulmonary arterial pressure)
3. SAP_mean (mean systemic arterial pressure / MAP)
4. LVEF
5. LVEDP

## Rasional Ilmiah per Target

### 1) Qp/Qs (Shunt Flow Ratio)
Qp/Qs adalah descriptor utama derajat shunt kiri-ke-kanan pada VSD dan menjadi pengikat langsung antara sirkuit sistemik-pulmoner melalui parameter shunt (contoh: R.vsd).

Persamaan ringkas:
$$
\frac{Q_p}{Q_s} = 1 + \frac{Q_{\text{shunt}}}{Q_s}, \qquad
Q_{\text{shunt}} \approx \frac{\overline{P}_{LV} - \overline{P}_{RV}}{R_{vsd}}
$$

Implikasi kalibrasi: error Qp/Qs akan memengaruhi hampir semua metrik tekanan dan volume hilir.

### 2) PAP_mean
PAP_mean merepresentasikan beban afterload RV dan menjadi target penting untuk parameter pulmoner (resistance/compliance).

Persamaan ringkas:
$$
\overline{P}_{PA} \approx Q_p \cdot R_{PAR} + \overline{P}_{LA}
$$

Implikasi kalibrasi: sangat informatif untuk menambatkan sisi pulmoner ketika ada shunt bermakna.

### 3) SAP_mean (MAP)
SAP_mean adalah constraint global sirkulasi sistemik yang stabil secara klinis, dan mengikat parameter resistansi sistemik.

Persamaan ringkas:
$$
\overline{P}_{SA} \approx Q_s \cdot R_{SAR} + \overline{P}_{RA}
$$

Implikasi kalibrasi: penting untuk menjaga model tidak drift ke solusi shunt-heavy namun perfusi sistemik tidak realistis.

### 4) LVEF
LVEF mengunci fungsi sistolik LV (elastance aktif), serta dipengaruhi preload dan afterload; sangat berguna untuk mencegah solusi numerik yang cocok tekanan tapi gagal fungsi pompa.

Persamaan ringkas:
$$
EF = \frac{V_{ED} - V_{ES}}{V_{ED}}
$$

### 5) LVEDP
LVEDP menambatkan fungsi diastolik LV dan menghubungkan left-heart filling ke tekanan pulmoner hulu.

Persamaan ringkas:
$$
LVEDP \approx E_{B,LV}\,(V_{ED} - V_{0,LV})
$$

Implikasi kalibrasi: metrik jembatan yang kuat antara dinamika diastolik LV dan hemodinamika paru.

## Tabel Implementasi untuk Pipeline
| Target | Tujuan fisiologis utama | Parameter model paling terkait |
|---|---|---|
| Qp/Qs | Besaran shunt L->R | R.vsd, R.PAR, R.SAR, C.PAR, C.SAR |
| PAP_mean | Beban pulmoner/RV afterload | R.PAR, C.PAR, R.PCOX, R.PVEN |
| SAP_mean | Perfusi sistemik | R.SAR, C.SAR, R.SC, C.SVEN |
| LVEF | Fungsi sistolik LV | E.LV.EA, V0.LV, afterload sistemik |
| LVEDP | Fungsi diastolik/filling LV | E.LV.EB, V0.LV, preload venous return |

## Toleransi Error yang Direkomendasikan
Untuk tujuan akhir, targetkan <5% pada top-5. Namun gunakan toleransi operasional sesuai modalitas data:
1. Catheter-derived (Qp/Qs Fick, PAP_mean, LVEDP): ideal <=5%, operasional <=8%
2. Non-invasive pressure (SAP_mean cuff/line): ideal <=5%
3. Echo-derived function/volume (LVEF, volume turunan): ideal <=5%, operasional <=10% jika kualitas data terbatas

## Status Verifikasi Referensi
Daftar referensi yang Anda berikan secara tema sudah konsisten dengan kerangka mekanistic modeling: shunt-physiology, pressure-flow coupling, sensitivity-identifiability-calibration-validation.

Namun untuk klaim final naskah, tetap lakukan verifikasi bibliografi akhir (judul, DOI, tahun, volume/issue) langsung dari DOI/PubMed publisher page agar tidak ada mismatch metadata.

Checklist final sebelum dipakai di thesis:
1. Verifikasi DOI dan metadata lengkap semua R1-R14
2. Tandai mana data berbasis cath vs echo pada tabel target
3. Cantumkan alasan jika satu metrik memakai toleransi >5% (uncertainty modality)
4. Sinkronkan top-5 ini ke objective function dan tabel validasi otomatis

## Referensi yang Dipetakan (sesuai daftar Anda)
1. Shimizu & Shishido, 2018
2. Ferrero et al., 2024
3. Ahmed et al., 2007
4. Colebank & Chesler, 2022
5. Kheyfets et al., 2023
6. Humbert et al. / JAHA primer, 2023
7. Sensitivity/optimization LPN, 2025 (PubMed 40598800)
8. Saxton et al., 2024
9. Bjordalsbakke et al., 2022
10. Haghebaert et al., 2025
11. Luo et al., 2011
12. Schiavazzi et al., 2017
13. Hanna et al., 2025
14. Colebank et al., 2024 guidelines
