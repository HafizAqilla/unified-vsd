# Pembaruan Progress — Model VSD Lumped Parameter
**Muhammad Hafiz Aqilla S & Aisyah Dzakiyah Muthmainnah**
**April 2026**

---

## Ringkasan Singkat

Sejak dokumen framework sebelumnya, terdapat beberapa perubahan besar yang telah diimplementasikan. Model yang sebelumnya hanya diuji dengan **data dummy** kini telah dijalankan menggunakan **data klinis pasien nyata (Pasien Reyna, usia 3 tahun 2 bulan, 13.4 kg)**. Selain itu, metode analisis sensitivitas diganti dari Sobol langsung ke pendekatan berbasis **PCE (Polynomial Chaos Expansion) menggunakan UQLab**, yang jauh lebih efisien. Proses kalibrasi juga diperkuat dengan strategi multi-start.

---

## A. Perbandingan Hasil: Sebelum vs Sekarang

### Hasil lama (data dummy)

Pada laporan sebelumnya, model dijalankan dengan data dummy. Hasilnya banyak metrik yang meleset jauh:

| Metrik | Error Lama |
|--------|-----------|
| PAP_min | **−62.8%** |
| PAP_max | **−34.4%** |
| PAP_mean | **−38.4%** |
| PVR | **−57.1%** |
| LVEDV | **−43.7%** |
| LVESV | **−37.8%** |
| RVEDV | **−36.1%** |
| RAP_mean | **+35.1%** |
| SAP_mean | 1.6% ✓ |
| QpQs | −1.1% ✓ |

Hanya 2 dari 13 metrik yang masuk kategori "sangat baik" (error < 5%).

---

### Hasil terbaru (data klinis pasien Reyna, iterasi terbaik)

Iterasi kalibrasi terbaik yang dicapai menghasilkan:

| Metrik | Unit | Klinis | Model | Error | Status |
|--------|------|--------|-------|-------|--------|
| PAP_min | mmHg | 10 | 9.89 | −1.1% | ✅ < 5% |
| PAP_max | mmHg | 20 | 20.06 | +0.3% | ✅ < 5% |
| PAP_mean | mmHg | 15 | 14.67 | −2.2% | ✅ < 5% |
| SAP_mean | mmHg | 71.3 | 70.86 | −0.6% | ✅ < 5% |
| QpQs | — | 1.194 | 1.18 | −1.2% | ✅ < 5% |
| CO_Lmin | L/min | 3.423 | 3.427 | +0.1% | ✅ < 5% |
| LVEF | — | 0.528 | 0.583 | +10.4% | ⚠ perlu perbaikan |
| SVR | WU | 19.37 | 22.28 | +15.0% | ⚠ perlu perbaikan |
| RAP_mean | mmHg | 5 | 5.76 | +15.2% | ⚠ batas terima |
| LVEDV | mL | 41 | 49.4 | +20.5% | ⚠ keterbatasan data |
| RVEDV | mL | 30.5 | 36.6 | +20.1% | ⚠ keterbatasan data |

**RMSE Terkalibrasi: 0.1084** (perbaikan: **+50.8%** dibanding baseline).

**6 dari 22 metrik masuk kategori sangat baik (error < 5%)**, yaitu seluruh tekanan paru (PAP), tekanan sistemik rata-rata (SAP_mean), rasio aliran (QpQs), dan cardiac output (CO).

---

### Hasil run terbaru (dengan GSA aktif + bobot baru)

Iterasi terbaru setelah pengaktifan GSA dan penyesuaian bobot kalibrasi:

| Metrik | Error | Status |
|--------|-------|--------|
| PAP_min | −3.1% | ✅ < 5% |
| PAP_max | +0.7% | ✅ < 5% |
| SAP_mean | −3.0% | ✅ < 5% |
| LVEF | **−2.1%** | ✅ < 5% |
| SVR | +6.0% | ↑ membaik dari +15% |
| PAP_mean | −5.2% | ≈ batas terima |
| QpQs | −15.1% | ↓ perlu perhatian |

LVEF kini masuk kategori < 5%, dan SVR turun drastis dari +15% ke +6%. QpQs mengalami regresi sementara karena perubahan distribusi bobot — akan diperbaiki pada run berikutnya.

---

## B. Apa yang Berubah dalam Kode

### 1. Data klinis nyata

Model sekarang menggunakan data kateterisasi dan ekokardiografi pasien Reyna secara langsung. Sebelumnya menggunakan data dummy yang tidak mencerminkan fisiologi pasien sebenarnya.

### 2. Analisis Sensitivitas: dari Sobol Langsung ke PCE via UQLab

**Sebelumnya:** analisis sensitivitas menggunakan metode Sobol langsung (Saltelli sampling), yang memerlukan `N × (d + 2)` simulasi ODE penuh.

- Dengan N=16 dan d=19: hanya **336 simulasi** → hasil tidak stabil, rawan salah pilih parameter aktif
- Untuk N=256: butuh **5.376 simulasi** → sangat lambat (3–6 jam hanya untuk GSA)

**Sekarang:** diganti dengan sobol **PCE (Polynomial Chaos Expansion) berbasis UQLab**.

- Hanya butuh **200 simulasi** untuk melatih surrogate polynomial
- Setelah surrogate selesai, indeks Sobol dihitung **secara analitik** dari koefisien polinomial — tanpa simulasi tambahan sama sekali
- Total: 200 simulasi → **840 nilai indeks Sobol** (21 parameter × 20 metrik × S1+ST)
- Waktu GSA turun dari **3–6 jam** menjadi **20–40 menit**

Dilakukan **dua kali** dalam satu pipeline:
- Sebelum kalibrasi: untuk menentukan parameter mana yang paling berpengaruh
- Setelah kalibrasi: untuk melihat bagaimana struktur sensitivitas berubah setelah parameter dikalibrasi ke data pasien

### 3. Kalibrasi multi-start (Stage 2)

Sebelumnya kalibrasi hanya berjalan dari satu titik awal (Stage 1). Sekarang ditambahkan **Stage 2: 5 restart dari titik awal yang bervariasi** (0.3×, 0.5×, 1.8×, 2.5×, dan acak). Solusi terbaik dari semua restart yang diambil, sehingga risiko terjebak di minimum lokal berkurang.

### 4. Penyesuaian bobot kalibrasi berbasis prioritas klinis

Bobot pada fungsi objektif disesuaikan berdasarkan kepercayaan data:

| Kelompok | Metrik | Bobot |
|----------|--------|-------|
| Prioritas utama | CO_Lmin, QpQs, SAP_mean | 5.0 – 6.0 |
| Tekanan paru | PAP_mean, PAP_max | 3.5 – 4.0 |
| Tekanan atrium | RAP_mean | 5.0 |
| Volume | LVEDV, RVEDV | 2.0 (diturunkan) |

Volume LVEDV/RVEDV sengaja diberi bobot rendah karena metode pengukuran klinisnya (Teichholz dari ekokardiografi M-mode) memiliki ketidakpastian ±15–20% pada pasien VSD dengan beban volume berlebih. Ini keterbatasan data, bukan keterbatasan model.

---

## C. Apa yang Sudah Terpenuhi dari Improvement Plan Sebelumnya

| Rencana Perbaikan (dokumen lama) | Status |
|----------------------------------|--------|
| Naikkan sampel Sobol ke N≥128 | ✅ Diselesaikan — diganti PCE 200 sampel, lebih efisien |
| Bangun ulang parameter aktif dari GSA yang lebih stabil | ✅ PCE LARS lebih stabil dari Sobol langsung N=16 |
| Bobot adaptif untuk metrik error > 20% | ✅ Bobot disesuaikan per kelompok klinis |
| Multi-start L-BFGS-B dari titik awal bervariasi | ✅ Stage 2 dengan 5 restart |
| Warm-start dari solusi terbaik sebelumnya | ✅ Stage 1 → Stage 2 menggunakan solusi Stage 1 sebagai referensi |

---

## D. Keterbatasan yang Terdokumentasi

**LVEDV/RVEDV masih ~20% error** — ini bukan kegagalan model, melainkan inkonsistensi internal data klinis:
- Cardiac output dari kateter (3.42 L/min) dengan HR=119 bpm → stroke volume = 28.8 mL
- Volume dari Teichholz: LVEDV=41 mL, LVESV=19.3 mL → stroke volume = 21.7 mL → CO hanya 2.58 L/min

Kedua angka tidak bisa dipenuhi secara bersamaan. Model mengutamakan CO dari kateter (kepercayaan lebih tinggi), sehingga volume bergeser. Hal ini telah didokumentasikan dalam kode dan akan dijelaskan sebagai keterbatasan pengukuran dalam laporan akhir.

---

## E. Langkah Selanjutnya

1. Jalankan pipeline lengkap untuk **Pasien Keisya** menggunakan konfigurasi yang sama
2. Bandingkan hasil GSA pre dan post kalibrasi untuk kedua pasien (kontribusi utama paper)
3. Finalisasi figur: pressure traces, PV loop, bar chart Sobol indices per metrik
4. Tulis bagian Methods dengan penjelasan PCE-Sobol dan justifikasi parameter kalibrasi
