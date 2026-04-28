# Penjelasan Matematis: Multi-Start Optimasi & PCE-Sobol

---

## B.3 Optimasi Multi-Start (Tahap 2)

Kalibrasi diselesaikan sebagai masalah optimasi terbatas:

$$\mathbf{x}^* = \arg\min_{\mathbf{x} \in [\mathbf{lb},\, \mathbf{ub}]} \; J(\mathbf{x}), \qquad J(\mathbf{x}) = \sum_{k} w_k \left(\frac{y_k(\mathbf{x}) - y_k^{\text{clin}}}{y_k^{\text{clin}}}\right)^2$$

di mana $\mathbf{x} \in \mathbb{R}^{19}$ adalah vektor parameter, $y_k(\mathbf{x})$ adalah keluaran model, dan $w_k$ adalah bobot klinis per metrik. Karena $J$ non-konveks, solver L-BFGS-B berbasis gradien dapat konvergen ke minimum lokal yang berbeda tergantung titik awal $\mathbf{x}_0$. Untuk mengatasinya, Stage 2 menjalankan lima restart independen dari titik awal:

$$\mathbf{x}_0^{(r)} = \text{clip}\!\left(s_r \cdot \mathbf{x}_0,\ \mathbf{lb},\ \mathbf{ub}\right), \quad s_r \in \{0.3,\ 0.5,\ 1.8,\ 2.5,\ \mathbf{u}\}, \quad \mathbf{u} \sim \mathcal{U}[\mathbf{lb}, \mathbf{ub}]$$

Solusi akhir diambil dari restart dengan nilai objektif terkecil: $\mathbf{x}^* = \arg\min_r J(\mathbf{x}^{(r)})$.

---

## C. PCE-Sobol via UQLab

Model ODE $\mathcal{M}: \mathbf{x} \mapsto y$ diasumsikan dapat diaproksimasi oleh ekspansi polinomial ber-orthonormal terhadap distribusi input $p(\mathbf{x})$:

$$y \approx \hat{y}(\mathbf{x}) = \sum_{\boldsymbol{\alpha} \in \mathcal{A}} c_{\boldsymbol{\alpha}}\, \Psi_{\boldsymbol{\alpha}}(\mathbf{x})$$

di mana $\Psi_{\boldsymbol{\alpha}}$ adalah basis polinomial multi-dimensi (Legendre untuk input Uniform), $c_{\boldsymbol{\alpha}}$ adalah koefisien yang dicari, dan $\mathcal{A}$ adalah himpunan multi-indeks yang dipilih oleh LARS. Indeks Sobol kemudian diperoleh **secara analitik** dari koefisien tersebut tanpa simulasi tambahan:

$$S_i = \frac{\sum_{\boldsymbol{\alpha} \in \mathcal{A}_i} c_{\boldsymbol{\alpha}}^2}{\sum_{\boldsymbol{\alpha} \in \mathcal{A}} c_{\boldsymbol{\alpha}}^2}, \qquad S_i^T = \frac{\sum_{\boldsymbol{\alpha} \in \mathcal{A}_i^T} c_{\boldsymbol{\alpha}}^2}{\sum_{\boldsymbol{\alpha} \in \mathcal{A}} c_{\boldsymbol{\alpha}}^2}$$

di mana $\mathcal{A}_i$ adalah himpunan multi-indeks yang **hanya** melibatkan parameter $x_i$, dan $\mathcal{A}_i^T$ adalah himpunan yang melibatkan $x_i$ dalam kombinasi apapun. Hasilnya setara dengan estimator Saltelli–Jansen pada Monte Carlo, namun hanya membutuhkan $N=200$ evaluasi ODE untuk training, bukan $N \times (d+2) = 4{,}620$ evaluasi pada metode langsung.
