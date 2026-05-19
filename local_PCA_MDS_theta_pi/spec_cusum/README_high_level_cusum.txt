1. SNP local PCA / MDS / lostruct
   → find candidate inversion-like regions

2. Candidate-level clustering
   → assign fish into bands / karyotype-like groups
   → e.g. band1, band2, band3

3. θπ matrix for the same region
   → per-sample θπ across windows

4. Group-wise θπ CUSUM
   → test whether each group has a sharp diversity changepoint near the candidate edges

5. Boundary consensus
   → compare CUSUM changepoints across groups
   → compare with D17, GHSL, SV, nSites