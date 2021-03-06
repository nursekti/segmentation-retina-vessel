# segmentation-retina-vessel

Metode Segmentasi Berbasis Thresholding dan Clustering
 
Metode ini merupakan metode paling optimal untuk melakukan segementasi pembuluh darah yang didapatkan praktikan hingga saat ini. Keoptimalan tersebut bisa ditunjukkan hasil citra yang ditunjukkan dari setiap langkah. Setiap langkah dalam metode tersebut memperkuat hasil citra yang dihasilkan pada langkah sebelumnya dan akhirnya mengarah pada hasil akhir yang optimal. Penjelasan rinci setiap tahapan/langkah adalah sebagai berikut.
	Input citra fundus retina dan evaluasi pada kanal hijau
Citra fundus retina yang digunakan sebagai input memiliki ekstensi .tif dan diambil dari database DRIVE. Citra diolah dalam kanal hijau karena kontrasnya yang paling tinggi jika dibandingkan dengan kedua kanal lainnya.
  
Gambar 2.1 Citra fundus asli dan citra fundus pada kanal hijau
	Preprocessing citra dengan metode CLAHE
CLAHE merupakan singkatan dari Contrast Limited Adaptive Histogram Equalization. Metode CLAHE ini menerapkan algoritma yang membagi citra ke dalam daerah-daerah kontekstual dan mengaplikasikan histogram equalization ke dalam setiap daerah tersebut. CLAHE menghitung histogram lokal pada setiap piksel citra retina dan melakukan histogram clipping, renormalisasi histogram, dan pemetaan piksel keluaran dengan intensitas yang sebanding dengan tingkatnya dalam histogram. Misalkan hi adalah histogram citra dan (m x m) adalah daerah kontekstual, maka tingkat rp untuk sebuah pixel dengan intensitas p dihitung menggunakan formula,
r_(p )=(∑_(i=0)^p▒〖min⁡(β,h_i )+(p+1)  x ((∑_(k=0)^N▒max⁡(0,h_k-β) )/(m x m)) 〗)
dengan clip limit β menentukan batas kontras peningkatan kualitas citra dan ∑_(i=0)^p▒min⁡(β,h_i )   mendeskripsikan tingkatnya di clipped histogram. Setiap daerah akan memiliki jumlah potongan piksel yang berbeda. Hal tersebut berguna untuk mendistribusikan kembali bagian histogram yang melebihi batas klip β secara merata di antara semua histogram yang ada untuk menormalisasi perhitungan tingkat di wilayah yang berbeda. Normalisasi tersebut diwakili oleh ∑_(i=0)^p▒((∑_(k=0)^N▒max⁡(0,h_k-β) )/(m x m))  dengan  h_k adalah histogram di daerah yang berbeda. Tingkat intensitas iin pada koordinat (x, y) dihitung dan diskalakan sehingga menghasilkan variabel r dengan 0,0 ≤ r ≤ 1,0. Keluaran level intensitas iout dalam beberapa interval gray scale di antara i1 dan i2 yang mengikuti persamaan
i_out=i_1+r x (i_2-i_1 )
Dalam kode yang diaplikasikan dengan metode ini, citra dibagi menjadi 75 baris dan 75 kolom (‘NumTiles’, [75 75]). Selain itu, clip limit yang digunakan adalah 0,09 sehingga menghasilkan kontras antara pembuluh darah yang akan disegmentasi dengan latar citranya.
 
Gambar 2.2 Citra Hasil Metode CLAHE
	Meningkatkan kualitas citra dengan filter average dan median
Citra hasil algoritma CLAHE masih menghasilkan noise sehingga harus dilakukan filterring. Dalam kode ini digunakan filter average yang mencari rata-rata dari piksel dengan piksel sekitarnya dan filter median yang mencari nilai tengah dari piksel terhadap piksel di sekitarnya juga. Filter average yang diaplikasikan menggunakan matriks 3x3 sehingga setiap piksel dirata-ratakan dengan delapan penjuru tetangganya. Selanjutnya, noise yang tersisa akan dihilangkan dengan filter median karena banyak noise-nya berupa salt & pepper. Selain itu, citra juga dilakukan peningkatan kualitas dengan adjustment yang memetakkan nilai intensitas suatu pixel ke nilai tertentu sehingga kontras pembuluh darah dengan latarnya menjadi semakin terlihat.
 
Gambar 2.3 Citra Hasil Filtering dan Peningkatan Kualitas
	Segmentasi dengan metode ISODATA Thresholding
Setelah dikontraskan dengan latarnya, dilakukan filtering serta ditingkatkan kualitasnya, segmentasi pembuluh darah dilakukan dengan metode ISODATA thresholding. Metode tersebut menghitung threshold global citra dengan melakukan perhitungan secara iteratif pada metode ISODATA. Awalnya, histogram citra disegmentasi menjadi dua bagian dengan threshold awal setengah dari rentang dinamis maksimum. Sampel pertama dihitung dari rata-rata piksel dari latar depan, sedangkan sampel kedua dihitung dari rata-rata piksel dari latar belakang. Threshold selanjutnya dihitung dari rata-rata kedua sampel tersebut. Perhitungan terus dilakukan hingga nilai threshold tidak berubah lagi. Nilai akhir threshold digunakan untuk membuat citra menjadi citra biner dengan batasnya adalah nilai threshold tersebut. 
 
Gambar 2.4 Citra Hasil Thresholding dengan Metode ISODATA
	Melakukan background substraction
Pada tahap ini citra dikurangkan dengan latarnya sehingga menghilangkan noise-noise terutama yang terletak di luar gambar bola mata. Metode ini digunakan agar objek yang disegmentasi bebas dari noise atau setidaknya meminimalisr kesalahan segementasi. 
 
Gambar 2.5 Citra Hasil Background Substraction
	Menghilangkan noise dengan operasi morfologi
Operasi morfologi yang digunakan adalah operasi opening. Operasi ini menghilangkan sisa-sisa noise dan kesalahan-kesalahan segementasi pada tahap sebelumnya. Pada kode, operasi opening ini digunakan untuk menghilangkan piksel-piksel dengan ukuran di bawah 150.
 
Gambar 2.6 Citra Hasil Akhir Segmentasi
