作成したプログラムのメモ
template ディレクトリに、プログラム作成に便利なテンプレートを配置(予定)
- DeltaF_fesdat.py
- deltaf.py
  fes.dat ファイルから状態間のΔF を導出
  前者は単一ファイルから、後者は複数ファイルから計算することができる。
- irrad_ortho_md.py
  照射石英作成用プログラム。筥崎さんのプログラムを加工
- ortho_crystal.py
  結晶化度合い導出プログラム。筥崎さんのプログラムを加工
- forDEBUG (ディレクトリ)
  XRD 計算が正しくできているか、デバックするときに使用したプログラム
  - Lorch.py
  - XRD.py
  - XRD2dim.py
  - numQ2dim.py
  - coeff_plot.py
- log2disp.py
  lammps log ファイル読み取りプログラム。ファイル形式に合わせて書き換えましょう
- timestep_reset.py
  2e9 timesteps を超えるLAMMPS 計算実行時、timestep を一度リセットする必要がある
  リセットした時間を読み取り、正しい時間に直すプログラム
- bin2trj.py
  照射石英作成時出力されるバイナリファイルをトラジェクトリファイルに変換
- grplot.py
- sq_plot.py
- num_sq_plot.py
  複数の動径分布関数、構造因子をプロット
- make2dimDat.py
  CV 2つ使用するメタダイナミクスの実行ファイル生成プログラム
- readcolvar.py
  メタダイナミクス実行結果から、集団変数プロットをプロットするプログラム
  単一または複数CVを自動で読み取りプロットします
- volumeLog.py
  ログファイルから体積や密度を読み取るプログラム？
- cif2cell.py
  ciffile から data ファイルを出力する、LAMMPS 配布プログラム
  これ使用するよりも、irrad_ortho_md.py のほうが楽かもしれない
- hills_plot1.py
- hills_plot.py
  単一および2つのCV で実行したメタダイナミクスについて、FES をプロットする。
  plumed sum_hills -hills HILLS --mintozero
  実行し、出力されるfes.dat を渡すとFES が出力される。
- morse_ev2kcalmol.py
  論文の力場をLAMMPS で使用する際、力の単位が違うので修正するために使用したプログラム
- remove_duplicate_data.py
  *<使用しない>*
  +ito stepjob 実行時、途中終了して重複した時間ステップを取り除くプログラム+
  時間ステップを取り除いても、サンプリングが狂うのでこのプログラムを使用して強引に時間ステップの辻褄を合わせることは禁止とする。
- xyz2cif.py
  xyz のトラジェクトリファイルからciffile を生成する
- num_hills_plot.py
  複数の fes.dat プロットできる。
- sumhills.py
  plumed sum_hills -hills HILLS --mintozero --stride
  stride オプションを使用してfes.dat を作成した後、FES の時間平均および経時変化をプロットする。
