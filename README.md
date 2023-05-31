作成したプログラムのメモ
template ディレクトリに、プログラム作成に便利なテンプレートを配置(予定)
** PLUMED 計算準備関係のプログラム **
- morse_ev2kcalmol.py
  morse 力場のパラメータを LAMMPS unit real に変換するプログラム
- buckParameter.py
  SHIK, buckingham 力場のパラメータを LAMMPS unit real に変換するプログラム
- irrad_ortho_md.py
  照射石英作成用プログラム。直方体に対応させた。
- timestep_reset.py
  2e9 timesteps を超えるLAMMPS 計算実行時、timestep を一度リセットする必要がある
  リセットした時間を読み取り、正しい時間に直すプログラム
- bin2trj.py
  照射石英作成時出力されるバイナリファイルをトラジェクトリファイルに変換
- metaDprep.py
  メタダイナミクスに使用するインプットファイルを作成するプログラム
  我々の環境ですぐに計算実行できる状態になるように、ファイルを書込み準備する。
** PLUMED 計算結果後処理関係のプログラム **
- DeltaF_fesdat.py
  メタダイナミクス計算の後処理に使用
  COLVAR ファイル読み取って、2つの準安定状態間の自由エネルギー差を導出
- deltaf.py
  fes.dat ファイルから状態間のΔF を導出
- ortho_crystal.py
  結晶化度合い導出プログラム。
- forDEBUG (ディレクトリ)
  XRD 計算が正しくできているか、デバックするときに使用したプログラム
  - Lorch.py
    Lorch function をプロットして関数形を確認した。
  - XRD.py, NaXRD.py
    XRDピークをプロットするプログラム。
  - XRD2dim.py
    2つの原子種類でXRDピークを算出する。
    現在は, Si-O だけに対応している。
  - numQ2dim.py
    plumed で XRDピーク を正しく導出できたか確認するために、計算のインプットファイルを作成するプログラム
  - coeff_plot.py
    原子散乱因子プロットするプログラム
  - compare_TrjCOLVAR.py
    LAMMPS, PLUMED それぞれの計算結果を入力して、XRD ピークを同時にプロットする。
    PLUMED の計算結果が正しいか確認する。
- log2disp.py
  lammps log ファイル読み取りプログラム。ファイル形式に合わせて書き換えましょう
- grplot.py
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
- remove_duplicate_data.py
  *<使用しない>*
  +ito stepjob 実行時、途中終了して重複した時間ステップを取り除くプログラム+
  時間ステップを取り除くとサンプリングが狂う。
  このプログラムを使用して強引に時間ステップの辻褄を合わせることは禁止とする。
  ito で計算を回すことを止めるべき。
- xyz2cif.py
  xyz のトラジェクトリファイルからciffile を生成する
- num_hills_plot.py
  複数の fes.dat プロットできる。
- sumhills.py
  plumed sum_hills -hills HILLS --mintozero --stride
  stride オプションを使用してfes.dat を作成した後、FES の時間平均および経時変化をプロットする。
