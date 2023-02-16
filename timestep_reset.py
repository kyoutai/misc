#!/usr/bin/env python
import numpy as np
import argparse

# reset_timestep オプションを適用し、時間をリセットした計算結果の後処理
# timestep が2e9を超えるとLAMMPSで計算不可能になるため、reset_timestep を適用
# して一旦数μ秒の時間をリセットしなければならない。
# リセットされたファイルの時間スケールを直すプログラム。
# # develop log # #
# 23/01/04 リセット1回(計算ステップ数4e9まで)に対応

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file')
par.add_argument('-o', '--outfile', help="outfile name", default="COLVAR_RESET")
par.add_argument('--reset', help="ステップ数リセット回数",
                 default=1, type=int)
par.add_argument('--option', help="インプットファイルの種類を指定",
                 choices=["colvar", "hills"], default="colvar")
args = par.parse_args()

if args.option == "colvar":
    w = open(args.outfile, "w")
    with open(args.file) as f:
        head = f.readline()
    w.write(head)
    if args.reset == 1:
        fr = np.loadtxt(args.file)
        maxidx = np.argmax(fr[:, 0])
        maxi = fr[maxidx, 0]
        fr[maxidx+1:, 0] = fr[maxidx+1:, 0] + maxi
        # for i in range(10):
        fmt = "{:8e}" + " {:9f}" * (fr.shape[1]-1) + "\n"
        print("number of data :", fr.shape[1]-1)
        for i in range(fr.shape[0]):
            # w.write("{:8e} {:9f}\n".format(fr[i, 0], fr[i, 1]))
            w.write(fmt.format(*list(fr[i])))
        exit("{} was created.".format(args.outfile))
    exit("invarid --reset argument")

elif args.option == "hills":
    w = open(args.outfile, "w")
    f = open(args.file)
    for i in range(3):
        head = f.readline()
        w.write(head)
    f.close()
    if args.reset == 1:
        fr = np.loadtxt(args.file)
        maxidx = np.argmax(fr[:, 0])
        maxi = fr[maxidx, 0]
        fr[maxidx+1:, 0] = fr[maxidx+1:, 0] + maxi
        # for i in range(10):
        for i in range(fr.shape[0]):
            w.write("{:8e} {:18.14f} {:3d} {:18.14f} {:3d}\n".format(
                fr[i, 0], fr[i, 1], int(fr[i, 2]), fr[i, 3], int(fr[i, 4])))
        exit("{} was created.".format(args.outfile))
    exit("invarid --reset argument")
exit("invarid --option argument")
