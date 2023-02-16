#!/usr/bin/env python
import copy
import argparse

# ito stepjob において生じた問題：計算途中で時間切れになるとrst 生成されない
# 出力されるデータに不具合が生じ計算が重複して、熱力学統計データが崩壊する。
# 例：10ns から20ns までの計算が終了せず、20ns スタートのrst ファイルが生成されないため、永遠と10ns スタートの計算を繰り返す
# これを防ぐために、in ファイルで実行時間を短か目に指定しなければならない。
# それでも途中終了してしまう場合があるため、途中終了した計算結果を整形するプログラムをここに記す。

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file')
par.add_argument('-o', "--outfile",  help='output file', default="hoge")
par.add_argument('--mode', help="spieces of input file", default="colvar",
                 choices=["colvar", "hills"])
args = par.parse_args()


def COLVARMode():
    w = open(args.outfile, "w")
    buf = ""
    with open(args.file) as f:
        # First: process INITIALIZE data
        # read first 2lines
        comment = f.readline()
        line = f.readline()
        # begin Conditional branch
        while (line[0] != comment[0]):     # "#! FIELDS" -> break
            buf += line
            line = f.readline()
        w.write(comment + buf)
        w.write(line)
        buf = ""
        # Second: process RUN data
        while line:
            line = f.readline()
            buf += line
            # store head
            while (line[0] != comment[0]):     # "#! FIELDS" -> break
                buf += line
                tail = copy.deepcopy(line)  # deep copy
                line = f.readline()
                if not line:    # EOF
                    w.write(comment + buf)
                    break
            if not line:        # EOF
                break
            # compare tail & head
            line = f.readline()
            head = copy.deepcopy(line)  # deep copy
            head, tail = float(head.split()[0]), float(tail.split()[0])
            if ((head - tail) == 0):
                w.write(comment + buf)
                w.write(line)
            buf = ""


def HILLSMode():
    # #! FIELDS time COORD sigma_COORD height biasf
    # #! SET multivariate false
    # #! SET kerneltype gaussian
    # 200  64.81873359429379  2  14.69387755102041     50
    w = open(args.outfile, "w")
    buf, comment = "", ""
    with open(args.file) as f:
        # First: process INITIALIZE data
        # read first 3lines
        line = f.readline()
        reference = copy.deepcopy(line)  # deep copy
        comment += line
        line = f.readline()
        comment += line
        line = f.readline()
        comment += line
        # begin Conditional branch
        while (line[:8] != reference[:8]):     # "#! FIELDS" -> break
            buf += line
            line = f.readline()
        w.write(comment + buf)
        w.write(line)
        buf = ""
        line = f.readline()
        line = f.readline()
        # Second: process RUN data
        while line:
            line = f.readline()
            buf += line
            # store head
            while (line[:8] != reference[:8]):     # "#! FIELDS" -> break
                buf += line
                tail = copy.deepcopy(line)  # deep copy
                line = f.readline()
                if not line:    # EOF
                    w.write(comment + buf)
                    break
            if not line:
                break
            # compare tail & head
            line = f.readline()
            line = f.readline()
            line = f.readline()
            head = copy.deepcopy(line)  # deep copy
            head, tail = float(head.split()[0]), float(tail.split()[0])
            if ((head - tail) == 400.0):
                w.write(comment + buf)
                w.write(line)
            buf = ""


if args.mode == "colvar":
    COLVARMode()
if args.mode == "hills":
    HILLSMode()
