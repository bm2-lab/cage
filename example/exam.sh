#! /bin/bash
#This script is used to test each command in cage

case $1 in
    'sg')
	python ../cage.py sg -s data/nmeth3015/nmeth3015.sgfq -o result/sg -g hg19
	;;
    'prep_se')
	python ../cage.py prep -s data/nmeth3015/nmeth3015.sg -f data/nmeth3015/nmeth3015_partial.fq -o result/prep/se -g hg19
	;;
    'prep_pe')
	python ../cage.py prep -s data/nbt2800/nbt2800.sg -f data/nbt2800/nbt2800_partial_1.fq -r data/nbt2800/nbt2800_partial_2.fq -o result/prep/pe -g mm9
	;;
    'mh')
	python ../cage.py mh -i data/nmeth3015/nmeth3015.samind -o result/mh -g hg19
	;;
    'indel')
	python ../cage.py indel -i data/nmeth3015/nmeth3015.samind -s data/nmeth3015/nmeth3015.sg -o result/indel -g hg19
	;;
    'fs_las')
	python ../cage.py fs -i data/nmeth3015/nmeth3015.iost -s data/nmeth3015/nmeth3015.sg -o result/fs/las/non-auto -g hg19 -u 35 -w 32 -m lasso
	;;
    'fs_las_a')
	python ../cage.py fs -i data/nmeth3015/nmeth3015.iost -s data/nmeth3015/nmeth3015.sg -o result/fs/las/auto -g hg19 -a -m lasso
	;;
    'fs_log')
	python ../cage.py fs -i data/nonribo/nonribo.st -s data/nonribo/nonribo.sg -o result/fs/log/non-auto -g hg19 -u 35 -w 32 -m logit
	;;
    'fs_log_a')
	python ../cage.py fs -i data/nonribo/nonribo.st -s data/nonribo/nonribo.sg -o result/fs/log/auto -g hg19 -a -m logit --init-radius 30 -r 10 --step 5
	;;
    'mt')
	python ../cage.py mt -i data/mt/n2800.iost data/mt/n3235.iost data/mt/nm3015.iost -s data/mt/n2800.iosg data/mt/n3235.iosg data/mt/nm3015.iosg -g mm9 mm9 hg19 -o result/mt
	;;
    'eval')
	python ../cage.py eval -c chr3 -b 195609062 -e 195609284 -m data/nmeth3015/nmeth3015.pkl -o result/eval -g hg19
	;;
    'eval_sg')
	python ../cage.py eval -s data/nmeth3015/nmeth3015.sg -m data/nmeth3015/nmeth3015.pkl -o result/eval_sg -g hg19
	;;
    'vis')
	python ../cage.py vis -f data/nmeth3015/nmeth3015_fesrep.xml -o result/vis
	;;
esac
