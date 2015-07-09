import subprocess
import pandas as pd

def AnnotateSg(str_f_sgsam, str_refgem, str_of_sgbed, str_of_sgmap):
    str_proc = r'''
ref=%s
sgsam=%s
sgbed=%s
sgmap=%s

samtools view -bT ${FASTADB}/${ref}.fa ${sgsam} | \
    bamToBed -i - > ${sgbed}

map_process() {
    cut -f1-4,6,10 "$1" | \
	awk 'BEGIN {FS = "\t"; OFS = "\t"; tmp = "";}
             {tmp = $6; $6 = $4; $4 = tmp; print $0;}' | \
	uniq -f5
}

map_l1() {
    bedtools intersect -a "$1" -b ${BEDDB}/${ref}ref.bed -wb | \
	map_process -
}

map_l2() {
    bedtools intersect -a "$1" -b ${BEDDB}/${ref}ref.bed -v | \
	bedtools intersect -a - -b ${BEDDB}/${ref}ucsc.bed -wb | \
	map_process -
}

map_l3() {
    bedtools intersect -a "$1" -b ${BEDDB}/${ref}ref.bed -v | \
	bedtools intersect -a - -b ${BEDDB}/${ref}ucsc.bed -v | \
	bedtools intersect -a - -b ${BEDDB}/${ref}gencode.bed -wb | \
	map_process -
}

map_l3b() {
    bedtools intersect -a "$1" -b ${BEDDB}/${ref}ref.bed -v | \
	bedtools intersect -a - -b ${BEDDB}/${ref}ucsc.bed -v | \
	awk 'BEGIN {FS = "\t"; OFS = "\t";}
             {$7 = $1; $8 = $2; $9 = $3; 
               $10 = "None"; $11 = 0; $12 = $6;print $0;}' | \
	map_process -
}

map_l4b() {
    bedtools intersect -a "$1" -b ${BEDDB}/${ref}ref.bed -v | \
	bedtools intersect -a - -b ${BEDDB}/${ref}ucsc.bed -v | \
	bedtools intersect -a - -b ${BEDDB}/${ref}gencode.bed -v | \
	awk 'BEGIN {FS = "\t"; OFS = "\t";}
             {$7 = $1; $8 = $2; $9 = $3; 
               $10 = "None"; $11 = 0; $12 = $6;print $0;}' | \
	map_process -
}

sgtest1=$(bedtools intersect -a ${sgbed} -b ${BEDDB}/${ref}ref.bed -v)
if [ "$sgtest1" = "" ]
then
    map_l1 ${sgbed} > ${sgmap}
    exit 1
fi

sgtest2=$(bedtools intersect -a ${sgbed} -b ${BEDDB}/${ref}ref.bed -v | \
		 bedtools intersect -a - -b ${BEDDB}/${ref}ucsc.bed -v)
if [ "$sgtest2" = "" ]
then
    cat <(map_l1 ${sgbed}) <(map_l2 ${sgbed}) > ${sgmap}
    exit 2
fi

if [ ! -f ${BEDDB}/${ref}gencode.bed ]
then
    cat <(map_l1 ${sgbed}) <(map_l2 ${sgbed}) <(map_l3b ${sgbed}) > ${sgmap}
    exit 4
fi

sgtest3=$(bedtools intersect -a ${sgbed} -b ${BEDDB}/${ref}ref.bed -v | \
		 bedtools intersect -a - -b ${BEDDB}/${ref}ucsc.bed -v | \
		 bedtools intersect -a - -b ${BEDDB}/${ref}gencode.bed -v)
if [ "$sgtest3" = "" ]
then
    cat <(map_l1 ${sgbed}) <(map_l2 ${sgbed}) <(map_l3 ${sgbed}) > ${sgmap}
    exit 3
else
    cat <(map_l1 ${sgbed}) <(map_l2 ${sgbed}) <(map_l3 ${sgbed}) <(map_l4b ${sgbed}) > ${sgmap}
    exit 5
fi
    '''% (str_refgem, str_f_sgsam, str_of_sgbed, str_of_sgmap)
    int_status = subprocess.call(str_proc, shell=True, executable='/bin/bash')
    return int_status
    
def __SeqRc(str_seq):
    lst_rcseq = [c for c in str_seq[::-1]]
    for i in xrange(len(lst_rcseq)):
        if lst_rcseq[i] == 'A':
            lst_rcseq[i] = 'T'
        elif lst_rcseq[i] == 'T':
            lst_rcseq[i] = 'A'
        elif lst_rcseq[i] == 'C':
            lst_rcseq[i] = 'G'
        elif lst_rcseq[i] == 'G':
            lst_rcseq[i] = 'C'
    return ''.join(lst_rcseq)

def OrganizeSgsam(str_f_sgsam, str_of_psg):
    f_sgsam = open(str_f_sgsam, 'r')
    f_psg = open(str_of_psg, 'w')

    gn_sgsam = (str_line.strip().split('\t') for str_line in f_sgsam)

    idx_sgid = 0
    idx_flag = 1
    idx_chr = 2
    idx_beg = 3
    idx_seq = 9
    
    str_strand = ''
    int_beg = 0
    int_end = 0
    int_csite = 0
    str_seq = ''
    for lst_sgsam in gn_sgsam:
        if lst_sgsam[idx_flag] == '16':
            str_strand = '-'
        elif lst_sgsam[idx_flag] == '0':
            str_strand = '+'

        if str_strand == '+':
            str_seq = lst_sgsam[idx_seq]
            int_beg = int(lst_sgsam[idx_beg])
            int_end = int(lst_sgsam[idx_beg]) - 1 + len(str_seq)
            int_csite = int_end - 5
        elif str_strand == '-':
            str_seq = __SeqRc(lst_sgsam[idx_seq])
            int_end = int(lst_sgsam[idx_beg])
            int_beg = int(lst_sgsam[idx_beg]) - 1 + len(str_seq)
            int_csite = int_end + 5
        lst_psg = [lst_sgsam[idx_sgid], lst_sgsam[idx_chr]]
        lst_psg.append(str_strand)
        lst_psg.append(str(int_beg))
        lst_psg.append(str(int_end))
        lst_psg.append(str_seq)
        lst_psg.append(str(int_csite))
        f_psg.write('\t'.join(lst_psg) + '\n')
    f_psg.close()

def MergeSg(str_f_psg, str_f_sgmap, str_of_sg):
    dfm_psg = pd.read_csv(str_f_psg, sep='\t', header=None)
    dfm_sgmap = pd.read_csv(str_f_sgmap, sep='\t', header=None)

    idx_sgid_psg = 0
    idx_sgid_sgmap = 5

    dfm_sg = pd.merge(dfm_psg, dfm_sgmap, left_on=idx_sgid_psg, right_on=idx_sgid_sgmap)
    dfm_sg.columns=range(len(dfm_sg.columns))

    idx_sgid = 0
    idx_chr = 2
    idx_strand = 3
    idx_gene = 11
    idx_beg = 4
    idx_end = 5
    idx_seq = 6
    idx_csite = 7
    dfm_sg = dfm_sg[[idx_sgid, idx_chr, idx_strand, idx_gene, idx_beg, idx_end, idx_seq, idx_csite]]

    dfm_sg.to_csv(str_of_sg, sep='\t', header=None, index=None)
