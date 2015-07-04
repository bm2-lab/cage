from __future__ import division
import re
from operator import itemgetter
from orderedset import OrderedSet
import pandas as pd

def __ReadSgtable(str_f_sg):
    f_sg = open(str_f_sg, 'r')
    gn_sg = (str_line.strip().split('\t') for str_line in f_sg)
    set_sg = {lst_sg[0] for lst_sg in gn_sg}
    f_sg.seek(0)
    gn_sg = (str_line.strip().split('\t') for str_line in f_sg)

    dict_chr = {}
    dict_sg = dict.fromkeys(set_sg)
    
    idx_sgid = 0
    idx_chr = 1
    idx_csite = 7
    for lst_sg in gn_sg:
        dict_sg[lst_sg[idx_sgid]] = [lst_sg[idx_sgid]] + lst_sg[2:]
        if lst_sg[1] in dict_chr.keys():
            dict_chr[lst_sg[1]].append({'sg':lst_sg[0], 'csite':int(lst_sg[-1])})
        else:
            dict_chr[lst_sg[1]] = [{'sg':lst_sg[0], 'csite':int(lst_sg[-1])}]
    for str_chr in dict_chr.keys():
        dict_chr[str_chr] = sorted(dict_chr[str_chr], key=itemgetter('csite'))
    return (dict_chr, dict_sg)

def __ParseCigar(str_cigar, re_pat):
    lst_beg = [re_obj.start() for re_obj in re_pat.finditer(str_cigar)]
    lst_end = [re_obj.end() for re_obj in re_pat.finditer(str_cigar)]
    lst_cigar = (''.join([str_cigar[i-1] for i in lst_end]), [int(str_cigar[i:j-1]) for (i, j) in zip(lst_beg, lst_end)])
    return lst_cigar

def __TestOverlap(int_beg, int_end, int_csite, int_loffset=20, int_roffset=20):
    if (int_csite + int_roffset) < int_beg:
        return -1
    elif (int_csite - int_loffset) > int_end:
        return 1
    else:
        return 0

def __LowerSg(i, j, func_comp):
    while i <= j:
        k = (i + j) // 2
        if (k == i and func_comp(k) == 0) or (func_comp(k) == 0 and func_comp(k-1) < 0):
            return k
        elif func_comp(k) >= 0:
            j = k - 1
        else:
            i = k + 1
    return -1

def __UpperSg(i, j, func_comp):
    while i <= j:
        k = (i + j) // 2
        if (k == j and func_comp(k) == 0) or (func_comp(k) == 0 and func_comp(k+1) > 0):
            return k
        elif func_comp(k) > 0:
            j = k - 1
        else:
            i = k + 1
    return -1

def __NavigateSg(int_beg, int_end, lst_dictchr, int_loffset=20, int_roffset=20):
    lst_sg = []
    i = 0
    j = len(lst_dictchr) - 1
    func_comp = lambda k: __TestOverlap(int_beg, int_end, lst_dictchr[k]['csite'], int_loffset, int_roffset)
    int_low = __LowerSg(i, j, func_comp)
    int_ups = __UpperSg(i, j, func_comp) + 1
    if int_low != -1 and int_ups != -1:
        lst_sg = [d['sg'] for d in lst_dictchr[int_low:int_ups]]
    return lst_sg

def __ProcessInd(lst_cur, dict_chr, dict_sg, re_cigar, re_ind, 
                  int_loffset=20, int_roffset=20):
    idx_idx = 0
    idx_chr = 2
    idx_beg = 3
    idx_cigar = 5
    
    lst_cigar = __ParseCigar(lst_cur[idx_cigar], re_cigar)
    int_sseqlen = sum([lst_cigar[1][i] for i in range(len(lst_cigar[0])) if lst_cigar[0][i] != 'I'])
    lst_bofstidx = [i.start() for i in re_ind.finditer(lst_cigar[0])]
    int_bofstidx = lst_bofstidx[0]
    int_bofst = sum(lst_cigar[1][0:int_bofstidx])
    int_eofstidx = [i.end() for i in re_ind.finditer(lst_cigar[0])][-1]
    int_eofst = sum(lst_cigar[1][int_eofstidx:])
    int_send = int(lst_cur[idx_beg]) - 1 + int_sseqlen
    int_idbeg = int(lst_cur[idx_beg]) - 1 + int_bofst
    int_idend = int_send + 1 - int_eofst
    int_idlen = abs(sum([lst_cigar[1][i] for i in range(len(lst_cigar[0])) if lst_cigar[0][i] == 'D']) - sum([lst_cigar[1][i] for i in range(len(lst_cigar[0])) if lst_cigar[0][i] == 'I']))
    str_idtype = ''
    str_fmstatus = ''

    if lst_cur[idx_chr] not in dict_chr.keys():
        raise KeyError
    
    lst_sg = __NavigateSg(int_idbeg, int_idend, dict_chr[lst_cur[idx_chr]], int_loffset, int_roffset)
    if not lst_sg:
        raise ValueError
    int_factor = len(lst_sg)
        
    if len(lst_bofstidx) == 1 and lst_cigar[0][int_bofstidx] == 'D':
        str_idtype = 'sdel'
    else:
        str_idtype = 'other'
    if int_idlen % 3 == 0:
        str_fmstatus = 'INF'
    else:
        str_fmstatus = 'OTF'
                
    str_batch = lst_cur[idx_idx].split('.')[0]
    for str_sg in lst_sg:
        lst_samind = [lst_cur[idx_idx], str_batch, lst_cur[idx_chr], lst_cur[idx_beg]]
        lst_samind.append(str(int_send))
        lst_samind.extend(dict_sg[str_sg])
        lst_samind.append(lst_cur[idx_cigar])
        lst_samind.append(str_idtype)
        lst_samind.append(str(int_idbeg))
        lst_samind.append(str(int_idend))
        lst_samind.append(str(int_idlen))
        lst_samind.append(str_fmstatus)
        lst_samind.append(int_factor)
        yield lst_samind

def __ProcessNonInd(lst_cur, dict_chr, dict_sg, re_cigar,
                    int_loffset=0, int_roffset=0):
    idx_idx = 0
    idx_chr = 2
    idx_beg = 3
    idx_cigar = 5

    lst_cigar = __ParseCigar(lst_cur[idx_cigar], re_cigar)
    int_sseqlen = sum([lst_cigar[1][i] for i in range(len(lst_cigar[0])) if lst_cigar[0][i] != 'I'])
    int_sbeg = int(lst_cur[idx_beg])
    int_send = int_sbeg - 1 + int_sseqlen
    int_idbeg = 0
    int_idend = 0
    int_idlen = 0
    str_idtype = 'none'
    str_fmstatus = 'None'
    
    if lst_cur[idx_chr] not in dict_chr.keys():
        raise KeyError
    
    lst_sg = __NavigateSg(int_sbeg, int_send, dict_chr[lst_cur[idx_chr]], int_loffset, int_roffset)
    if not lst_sg:
        raise ValueError
    int_factor = len(lst_sg)

    str_batch = lst_cur[idx_idx].split('.')[0]
    for str_sg in lst_sg:
        lst_samind = [lst_cur[idx_idx], str_batch, lst_cur[idx_chr], lst_cur[idx_beg]]
        lst_samind.append(str(int_send))
        lst_samind.extend(dict_sg[str_sg])
        lst_samind.append(lst_cur[idx_cigar])
        lst_samind.append(str_idtype)
        lst_samind.append(str(int_idbeg))
        lst_samind.append(str(int_idend))
        lst_samind.append(str(int_idlen))
        lst_samind.append(str_fmstatus)
        lst_samind.append(int_factor)
        yield lst_samind
    

def __IntegrateSg(str_f_sam, str_f_sg, bln_pair=True, int_idloffset=20,
                int_idroffset=20, int_loffset=0, int_roffset=0):
    dict_chr, dict_sg = __ReadSgtable(str_f_sg)
    
    f_sam = open(str_f_sam, 'r')

    gn_sam = (str_line.strip().split('\t') for str_line in f_sam)
    re_cigar = re.compile(r'\d+[A-Z]')
    re_ind = re.compile(r'[DI]')

    idx_cigar = 5
    lst_cur = []
    bln_ind = True
    if bln_pair == True:
        lst_fwd = []
        lst_rev = []
        while True:
            try:
                lst_fwd = next(gn_sam)
                lst_rev = next(gn_sam)
                bln_fwdind = True if re_ind.search(lst_fwd[idx_cigar]) else False
                bln_revind = True if re_ind.search(lst_rev[idx_cigar]) else False
                lst_cur = lst_rev if bln_revind else lst_fwd
                bln_ind = True if (bln_fwdind or bln_revind) else False
                gn_samind = None
                if bln_ind == True:
                    try:
                        gn_samind = __ProcessInd(lst_cur,dict_chr, dict_sg, re_cigar, re_ind, int_idloffset, int_idroffset)
                        for lst_samind in gn_samind:
                            yield lst_samind
                    except KeyError as e:
                        continue
                    except ValueError as e:
                        continue
                else:
                    try:
                        gn_samind = __ProcessNonInd(lst_cur, dict_chr, dict_sg, re_cigar, int_loffset, int_roffset)
                        for lst_samind in gn_samind:
                            yield lst_samind
                    except KeyError as e:
                        continue
                    except ValueError as e:
                        continue
                
            except StopIteration:
                break
    else:
        for lst_sam in gn_sam:
            lst_cur = lst_sam
            bln_ind = True if re_ind.search(lst_cur[idx_cigar]) else False
            if bln_ind == True:
                try:
                    gn_samind = __ProcessInd(lst_cur,dict_chr, dict_sg, re_cigar, re_ind, int_idloffset, int_idroffset)
                    for lst_samind in gn_samind:
                        yield lst_samind
                except KeyError as e:
                    continue
                except ValueError as e:
                    continue
            else:
                try:
                    gn_samind = __ProcessNonInd(lst_cur, dict_chr, dict_sg, re_cigar, int_loffset, int_roffset)
                    for lst_samind in gn_samind:
                        yield lst_samind
                except KeyError as e:
                    continue
                except ValueError as e:
                    continue
            
def __Dedupe(items, key=None):
    seen = OrderedSet()
    num_seen = list()
    gn_item = (item for item in items)
    while True:
        try:
            item = gn_item.next()
        except Exception as e:
            yield (None, num_seen)
            break
        else:
            val = item if key is None else key(item)
            if val not in seen:
                yield (item, None)
                seen.add(val)
                num_seen.append(1)
            else:
                num_seen[seen.index(val)] += 1

def ProcessSam(str_f_sam, str_f_sg, str_of_samind, bln_pair=True,
               int_idloffset=20, int_idroffset=20, int_loffset=0, int_roffset=0):
    idx_chr = 2
    idx_sgid = 5
    idx_idtype = 13
    idx_idbeg = 14
    idx_idlen = 16
    idx_factor = 18

    gn_lst_psamind = __IntegrateSg(str_f_sam, str_f_sg, bln_pair, int_idloffset, int_idroffset, int_loffset, int_roffset)
    gn_lst_samind = __Dedupe(gn_lst_psamind, key=lambda lst_psamind: (lst_psamind[idx_chr], lst_psamind[idx_sgid], lst_psamind[idx_idtype], lst_psamind[idx_idbeg], lst_psamind[idx_idlen], lst_psamind[idx_factor]))

    dfm_samind = pd.DataFrame()
    for (lst_samind, lst_num) in gn_lst_samind:
        if lst_num == None:
            dfm_samind = dfm_samind.append([lst_samind])
        else:
            dfm_samind[len(dfm_samind.columns)] = lst_num
    didx_count = len(dfm_samind.columns) - 1
    didx_factor = len(dfm_samind.columns) - 2
    dfm_samind[len(dfm_samind.columns)] = dfm_samind[didx_count] / dfm_samind[didx_factor]
    dfm_samind.to_csv(str_of_samind, sep='\t', header=None, index=None)
