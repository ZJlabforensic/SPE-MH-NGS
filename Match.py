# -*- coding: utf-8 -*-
# @Time    : 2022
# @Author  : WangHaoyu
# @E-mail  : wanghy0707@gmail.com
# @Github  :
# @desc    :


import numpy as np
import Toolsbox


class MHMatch:
    """
    用于基于微单倍型字典，对测序结果进行MH匹配，并进行reads计数
    //2021.10.19 PS.新增对引物和UMI的match相关函数
    //2021.11.30 输出的Match文件顺序为：MH_ID、primer、Allele、UMI、count
    """

    def __init__(self, dict_chr_info, ls_id, dict_chr_refalt, dict_mh_primer):
        """
        dict_chr_refalt: 由PreVcf得到的结果（嵌套字典，详细格式见PreProcessing_PreVcf）
        dict_chr_info: 由PreInfo得到的结果（嵌套字典，详细格式见PreProcessing_PreInfo）
        ls_id: 由PreInfo得到的结果（储存所有MH ID的列表）
        dict_mh_primer: 由PrePrimer得到的结果（储存各MH ID对应引物列表的字典）
        """
        self.dict_chr_info = dict_chr_info
        self.ls_id = ls_id
        self.dict_chr_refalt = dict_chr_refalt
        self.dict_mh_primer = dict_mh_primer

    """
    //2021.11.15 新增：将结果储存相关函数从match函数中独立
                 对测序结果进行储存（用于汇总和count）
                 1.输出数据：以ID-Primer-Allele-UMI-Count的嵌套字典形式
                 2.用于purity过滤的数据：ID-UMI-Allele-Count的嵌套字典形式，用于对每个MH下的UMI进行purity的计算并过滤
                 3.错配数据（mismatch）：储存MH等位基因提取时，存在未与参考分型对应SNP的MH相关信息，即原MH文件
                 4.错配数据（no_primer）：即Primer和UMI匹配后为None的数据
    """
    def _save_match(self, dict_match_results, mh_id, primer, umi, seq_extract_mh):
        """
        //2021.11.15 新增：以嵌套字典形式储存匹配结果（即1），无返回值
        """
        if primer not in dict_match_results[mh_id].keys():  # 如果该MH下的primer序列未出现过
            dict_match_results[mh_id][primer] = {seq_extract_mh: {umi: 1}}  # 创建primer字典（包括umi等）

        else:  # 如果primer序列出现过，则已有了对应primer的子字典

            if seq_extract_mh not in dict_match_results[mh_id][primer].keys():  # 如果primer下allele未出现过
                dict_match_results[mh_id][primer][seq_extract_mh] = {umi: 1}

            else:
                if umi not in dict_match_results[mh_id][primer][seq_extract_mh].keys():  # 如果等位基因未出现过
                    dict_match_results[mh_id][primer][seq_extract_mh][umi] = 1
                else:
                    dict_match_results[mh_id][primer][seq_extract_mh][umi] += 1

    def _save_for_purity(self, dict_for_purity, mh_id, umi, seq_extract_mh):
        """
        //2021.11.15 新增：以嵌套字典形式储存用于purity计算和过滤的结果（即2），无返回值
        //2022.03.15 标签格式为：ID、UMI、Allele、count
        """
        if umi not in dict_for_purity[mh_id].keys():  # 如果该MH下的UMI未出现过
            dict_for_purity[mh_id][umi] = {seq_extract_mh: 1}

        else:  # 如果UMI已出现过，则已有了对应UMI的子字典

            if seq_extract_mh not in dict_for_purity[mh_id][umi].keys():
                dict_for_purity[mh_id][umi] = {seq_extract_mh: 1}
            else:
                dict_for_purity[mh_id][umi][seq_extract_mh] = dict_for_purity[mh_id][umi][seq_extract_mh] + 1

    def _save_mismatch(self, ls_mismatch, mh_id, primer, umi, seq_extract_mh, seq_extract_mh_z):
        """
        //2021.11.15 新增：以列表形式储存mismatch结果（即3），无返回值
        """
        if "Z" in seq_extract_mh_z:  # 如果存在非参考基因型的碱基
            ls_mismatch.append([mh_id, primer, umi, seq_extract_mh, seq_extract_mh_z])
        else:
            pass

    def _save_no_primer(self, ls_no_primer, primer, sequence_single):
        """
        //2021.11.15 新增：以列表形式储存no_primer结果（即4），无返回值
        """
        if primer == "None":
            ls_no_primer.append(sequence_single)
        else:
            pass

    """//2021.11.30 序列信息提取：MHID、引物、UMI、匹配MH等位基因（无Z）、匹配MH等位基因（含Z）"""
    """提取MH等位基因"""
    def _single_match_mhpos(self, sequence_startpos, dict_sequence_info):
        """
        对单条序列进行MH匹配：与序列所在染色体上的某条MH进行匹配
        输出结果为匹配到的MH的起始POS（用于匹配dict_sequence_info中的项）
        //2022.03.14 复核：为序列匹配MH无误（基于序列起始pos和同染色体上所有MH的起始pos的相减结果，取最小值）
        """
        array_mh_startpos_all = np.array(list(dict_sequence_info.keys()))  # 提取info字典中该染色体上所有MH的起始pos并储存至数组
        matched_mh_startpos = array_mh_startpos_all[np.argmin(abs(array_mh_startpos_all - int(sequence_startpos)))]
        return matched_mh_startpos

    def _single_match_extract_f(self, sequence_end, sequence_seq_m, array_mh_pos, dict_sequence_refalt):
        """
        对dif_mh_seq>=0的情况（测序序列在覆盖MH之前）进行sequence等位基因提取
        sequence_end: 测序序列末尾POS
        sequence_seq_m: 测序序列（仅保留M+填充缺失）
        array_mh_pos: 转化为数组的、该条序列对应的MH的组成SNP POS
        dict_sequence_refalt: 该条染色体上所有MH组成SNP的参考碱基字典
        //2022.03.15 复核：提取碱基是从序列末尾POS开始数，测序序列在覆盖MH之前的碱基提取和check无问题
        """
        dif_seqend_mh = sequence_end - array_mh_pos  # 计算MH的每个POS是否在sequence范围内

        seq_extract_mh = ""
        seq_extract_mh_z = ""
        for dif, snppos in zip(dif_seqend_mh, array_mh_pos):  # 对单个组成SNP进行提取和检查（是否属于参考碱基范围）

            if dif >= 0:  # POS在sequence范围内
                base = sequence_seq_m[-1-dif]  # 提取该碱基
                seq_extract_mh += base  # 将该碱基纳入sequence提取的等位基因

                if base in dict_sequence_refalt[snppos]:  # 如果该碱基存在于VCF中，非新碱基，则直接添加
                    seq_extract_mh_z += base
                elif base not in dict_sequence_refalt[snppos]:  # 如果该碱基不存在于VCF的refalt中，则将其替换为"Z"
                    seq_extract_mh_z += "Z"

            elif dif < 0:  # POS在sequence范围外，则使用"N"作为占位符替代碱基
                seq_extract_mh += "N"
                seq_extract_mh_z += "N"

        return seq_extract_mh, seq_extract_mh_z

    def _single_match_extract_b(self, sequence_start, sequence_end, sequence_seq_m, array_mh_pos, dict_sequence_refalt):
        """
        对dif_mh_seq<0的情况（测序序列在覆盖MH之后）进行sequence等位基因提取
        //2022.03.15 复核：提取碱基是从序列末尾POS开始数，测序序列在覆盖MH之后的碱基提取和check无问题
        """
        dif_mh_seqend = array_mh_pos - sequence_end  # 计算POS是否在sequence范围内
        dif_mh_seqstart = array_mh_pos - sequence_start  # 计算POS是否在sequence范围内

        seq_extract_mh = ""
        seq_extract_mh_z = ""
        for dif1, dif2, snppos in zip(dif_mh_seqend, dif_mh_seqstart, array_mh_pos):

            if dif1 <= 0 and dif2 >= 0:  # POS在sequence范围内：既大于等于序列起始，也小于等于序列终止
                base = sequence_seq_m[-1+dif1]  # 提取该碱基
                seq_extract_mh += base  # 将该碱基纳入sequence提取的等位基因

                if base in dict_sequence_refalt[snppos]:  # 如果该碱基存在于VCF中，非新碱基，则直接添加
                    seq_extract_mh_z += base
                elif base not in dict_sequence_refalt[snppos]:  # 如果该碱基不存在于VCF的refalt中，则将其替换为"Z"
                    seq_extract_mh_z += "Z"

            elif dif1 > 0 or dif2 < 0:  # POS在sequence范围外：小于序列起始，或大于序列终止
                seq_extract_mh += "N"
                seq_extract_mh_z += "N"

        return seq_extract_mh, seq_extract_mh_z
    
    """对测序序列的测序方向进行判断（0代表正向，1代表反向），并通过判断测序方向获取UMI提取位置（序列起始或末尾）"""
    def _single_match_extract_umi_primer_check(self, sequence_seq_m, primer):
        """
        对序列是否对应某引物进行判断（仅对端部，即最起始或最末尾进行判断）
        //2021.10.19 新增对引物和UMI的match相关函数
        //2021.11.02 明确了windows_check的大小会影响引物的匹配情况（过大可导致primer_check_0和primer_check_1均大于0）
                     需进一步明确fault_tolerance_rate的设置对引物匹配情况的影响
        //2022.03.14 目前179 MH体系中最短的引物为20bp，最长为44bp
        //2022.03.14 复核：调整引物检测所使用引物长度（原本为引物的后15个片段），改为引物长度-3（允许三个错配）。检测的引
                     物长度均以此长度为准进行检测
        //2022.03.15 依据dict_mh_primer的格式修改，输入的primer变更为[0, primer]或[1, primer]，据此，更改此处的引物
                     check逻辑：如果为0引物，则只进行0引物检测；如果是1引物，则只进行1引物检测。返回结果为引物种类（0或
                     1）或None
        """
        """引物检测参数"""
        len_primer = len(primer[1]) - 3  # 获取用于检测的引物序列的长度（引物长度-3）
        windows_check = len_primer  # 滑窗检测次数（对齐sequence_m算1次，之后每滑动一次即追加一次，滑窗次数与检测引物序列长度保持一致）
        fault_tolerance_rate = 0.2  # 引物检测的容错率
        primer_fault_tolerance_rate = int(len_primer * fault_tolerance_rate)  # 该引物在检测时允许不匹配的碱基数

        """引物check，按照传入引物的是0引物还是1引物进行分别匹配"""
        if primer[0] == "0":  # 若当前输入的为0引物，则只按0引物进行匹配

            primer_f_b = [primer[1][:len_primer], primer[1][-len_primer:]]  # 引物实际检测序列，分别为取前xx碱基或后xx碱基，提高匹配率
            primer_check_0 = 0  # 对0引物检测结果进行记录
            for primer in primer_f_b:  # 对前xx碱基和后xx碱基都进行匹配

                for i in range(windows_check):  # 在检测次数范围内循环检测

                    count = 0  # 记录单次检测中不匹配的碱基数量
                    for base_seq, base_primer in zip(sequence_seq_m[i:len_primer+i], primer):

                        if base_seq != base_primer:
                            count += 1
                        else:
                            continue

                    if count <= primer_fault_tolerance_rate:  # 如果不匹配碱基数目在允许范围内
                        primer_check_0 += 1
                    else:
                        continue

            if primer_check_0 > 0:  # 如果0引物匹配成功，输出0
                return "0"
            else:  # 如果0引物匹配不成功，则输出None
                return None

        elif primer[0] == "1":  # 如果当前输入的是1引物

            primer_f_b = [Toolsbox.FormatTools.convert_base(primer[1][::-1])[:len_primer],
                          Toolsbox.FormatTools.convert_base(primer[1][::-1])[-len_primer:]]  # 引物实际检测序列，分别为取前xx碱基或后xx碱基，提高匹配率
            primer_check_1 = 0  # 对1引物检测进行记录
            for primer in primer_f_b:

                for i in range(windows_check):

                    if i == 0:  # 如果是第一次check，索引有所不同

                        count = 0
                        for base_seq, base_primer in zip(sequence_seq_m[-len_primer:], primer):

                            if base_seq != base_primer:
                                count += 1
                            else:
                                continue

                        if count <= primer_fault_tolerance_rate:  # 如果不匹配碱基数目在允许范围内
                            primer_check_1 += 1
                        else:
                            continue

                    elif i > 0:  # 非第一次check

                        count = 0
                        for base_seq, base_primer in zip(sequence_seq_m[-len_primer-i:-i], primer):

                            if base_seq != base_primer:
                                count += 1
                            else:
                                continue

                        if count <= primer_fault_tolerance_rate:  # 如果不匹配碱基数目在允许范围内
                            primer_check_1 += 1
                        else:
                            continue

            if primer_check_1 > 0:  # 如果1引物匹配成功，输出0
                return "1"
            else:  # 如果1引物匹配不成功，则输出None
                return None

    def _single_match_extract_umi(self, sequence_seq, sequence_seq_m, ls_primer):
        """
        依据MH的引物序列列表，提取该条测序序列所对应的UMI序列及引物序列
        sequence: 完整的测序序列
        ls_primer: 该测序序列对应的MH的测序引物列表
        //2021.10.19 新增：对引物和UMI的match相关函数
        //2022.03.15 复核：提取UMI的位置无误
        """
        for primer in ls_primer:  # 对单条引物序列（包含0、1引物信息的列表）进行判断
            
            primer_test_result = self._single_match_extract_umi_primer_check(sequence_seq_m, primer)
            if primer_test_result:  # 如果匹配上，则输出引物序列及对应的UMI

                if primer_test_result == "0":  # 如果是0引物
                    return primer[1], sequence_seq[-12:]
                elif primer_test_result == "1":  # 如果是1引物
                    return primer[1], sequence_seq[:12]

            else:  # 如果均未匹配上，则匹配另一条引物（如果存在），或输出None
                continue
    
    """提取单条测序结果的ID、primer、UMI、Allele、wrongAllele"""
    def _single_match_extract(self, sequence_startpos, sequence_seq_m, sequence_seq, sequence_len, matched_mh_startpos,
                              dict_sequence_refalt, dict_sequence_info):
        """
        对单条序列进行MH匹配：依据匹配MH的信息，提取sequence中的MH等位基因
        输出结果为该sequence中的等位基因
        //2022.03.15 复核：单条测序结果的ID、primer、UMI、allele、wrongallele提取无误
        """
        matched_mh_value = dict_sequence_info[matched_mh_startpos]  # 输出匹配到MH的value（即包含ID、包含SNP数目、SNPpos的列表）
        
        """进行引物匹配和UMI提取，依据提取结果进行赋值"""
        ls_primer = self.dict_mh_primer[matched_mh_value[0]]  # 提取匹配到MH的MH ID所对应的引物序列，提取到的是一个二维列表
        _ = self._single_match_extract_umi(sequence_seq, sequence_seq_m, ls_primer)  # 提取匹配的引物、UMI序列

        if _:  # 如果引物匹配成功、UMI提取成功，则分别赋值
            primer, umi = _
        else:  #
            primer, umi = "None", "None"
        
        """进行sequence中MH等位基因的提取（与primer和UMI无关）"""
        array_mh_pos = np.array(matched_mh_value[2:])  # 将MH的组成SNP POS列表转化为数组
        sequence_end = sequence_startpos + int(sequence_len) - 1  # 获取sequence的末尾POS
        dif_mh_seq = matched_mh_startpos - sequence_startpos  # 计算MH起始与sequence起始之间的距离

        if dif_mh_seq >= 0:  # 如果测序序列起始在MH起始之前
            seq_extract_mh, seq_extract_mh_z = self._single_match_extract_f(sequence_end, sequence_seq_m,
                                                                            array_mh_pos, dict_sequence_refalt)
            return matched_mh_value[0], primer, umi, seq_extract_mh, seq_extract_mh_z

        elif dif_mh_seq < 0:  # 如果测序序列起始在MH起始之后
            seq_extract_mh, seq_extract_mh_z = self._single_match_extract_b(sequence_startpos, sequence_end, sequence_seq_m,
                                                                            array_mh_pos, dict_sequence_refalt)
            return matched_mh_value[0], primer, umi, seq_extract_mh, seq_extract_mh_z

    def single_match(self, sequence_single):
        """
        对单条序列进行MH匹配
        输出结果为该sequence匹配的等位基因（包括错配替换前/后的情况（如果存在错配））、该等位基因所属MH的ID（用于字典匹配）
        //2022.03.15 复核：单条测序结果的ID、primer、UMI、allele、wrongallele提取无误
        """
        sequence_chr = sequence_single[0]  # 提取序列的染色体信息
        sequence_startpos = int(sequence_single[1])  # 提取序列的起始POS信息
        sequence_len = sequence_single[2]  # 提取序列的长度
        sequence_seq_m = sequence_single[3]  # 提取序列的序列信息（仅包含M）
        sequence_seq = sequence_single[4]  # 提取序列的序列信息（完整序列）

        dict_sequence_refalt = self.dict_chr_refalt[sequence_chr]  # 依据染色体信息提取dict_refalt中的子字典
        dict_sequence_info = self.dict_chr_info[sequence_chr]  # 依据染色体信息提取dict_info的子字典

        matched_mh_startpos = self._single_match_mhpos(sequence_startpos, dict_sequence_info)  # 输出sequence匹配到MH的起始POS

        return self._single_match_extract(sequence_startpos, sequence_seq_m, sequence_seq, sequence_len, matched_mh_startpos,
                                          dict_sequence_refalt, dict_sequence_info)

    """
    //2021.11.30 对匹配后的结果进行过滤（依据purity、依据UMI绝对值、依据allele的比例）
                 匹配结果的过滤调整为可选择过滤条件
    """
    """基于dict_for_purity的结果，对特定MH下UMI的purity进行计算"""     
    def _purity_calc(self, dict_allele_count):
        """
        //2021.11.15 新增：对包含多个Allele的UMI进行purity计算，并输出count占比最高的等位基因名称
        //2022.03.15 修正获取单个UMI下count占比最高等位基因的代码
        """
        ls_allele = []
        ls_allele_count = []
        for key, value in dict_allele_count.items():  # 对单个UMI下的每个Allele（key）-Count（value）操作
            ls_allele.append(key)
            ls_allele_count.append(value)
        
        purity = max(ls_allele_count) / sum(ls_allele_count)  # 计算purity
        max_allele = ls_allele[ls_allele_count.index(max(ls_allele_count))]  # 获取单个UMI下Count占比最高的Allele
        return purity, max_allele
    
    def _umi_judgement(self, dict_for_purity, id, purity_filter):
        """
        //2021.11.15 新增：对特定MH下所有UMI中包含Allele数目进行判断及purity计算，基于purity计算结果输出用于purity过滤的内容
        //           purity < 0.8：输出UMI（特定MH下包含该UMI的结果都将删除）
        //           0.8 <= putity < 1：输出{UMI:Allele}（特定MH下，该UMI对应的特定Allele将保留）
        """
        dict_spe_mh_purity = dict_for_purity[id]  # 提取特定MH的purity字典：{UMI: {allele1: count, allele2:count, ...}, ...}
        
        ls_umi_excluded = []  # 储存purity < 0.8的结果
        dict_umi_max_allele = {}  # 储存0.8 <= purity < 1的结果
        for umi, dict_allele_count in dict_spe_mh_purity.items():  # 对单个UMI进行操作：UMI（key）-{Allele: Count}（value）
            
            if len(dict_allele_count) == 1:  # 如果单个UMI下仅有一个Allele，则purity必为1，则跳过
                continue
            
            else:  # 单个UMI下有多个Allele
                purity, max_allele = self._purity_calc(dict_allele_count)
                
                if purity < purity_filter:  # <设定值，输出UMI（用于排除）
                    ls_umi_excluded.append(umi)
                elif purity_filter <= purity < 1:  # 设定值 <= < 1：输出{umi：max_allele}
                    dict_umi_max_allele.update({umi: max_allele})
        
        return ls_umi_excluded, dict_umi_max_allele

    def _allele_count_filter(self, ls_spe_mh_filter_match_results, allele_num_filter):
        """
        //2022.11.17 基于每个等位基因的绝对值进行过滤
        """
        """提取等位基因:计数字典"""
        dict_allele_count = {}
        count_all = 0  # 等位基因count数之和
        for single_match in ls_spe_mh_filter_match_results:

            if single_match[2] not in dict_allele_count.keys():
                dict_allele_count[single_match[2]] = single_match[4]
                count_all += single_match[4]
            else:
                dict_allele_count[single_match[2]] += single_match[4]
                count_all += single_match[4]

        """过滤"""
        ls_allele_retain = [key for key, value in dict_allele_count.items() if value > allele_num_filter]
        ls_retain_results = [line for line in ls_spe_mh_filter_match_results if line[2] in ls_allele_retain]
        return ls_retain_results

    def _allele_count_filter_ratio(self, ls_spe_mh_filter_match_results, allele_proportion_filter):
        """
        //2021.11.30 基于每个等位基因占所有等位基因的比例进行过滤
        """
        """提取等位基因:计数字典"""
        dict_allele_count = {}
        count_all = 0  # 等位基因count数之和
        for single_match in ls_spe_mh_filter_match_results:

            if single_match[2] not in dict_allele_count.keys():
                dict_allele_count[single_match[2]] = single_match[4]
                count_all += single_match[4]
            else:
                dict_allele_count[single_match[2]] += single_match[4]
                count_all += single_match[4]

        """确定单MH下要保留的等位基因类型，并进行过滤"""
        threshold = int(count_all * allele_proportion_filter + 0.5)  # 计算allele的阈值（+0.5是为了四舍五入）
        ls_allele_retain = [key for key, value in dict_allele_count.items() if value > threshold]  # 提取count数高于阈值的等位基因
        ls_retain_results = [line for line in ls_spe_mh_filter_match_results if line[2] in ls_allele_retain]
        return ls_retain_results

    """基于purity和UMI-count对dict_match_results进行过滤"""
    def match_filter(self, dict_match_results, dict_for_purity,
                     purity_filter, umi_count_filter, allele_num_filter, allele_proportion_filter):
        """
        //2021.11.15 新增：基于purity和UMI-count对dict_match_results进行过滤。过滤后输出的过滤后match文件为列表形式
        //2021.11.30 可在函数输入值中设置使用的过滤条件
        //2022.03.15 修改：对于purity条件和UMI-count条件，若同时设置该两类过滤，则单条结果必须满足两类条件才可保存；若
                     只设置单类过滤，则单条结果只需要满足单类条件即可保存。在基于上述两个条件的过滤后，再进行等位基因比例
                     过滤。若无purity过滤和UMI-count条件，则直接将索引的单条match数据储存至过滤结果中进行下一步
        //2022.11.17 新增对等位基因绝对值的过滤
        """
        ls_filter_match_results = []  # 储存过滤后的结果
        for id in self.ls_id:  # 对每一个MH进行操作
            
            ls_spe_mh_results = Toolsbox.FormatTools._match_dict_to_list(dict_match_results, id)  # 提取指定MH信息（列表）

            if purity_filter:  # 如果需要进行purity过滤
                ls_umi_excluded, dict_umi_max_allele = self._umi_judgement(dict_for_purity, id, purity_filter)  # 输出purity排除的信息
            else:
                pass

            ls_spe_mh_filter_match_results = []
            for single_match in ls_spe_mh_results:  # 对单条match数据进行操作

                """先依据purity进行过滤"""
                if purity_filter:

                    if single_match[3] in ls_umi_excluded:  # 如果UMI在排除列表中，则跳过并结束此轮循环
                        continue

                    if single_match[3] in dict_umi_max_allele.keys():

                        if dict_umi_max_allele[single_match[3]] == single_match[2]:  # 如果UMI和Allele均匹配，则储存并结束此轮循环

                            if umi_count_filter:  # 如果有UMI过滤，则对purity过滤后的结果进行UMI判断，符合则储存

                                if single_match[4] >= umi_count_filter:  # 如果count>设置值，则储存并结束此轮循环
                                    ls_spe_mh_filter_match_results.append(single_match)
                                    continue

                                else:  # 未达到设置值，则跳过循环
                                    continue

                            else:  # 如果没有UMI过滤，则直接储存并跳出循环
                                ls_spe_mh_filter_match_results.append(single_match)
                                continue

                else:
                    pass  # 不跳出对单条match进行操作的循环，继续执行下面的代码
            
                """依据UMI-count进行过滤（在不存在purity过滤的情况下，直接依据UMI-count进行过滤）"""
                if umi_count_filter:

                    if single_match[4] >= umi_count_filter:  # 如果count>设置值，则储存并结束此轮循环
                        ls_spe_mh_filter_match_results.append(single_match)
                        continue
                    else:
                        continue

                else:
                    pass

                """如果不进行purity过滤和UMI-count过滤，则直接储存结果（用于等位基因比例过滤）"""
                if not purity_filter and not umi_count_filter:
                    ls_spe_mh_filter_match_results.append(single_match)

            """依据等位基因绝对值或所占比例进行过滤（先绝对值，后比例）"""
            if allele_num_filter:
                ls_spe_mh_filter_match_results = self._allele_count_filter(ls_spe_mh_filter_match_results,
                                                                           allele_num_filter)
            else:
                pass

            if allele_proportion_filter:
                ls_spe_mh_filter_match_results = self._allele_count_filter_ratio(ls_spe_mh_filter_match_results,
                                                                                 allele_proportion_filter)
            else:
                pass

            ls_filter_match_results.extend(sorted(ls_spe_mh_filter_match_results, key=lambda x: x[4], reverse=True))
        return ls_filter_match_results

    """主函数"""
    def match(self, sam_filter_filepath=None, sam_match_filepath=None, sam_mismatch_filepath=None,
              filter=False, purity_filter=False, umi_count_filter=False,
              allele_num_filter=False, allele_proportion_filter=False):
        """
        sam_match_filepath: 储存匹配结果文件的路径
        sam_mismatch_filepath: 储存匹配错误结果文件、引物文件为None的filter行的路径
        filter: 是否对匹配结果进行过滤，只有当filter为True时才执行后续过滤
        purity_filter: 基于purity进行过滤，如果purity值为False则不进行（可选过滤），可输入浮点数
        umi_count_filter: 基于单一等位基因下单个UMI下的count数进行过滤，如果为False则不进行（可选过滤），可输入整数
        allele_proportion_filter: 基于单一MH下单个allele的count数占该MH的比例进行过滤，如果为False则不进行（可选过滤），
                                  可输入浮点数
        用于抓取匹配结果，进行计数并输出（嵌套字典）：
        dict_match_results: {mhid1: {primer1: {allele1: {umi1: count, umi2: count, ...}, allele2: ...}, primer2: ...},
                             mhid2: {primer1: {allele1: {umi1: count, umi2: count, ...}, allele2: ...}, primer2: ...}, ...}
        ls_mismatch: [mhid1, primer, umi, mismatch(origin), mismatch(replaced)],
                      mhid2, primer, umi, mismatch(origin), mismatch(replaced)], ...]
        //2021.11.01 新增：若引物匹配为结果为None，则输出该序列
        //           思路：新增一个判断，当输出的primer为None的时候，将sequence_single单独保存一个文件
        //2021.11.15 将文件储存相关函数从match函数中独立
        //2021.11.30 将Match文件的过滤更改为可选
        //2022.03.01 将NNNNATNNN类似的等位基因剔除
        //2022.11.17 新增通过等位基因绝对值的过滤
        """
        if sam_filter_filepath is None:
            raise ValueError('"sam_filter_filepath" is None')
        if sam_match_filepath is None:
            raise ValueError('"sam_match_filepath" is None')
        if sam_mismatch_filepath is None:
            raise ValueError('"sam_mismatch_filepath" is None')

        """列出储存去冗余及过滤后的sam文件名"""
        sam_filter_filename = Toolsbox.FileTools.ls_directory(sam_filter_filepath)

        """对每一个sam文件进行操作"""
        for filename in sam_filter_filename:
            sam_filter_file = Toolsbox.FileTools.open_vcf_file(sam_filter_filepath + filename)  # 打开sam文件

            dict_match_results = {id: {} for id in self.ls_id}  # 储存匹配结果（即1）
            dict_for_purity = {id: {} for id in self.ls_id}  # 储存用于purity计算的结果（即2）
            ls_mismatch = []  # 储存有碱基不匹配情况的等位基因（即3）
            ls_no_primer = []  # 储存primer值为None的filter行（即4）

            for sequence_single in sam_filter_file:  # 对一条测序信息进行操作，并判断该序列是否可以被纳入结果中
                mh_id, primer, umi, seq_extract_mh, seq_extract_mh_z = self.single_match(sequence_single)

                """储存结果
                满足存在引物, seq_extract_mh中无缺失（D）、seq_extract_mh非两头都是N，则认为是可使用序列
                """
                if primer != "None" and "D" not in seq_extract_mh:

                    if seq_extract_mh.startswith("N") and seq_extract_mh.endswith("N"):
                        continue

                    else:
                        self._save_match(dict_match_results, mh_id, primer, umi, seq_extract_mh)  # 储存match结果
                        self._save_for_purity(dict_for_purity, mh_id, umi, seq_extract_mh)  # 储存purity计算和过滤相关结果
                        self._save_mismatch(ls_mismatch, mh_id, primer, umi, seq_extract_mh, seq_extract_mh_z)  # 储存mismatch结果

                else:
                    continue

                """储存结果（primer == "None"）"""
                self._save_no_primer(ls_no_primer, primer, sequence_single)  # 储存no_primer结果

            """过滤"""
            if filter:  # 如果要进行过滤
                ls_filter_match_results = self.match_filter(dict_match_results, dict_for_purity,
                                                            purity_filter, umi_count_filter,
                                                            allele_num_filter, allele_proportion_filter)
                Toolsbox.FileTools.save_file(Toolsbox.FormatTools.str_ls(ls_filter_match_results),
                                             sam_match_filepath+"_Filter_"+filename)
            else:  # 如果不进行过滤
                Toolsbox.FileTools.save_file_match(dict_match_results, sam_match_filepath+"_No_Filter_"+filename)

            """储存文件(mismatch和no_primer)"""
            Toolsbox.FileTools.save_file(ls_mismatch, sam_mismatch_filepath+filename)  # 储存错误匹配文件
            Toolsbox.FileTools.save_file(ls_no_primer, sam_mismatch_filepath+"_No_Primer_"+filename)  # 储存引物为None的文件
            print("{}已匹配完成".format(filename))