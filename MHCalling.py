# -*- coding: utf-8 -*-
# @Time    : 2022
# @Author  : WangHaoyu
# @E-mail  : wanghy0707@gmail.com
# @Github  :
# @desc    :


import Toolsbox


class Analysis_call:
    """
    将样本进行MH calling得到的含N等位基因转化为无N等位基因，并同步修改MH的ID
    """
    @classmethod
    def _get_ref_info_snp(cls, info_ref_filepath):
        """
        获取INFO参考文件的ID和SNP组成信息
        eg. [[ID1, snpstart_snp_end(str)],
            [ID2, snp_start_snp_end(str)], ...]
        """
        info_ref = Toolsbox.FileTools.open_ls_file(info_ref_filepath)

        ls_ref_info_snp = []
        for line in info_ref:  # 对每一行进行操作
            id = line[0]  # ID
            snpstart = line[11]  # 起始SNP
            snpend = line[11+(int(line[8])-1)*5]  # 终止SNP
            ls_ref_info_snp.append([id, snpstart+"_"+snpend])
        return ls_ref_info_snp

    def __init__(self, info_ref_filepath):
        """
        info_ref_filepath: 进行SingleCalling或MixtureCalling时用于对含N MH进行降级的INFO参考文件路径
        """
        self.info_ref = Toolsbox.FileTools.open_ls_file(info_ref_filepath)
        self.ls_info_ref_snp = self._get_ref_info_snp(info_ref_filepath)  # 直接获取INFO参考文件的ID和SNP组成信息

    def _extract_mh_calling_snp_ls(self, chr, id, dict_chr_info):
        """
        提取单一样本MH Calling结果中每一个MH的SNP列表
        """
        for values in dict_chr_info[chr].values():

            if id != values[0]:
                pass
            elif id == values[0]:  # 如果匹配成功则返回SNP列表
                return values[2:]

    def _extract_mh_calling(self, mh_geno_origin, dict_chr_info):
        """
        提取样本的单一MH Calling结果中的信息：chr、ID、allele
        单一MH Calling结果实例：[ID, allele1, count1, allele2, count2]
        """
        id = mh_geno_origin[0]  # ID号
        chr = str(int(id[6:8]))  # 染色体号
        ls_snppos = self._extract_mh_calling_snp_ls(chr, id, dict_chr_info)  # MH的组成SNP列表
        ls_allele_tmp = mh_geno_origin[1:][::2]  # 等位基因列表
        ls_allele = [allele for allele in ls_allele_tmp if "puzzle" not in allele and "Down" not in allele]
        ls_count = mh_geno_origin[1:][1::2]  # 每个等位基因对应的count列表
        return ls_snppos, ls_allele, ls_count

    def _get_new_allele_snp(self, ls_snppos, ls_allele):
        """
        对等位基因中的N进行判断，同时获取新MH的SNP起始和终止信息（posstart_pos_end(str)）、新等位基因列表
        """
        ls_new_snppos = []
        ls_new_allele = [allele.replace("N", "") for allele in ls_allele]  # 将N去掉，就是新MH

        """获取新MH的SNP pos信息"""
        for snppos, base in zip(ls_snppos, ls_allele[0]):
            if base == "N":
                continue
            else:
                ls_new_snppos.append(str(snppos))  # 提取新pos时，直接将其转化为字符串

        new_snp = ls_new_snppos[0] + "_" + ls_new_snppos[-1]
        return new_snp, ls_new_allele

    def _get_new_id(self, new_snp):
        """
        依据新的SNP信息，获取新MH的对应ID
        """
        for line in self.ls_info_ref_snp:
            if new_snp == line[1]:
                return line[0]
            else:
                continue

    def _get_new_geno(self, new_id, ls_new_allele, ls_count):
        """
        获取新MH的分型信息
        """
        mh_geno = []
        for allele, count in zip(ls_new_allele, ls_count):
            mh_geno.append(allele)
            mh_geno.append(count)

        if new_id:  # 如果新id存在
            mh_geno = [new_id] + mh_geno
        else:  # 如果不存在新id（new_snp为单一SNP）
            mh_geno = ["puzzle_results_!_please_check_origin_results"]
        return mh_geno

    def _get_new_info(self, new_id):
        """
        获取新的MH对应的INFO信息（完整信息），输出为列表
        """
        if new_id:  # 如果存在新id
            for line in self.info_ref:
                if line[0] == new_id:
                    return line
                else:
                    continue

        else:  # 如果不存在新id
            return ["None"]

    def analysis_calling(self, mh_geno_origin, dict_chr_info):
        """
        将样本进行MH calling得到的含N等位基因转化为无N等位基因，并同步修改MH的ID的主函数，输出两行：
        1.新MH的分型信息
        2.新MH的INFO信息
        //2022.03.04 新增对原始结果的判断，只有当非no_match或puzzle时才会进行下一步
        """
        if len(mh_geno_origin) == 2:
            return ["no_match/puzzle"], ["None"]
        else:
            ls_snppos, ls_allele, ls_count = self._extract_mh_calling(mh_geno_origin, dict_chr_info)
            new_snp, ls_new_allele = self._get_new_allele_snp(ls_snppos, ls_allele)  # 获取新MH的SNP起始终止信息、新等位基因列表
            new_id = self._get_new_id(new_snp)  # 可能为None

            """输出文件"""
            mh_geno_new = self._get_new_geno(new_id, ls_new_allele, ls_count)
            mh_info_new = self._get_new_info(new_id)
            return mh_geno_new, mh_info_new


class SingleCall:
    """
    对单个体样本进行MH calling
    """

    def __init__(self, match_filepath, info_ref_filepath, dict_chr_info, ls_id, dict_chr_relalt):
        """
        match_filepath: 由Match.match得到的Match文件储存的路径
        info_ref_filepath: 进行SingleCalling或MixtureCalling时用于对含N MH进行降级的INFO参考文件路径
        dict_chr_info: 由PreInfo得到的结果（嵌套字典，详细格式见PreProcessing_PreInfo）
        ls_id: 储存MH ID的列表
        dict_chr_relalt: 由PreVcf得到的结果（嵌套字典，详细格式见PreProcessing_PreVcf）
        """
        self.match_filepath = match_filepath
        self.info_ref_filepath = info_ref_filepath
        self.dict_chr_info = dict_chr_info
        self.ls_id = ls_id
        self.dict_chr_relalt = dict_chr_relalt

    """
    提取特定MH的列表
    """
    def _extract_mh(self, match, id):
        """
        按照match文件中的原始顺序提取单个MH的内容
        //2021.12.17 暂且保留使用列表格式，防止需要中间输出tmp文件
        //2022.02.09 未考虑到未匹配到MH match的情况，此时输出空列表
        """
        ls_spe_mh = []
        for line in match:

            if line[0] == id:
                ls_spe_mh.append(line[:4]+[int(line[4])])  # 保留所有信息，将count转换为int类型
            else:
                continue

        if ls_spe_mh:
            len_mh_allele = len(ls_spe_mh[0][2])  # 获取当前操作MH的等位基因包含的SNP数
        else:
            len_mh_allele = 0
        return ls_spe_mh, len_mh_allele

    """
    提取该MH中各SNP的参考分型列表
    """
    def _ls_spe_mh_info(self, id, mh_chr):
        """
        提取特定MH的INFO信息
            ['CHBCHR01_0674112', 12, 247032240, 247032247, 247032249, ...]
        """
        for ls_mh_info in self.dict_chr_info[mh_chr].values():  # 对于每一条INFO信息进行操作
            if ls_mh_info[0] == id:
                return ls_mh_info
            else:
                continue

    def _ls_spe_mh_refalt(self, id):
        """
        创建特定MH中各SNP的参考分型列表（一维列表）
            ["TA", "GC", "AC", ...]
        """
        mh_chr = str(int(id[6:8]))  # 获取当前特定MH的染色体号
        ls_spe_mh_info = self._ls_spe_mh_info(id, mh_chr)  # 提取特定MH的INFO信息列表

        dict_spe_mh_relalt = self.dict_chr_relalt[mh_chr]
        ls_spe_mh_refalt = []
        for snp in ls_spe_mh_info[2:]:  # 按顺序提取特定MH的SNP参考分型列表
            ls_spe_mh_refalt.append(dict_spe_mh_relalt[snp])
        return ls_spe_mh_refalt

    """
    依据SNP参考分型对结果进行过滤
    """
    def _dict_base_count(self, len_mh_allele):
        """
        依据组成MH的SNP数目创建用于SNP计数的字典列表，每一个子字典对应一个SNP
        """
        return [{"A": 0, "T": 0, "C": 0, "G": 0} for i in range(len_mh_allele)]

    def _snp_count(self, ls_spe_mh, len_mh_allele):
        """
        对特定MH中所有等位基因的每个组成的SNP碱基类型进行计数，得到特定MH所有有效reads中的SNP计数结果，每一个子字典对应一个SNP
            [{"A": 200, "T": 150, "C": 5, "G": 0},  # the first snp
            {"A": 30, "T": 1, "C": 400, "G": 6}, ...]  # the second, ... snp
        //2022.02.09 当组成SNP的碱基为N时，需进行跳过
        //2022.05.05 当组成SNP的碱基为D时，需进行跳过
        """
        ls_dict_snp_count = self._dict_base_count(len_mh_allele)  # 创建碱基计数字典
        for line in ls_spe_mh:  # 对每一个计数结果进行操作

            for snp, dict_snp_count in zip(line[2], ls_dict_snp_count):  # 对单个SNP进行计数

                if snp == "N" or snp == "D":
                    continue
                else:
                    dict_snp_count[snp] += line[-1]
        return ls_dict_snp_count

    def _dict_base_count_filter(self, ls_spe_mh, len_mh_allele, id):
        """
        【已取消使用】
        基于特定MH中各SNP中的参考分型列表，对SNP碱基计数字典进行过滤，仅保留符合参考分型的结果
        //2022.02.09 字典中去除键值对后对for循环会产生影响，更改了去除键值对的方式
        //2022.02.14 输出字典更正为ls_dict_snp_count_pop
        //2022.03.01 依据参考分型列表对SNP碱基计数字典进行过滤会导致将样本中真正的SNP突变过滤掉从而得到错误的纯/杂合结果，因此这
                     一步中不再进行而直接输出ls_dict_snp_count
        """
        ls_spe_mh_refalt = self._ls_spe_mh_refalt(id)
        ls_dict_snp_count = self._snp_count(ls_spe_mh, len_mh_allele)

        ls_dict_snp_count_pop = []
        for spe_mh_refalt, dict_snp_count in zip(ls_spe_mh_refalt, ls_dict_snp_count):

            ls_popkey = []
            for key in dict_snp_count.keys():

                if key not in spe_mh_refalt:
                    ls_popkey.append(key)
                else:
                    pass

            for popkey in ls_popkey:
                dict_snp_count.pop(popkey)

            ls_dict_snp_count_pop.append(dict_snp_count)
        return ls_dict_snp_count_pop

    """
    基于过滤后的SNP字典，选取top2/top1作为过滤比对碱基，并对match文件进行过滤
    //2022.02.14 去除了mh_filter中基于过滤对比碱基进一步对ls_spe_mh进行过滤的操作，但mh_filter函数仍暂时保留
    """
    def _extract_top(self, ls_dict_snp_count, peak_height_ratio):
        """
        从过滤后的SNP计数字典中，提取top1/top2的碱基
        list_dict_snp_count: 对MH每一个SNP进行ATCG碱基记数的字典，该字典已经过过滤，仅保留符合SNP参考分型的记数结果
            {"A": 0, "T": 454}
        peak_height_ratio: 用于判断该SNP上等位基因是纯合子还是杂合子的比例，小于该比例为纯合子，大于该比例则为杂合子
        //2022.02.10 匹配的MH match结果在某一SNP上全为N时，SNP计数字典中各项为0，无法被除，需跳过此类情况
        //2022.02.14 由于取消了mh_filter功能，将此处ls_retain_base相关的代码全部取消（包括不再输出结果），仅保留纯杂合判断
        """
        ls_he_ho = []
        for dict_snp_count in ls_dict_snp_count:  # 按顺序对特定MH中的每一个SNP进行操作

            sorted_dict_key = sorted(dict_snp_count, key=dict_snp_count.get, reverse=True)  # 按照value由大到小对字典的key进行排序

            if dict_snp_count[sorted_dict_key[0]] == 0:  # 如果在某一SNP位置上全是N时
                ls_he_ho.append("ho")
            else:
                if dict_snp_count[sorted_dict_key[1]] / dict_snp_count[sorted_dict_key[0]] < peak_height_ratio:  # 判定纯合子
                    ls_he_ho.append("ho")
                else:  # 判定杂合子
                    ls_he_ho.append("he")
        return ls_he_ho

    def mh_filter(self, ls_spe_mh, ls_dict_snp_count, peak_height_ratio):
        """
        【已取消使用】
        基于过滤比对碱基，对ls_spe_mh进行过滤，去除碱基不包含在过滤对比碱基中的序列
        //2022.02.14 取消了mh_filter功能，此函数不再使用，但仍保留于此
        """
        ls_retain_base, ls_he_ho = self._extract_top(ls_dict_snp_count, peak_height_ratio)  # 创建保留SNP碱基列表

        ls_retain_spe_mh = []  # 保存过滤后的match结果文件
        for line in ls_spe_mh:

            count = 0  # 对每一行重置count
            for base_snp, ratain_base in zip(line[2], ls_retain_base):
                if base_snp != "N":

                    if base_snp not in ratain_base:  # 出现一次SNP的不匹配情况，count+1
                        count += 1
                    else:
                        continue

                else:
                    continue

            if count == 0:
                ls_retain_spe_mh.append(line)
            else:
                pass
        return ls_retain_spe_mh, ls_he_ho

    """
    基于特定MH经过SNP过滤后的信息，提取等位基因分型
    """
    def _divide_spe_mh(self, ls_spe_snp_filter_mh):
        """
        将特定MH的信息按等位基因是否完整进行分割，以字典形式储存：allele-count
        //2022.02.11 原字典储存只储存了最后一项，所有字典均为单项字典。未合计相同等位基因的count数，先进行了修改
        //2022.05.04 补充判断条件，N存在或D存在时均保留到不完整列表
        """
        dict_integrate = {}
        dict_n = {}
        for line in ls_spe_snp_filter_mh:

            if "N" not in line[2]:  # 判断为完整则将allele-count储存

                if line[2] not in dict_integrate.keys():  # 如果是第一次储存
                    dict_integrate[line[2]] = line[-1]
                else:  # 如果不是第一次储存
                    dict_integrate[line[2]] = dict_integrate[line[2]] + line[-1]

            elif "N" in line[2] or "D" in line[2]:  # 判断为不完整时则保留至不完整列表

                if line[2] not in dict_n.keys():  # 如果是第一次储存
                    dict_n[line[2]] = line[-1]
                else:  # 如果不是第一次储存
                    dict_n[line[2]] = dict_n[line[2]] + line[-1]
        return dict_integrate, dict_n

    def _dicide_spe_mh_n(self, dict_n):
        """
        将特定MH的不完整等位基因按照N在前还是N在后区分为两个子字典，同时判断使用哪一个字典
        """
        dict_n_f = {key: value for key, value in dict_n.items() if key.startswith("N")}  # N在前
        dict_n_b = {key: value for key, value in dict_n.items() if key.endswith("N")}  # N在后

        """2022.10.12 新增：暂时按照f字典或b字典直接进行输出，不再进行其他判断"""
        judge = "f"  # 若为f则输出f，若为b则输出b
        if judge == "f":
            return dict_n_f
        elif judge == "b":
            return dict_n_b

        if not dict_n_f:  # 如果一个字典为空，则直接返回另一个字典
            return dict_n_b
        elif not dict_n_b:
            return dict_n_f
        else:  # 如果两个字典均不为空，则首先判断含N量及其对应的等位基因数，若相同时，则进一步判断count数

            def _n_count(dict_n):
                """
                判断dict_n（分为N在前和N在后）中不同含N量下包含的等位基因类别数，并依据含N量输出top2的等位基因类别数及所有等位基因
                的count总数
                    [[最小含N量1, 具最小含N量等位基因1, 具最小含N量等位基因count1],
                    [最小含N量2, 具最小含N量等位基因2, 具最小含N量等位基因count2]]
                //2022.02.10 未考虑含N量相同的情况，以及字典中只有一个等位基因的情况
                    [最小含N量, 具最小含N量等位基因, 具最小含N量等位基因count]
                //2022.03.04 简化了含N量情况等于1时的判断
                """
                dict_n_type_count = {}
                for key in dict_n.keys():
                    n_count = key.count("N")

                    if n_count not in dict_n_type_count.keys():
                        dict_n_type_count[n_count] = [1, dict_n[key]]
                    else:
                        dict_n_type_count[n_count] = [dict_n_type_count[n_count][0] + 1,
                                                      dict_n_type_count[n_count][1] + dict_n[key]]

                sorted_dict_key = sorted(dict_n_type_count.keys())  # 将key按照key（含N量）由小到大排序

                if len(sorted_dict_key) > 1:  # 如果含N量情况大于1，则正常输出
                    ls_top_2 = [[sorted_dict_key[0]] + dict_n_type_count[sorted_dict_key[0]],
                                [sorted_dict_key[1]] + dict_n_type_count[sorted_dict_key[1]]]  # 等位基因含N量及该含N量下等位基因数

                elif len(sorted_dict_key) == 1:  # 如果含N量情况等于1，有多个/1个含N量相同的等位基因
                    ls_top_2 = [[sorted_dict_key[0]] + dict_n_type_count[sorted_dict_key[0]]]
                return ls_top_2

            def _f_b(ls_top_2_n_f, ls_top_2_n_b):
                """
                判断输出的是f字典还是b字典
                先比较f字典和b字典的top1等位基因类别数，用于决定f字典和b字典的比较是使用top1等位基因还是top2等位基因
                之后比较等位基因的含N量，输出少的一方
                若等位基因的含N量相同，则比较含N量相同的所有等位基因的累积count数，输出多的一方
                如果count数相同，则使用f字典
                //2022.03.04 新增当ls_top_2_f和ls_top_2_b均只包含一项等位基因时的判断
                """
                if ls_top_2_n_f[0][1] >= 2 and ls_top_2_n_b[0][1] >= 2:  # 如果两个字典的top1等位基因类别数均大于2，则直接比较top1

                    if ls_top_2_n_f[0][0] < ls_top_2_n_b[0][0]:  # 比较两个字典top1的含N量，输出少的一方
                        return "f"
                    elif ls_top_2_n_f[0][0] > ls_top_2_n_b[0][0]:
                        return "b"
                    elif ls_top_2_n_f[0][0] == ls_top_2_n_b[0][0]:  # 如果二者相同，则进一步判断count数

                        if ls_top_2_n_f[0][2] < ls_top_2_n_b[0][2]:
                            return "b"
                        elif ls_top_2_n_f[0][2] > ls_top_2_n_b[0][2]:
                            return "f"
                        elif ls_top_2_n_f[0][2] == ls_top_2_n_b[0][2]:  # 如果count相同，则输出f
                            return "f"

                elif ls_top_2_n_f[0][1] >= 2 and ls_top_2_n_b[0][1] == 1:  # 如果一个字典的top1等位基因类别数为1，则比较第二位

                    if len(ls_top_2_n_b) == 1:  # 如果top1等位基因类别数为1的字典只包含一项
                        return "f"

                    else:  # 如果top1等位基因类别数为1的字典只包含多项
                        if ls_top_2_n_f[0][0] < ls_top_2_n_b[1][0]:  # 比较两个字典top1和top2的含N量，输出少的一方
                            return "f"
                        elif ls_top_2_n_f[0][0] > ls_top_2_n_b[1][0]:
                            return "b"
                        elif ls_top_2_n_f[0][0] == ls_top_2_n_b[1][0]:  # 如果二者相同，则进一步判断count数

                            if ls_top_2_n_f[0][2] < ls_top_2_n_b[1][2]+ls_top_2_n_b[0][2]:
                                return "b"
                            elif ls_top_2_n_f[0][2] > ls_top_2_n_b[1][2]+ls_top_2_n_b[0][2]:
                                return "f"
                            elif ls_top_2_n_f[0][2] == ls_top_2_n_b[1][2]+ls_top_2_n_b[0][2]:  # 如果count相同，则输出f
                                return "f"

                elif ls_top_2_n_f[0][1] == 1 and ls_top_2_n_b[0][1] >= 2:

                    if len(ls_top_2_n_f) == 1:  # 如果top1等位基因类别数为1的字典只包含一项
                        return "b"

                    else:  # 如果top1等位基因类别数为1的字典只包含多项
                        if ls_top_2_n_f[1][0] < ls_top_2_n_b[0][0]:
                            return "f"
                        elif ls_top_2_n_f[1][0] > ls_top_2_n_b[0][0]:
                            return "b"
                        elif ls_top_2_n_f[1][0] == ls_top_2_n_b[0][0]:

                            if ls_top_2_n_f[1][2]+ls_top_2_n_f[0][2] < ls_top_2_n_b[0][2]:
                                return "b"
                            elif ls_top_2_n_f[1][2]+ls_top_2_n_f[0][2] > ls_top_2_n_b[0][2]:
                                return "f"
                            elif ls_top_2_n_f[1][2]+ls_top_2_n_f[0][2] == ls_top_2_n_b[0][2]:
                                return "f"

                elif ls_top_2_n_f[0][1] == 1 and ls_top_2_n_b[0][1] == 1:

                    if len(ls_top_2_n_f) == 1 and len(ls_top_2_n_b) > 1:  # 如果F类型只有1个，且B类型大于1个
                        return "b"

                    elif len(ls_top_2_n_f) > 1 and len(ls_top_2_n_b) == 1:  # 如果B类型只有1个，且F类型大于1个
                        return "f"

                    elif len(ls_top_2_n_f) == 1 and len(ls_top_2_n_b) == 1:
                        if ls_top_2_n_f[0][2] >= ls_top_2_n_b[0][2]:
                            return "f"
                        elif ls_top_2_n_f[0][2] < ls_top_2_n_b[0][2]:
                            return "b"

                    else:
                        if ls_top_2_n_f[1][0] < ls_top_2_n_b[1][0]:
                            return "f"
                        elif ls_top_2_n_f[1][0] > ls_top_2_n_b[1][0]:
                            return "b"
                        elif ls_top_2_n_f[1][0] == ls_top_2_n_b[1][0]:

                            if ls_top_2_n_f[1][2]+ls_top_2_n_f[0][2] < ls_top_2_n_b[1][2]+ls_top_2_n_b[0][2]:
                                return "b"
                            elif ls_top_2_n_f[1][2]+ls_top_2_n_f[0][2] > ls_top_2_n_b[1][2]+ls_top_2_n_b[0][2]:
                                return "f"
                            elif ls_top_2_n_f[1][2]+ls_top_2_n_f[0][2] == ls_top_2_n_b[1][2]+ls_top_2_n_b[0][2]:
                                return "f"

            ls_top_2_n_f = _n_count(dict_n_f)
            ls_top_2_n_b = _n_count(dict_n_b)

            choice = _f_b(ls_top_2_n_f, ls_top_2_n_b)  # 判断所使用的字典类型
            if choice == "f":
                return dict_n_f
            elif choice == "b":
                return dict_n_b

    def _he_ho(self, ls_he_ho):
        """
        依据ls_he_ho判断特定MH是纯合子还是杂合子
        """
        if "he" in ls_he_ho:
            return "he"
        else:
            return "ho"

    def _extract_geno_he_xn_xdown(self, dict_x, dict_n, n_ref, peak_height_ratio):
        """
        dict_x: 待降级的等位基因-count字典（可为完整等位基因字典，可为不完整等位基因字典）
        dict_n: 含N等位基因-count字典（已经过选择的f_n字典或b_n字典），用于更新降级后等位基因的count
        n_ref: 等位基因降级所参考的不完整等位基因
        peak_height_ratio:
        //2022.02.21 变更降级思路：在降级一次后可能仍然无法满足条件，因此需要进入迭代。首先明确需要进行降级的必然为he类型，输出等
        //           位基因必须>=2。传入待降级的字典，其可能为一项，可能为多项，若为第n次（n>1）降级，其可能包含多项内容。将待降
        //           级字典中的每一个key均执行降级操作，并将降级后的字典update至dict_n中更新count数，同时要在字典中去掉原key。
        //           判断dict_n中含N量最少的项数，若项数仍小于2，则进行下一轮迭代，若大于2，则输出含N量最少的、count数最多的两
        //           项。若count数相同，则全输出，并标明puzzle。
        """
        """分离待降级等位基因及count"""
        ls_x_key = list(dict_x)
        ls_x_value = list(dict_x.values())  # 将待降级字典转化为列表，两个列表分别包含key和value

        """将待降级等位基因降级"""
        ls_x_down_key = []  # 储存降级后的完整等位基因，及其修正后的count数（二维列表）
        for x_key in ls_x_key:  # 对每一个待降级的等位基因均依照参考不完整等位基因进行降级

            x_down = ""
            for base_n, base_x in zip(n_ref, x_key):

                if base_n != "N":
                    x_down = x_down + base_x
                elif base_n == "N":
                    x_down = x_down + "N"

            ls_x_down_key.append(x_down)

        """判断降级后的等位基因是否包含在不完整等位基因集合中，依据判断结果对原dict_n的count进行调整，同时删除旧key"""
        for x_down_key, x_value in zip(ls_x_down_key, ls_x_value):
            if x_down_key in dict_n.keys():  # 如果在不完整等位基因集合中存在
                dict_n[x_down_key] = dict_n[x_down_key] + x_value
            else:  # 如果在不完整等位基因集合中不存在
                dict_n[x_down_key] = x_value

        for popkey in ls_x_key:  # 当降级前的旧key存在于dict_n时（若降级前的旧key为完整等位基因，则不存在），删除旧key
            if popkey in dict_n.keys():
                dict_n.pop(popkey)
            else:
                continue

        """对降级后的dict_n进行含N量判断，看是否能输出结果，若不能，则继续进行降级循环"""
        dict_n_count = {}  # 储存不完整等位基因含N计数，key:value=等位基因: 含N量
        for key in dict_n.keys():
            dict_n_count[key] = key.count("N")
        sorted_dict_key_n_count = sorted(dict_n_count, key=dict_n_count.get, reverse=False)  # 按含N量由小到大对等位基因排序
        min_n = dict_n_count[sorted_dict_key_n_count[0]]  # 获取含N最少等位基因的含N量
        ls_min_n_allele = [key for key, value in dict_n_count.items() if value == min_n]  # 获取含N最少的等位基因列表
        dict_min_n = {key: dict_n[key] for key in ls_min_n_allele}  # 获取含N最少的dict_n字典，同时可作为下轮迭代的key

        """判断继续迭代的条件：最大含N量不等于最小含N量"""
        if dict_n_count[sorted_dict_key_n_count[0]] != dict_n_count[sorted_dict_key_n_count[-1]]:
            judge = "yes"
        else:
            judge = "no"

        if len(ls_min_n_allele) > 2:  # 若含N量最少的等位基因数超过2个，可经过判断后输出结果
            sorted_dict_key_count = sorted(dict_min_n, key=dict_min_n.get, reverse=True)  # 按count由大到小对等位基因排序

            """自第3位起，判断count是否与第2位相同，若相同，则输出相同key"""
            ls_same_to_2 = []  # 储存与第2位count数相同的等位基因的key
            for key in sorted_dict_key_count[2:]:  # 自第3位开始比较

                if dict_n[key] == dict_n[sorted_dict_key_count[1]]:  # 若与第2位相同，则储存
                    ls_same_to_2.append(key)
                else:
                    pass

            """依据上述比较结果，判断是否能进行正常杂合子等位基因输出"""
            ls_result = [sorted_dict_key_count[0], dict_n[sorted_dict_key_count[0]],
                         sorted_dict_key_count[1], dict_n[sorted_dict_key_count[1]]]  # 输出前两位
            if ls_same_to_2:  # 如果有多位与第2位的count数相同，则全部输出并标记
                for key in ls_same_to_2:
                    ls_result.extend([key, dict_n[key]])
                ls_result.append("Down_&_puzzle_he_in_nn_multi=2")

            else:  # 如果没有等位基因与第2位的count数相同，说明可正常输出杂合子的两个等位基因
                if ls_result[3]/ls_result[1] < peak_height_ratio or ls_result[3]/ls_result[1] > 1/peak_height_ratio:
                    ls_result.append("Down_&_puzzle_he_in_nn_phr<0.3")
                else:
                    ls_result.append("Down")
                    pass

            return ls_result

        elif len(ls_min_n_allele) == 2:  # 若含N量最少的等位基因数为2个，可直接输出
            sorted_dict_key_count = sorted(dict_min_n, key=dict_n_count.get, reverse=True)  # 按count由大到小对等位基因排序
            ls_result = [sorted_dict_key_count[0], dict_n[sorted_dict_key_count[0]],
                         sorted_dict_key_count[1], dict_n[sorted_dict_key_count[1]]]  # 输出两位

            if ls_result[3]/ls_result[1] < peak_height_ratio or ls_result[3]/ls_result[1] > 1/peak_height_ratio:
                ls_result.append("Down_&_puzzle_he_in_nn_phr<0.3")
            else:
                ls_result.append("Down")
            return ls_result

        elif len(ls_min_n_allele) == 1 and judge == "yes":  # 当含N量最少的等位基因数仍为1，且满足继续迭代条件时，继续迭代
            n_ref = sorted_dict_key_n_count[1]  # 含N量最少的等位基因序号为0，则含N量第二少的等位基因序号为1
            return self._extract_geno_he_xn_xdown(dict_min_n, dict_n, n_ref, peak_height_ratio)

        elif judge == "no":  # 当不满足继续迭代条件时，终止并输出
            return ["puzzle_he_in_nn_no_results_for_x_down"]

    def _extract_geno_he_ii(self, dict_integrate, peak_height_ratio):
        """
        从>2个完整等位基因中提取杂合子基因型：取多
        如果有多等位基因count与第2位相同，则全部输出并标记为"puzzle_he_ii_multi=2"
        对于所取的两个等位基因，还需要判断二者reads的比例关系，如果二者的比例低于0.3，则在输出时标记为"puzzle_he_ii_phr<0.3"
        """
        sorted_dict_key = sorted(dict_integrate, key=dict_integrate.get, reverse=True)  # 将key按value由大到小顺序

        """自第3位起，判断count是否与第2位相同，若相同，则输出相同key"""
        ls_same_to_2 = []  # 储存与第2位count数相同的等位基因的key
        for key in sorted_dict_key[2:]:  # 自第3位开始比较

            if dict_integrate[key] == dict_integrate[sorted_dict_key[1]]:  # 若与第2位相同，则储存
                ls_same_to_2.append(key)
            else:
                pass

        """依据上述比较结果，判断是否能进行正常杂合子等位基因输出"""
        ls_result = [sorted_dict_key[0], dict_integrate[sorted_dict_key[0]],
                     sorted_dict_key[1], dict_integrate[sorted_dict_key[1]]]  # 输出前两位（按照count大小顺序）
        if ls_same_to_2:  # 如果有多位与第2位的count数相同，则全部输出并标记
            for key in ls_same_to_2:
                ls_result.extend([key, dict_integrate[key]])
            ls_result.append("puzzle_he_ii_multi=2")

        else:  # 如果没有等位基因与第2位的count数相同，说明可正常输出杂合子的两个等位基因
            if ls_result[3]/ls_result[1] < peak_height_ratio:  # 如果比例小于0.3，则可能与判断的杂合子结果有出入，需标记puzzle
                ls_result.append("puzzle_he_ii_phr<0.3")
            else:
                pass
        return ls_result

    def _extract_geno_he_nn(self, dict_n, peak_height_ratio):
        """
        从不完整等位基因中提取杂合子基因型（most hard）
        注意：在nn类型中在选择或降级前还需要再考虑所选择的含N等位基因的N是在前还是在后，基于所选择的子集内进行选择。所以此处应该是不准的
        """
        """判断所使用的字典类型为f_N还是b_N"""
        dict_n = self._dicide_spe_mh_n(dict_n)

        """
        将不完整等位基因按含N数量，由小到大排序
        //2022.02.10 此处同样需要考虑含N量相同的情况，以及字典中只有一个等位基因的情况
        """
        dict_n_count = {}  # 储存不完整等位基因含N计数，key:value=等位基因: 含N量
        for key in dict_n.keys():
            dict_n_count[key] = key.count("N")
        sorted_dict_value = sorted(set(dict_n_count.values()), reverse=False)  # 将N计数字典value按value（含N量）由小到大排序
        min_n = sorted_dict_value[0]  # 获取含N最少等位基因的含N量
        ls_min_n_allele = [key for key, value in dict_n_count.items() if value == min_n]  # 获取含N最少的等位基因列表

        """
        判断具有最小含N量的含N等位基因数目
        =1：对具最小含N量的等位基因进行降级
        =2：直接输出，无需降级
        >2：检测每个等位基因的count，选择count最大的两个输出，无需降级
        对于所取的两个等位基因，还需要判断二者reads的比例关系，如果二者的比例低于0.3，则在输出时标记为"puzzle_he_nn_phr<0.3"
        """
        if len(ls_min_n_allele) == 1:  # =1时（只有具最小含N量的等位基因数为1时，才需要执行降级操作）
            dict_min_n = {ls_min_n_allele[0]: dict_n[ls_min_n_allele[0]]}  # 提取含N最少的等位基因-count子字典，dict_x

            """判断当前等位基因是否有能力执行降级操作：含N类型超过1种"""
            if len(sorted_dict_value) > 1:  # 如果N计数字典含N类型超过一种，有能力执行降级操作
                min2_n = sorted_dict_value[1]
                ls_min2_n_allele = [key for key, value in dict_n_count.items() if value == min2_n]  # 获取含N次少的等位基因列表，ref_n
                return self._extract_geno_he_xn_xdown(dict_min_n, dict_n, ls_min2_n_allele[0], peak_height_ratio)

            elif len(sorted_dict_value) == 1:  # 如果N计数字典含N类型仅有一种，无能力执行降级操作（这个判断好像是不必要的？）
                return ["puzzle_he_nn_n_count_type=1"]

        elif len(ls_min_n_allele) == 2:  # =2时直接输出
            ls_result = [ls_min_n_allele[0], dict_n[ls_min_n_allele[0]],
                         ls_min_n_allele[1], dict_n[ls_min_n_allele[1]]]

            if ls_result[3]/ls_result[1] < peak_height_ratio or ls_result[3]/ls_result[1] > 1/peak_height_ratio:
                ls_result.append("puzzle_he_nn_phr<0.3")
            else:
                pass
            return ls_result

        else:  # >2时需进行进一步判断
            dict_min_n = {key: dict_n[key] for key in ls_min_n_allele}  # 提取含N最少的等位基因-count子字典
            sorted_dict_min_n_key = sorted(dict_min_n, key=dict_min_n.get, reverse=True)  # 将key按value由大到小排序

            """自第3位起，判断count是否与第2位相同，若相同，则输出相同key"""
            ls_same_to_2 = []  # 储存与第2位count数相同的等位基因的key
            for key in sorted_dict_min_n_key[2:]:  # 自第3位开始比较

                if dict_n[key] == dict_n[sorted_dict_min_n_key[1]]:  # 若与第2位相同，则储存
                    ls_same_to_2.append(key)
                else:
                    pass

            """依据上述比较结果，判断是否能进行正常杂合子等位基因输出"""
            ls_result = [sorted_dict_min_n_key[0], dict_n[sorted_dict_min_n_key[0]],
                         sorted_dict_min_n_key[1], dict_n[sorted_dict_min_n_key[1]]]  # 输出前两位
            if ls_same_to_2:  # 如果有多位与第2位的count数相同，则全部输出并标记
                for key in ls_same_to_2:
                    ls_result.extend([key, dict_n[key]])
                ls_result.append("puzzle_he_nn_multi=2")

            else:  # 如果没有等位基因与第2位的count数相同，说明可正常输出杂合子的两个等位基因
                if ls_result[3]/ls_result[1] < peak_height_ratio or ls_result[3]/ls_result[1] > 1/peak_height_ratio:
                    ls_result.append("puzzle_he_nn_phr<0.3")
                else:
                    pass
            return ls_result

    def _extract_geno_he_in(self, dict_integrate, dict_n, peak_height_ratio):
        """
        从完整和一个不完整等位基因中提取在杂合子基因型（依据不完整等位基因，需对完整等位基因进行降级处理）
        """
        """判断所使用的字典类型为f_N还是b_N"""
        dict_n = self._dicide_spe_mh_n(dict_n)

        """将不完整等位基因按含N数量，由小到大排序，获取ref_n"""
        dict_n_count = {}  # 储存不完整等位基因含N计数，key:value=等位基因: 含N量
        for key in dict_n.keys():
            dict_n_count[key] = key.count("N")
        sorted_dict_key = sorted(dict_n_count, key=dict_n_count.get, reverse=False)  # 按含N量由小到大对等位基因排序
        ref_n = sorted_dict_key[0]

        """输出降级结果"""
        return self._extract_geno_he_xn_xdown(dict_integrate, dict_n, ref_n, peak_height_ratio)

    def _extract_geno_he(self, dict_integrate, dict_n, peak_height_ratio):
        """
        提取杂合子基因型（以列表形式），以count多少顺序进行排列
        """
        if dict_integrate:  # 如果存在完整等位基因

            if len(dict_integrate) == 2:  # 只有2个时，直接输出
                sorted_dict_key = sorted(dict_integrate, key=dict_integrate.get, reverse=True)
                return [sorted_dict_key[0], dict_integrate[sorted_dict_key[0]],
                        sorted_dict_key[1], dict_integrate[sorted_dict_key[1]]]

            elif len(dict_integrate) > 2:  # 超过2个时，走ii类型
                return self._extract_geno_he_ii(dict_integrate, peak_height_ratio)

            elif len(dict_integrate) == 1:  # 只有1个时，走in类型

                if dict_n:  # 如果存在不完整等位基因，可进行降级calling
                    return self._extract_geno_he_in(dict_integrate, dict_n, peak_height_ratio)
                else:  # 如果不存在不完整等位基因，且杂合子又仅有1个完整等位基因，则输出puzzle
                    return ["puzzle_he_only_one_integrate_allele_no_n_allele"]

        else:  # 无完整等位基因，走nn类型

            if dict_n:  # 存在不完整等位基因
                return self._extract_geno_he_nn(dict_n, peak_height_ratio)
            else:
                return ["puzzle_he_no_integrate_allele_and_n_allele"]

    def _extract_geno_ho_n(self, dict_n):
        """
        从不完整等位基因中提取纯合子基因型：有长取长，同长取多
        如果存在两个等位基因的包含的N相同，且reads数相同，则均输出并标记为"puzzle_ho"
        """
        """将不完整等位基因按含N数量，由小到大排序"""
        dict_n_count = {}  # 储存不完整等位基因含N计数
        for key in dict_n.keys():
            dict_n_count[key] = key.count("N")
        sorted_dict_key = sorted(dict_n_count, key=dict_n_count.get, reverse=False)  # 将N计数字典key按value由小到大顺序
        min_n = dict_n_count[sorted_dict_key[0]]  # 获取含N最少等位基因的含N量
        ls_min_n_allele = [key for key, value in dict_n_count.items() if value == min_n]  # 获取含N最少的等位基因列表

        """
        判断是否有>1个等位基因具有相同的最小含N量
        若否，则返回含N最少的等位基因；
        若是，则返回具有相同长度等位基因中count最多的，若最多的count数相同，则附加输出标签"puzzle_ho_multi_count="
        """
        if len(ls_min_n_allele) == 1:  # 若只有一个等位基因具有最小N数，则直接返回等位基因和count
            return [ls_min_n_allele[0], dict_n[ls_min_n_allele[0]]]

        else:  # 如果有超过1个等位基因具有最小N数，则进一步判断他们的count，取最大者
            ls_min_n_allele_count = [dict_n[key] for key in ls_min_n_allele]  # 按顺序提取含N最少的等位基因对应的count数
            ls_min_n_allele_count_max_index = [i for i, count in enumerate(ls_min_n_allele_count)
                                               if count == max(ls_min_n_allele_count)]  # 获取count最大值的索引

            if len(ls_min_n_allele_count_max_index) == 1:  # 如果count最大值只有一个，输出该值对应的等位基因和count
                return [ls_min_n_allele[ls_min_n_allele_count_max_index[0]],
                        dict_n[ls_min_n_allele[ls_min_n_allele_count_max_index[0]]]]

            else:  # 如果count最大值有多个，则输出每一个等位基因、count，并标记"puzzle"
                ls_puzzle = []
                for index in ls_min_n_allele_count_max_index:
                    ls_puzzle.extend([ls_min_n_allele[index], dict_n[ls_min_n_allele[index]]])
                ls_puzzle.append("puzzle_ho_multi_max_count>=2")
                return ls_puzzle

    def _extract_geno_ho(self, dict_integrate, dict_n):
        """
        提取纯合子基因型（以列表形式）
        """
        if dict_integrate:  # 如果存在完整等位基因
            sorted_dict_key = sorted(dict_integrate, key=dict_integrate.get, reverse=True)
            return [sorted_dict_key[0], dict_integrate[sorted_dict_key[0]]]

        else:  # 如果不存在完整等位基因

            if dict_n:  # 存在不完整等位基因
                return self._extract_geno_ho_n(dict_n)
            else:
                return ["puzzle_ho_no_integrate_allele_and_n_allele"]

    def _extract_geno(self, ls_spe_snp_filter_mh, ls_he_ho, peak_height_ratio):
        """
        提取特定MH的基因型
        MH基因型表现方式：
        Ho: ID allele count
        He: ID allele1 count1 allele2 count2
        """
        dict_integrate, dict_n = self._divide_spe_mh(ls_spe_snp_filter_mh)  # 将特定MH按等位基因完整与否分割
        geno_type = self._he_ho(ls_he_ho)  # 判断特定MH是纯合子还是杂合子

        if geno_type == "ho":  # 提取纯合子基因型
            return self._extract_geno_ho(dict_integrate, dict_n)
        else:  # 提取杂合子基因型
            return self._extract_geno_he(dict_integrate, dict_n, peak_height_ratio)

    def _singlecall(self, filename, peak_height_ratio, save, save_path):
        """
        单核运行的主函数
        //2022.02.09 未考虑到未匹配到MH match的情况，此时输出空列表
        """
        match_file = Toolsbox.FileTools.open_ls_file(self.match_filepath + filename)
        a_call = Analysis_call(self.info_ref_filepath)  # 降级器！

        ls_geno_origin_results = []  # 储存单个样本文件的结果（原始结果）
        ls_geno_new_results = []  # 储存单个样本文件的结果（将原始MH进行降级后的结果）
        ls_info_new_results = []  # 储存降级后MH对应的INFO文件
        for id in self.ls_id:  # 对每一个MH进行操作
            ls_spe_mh, len_mh_allele = self._extract_mh(match_file, id)  # 提取特定的MH的match信息以及allele的长度

            if ls_spe_mh:  # 如果特定MH的match信息存在，则进行正常流程的过滤
                """
                依据SNP信息进行match过滤
                //2022.03.01 取消了依据参考分型进行SNP过滤的步骤，此处原本的_dict_base_count_filter函数修改为_snp_count
                """
                ls_dict_snp_count = self._snp_count(ls_spe_mh, len_mh_allele)  # 获取过滤后的特定MH SNP计数字典
                ls_he_ho = self._extract_top(ls_dict_snp_count, peak_height_ratio)  # 依据SNP进行过滤

                """
                输出特定MH的基因型信息（2种）：1.原始结果；2.将原始MH进行降级后的结果
                同时输出降级后结果对应的INFO文件
                //2022.02.14 由于取消了mh_filter，此处_extract_geno中改为直接使用ls_spe_mh
                """
                mh_geno_origin = [id] + self._extract_geno(ls_spe_mh, ls_he_ho, peak_height_ratio)
                mh_geno_new, mh_info_new = a_call.analysis_calling(mh_geno_origin, self.dict_chr_info)

            else:  # 如果特定MH的match信息不存在，则输出空
                mh_geno_origin = [id, "no_match"]
                mh_geno_new = ["no_match"]
                mh_info_new = ["no_match"]

            ls_geno_origin_results.append(mh_geno_origin)
            ls_geno_new_results.append(mh_geno_new)
            ls_info_new_results.append(mh_info_new)

        """文件储存"""
        if save:
            geno_origin_name = filename.replace(".sam", ".txt")
            geno_new_name = filename.replace(".sam", "_new.txt")
            info_new_name = filename.replace(".sam", "_new.info")

            Toolsbox.FileTools.save_file(Toolsbox.FormatTools.str_ls(ls_geno_origin_results), save_path+geno_origin_name)
            Toolsbox.FileTools.save_file(Toolsbox.FormatTools.str_ls(ls_geno_new_results), save_path+geno_new_name)
            Toolsbox.FileTools.save_file(ls_info_new_results, save_path+info_new_name)
        else:
            pass
        print("{}文件MH calling已结束".format(filename))

    def singlecall(self, peak_height_ratio=0.3, save=True, save_path=""):
        """
        主函数puzzle_no_n_allele
        """
        if not isinstance(peak_height_ratio, float):
            raise ValueError('"peak_height_ratio" must be float type')
        if not isinstance(save, bool):
            raise ValueError('"save" must be true or false')
        if not isinstance(save_path, str):
            raise ValueError('"save_path" must be string')

        if peak_height_ratio >= 1 or peak_height_ratio <= 0:
            raise ValueError('"peak_height_ratio" must between 0 and 1')

        ls_filename = Toolsbox.FileTools.ls_directory(self.match_filepath)

        for filename in ls_filename:  # 对单个match文件进行操作
            print("{}文件开始进行MH calling".format(filename))
            self._singlecall(filename, peak_height_ratio, save, save_path)
