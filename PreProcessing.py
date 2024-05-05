# -*- coding: utf-8 -*-
# @Time    : 2022
# @Author  : WangHaoyu
# @E-mail  : wanghy0707@gmail.com
# @Github  :
# @desc    :


import re
import Toolsbox


class Extract_INFO:
    """
    提取INFO文件（精简后）中的信息
    提取的信息均为列表或列表字典嵌套形式
    """

    def __init__(self, INFO_file):
        """
        INFO_file: 导入的INFO文件（精简后）
        """
        self.INFO_file = INFO_file

    def extract_chr(self):
        """
        获取INFO文件中所有MH的染色体信息（用于对SNP进行再编号）
        """
        return [mh[0][3:8] for mh in self.INFO_file]

    def extract_id(self):
        """
        获取INFO文件中所有MH的ID信息
        """
        return [mh[0] for mh in self.INFO_file]

    def extract_bp(self):
        """
        获取INFO文件中所有MH的bp信息（int）
        """
        return [int(mh[1]) for mh in self.INFO_file]

    def extract_ae(self):
        """
        获取INFO文件中所有MH的AE值信息（float）
        """
        return [float(mh[2]) for mh in self.INFO_file]

    def extract_het(self):
        """
        获取INFO文件中所有MH的杂合度信息（float）
        """
        return [float(mh[3]) for mh in self.INFO_file]

    def extract_dp(self):
        """
        获取INFO文件中所有MH的个人识别概率信息（float）
        """
        return [float(mh[4]) for mh in self.INFO_file]

    def extract_snp_num(self):
        """
        获取INFO文件中所有MH的SNP数目信息（int）
        """
        return [int(mh[5]) for mh in self.INFO_file]

    def extract_allele_num(self):
        """
        获取INFO文件中所有MH的等位基因数目信息（int）
        """
        return [int(mh[6]) for mh in self.INFO_file]

    def extract_snp(self):
        """
        获取INFO文件中所有MH的SNP位置信息（二维列表形式，返回列表）
        [[mh1snppos1, mh1snppos2, ...],
        [mh2snppos1, mh2snppos2, ...]]
        """
        ls_snp_num = self.extract_snp_num()
        return [mh[7:7 + snp_num] for mh, snp_num in zip(self.INFO_file, ls_snp_num)]


class PreVcf:
    """
    用于对vcf进行预处理：保留染色体名称、SNP POS、参考等位基因分型
    //2021.09.10
    //pre_vcf = PreVcf("C:/Users/Hao_Yu/Desktop/CHB/")
    //dict_chr_refalt = pre_vcf.vcf_filter()
    //>>> {'1': {4658311(int): 'TC', 4658339: 'TA', 4658346: 'TA', ...}, ...}
    """

    def __init__(self, vcf_filepath):
        """
        vcf_filepath: 储存vcf文件的文件夹路径
        """
        self.vcf_filepath = vcf_filepath

    def _vcf_ref_alt(self, ls_vcf_snp):
        """
        用于判断单条SNP共需要保留几项REF和ALT，并进行相应去尾处理：
        若为PASS，说明为二等位SNP，则去掉倒数两项；若为100，说明为三等位SNP，则去掉倒数一项；若不为前二者，说明为四等位SNP，全部保留
        """
        if ls_vcf_snp[-1] == "PASS":
            return "".join(ls_vcf_snp[:-2])
        elif ls_vcf_snp[-1] == "100":
            return "".join(ls_vcf_snp[:-1])
        else:
            return "".join(ls_vcf_snp)

    def vcf_filter(self):
        """
        主函数
        用于输出最终的结果（以嵌套字典形式）：
        {1: {pos1: refalt1, pos2: refalt2, ...},
        2: {pos1: refalt1, pos2: refalt2, ...}, ...}
        """
        dict_chr_refalt = {}  # 储存最终结果的字典

        """列出文件夹中所有vcf文件的名字"""
        vcf_filename = Toolsbox.FileTools.ls_directory(self.vcf_filepath)

        """对单个vcf文件进行过滤操作"""
        for filename in vcf_filename:
            chrname = filename.split(".")[0][3:]  # 提取当前vcf文件的染色体号
            vcf_file = Toolsbox.FileTools.open_vcf_file(self.vcf_filepath + filename)  # 以去掉注释行的形式打开vcf文件

            dict_pos_refalt = {}
            for snp in vcf_file:
                pos = int(snp[1])  # 提取pos时使用int的方式，方便后续sam与info文件计算后直接进行检索
                str_vcf_ref_alt = self._vcf_ref_alt(snp[3:7])  # 提取该SNP的REF和ALT分型
                dict_pos_refalt.update({pos: str_vcf_ref_alt})

            dict_chr_refalt.update({chrname: dict_pos_refalt})
        print("vcf已过滤")
        return dict_chr_refalt


class PreInfo:
    """
    用于对info文件进行预处理：保留染色体名称，MHID、包含的SNP数目、SNP POS列表
    //2021.09.10：
    //pre_info = PreInfo("C:/Users/Hao_Yu/Desktop/179.txt")
    //dict_chr_info = pre_info.info_filter()
    //>>> {'1': {247032240: ['CHBCHR01_0674112', 12, 247032240, 247032247, 247032249, ...], ...}, ...}
    """

    def __init__(self, info_filepath):
        """
        info_filepath: info文件的绝对路径
        """
        self.info_filepath = info_filepath

    def _info_snp(self, ls_snp_str):
        """
        将单个MH的SNP列表由str转化为int
        """
        ls_snp_int = [int(snp) for snp in ls_snp_str]
        return ls_snp_int

    def info_filter(self):
        """
        主函数
        用于输出最终的结果（以嵌套字典的形式）：
        {'1': {start_pos1: [id, snp_num, pos1, pos2,...], start_pos2: [], ...},
        '2': {start_pos1: [id, snp_num, pos1, pos2,...], ...}, ...}
        """
        dict_chr_info = {str(chrname): {} for chrname in range(1, 23)}  # 储存最终结果的字典

        """打开INFO文件对INFO文件进行操作"""
        info_file = Toolsbox.FileTools.open_ls_file(self.info_filepath)  # 打开info文件
        info_c = Extract_INFO(info_file)

        ls_id = info_c.extract_id()  # 提取INFO文件中MH的ID（str）
        ls_snp_num = info_c.extract_snp_num()  # 提取INFO文件中的SNP数目（int）
        ls_snp = info_c.extract_snp()

        for id, snp_num, snp in zip(ls_id, ls_snp_num, ls_snp):
            chrname = str(int(id[6:8]))  # 提取该MH的染色体信息（用于update字典）

            ls_snp_int = self._info_snp(snp)  # 将SNP列表转化为int形式
            snp_first = ls_snp_int[0]  # 提取SNP列表中的第一个SNP作为字典的key

            ls_id_snp_num = [id, snp_num]
            ls_id_snp_num.extend(ls_snp_int)  # 将同一个MH的信息整合在一个列表中

            dict_merge = {snp_first: ls_id_snp_num}
            dict_chr_info[chrname].update(dict_merge)  # 将单一MH的信息update到总字典中

        print("INFO已过滤")
        return dict_chr_info, ls_id


class PrePrimer:
    """
    用于对引物文件进行预处理：输出MH ID、引物序列信息
    //2021.10.19 新增对Primer文件的处理
    //2022.03.15 已知dict_mh_primer仅在match过程中使用，故将调整dict_mh_primer的格式，以便于修改match中的引物匹配规则
                 新修改后的输出格式为：
                 {MHID1: [[0, primer1], [1, primer]], ...
                 MHIDn: [[0, primer1]]}
    """

    def __init__(self, primer_filepath, ls_id):
        """
        primer_filepath: primer文件的绝对路径
        ls_id: 由info_filter函数生成的MH ID列表
        """
        self.primer_filepath = primer_filepath
        self.ls_id = ls_id

    def primer_load(self):
        """
        导入primer文件，以字典形式:
        {id1: [primer1], ...,
        idn: [primer1, primer2]}
        """
        dict_mh_primer = {id: [] for id in self.ls_id}  # 储存最终结果的字典

        """打开primer文件对primer文件进行操作"""
        primer_file = Toolsbox.FileTools.open_ls_file(self.primer_filepath)

        for primer in primer_file:  # 对单条primer信息进行操作
            dict_mh_primer[primer[0]].append(primer[1:3])
        print("primer文件已处理")
        return dict_mh_primer


class PreSam:
    """
    用于对sam文件进行预处理：保留染色体、起始POS、原始序列信息
    过滤后的信息将保留一份在本地
    //2022.03.14 复核：sam文件预过滤无问题
    """

    def __init__(self, sam_filepath, dict_chr_info):
        """
        sam_filepath: 储存sam文件的文件夹路径（对测序得到的多个样本文件均进行处理）
        dict_chr_info: 由PreInfo得到的结果（嵌套字典）
        """
        self.sam_filepath = sam_filepath
        self.dict_chr_info = dict_chr_info

    def _sam_del_redundancy_chr(self, chrname):
        """
        去冗余信息时对染色体信息进行转化
        //2021.09.10 此函数的转换内容视测序中对染色体的编号变化而变更
        """
        if "NC" in chrname:
            chrname_str = int(chrname.split(".")[0][-2:])

            if 1 <= chrname_str <= 22:
                return str(chrname_str)

    def _sam_del_redundancy_sequence(self, sequence_info, sequence):
        """
        依据序列信息对序列中的insert和deletion进行处理，输出sequence_m
        //2021.10.13 此函数最终输出结果仅保留了M部分的结果，去除了前S和后S
        //2021.10.19 补充考虑当序列信息为H起始时，原始序列仅包含M部分的情况，对于此类序列做删除处理
        //2022.03.14 复核：I与D的添加、删除问题。在原始序列中，包含标记为I的碱基，不包含标记为D的碱基。若要想测序序列碱基的pos与
                     参考序列的pos一致（此处pos指的是M部分的pos），则需要去除序列中I对应的部分，补全D对应的部分
        //2022.03.18 新增S在前或在后均进行输出的判断；新增纳入序列信息中包含H的序列提取
        """
        find = re.findall(r'[0-9]+|[A-Z]+', sequence_info)  # 匹配序列信息中的数字和字母

        num = [int(i) for i in find[::2]]  # 提取序列信息中的数字，储存在列表中（eg. [22, 17, 1, 112, 43]）
        type_str = find[1::2]  # 提取序列信息中的字母，储存在列表中（eg. [S, M, I, M, S]）
        len_num = len(num)

        if type_str[0] == "S" or type_str[-1] == "S":  # 如果序列以S开头或结尾（序列中有S）

            for i in range(len_num):

                if type_str[i] == "I":  # 如果类型是插入
                    sequence = sequence[:sum(num[:i])] + sequence[sum(num[:i + 1]):]
                    num[i] = 0
                elif type_str[i] == "D":  # 如果类型是缺失
                    sequence = sequence[:sum(num[:i])] + num[i] * "D" + sequence[sum(num[:i]):]

            if type_str[0] == "S" and type_str[-1] == "S":  # 去头去尾
                return sequence[num[0]:-num[-1]]
            elif type_str[0] == "S" and type_str[-1] != "S":  # 去头留尾
                return sequence[num[0]:]
            elif type_str[0] != "S" and type_str[-1] == "S":  # 留头去尾
                return sequence[0:-num[-1]]

        elif type_str[0] == "H" or type_str[-1] == "H":  # 如果序列与以H开头或结尾，则序列本身就是要输出的序列，只需去I补D即可

            if type_str[0] != "H":  # 如果头部无H，检索位置无需改动

                for i in range(len_num):

                    if type_str[i] == "I":  # 如果类型是插入
                        sequence = sequence[:sum(num[:i])] + sequence[sum(num[:i + 1]):]
                        num[i] = 0
                    elif type_str[i] == "D":  # 如果类型是缺失
                        sequence = sequence[:sum(num[:i])] + num[i] * "D" + sequence[sum(num[:i]):]
                return sequence

            elif type_str[0] == "H":

                for i in range(1, len_num):

                    if type_str[i] == "I":  # 如果类型是插入
                        sequence = sequence[:sum(num[1:i])] + sequence[sum(num[1:i + 1]):]
                        num[i] = 0
                    elif type_str[i] == "D":  # 如果类型是缺失
                        sequence = sequence[:sum(num[1:i])] + num[i] * "D" + sequence[sum(num[1:i]):]
                return sequence

        else:  # 如果序列以除开头、结尾
            return None

    def _sam_filter_single(self, sequence_single):
        """
        依据单条测序结果的起始POS和dict_chr_info中的信息，去除冗余的测序结果
        sequence_single: 去除冗余信息后的单条sequence结果，包括[chrname, startpos, len_sequence_m, sequence_m, sequence]
        //2021.09.11  panel中最大MH的长度可依据所使用panel的不同而进行更改
        //2022.03.14  复核：预过滤sam文件中判断序列是否覆盖MH无问题。同时简化了输入参数。此外注意序列的startpos对应sequence_m
        """
        ls_mh_startpos = self.dict_chr_info[sequence_single[0]].keys()  # 依据测序结果的染色体信息，提取该染色体上所有MH的起始POS

        for snppos in ls_mh_startpos:  # 对每一个MH的起始POS进行操作，直到有一个MH的起始POS能与该序列匹配才输出结果

            if snppos - sequence_single[1] > 0:  # 如果测序序列在MH的“前面”

                if snppos - sequence_single[1] > sequence_single[2]:  # 如果测序序列起始与MH的距离超过测序序列长度，则测序序列未测到当前MH
                    continue
                elif snppos - sequence_single[1] <= sequence_single[2]:  # 如果测序序列起始与MH的距离小于测序序列长度，则测序序列测到了当前MH（至少测到第一个SNP）
                    return sequence_single

            elif snppos - sequence_single[1] <= 0:  # 如果测序序列在MH的“后面”（包含测序序列起始与MH起始重叠的情况）

                if sequence_single[1] - snppos > 200:  # 如果测序序列起始与MH的距离超过panel中MH的最大长度，则测序序列未测到当前MH
                    continue
                elif sequence_single[1] - snppos <= 200:  # 如果测序序列起始与MH的距离小于panel中MH的最大长度，则测序序列测到了当前MH（至少测到最后一个SNP）
                    return sequence_single

    def _sam_filter(self, sam_file):
        """
        依据每条测序结果的起始POS和dict_chr_info中的信息，去除冗余的测序结果
        sam_file: 导入的sam文件
        """
        ls_sam_filter = []
        for sequence_single in sam_file:  # 对每一条测序结果进行操作
            chrname = self._sam_del_redundancy_chr(sequence_single[2])  # 提取染色体信息

            if chrname:
                startpos = int(sequence_single[3])  # 提取序列起始POS信息
                sequence_info = sequence_single[5]  # 提取序列的INFO信息（eg. 22S17M1I112M43S）
                sequence = sequence_single[9]  # 提取序列信息

                sequence_m = self._sam_del_redundancy_sequence(sequence_info, sequence)  # 将序列信息转化(仅保留M部分)
                if sequence_m:  # 如果序列能依据序列信息提取出sequence_m则继续操作
                    len_sequence_m = len(sequence_m)
                    ls_sequence_single = [chrname, startpos, len_sequence_m, sequence_m, sequence]  # 去冗余信息后的单条测序结果

                    sequence_filter = self._sam_filter_single(ls_sequence_single)  # 对单一测序结果进行判断

                    if sequence_filter:  # 如果单一测序结果可与某一MH有重叠，则将其保留
                        ls_sam_filter.append(sequence_filter)
                    else:
                        continue
                else:
                    continue
            else:
                continue

        return ls_sam_filter

    def sam_filter(self, sam_filter_filepath):
        """
        主函数
        用于输出最终的去冗余、过滤结果（以列表的形式）：
        [[chr, start_pos, len_sequence, sequence],
        [chr, start_pos, len_sequence, sequence], ...]
        """

        """列出文件夹中所有sam文件的名字"""
        sam_filename = Toolsbox.FileTools.ls_directory(self.sam_filepath)
        print(sam_filename)

        """对每一个sam文件进行操作"""
        for filename in sam_filename:
            sam_file = Toolsbox.FileTools.open_ls_file(self.sam_filepath + filename)  # 打开sam文件

            sam_filter = self._sam_filter(sam_file)

            """储存临时结果文件"""
            save_path = sam_filter_filepath + "{}".format(filename)  # 储存临时结果文件的路径（仍以原文件名命名）
            Toolsbox.FileTools.save_file(Toolsbox.FormatTools.str_ls(sam_filter), save_path)
            print("{}文件已过滤".format(filename))