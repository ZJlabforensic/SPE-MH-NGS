# -*- coding: utf-8 -*-
# @Time    : 2021//
# @Author  : WangHaoyu
# @E-mail  : wanghy0707@gmail.com
# @Github  :
# @desc    :


import os
import collections
import copy


class FileTools:
    """
    文件操作类工具
    """

    @classmethod
    def ls_directory(cls, path):
        """
        获取文件夹中所有文件的名称
        """
        return os.listdir(path)

    @classmethod
    def open_ls_file(cls, path, separate_sym="\t", type="str"):
        """
        依据文件的行分隔符，以列表形式打开文件（一行为一个列表）
        """
        ls_file = []
        with open(path, "r", encoding="utf-8") as file:
            if type == "str":
                for line in file.readlines():
                    ls_file.append(list(line.strip("\n").split(separate_sym)))
            elif type == "float":
                for line in file.readlines():
                    line_str = list(line.strip("\n").split(separate_sym))
                    line_float = [float(i) for i in line_str]
                    ls_file.append(line_float)
        return ls_file

    @classmethod
    def open_dict_file(cls, path, separate_sym="\t"):
        """
        依据文件的行分隔符，以字典形式打开文件（行第一项为key，其后所有项为value）
        """
        dict_file = {}
        with open(path, "r") as file:
            for line in file.readlines():
                id, *_ = list(line.strip("\n").split(separate_sym))
                dict_file[id] = _
        return dict_file

    @classmethod
    def open_allele_file(cls, path, separate_sym="\t"):
        """
        依据文件的行分隔符，以列表形式打开文件（一行为一个列表），专用于打开Allele文件
        """
        ls_allele = []
        with open(path, "r") as file:
            for line in file.readlines():
                ls_allele.append(list(line.strip("\n").split(separate_sym))[1:])
        return ls_allele

    @classmethod
    def open_vcf_file(cls, path, separate_sym="\t"):
        """
        打开vcf文件（不包括注释行）
        """
        ls_vcf = []
        with open(path, "r") as file:
            for line in file.readlines():
                if line.startswith("#"):
                    pass
                else:
                    ls_vcf.append(list(line.strip("\n").split(separate_sym)))
        return ls_vcf

    @classmethod
    def save_file(cls, output_file, output_path):
        """
        以"\t"为分隔符，储存生成的文件（二维列表）
        """
        str_file = "\n".join(["\t".join(ls) for ls in output_file])

        with open(output_path, "w") as file:
            file.write(str_file)
            file.close()

    @classmethod
    def save_file_match(cls, output_file, output_path):
        """
        储存不进行过滤匹配结果：标签顺序为：MHID, primer, allele, UMI
        """
        ls_mh_match_all_sorted = []
        for mh in output_file.keys():  # 对每一个MH进行操作

            for primer in output_file[mh].keys():  # 对MH下的每一个primer进行操作

                for allele in output_file[mh][primer].keys():  # 对primer下的每一个umi进行操作

                    ls_allele_umi_count = []
                    for umi, count in output_file[mh][primer][allele].items():  # 对umi下的等位基因及计数进行操作

                        line = [mh, primer, allele, umi, count]
                        ls_allele_umi_count.append(line)

                    ls_allele_umi_count_sorted = sorted(ls_allele_umi_count, key=lambda x: x[4], reverse=True)
                    ls_mh_match_all_sorted.extend(ls_allele_umi_count_sorted)

        ls_mh_match_all_sorted_str = FormatTools.str_ls(ls_mh_match_all_sorted)
        cls.save_file(ls_mh_match_all_sorted_str, output_path)


class FormatTools:
    """
    格式处理相关小工具
    """

    """一般列表操作"""
    @classmethod
    def merge_ls(cls, ls1, ls2, type="non"):
        """
        使用extend方法合并两个维度相同的列表
        """
        if type == "non":
            for i, j in zip(ls1, ls2):
                i.extend(j)

        elif type == "return":
            ls3 = copy.deepcopy(ls1)
            for i, j in zip(ls3, ls2):
                i.extend(j)
            return ls3

    @classmethod
    def set_ls(cls, ls):
        """
        对二维列表中的每一项进行去重
        """
        return [list(set(i)) for i in ls]

    @classmethod
    def str_ls(cls, ls):
        """
        将二维列表中的非string项（主要针对int和float）转化为string项，在储存文件时使用
        """
        ls_str = []
        for line in ls:
            new_line = [str(i) for i in line]
            ls_str.append(new_line)
        return ls_str

    @classmethod
    def difference_ls(cls, ls1, ls2):
        """
        获取两个列表的差集：ls1中有而ls2没有的元素
        """
        return [i for i in ls1 if i not in ls2]

    @classmethod
    def del_one_ls(cls, ls, index, type="merge"):
        """
        获取从列表中去掉指定index处的项后剩余的列表
        """
        if index > len(ls) - 1:
            raise ValueError('"index" out of range')

        if type == "merge":
            return ls[:index] + ls[index+1:]
        elif type == "divide":
            return ls[:index], ls[index+1:]

    @classmethod
    def flatten_ls(cls, ls):
        """
        将多维列表（列表中可包含不同维度、可包含列表、元组等可迭代对象）降维为一维列表
        """
        ls_flatten = []
        def fly(ls):
            for item in ls:
                if isinstance(item, collections.abc.Iterable) and not isinstance(item, (str, bytes)):
                    fly(item)
                else:
                    ls_flatten.append(item)
        fly(ls)
        return ls_flatten

    @classmethod
    def convert_base(cls, string):
        """
        将序列中的ATCG碱基分别进行转换：A→T，T→A，C→G，G→C
        """
        string_tmp = string.replace("A", "a").replace("T", "t").replace("C", "c").replace("G", "g")
        string = string_tmp.replace("a", "T").replace("t", "A").replace("c", "G").replace("g", "C")
        return string
    
    """
    本套脚本内部使用的功能（前缀加_）
    """
    @classmethod
    def _match_dict_to_list(cls, dict_match_results, id):
        """
        将汇总为字典形式的match文件中特定MHID的内容转化为列表并单独提取出来
        dict_match_results: 储存所有匹配结果的字典：ID-Primer-Allele-UMI-Count
        id: 特定MH的ID
        """
        ls_spe_mh_match = []  # 储存特定MHID下列表的结果
        
        dict_spe_mh_match = dict_match_results[id]  # 提取特定MHID对应的字典
        for primer in dict_spe_mh_match.keys():  # 对MH下的每一个primer进行操作
            
            for allele in dict_spe_mh_match[primer].keys():  # 对primer下的每一个umi进行操作
                    
                for umi, count in dict_spe_mh_match[primer][allele].items():  # 对umi下的等位基因及计数进行操作
                    
                    line = [id, primer, allele, umi, count]
                    ls_spe_mh_match.append(line)
                    
        return ls_spe_mh_match