# -*- coding: utf-8 -*-
# @Time    : 2022
# @Author  : WangHaoyu
# @E-mail  : wanghy0707@gmail.com
# @Github  :
# @desc    :


import PreProcessing
import Match
import MHCalling


if __name__ == "__main__":
    """
    文件路径设置
    """
    """【输入文件】"""
    """输入文件储存路径"""
    info_filepath = "C:/Users/Hao_Yu/Desktop/use/179_n.txt"  # INFO文件的路径
    primer_filepath = "C:/Users/Hao_Yu/Desktop/use/primer_n.txt"  # primer文件的路径
    info_ref_filepath = "C:/Users/Hao_Yu/Desktop/use/info_ref_n.txt"  # 在SingleCalling降级时的INFO参考库路径（其格式非精简后的INFO，文件来自小潘师兄）
    vcf_filepath = "C:/Users/Hao_Yu/Desktop/use/CHB/"  # 储存群体VCF文件的文件夹路径（用于进行SNP参考分型核对，在更改群体时需注意）
    sam_filepath = "C:/Users/Hao_Yu/Desktop/use/samfile/"  # 储存所有样本SAM文件的文件夹路径

    """【输出文件】"""
    """filter和match文件储存路径"""
    sam_filter_filepath = "C:/Users/Hao_Yu/Desktop/use/sam_filter/"  # 储存过滤后SAM文件的文件夹路径
    sam_match_filepath = "C:/Users/Hao_Yu/Desktop/use/sam_match/"  # 储存匹配结果文件的路径
    sam_mismatch_filepath = "C:/Users/Hao_Yu/Desktop/use/sam_mismatch/"  # 储存匹配错误结果文件的路径

    """MH Calling相关储存路径"""
    mh_calling_single_filepath = "C:/Users/Hao_Yu/Desktop/use/result/MH calling (single)/"  # 储存mh calling(single)结果

    """
    具体函数
    """
    """文件预处理（必须）"""
    dict_chr_info, ls_id = PreProcessing.PreInfo(info_filepath).info_filter()
    dict_chr_relalt = PreProcessing.PreVcf(vcf_filepath).vcf_filter()
    dict_mh_primer = PreProcessing.PrePrimer(primer_filepath,
                                             ls_id).primer_load()

    """Sam文件过滤"""
    PreProcessing.PreSam(sam_filepath,
                         dict_chr_info).sam_filter(sam_filter_filepath)  # 由于SAM文件预处理会储存处理后的文件，故不必赋值

    """MH匹配"""
    Match.MHMatch(dict_chr_info,
                  ls_id,
                  dict_chr_relalt,
                  dict_mh_primer).match(sam_filter_filepath,
                                        sam_match_filepath,
                                        sam_mismatch_filepath,
                                        filter=False,
                                        purity_filter=0.9,
                                        umi_count_filter=10,
                                        allele_num_filter=10,
                                        allele_proportion_filter=0.01)

    """MH Calling"""
    MHCalling.SingleCall(sam_match_filepath,
                         info_ref_filepath,
                         dict_chr_info,
                         ls_id,
                         dict_chr_relalt).singlecall(peak_height_ratio=0.3,
                                                     save=True,
                                                     save_path=mh_calling_single_filepath)