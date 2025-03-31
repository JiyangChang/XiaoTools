#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re
import os
import argparse

class Maxlen:
    def __init__(self, gff, outgff):
        self.gff = gff
        self.outgff = outgff
    
    def _getmax(self, genelen):
        """获取每个基因中最长的mRNA ID
        
        Args:
            genelen (dict): 包含基因和对应mRNA长度的字典
            
        Returns:
            list: 包含最长mRNA ID的列表
        """
        mrnaids = []
        for geneid in genelen:
            # 按长度降序排序mRNA
            newidlen = sorted(genelen[geneid].items(), key=lambda x: x[1], reverse=True)
            # 取最长的mRNA ID
            mrnaids.append(newidlen[0][0])  # (id, len)
        return mrnaids
    
    def _readgff(self):
        """读取GFF文件并提取相关信息
        
        Returns:
            tuple: 包含三个字典
                - gene_mrna_len: 基因到mRNA长度的映射
                - mrna_info: mRNA到其GFF行信息的映射
                - rna_proid: mRNA ID映射(保留但未使用)
                - mrna_gene: mRNA到基因的映射
        """
        gene_mrna_len = {}
        mrna_info = {}
        rna_proid = {}
        mrna_gene = {}
        
        with open(self.gff, 'r') as con:
            for line in con:
                if line.startswith("#") or not line.strip():
                    continue
                
                lines = line.strip().split("\t")
                if len(lines) < 9:
                    continue
                
                if lines[2] == "mRNA":
                    # 提取mRNA ID和父基因ID
                    idinfo = re.findall(r'ID=([^;]+);.*Parent=([^;]+)', lines[8])
                    if not idinfo:
                        continue
                    
                    mrna_id, gene_id = idinfo[0]
                    mrna_gene[mrna_id] = gene_id
                    
                    # 初始化基因到mRNA长度的映射
                    if gene_id not in gene_mrna_len:
                        gene_mrna_len[gene_id] = {}
                    gene_mrna_len[gene_id][mrna_id] = 0
                    
                    # 保存mRNA行信息
                    mrna_info[mrna_id] = line
                
                elif lines[2] == "CDS":
                    # 提取CDS的父mRNA ID
                    idinfocds = re.findall(r'Parent=([^;]+)', lines[8])
                    if not idinfocds:
                        continue
                    
                    mrna_id = idinfocds[0]
                    rna_proid[mrna_id] = mrna_id
                    
                    if mrna_id in mrna_info:
                        # 追加CDS行信息
                        mrna_info[mrna_id] += line
                        
                        # 计算CDS长度并累加
                        gene_id = mrna_gene[mrna_id]
                        cdslen = int(lines[4]) - int(lines[3]) + 1
                        gene_mrna_len[gene_id][mrna_id] += cdslen
                    else:
                        sys.stderr.write(f"Warning: {mrna_id} has no mRNA info!\n")
        
        return gene_mrna_len, mrna_info, rna_proid
    
    def main(self):
        """主处理函数，执行整个流程"""
        try:
            with open(self.outgff, 'w') as out:
                g_m_l, m_i, r_p = self._readgff()
                maxmrnaidlen = self._getmax(g_m_l)
                
                # 写入最长的mRNA及其CDS信息
                for mrna_id in maxmrnaidlen:
                    if mrna_id in m_i:
                        out.write(m_i[mrna_id])
                    else:
                        sys.stderr.write(f"Warning: {mrna_id} not found in mRNA info\n")
        except IOError as e:
            sys.stderr.write(f"Error: {str(e)}\n")
            sys.exit(1)

def main():
    # 使用argparse处理命令行参数
    parser = argparse.ArgumentParser(
        description='Select the longest mRNA (based on CDS length) for each gene from GFF file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('gff', help='Input GFF file')
    parser.add_argument('outgff', help='Output GFF file containing only the longest mRNAs')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.isfile(args.gff):
        sys.stderr.write(f"Error: Input file {args.gff} does not exist\n")
        sys.exit(1)
    
    # 执行处理
    processor = Maxlen(args.gff, args.outgff)
    processor.main()

if __name__ == "__main__":
    main()
