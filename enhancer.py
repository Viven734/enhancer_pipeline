import argparse

parser = argparse.ArgumentParser(description="RNA-seq analysis pipeline for non-model organisms.")
parser.add_argument("-o","--output_dir", required = True, type = str, help = "")
parser.add_argument("-c","--cancer_erna", required = True, type = str, help = "")
parser.add_argument("-n","--normol_erna", required = True, type = str, help = "")
parser.add_argument("-cl","--clinical", required = True, type = str, help = "")
parser.add_argument("-e","--count", required = True, type = str, help = "")
args = parser.parse_args()

print("checking modules")
import glob,os,sys,time,numpy
import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri, default_converter
from rpy2.robjects.conversion import localconverter
import pandas as pd
print("imports done")

Rscript_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'Rscript')

os.chdir(args.output_dir)

print("Begin ifferentially analysis...")
#1.Differentially analysis
out_dir = os.path.join(args.output_dir,'01_Differentially_analysis')
if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
if os.path.isfile((out_dir+'/all_sample_eRNA.txt')) == False:
    # 读取两个CSV文件，假设我们以第一列作为行名
    df1 = pd.read_csv(args.normol_erna, index_col=0)
    df2 = pd.read_csv(args.cancer_erna, index_col=0)
    # 找出两个数据框中都存在的列名
    common_columns = df1.columns.intersection(df2.columns)

    # 在第二个数据框中删除这些列
    df2 = df2.drop(columns=common_columns)
    # 合并两个CSV文件
    df = pd.merge(df1, df2, left_index=True, right_index=True)
    # 将结果保存到新的TXT文件，使用制表符作为分隔符
    df.to_csv(out_dir+'/all_sample_eRNA.txt', sep='\t')
if os.path.isfile((out_dir+'/filter_count_eRNA.txt')) == False:
    numpy2ri.activate()
    with localconverter(default_converter + numpy2ri.converter):
        robjects.r.source(Rscript_path+'/2_differentially_analysis.R')
        robjects.r.differentially_analysis(out_dir)
print("Differentially analysis has finished")

print("Begin prognostic analysis...")
int_dir = os.path.join(args.output_dir,'01_Differentially_analysis')
out_dir = os.path.join(args.output_dir,'02_Prognostic_analysis_RPKM_(7vs3)')
if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
#cl.txt
if os.path.isfile((out_dir+'/cl.txt')) == False:
    df = pd.read_excel(args.clinical)
    # 只保留前三列
    df = df.iloc[:, 1:4]
    df.columns.values[1]='OS'
    df['OS'] = df['OS'].replace(0,1)
    # 输出成制表符分隔的txt文件
    df.to_csv(out_dir+'/cl.txt', sep='\t', index=False)
#exp.txt
# 读取txt文件和xls文件
df_txt = pd.read_csv(int_dir+"/filter_count_eRNA.txt", sep='\t', index_col=0)
df_ids = pd.read_csv(int_dir+"/DIFF_all.txt", sep='\t')
df2 = pd.read_csv(out_dir+'/cl.txt', sep='\t', index_col=0)
if os.path.isfile((out_dir+'/exp.txt')) == False:
    # 只保留txt文件中行名与id列值相同的行
    df_txt = df_txt[df_txt.index.isin(df_ids['id'])]
    common_cols = list(set(df_txt.columns) & set(df2.index))
    df_txt = df_txt[common_cols]
    # 输出结果
    df_txt.to_csv(out_dir+'/exp.txt', sep='\t')
if os.path.isfile((out_dir+'/test_cl.csv')) == False:
    numpy2ri.activate()
    with localconverter(default_converter + numpy2ri.converter):
        robjects.r.source(Rscript_path+'/3_training.R')
        robjects.r.training_and_testing(out_dir)
if os.path.isfile((out_dir+'/diff_unicox.csv')) == False:
    numpy2ri.activate()
    with localconverter(default_converter + numpy2ri.converter):
        robjects.r.source(Rscript_path+'/4_cox.R')
        robjects.r.cox(out_dir)

int_dir = out_dir
out_dir1=os.path.join(out_dir,'train')
if os.path.isdir(out_dir1) == False: os.mkdir(out_dir1)
out_dir2=os.path.join(out_dir,'test')
if os.path.isdir(out_dir2) == False: os.mkdir(out_dir2)
out_dir3=os.path.join(out_dir,'combined')
if os.path.isdir(out_dir3) == False: os.mkdir(out_dir3)
if os.path.isfile((out_dir1+"/train_lasso_min_eRNAs_unicox_p0.001.txt")) == False:
    numpy2ri.activate()
    with localconverter(default_converter + numpy2ri.converter):
        robjects.r.source(Rscript_path+'/5_lasso_multicox.R')
        robjects.r.lasso(int_dir,out_dir1,out_dir2,out_dir3)
#KM
KM1=os.path.join(out_dir1,'KM')
if os.path.isdir(KM1) == False: os.mkdir(KM1)
KM2=os.path.join(out_dir2,'KM')
if os.path.isdir(KM2) == False: os.mkdir(KM2)
KM3=os.path.join(out_dir3,'KM')
if os.path.isdir(KM3) == False: os.mkdir(KM3)
if os.path.isfile((KM1+"/Training_set.pdf")) == False:
    numpy2ri.activate()
    with localconverter(default_converter + numpy2ri.converter):
        robjects.r.source(Rscript_path+'/6_train_KM.R')
        robjects.r.KM(out_dir1,KM1)
if os.path.isfile((KM2+"/Testing_set.pdf")) == False:
    numpy2ri.activate()
    with localconverter(default_converter + numpy2ri.converter):
        robjects.r.source(Rscript_path+'/6_test_KM.R')
        robjects.r.KM(out_dir2,KM2)
if os.path.isfile((KM1+"/Combined_set.pdf")) == False:
    numpy2ri.activate()
    with localconverter(default_converter + numpy2ri.converter):
        robjects.r.source(Rscript_path+'/6_combined_KM.R')
        robjects.r.KM(out_dir3,KM3)
#ROC
ROC1=os.path.join(out_dir1,'ROC')
if os.path.isdir(ROC1) == False: os.mkdir(ROC1)
ROC2=os.path.join(out_dir2,'ROC')
if os.path.isdir(ROC2) == False: os.mkdir(ROC2)
ROC3=os.path.join(out_dir3,'ROC')
if os.path.isdir(ROC3) == False: os.mkdir(ROC3)
if os.path.isfile((ROC1+"/Train_ROC.pdf")) == False:
    numpy2ri.activate()
    with localconverter(default_converter + numpy2ri.converter):
        robjects.r.source(Rscript_path+'/7_train_ROC.R')
        robjects.r.ROC(out_dir1,ROC1)
if os.path.isfile((ROC2+"/Test_ROC.pdf")) == False:
    numpy2ri.activate()
    with localconverter(default_converter + numpy2ri.converter):
        robjects.r.source(Rscript_path+'/7_test_ROC.R')
        robjects.r.ROC(out_dir2,ROC2)
if os.path.isfile((ROC1+"/Combined_ROC.pdf")) == False:
    numpy2ri.activate()
    with localconverter(default_converter + numpy2ri.converter):
        robjects.r.source(Rscript_path+'/7_combined_ROC.R')
        robjects.r.ROC(out_dir3,ROC3)
print("Prognostic analysis has finished")
print("Begin GSEA...")
out_dir = os.path.join(args.output_dir,'03_GSEA')
if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
if os.path.isfile((out_dir+"/filter_count_gene_exp.txt")) == False:
    numpy2ri.activate()
    with localconverter(default_converter + numpy2ri.converter):
        robjects.r.source(Rscript_path+'/8_GESA.R')
        robjects.r.GESA(args.count,out_dir)
print("GSEA has finished")
print("Begin Drug prediction...")
int_dir=out_dir
out_dir = os.path.join(args.output_dir,'04_Drug_prediction')
if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
currentdir = os.path.dirname(os.path.abspath(__file__))
df_txt = pd.read_csv(int_dir+"/filter_count_gene_exp.txt", sep='\t', index_col=0)
df_ids = pd.read_csv(int_dir+"/diff_gene.txt", sep='\t')
if os.path.isfile((out_dir+'/exp.txt')) == False:
    # 只保留txt文件中行名与id列值相同的行
    df_txt = df_txt[df_txt.index.isin(df_ids['id'])]
    # 输出结果
    df_txt.to_csv(out_dir+'/diff_gene_exp.txt', sep='\t')
# if os.path.isfile((out_dir+"/filter_count_gene_exp.txt")) == False:
    numpy2ri.activate()
    with localconverter(default_converter + numpy2ri.converter):
        robjects.r.source(Rscript_path+'/9_drug.R')
        robjects.r.drug(currentdir,out_dir)

