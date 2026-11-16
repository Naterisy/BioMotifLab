"""Run with: streamlit run app.py"""
import time
import pandas as pd
import streamlit as st
from pathlib import Path
from biomotiflab.io import read_fasta
from biomotiflab.core import discover, reverse_complement
from biomotiflab.stats import calibrate
from biomotiflab.report import logo_figure, json_text, markdown_report

st.set_page_config(page_title='BioMotifLab',page_icon='🧬',layout='wide')
st.title('BioMotifLab · DNA 基序研究工作台')
st.caption('双链搜索 · 可解释 PWM · 独立留出检验 · 可导出实验记录')
st.info('OOPS 假设每条序列恰有一个基序。结果用于研究探索；随机背景显著性不等同于真实生物结合。')
with st.sidebar:
    st.header('分析参数')
    methods = st.multiselect('对比算法',['uniform','global','adaptive','kmer'],default=['global','adaptive'])
    width = st.slider('基序宽度',6,24,10)
    starts = st.slider('初始化次数',2,16,8)
    iterations = st.slider('最大 EM 迭代',10,100,60)
    shrinkage = st.number_input('收缩强度 τ',min_value=0.,max_value=10000.,value=40.)
    seed = st.number_input('随机种子',min_value=0,value=42)
    shuffles = st.select_slider('留出置换次数',options=[49,99,199,499],value=199)
train_file = st.file_uploader('上传发现集 FASTA（ACGT；N 窗口会排除）',type=['fa','fasta','fna','txt','s'])
test_file = st.file_uploader('可选：独立测试集 FASTA，用于固定基序置换检验',type=['fa','fasta','fna','txt','s'])
example = st.checkbox('使用内置可复现合成示例',value=True)
st.caption('限制：≤200条序列，每条≤2000 bp，总计≤100000 bp；公开演示时请谨慎上传私人数据。')
if st.button('执行真实分析',type='primary'):
    try:
        if not methods: raise ValueError('请至少选择一个算法。')
        if train_file is not None:
            seqs = read_fasta(train_file.getvalue().decode('utf-8'),is_text=True)
        elif example:
            seqs = read_fasta(str(Path(__file__).parent/'data/example/train.fasta'))
        else: raise ValueError('请上传 FASTA 或使用示例。')
        test = read_fasta(test_file.getvalue().decode('utf-8'),is_text=True) if test_file else []
        for batch in [seqs,test]:
            if batch and (len(batch)>200 or max(map(len,batch))>2000 or sum(map(len,batch))>100000):
                raise ValueError('数据超出交互计算限制，请使用命令行。')
        if {min(s, reverse_complement(s)) for s in seqs} & {min(s, reverse_complement(s)) for s in test}: raise ValueError('训练与测试存在相同或反向互补序列，不能作为独立留出检验。')
        results = []
        with st.spinner('正在搜索双链窗口并优化 PWM…'):
            for method in methods:
                r = discover(seqs,width,method=method,seed=int(seed),n_starts=starts,max_iter=iterations,shrinkage=shrinkage)
                r['n_train'] = len(seqs)
                r['configuration'] = dict(width=width,seed=int(seed),n_starts=starts,max_iter=iterations,shrinkage=shrinkage,shuffles=shuffles)
                if test:
                    r['heldout'] = calibrate(test,r['pwm'],method=method if method!='kmer' else 'global',n_shuffles=shuffles,seed=int(seed),shrinkage=shrinkage,global_background=r['global_background'])
                results.append(r)
        st.session_state['results'] = results
    except (ValueError,OSError,UnicodeError) as exc:
        st.error(str(exc))
        st.session_state.pop('results',None)
if 'results' in st.session_state:
    results = st.session_state['results']
    st.subheader('算法比较')
    st.dataframe(pd.DataFrame([dict(算法=r['method'],共识=r['consensus'],耗时秒=r['runtime_seconds'],迭代=r['iterations']) for r in results]),hide_index=True)
    st.caption('不同背景模型的训练似然不宜直接视为准确率；有植入真值时才计算定位指标。')
    for r in results:
        with st.expander(r['method']+' · '+r['consensus'],expanded=True):
            st.pyplot(logo_figure(r['pwm']))
            st.dataframe(pd.DataFrame(r['pwm'],columns=list('ACGT')))
            st.dataframe(pd.DataFrame(r['sites']),hide_index=True)
            if 'heldout' in r:
                st.write('独立测试：序列级最大双链窗口统计量，BH 校正同一算法的测试序列。')
                st.dataframe(pd.DataFrame(r['heldout']),hide_index=True)
                st.caption('跨算法、跨宽度反复筛选尚未校正；不同序列可能相关，BH 的适用假设需单独审查。')
            st.download_button('下载 JSON',json_text(r),file_name=r['method']+'.json',mime='application/json',key=r['method']+'json')
            st.download_button('下载分析报告',markdown_report(r),file_name=r['method']+'.md',mime='text/markdown',key=r['method']+'md')
