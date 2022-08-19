# -*- coding: utf-8 -*-
"""
Created on Fri May  3 23:08:35 2019

@author: djric
"""

import pandas as pd

work_book = pd.ExcelFile(r'..\data_gather_non_normalized.xlsx')

chemical_index = {}
for sheet in work_book.sheet_names:
    df = work_book.parse(sheet)
    chemical = df.columns[0]
    chemical_index[sheet] = chemical
    new_df = extract_wanted_data(df)
    path = sheet + '.csv'
    new_df.to_csv(path, index = False)
    
    
def extract_wanted_data(df):
    start_row = df.index[df.iloc[:,2] == 1].tolist()[0]
    
    new_df = df.iloc[start_row:, :]

    time = new_df.iloc[:,0]
    new_df.iloc[:,0] = time - time.iloc[0]
    
    new_cols = ['C' + str(i + 1) for i in range(11)]
    new_cols = ['Time'] + new_cols
    
    new_cols_dic = dict(zip(list(df.columns), new_cols))
    
    new_df.rename(columns = new_cols_dic, inplace = True)
    
    return new_df