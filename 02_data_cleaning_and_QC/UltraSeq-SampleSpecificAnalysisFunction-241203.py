import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import math
import seaborn as sns
import regex
import scipy.stats
from scipy.stats import rankdata


def Calculate_Sample_Specific_Relative_Normalized_Metrics(input_df1, input_df2, percentile_list, group_trait='gRNA'):
    """
    Calculates sample-specific relative normalized metrics.
    """
    if input_df2 is not None:
        temp_df = input_df1.groupby(
            [group_trait, 'Sample_ID', 'Type'], as_index=False
        ).apply(Cal_Tumor_Size_simple, percentile_list)
        
        temp_df2 = input_df2.groupby(
            [group_trait, 'Type'], as_index=False
        ).apply(Cal_Tumor_Size_Cas9_negative)

        normalized_metrics_list = []
        for sample_id, group in temp_df.groupby('Sample_ID'):
            normalized_metrics = Generate_Normalized_Metrics(group, temp_df2, ['TTN', 'TTB'], group_trait)
            normalized_metrics['Sample_ID']=sample_id
            normalized_metrics_list.append(normalized_metrics)

        temp_out = pd.concat(normalized_metrics_list, ignore_index=True)
        temp_df = temp_df.merge(temp_out, on=[group_trait, 'Sample_ID'])
    else:
        temp_df = input_df1.groupby(
            [group_trait, 'Sample_ID', 'Type'], as_index=False
        ).apply(Cal_Tumor_Size_simple, percentile_list, 'size')

    # Add cohort-specific relative metrics
    new_list = []
    for sample_id, group in temp_df.groupby('Sample_ID'):
        Add_Corhort_Specific_Relative_Metrics(group,group_trait)
        new_list.append(group)

    output_df = pd.concat(new_list, ignore_index=True)
    return output_df

def Generate_Normalized_Metrics(input_df1, input_df2, trait_list,group_trait='gRNA'):
    """
    This function normalizes input_df1 using metrics defined by trait_list based on input_df2.
    input_df1 is the experimental group, input_df2 is the control group.
    
    Parameters:
    input_df1 (pd.DataFrame): Experimental group data.
    input_df2 (pd.DataFrame): Control group data.
    trait_list (list): List of traits/metrics to normalize.
    
    Returns:
    pd.DataFrame: DataFrame with normalized metrics.
    """
    # Calculate the sum for each trait in both dataframes
    dict1 = {trait: input_df1[trait].sum() for trait in trait_list}
    dict2 = {trait: input_df2[trait].sum() for trait in trait_list}

    # Set the index to 'gRNA' for both dataframes
    temp1 = input_df1.set_index(group_trait)
    temp2 = input_df2.set_index(group_trait)
    
    # Ensure temp2 only contains rows that are also in temp1
    temp2 = temp2.loc[temp1.index]
    
    # Initialize the output DataFrame
    temp_output_df = pd.DataFrame({group_trait: temp1.index.values})
    
    # Normalize the metrics
    for trait in trait_list:
        normalized_trait = trait + '_normalized'
        temp_output_df[normalized_trait] = np.array(temp1[trait].to_list())/np.array(temp2[trait].to_list())*dict2[trait]/dict1[trait]
    
    return temp_output_df

def Add_Corhort_Specific_Relative_Metrics(input_df,group_trait='gRNA'):
    # Add relative metrics for LN_mean, GEO_mean etc using the median of inert
    temp_sub = input_df[input_df['Type']=='Inert']
    for temp_cname in input_df.drop(columns=[group_trait,'Type','Sample_ID'],inplace = False).columns:
        temp_name = temp_cname+'_relative'
        input_df[temp_name] = input_df[temp_cname]/temp_sub[temp_cname].median()

################
def Cal_Tumor_Size_simple(x,input_percentile,mode='None'):
    d = {}
    temp_vect = x['Cell_number']
    if type (temp_vect) == 'int':
        temp_vect = [temp_vect]
    # measure size
    d['LN_mean'] = LN_Mean(temp_vect)
    d['Geo_mean'] = Geometric_Mean(temp_vect)
    Percentile_list = list(np.percentile(temp_vect,input_percentile))
    for c,y in enumerate(input_percentile):
        temp_name = str(y)+'_percentile'
        d[temp_name] = Percentile_list[c]
    # measure number and burden
    if mode != 'size':
        d['TTN'] = len(temp_vect) # this is total tumor number
        d['TTB'] = sum(temp_vect)
    return pd.Series(d, index=list(d.keys())) 

def Cal_Tumor_Size_Cas9_negative(x):
    d = {}
    temp_vect = x['Cell_number']
    if type (temp_vect) == 'int':
        temp_vect = [temp_vect]
    d['TTN'] = len(temp_vect) # this is total tumor number
    d['TTB'] = sum(temp_vect)
    return pd.Series(d, index=list(d.keys())) 

def LN_Mean(input_vector):
    log_vector = np.log(input_vector)
    temp_mean = log_vector.mean()
    temp_var = log_vector.var()
    if len(log_vector)==1:
        temp_var = 0 # if only one clonal
    return (math.exp(temp_mean + 0.5*temp_var))

# calculate the Geometric mean from a vector of number
def Geometric_Mean(input_vector):
    log_vector = np.log(input_vector)
    temp_mean = log_vector.mean()
    return (math.exp(temp_mean))


def fdr(p_vals):
    p = np.asfarray(p_vals) # make input as float array
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    p = p[by_descend] # sort pvalue from small to large
    ranked_p_values = rankdata(p,method ='max') # this max is very important, when identical, use largest
    fdr = p * len(p) / ranked_p_values
    fdr = np.minimum(1, np.minimum.accumulate(fdr))

    return fdr[by_orig]



