import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import nbinom
import re
import math
import arviz as az
from collections import OrderedDict
from collections import defaultdict
import json


###########
def get_common_columns(df1: pd.DataFrame, df2: pd.DataFrame) -> set:
    """
    Returns a set of column names that are common to both df1 and df2.
    """
    return set(df1.columns).intersection(df2.columns)

# ##########
# Plotting

#Shapley waterfall
def plot_shapley_waterfall(dict_shapley, title=None, show_ci=True, group_by_sign=True,
                                 total_label="Total baseline adjustment", figsize=(8, 6),
                                 annotate_values=False):
    """
    dict_shapley format:
      {"var": {"mean":..., "q2p5":..., "q97p5":...}, ...}
    """

    # --- Unpack ---
    rows = []
    for var, s in dict_shapley.items():
        mean = float(s["mean"])
        lo   = float(s["q2p5"])
        hi   = float(s["q97p5"])
        rows.append((var, mean, lo, hi))

    # --- Sort / group ---
    if group_by_sign:
        pos = [r for r in rows if r[1] >= 0]
        neg = [r for r in rows if r[1] < 0]
        pos.sort(key=lambda r: abs(r[1]), reverse=True)
        neg.sort(key=lambda r: abs(r[1]), reverse=True)
        rows = pos + neg
    else:
        rows.sort(key=lambda r: abs(r[1]), reverse=True)

    names = [r[0] for r in rows]
    means = np.array([r[1] for r in rows], dtype=float)
    lows  = np.array([r[2] for r in rows], dtype=float)
    highs = np.array([r[3] for r in rows], dtype=float)

    # --- Waterfall cumulative starts/ends ---
    starts = np.concatenate(([0.0], np.cumsum(means)[:-1]))
    ends   = starts + means
    total  = means.sum()

    # --- Plot ---
    fig, ax = plt.subplots(figsize=figsize)

    # sign mask
    is_pos = means >= 0

    # Draw bars (use default cycle colors via two separate calls)
    for i, (s, e, m, pos_flag) in enumerate(zip(starts, ends, means, is_pos)):
        left  = min(s, e)
        width = abs(e - s)

        ax.barh(i, width, left=left, height=0.7)

        # CI for the increment, anchored at start: [s+lo, s+hi]
        if show_ci:
            ci_lo = s + lows[i]
            ci_hi = s + highs[i]
            ax.plot([ci_lo, ci_hi], [i, i], linewidth=2)
            ax.plot([ci_lo, ci_lo], [i - 0.15, i + 0.15], linewidth=2)
            ax.plot([ci_hi, ci_hi], [i - 0.15, i + 0.15], linewidth=2)

        if annotate_values:
            ax.text(e, i, f" {m:+.3f}", va="center", ha="left")

    # Connectors
    for i in range(len(names) - 1):
        ax.plot([ends[i], ends[i]], [i + 0.35, i + 0.65], linewidth=1)

    # Total bar
    y_total = len(names)
    ax.barh(y_total, abs(total), left=min(0, total), height=0.7)
    if annotate_values:
        ax.text(total, y_total, f" {total:+.3f}", va="center", ha="left")

    # Labels
    ax.set_yticks(range(len(names) + 1))
    ax.set_yticklabels(names + [total_label])

    # Bold zero line
    ax.axvline(0, linewidth=2)

    plt.xticks(rotation=45)

    ax.set_xlabel("Contribution to baseline-adjustment component")
    if title:
        ax.set_title(title)

    plt.tight_layout()
    plt.show()


def posterior_plot_paper(
    samples,
    label,
    units=None,
    hdi_prob=0.95,
    decimals=2,
    figsize=(3.4, 2.2),
    dpi=300,
    font_base=9,
    legend_outside=True,
):
    # Optional: keep ArviZ style (if your ArviZ version has styles)
    # Comment out if it errors on your setup.
    try:
        az.style.use("arviz-whitegrid")
    except Exception:
        pass

    x = np.asarray(samples)
    if x.ndim == 1:
        x_plot = x[None, :]
    elif x.ndim == 2:
        x_plot = x
    else:
        # scalar plot: flatten all draws
        x_plot = x.reshape(1, -1)

    x_flat = x_plot.ravel()

    lo, hi = az.hdi(x_flat, hdi_prob=hdi_prob)
    med = np.median(x_flat)

    # Fonts: small + consistent
    plt.rcParams.update({
        "font.size": font_base,
        "axes.labelsize": font_base,
        "xtick.labelsize": font_base - 1,
        "ytick.labelsize": font_base - 1,
        "legend.fontsize": font_base - 1,
    })

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Draw the nice ArviZ posterior (but no point estimate annotation)
    az.plot_posterior(
        {"_": x_plot},
        hdi_prob=hdi_prob,
        point_estimate=None,     # avoids ArviZ adding median/mean text/line
        ax=ax
    )

    # Remove any ArviZ-added text blocks (varies by version)
    for t in list(ax.texts):
        t.remove()
    ax.set_title("")
    ax.set_ylabel("")

    # Thin, unobtrusive reference lines (paper-ish)
    vline_kw = dict(linewidth=0.8, alpha=0.9, zorder=5, solid_capstyle="butt")
    ax.axvline(med, **vline_kw, label=f"Median {med:.{decimals}f}")
    ax.axvline(lo, linestyle="--", **vline_kw,
               label=f"{int(hdi_prob*100)}% HDI [{lo:.{decimals}f}, {hi:.{decimals}f}]")
    ax.axvline(hi, linestyle="--", **vline_kw, label="_nolegend_")  # no duplicate entry

    # X label
    xlabel = label + (f" ({units})" if units else "")
    ax.set_xlabel(xlabel)

    # Legend placement that won't overlap the density
    if legend_outside:
        ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
        fig.tight_layout(rect=(0, 0, 0.80, 1))  # leave room on the right for legend
    else:
        ax.legend(frameon=False, loc="upper right")
        fig.tight_layout()

    # Keep the clean ArviZ vibe
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return fig, ax
###########
#Wrangling
def data_wrangling(df, target):

    df = df.rename(columns={"Koko_MRI":"Size_MRI_preop", "Koko_UA":"Size_US_preop", "Multifokaalinen":"Multifocal_preop", "Imusolmukesytologia":"Node_cytology_preop"})
    
    # Creating age groups
    if(target=="detailed"):
        bins = [0, 40, 50, 60, 70, 120]  # Define the bins
        labels = ['Under 40', '41-50', '51-60', '61-70', '71+']  # Define the group labels
        df['Age_group'] = pd.cut(df['Age'], bins=bins, labels=labels)
    elif(target in ["Helen_006"]):
        # Creating age groups around RCT median
        bins = [0, 50, 70, 120]  # Define the bins
        labels = ['Under 50', '50 - 70', 'Over 70']  # Define the group labels
    
        df['Age_group_med'] = pd.cut(df['Age'], bins=bins, labels=labels)
    elif(target in ["Neosphere"]):
        # Creating age groups around RCT median
        bins = [0, 50, 77, 120]  # Define the bins
        labels = ['Under 50', '50 - 77', 'Over 77']  # Define the group labels
    
        df['Age_group_med'] = pd.cut(df['Age'], bins=bins, labels=labels)
    elif(target in ["Tryphaena", "Kristine"]):
        # Creating age groups around RCT median
        bins = [0, 49, 120]  # Define the bins
        labels = ['Under 49', 'Over 49']  # Define the group labels
    
        df['Age_group_med'] = pd.cut(df['Age'], bins=bins, labels=labels)
    
    # Creating ER and PR groups
    ##Data imputation to empty values
    df.loc[(df["ER_pnb_procent_preop_max"].isnull()) & (df["ER_pnb_status_preop"]=="Estrogeenireseptori negatiivinen"), "ER_pnb_procent_preop_max"] = 0
    df.loc[(df["PR_pnb_procent_preop_max"].isnull()) & (df["PR_pnb_status_preop"]=="Proestrogeenireseptori negatiivinen"), "PR_pnb_procent_preop_max"] = 0
    
    ##For ER*PR interaction term, clear the mixed types
    df.loc[df["ER_pnb_status_preop"].str.contains("positiivinen") , "ER_pnb_status_preop"] = "Estrogeenireseptori positiivinen"
    df.loc[df["PR_pnb_status_preop"].str.contains("positiivinen") , "PR_pnb_status_preop"] = "Proestrogeenireseptori positiivinen"
    
    df["ER_PR_status_inter"] = df["ER_pnb_status_preop"] + ", " + df["PR_pnb_status_preop"]


    #Hormone status from RCTs

    cond_er_pos = (df["ER_pnb_status_preop"].str.contains("positiivinen"))
    cond_pr_pos = (df["PR_pnb_status_preop"].str.contains("positiivinen") )
    
    df['hormone_status'] = pd.Series(index=df.index)
    df.loc[cond_er_pos | cond_pr_pos, 'hormone_status'] = "Positive"
    df.loc[(~cond_er_pos) | (~cond_pr_pos), 'hormone_status'] = "Negative"
    
    if(target=="detailed"):
        bins = [0, 10, 75, 100]  # Define the bins
        labels = ['1-10', '11-75', '76-100']  # Define the group labels
        
        df['ER_group'] = pd.cut(df['ER_pnb_procent_preop_max'], bins=bins, labels=labels)
        df['ER_group'] = df['ER_group'].cat.add_categories("0")
        df.loc[df["ER_pnb_procent_preop_max"] == 0, "ER_group"] = "0"
        df['ER_group'] = df['ER_group'].cat.reorder_categories(["0"] + labels, ordered=True)
    
        df['PR_group'] = pd.cut(df['PR_pnb_procent_preop_max'], bins=bins, labels=labels)
        df['PR_group'] = df['PR_group'].cat.add_categories("0")
        df.loc[df["PR_pnb_procent_preop_max"] == 0, "PR_group"] = "0"
        df['PR_group'] = df['PR_group'].cat.reorder_categories(["0"] + labels, ordered=True)
    
    ##HER2 IHC groups
    #df["HER2_IHC_min_max"] = "Min " + df["HER2_IHC_pnb_grade_preop_min"] + ", Max " + df["HER2_IHC_pnb_grade_preop_max"]
    
    #Histology groups

    if(target=="NA"):
        conditions = [
            df['Histology_pnb_preop'].isin(['Carcinoma ductale', 'Carcinoma ductale, Carcinoma ductale in situ', 'Carcinoma ductale invasivum', 'Carcinoma ductale, Carcinoma inflammatoricum']),
            df['Histology_pnb_preop'].isin(['Carcinoma lobulare', 'Carcinoma lobulare, Carcinoma ductale in situ']),
            df['Histology_pnb_preop'].isin(['Carcinoma ductale, Carcinoma lobulare', 'Carcinoma ductale & Carcinoma lobulare', 'Carcinoma lobulare & Carcinoma ductale', 'Carcinoma ductale, Carcinoma lobulare, Carcinoma ductale in situ']),
            df['Histology_pnb_preop'].isin(['Carcinoma micropapillare', 'Carcinoma micropapillare, Carcinoma apocrinum', 'Carcinoma micropapillare,Carcinoma ductale', 'Carcinoma ductale, Carcinoma micropapillare'])
            ]
        choices = ['Carcinoma ductale', 'Carcinoma lobulare', 'Carcinoma ductale & Carcinoma lobulare', 'Carcinoma micropapillare', 'other']
        # Apply the conditions to create the new column
        df['histology'] = pd.Series(index=df.index)
        df.loc[conditions[0], 'histology'] = choices[0]
        df.loc[conditions[1], 'histology'] = choices[1]
        df.loc[conditions[2], 'histology'] = choices[2]
        df.loc[conditions[3], 'histology'] = choices[3]
        df.loc[((~conditions[0])&(~conditions[1])&(~conditions[2])&(~conditions[3])), 'histology'] = choices[4]
    elif(target in ["Helen_006", "detailed"]):
        conditions = [
            df['Histology_pnb_preop'].isin(['Carcinoma ductale', 'Carcinoma ductale, Carcinoma ductale in situ', 'Carcinoma ductale invasivum', 'Carcinoma ductale, Carcinoma inflammatoricum']),
            df['Histology_pnb_preop'].isin(['Carcinoma lobulare', 'Carcinoma lobulare, Carcinoma ductale in situ']),
            df['Histology_pnb_preop'].isin(['Carcinoma micropapillare', 'Carcinoma micropapillare, Carcinoma apocrinum', 'Carcinoma micropapillare,Carcinoma ductale', 'Carcinoma ductale, Carcinoma micropapillare','Carcinoma ductale, Carcinoma lobulare', 'Carcinoma ductale & Carcinoma lobulare', 'Carcinoma lobulare & Carcinoma ductale', 'Carcinoma ductale, Carcinoma lobulare, Carcinoma ductale in situ'])
            ]
        choices = ['Carcinoma ductale', 'Carcinoma lobulare', 'other']
        # Apply the conditions to create the new column
        df['histology'] = pd.Series(index=df.index)
        df.loc[conditions[0], 'histology'] = choices[0]
        df.loc[conditions[1], 'histology'] = choices[1]
        df.loc[((~conditions[0])&(~conditions[1])), 'histology'] = choices[2]
    else:
        conditions = [
            df['Histology_pnb_preop'].isin(['Carcinoma ductale', 'Carcinoma ductale, Carcinoma ductale in situ', 'Carcinoma ductale invasivum', 'Carcinoma ductale, Carcinoma inflammatoricum']),
            df['Histology_pnb_preop'].isin(['Carcinoma lobulare', 'Carcinoma lobulare, Carcinoma ductale in situ']),
            df['Histology_pnb_preop'].isin(['Carcinoma ductale, Carcinoma lobulare', 'Carcinoma ductale & Carcinoma lobulare', 'Carcinoma lobulare & Carcinoma ductale', 'Carcinoma ductale, Carcinoma lobulare, Carcinoma ductale in situ']),
            df['Histology_pnb_preop'].isin(['Carcinoma micropapillare', 'Carcinoma micropapillare, Carcinoma apocrinum', 'Carcinoma micropapillare,Carcinoma ductale', 'Carcinoma ductale, Carcinoma micropapillare'])
            ]
        choices = ['Carcinoma ductale', 'Carcinoma lobulare', 'Carcinoma ductale & Carcinoma lobulare', 'other']
        # Apply the conditions to create the new column
        df['histology'] = pd.Series(index=df.index)
        df.loc[conditions[0], 'histology'] = choices[0]
        df.loc[conditions[1], 'histology'] = choices[1]
        df.loc[conditions[2], 'histology'] = choices[2]
        df.loc[((~conditions[0])&(~conditions[1])&(~conditions[2])), 'histology'] = choices[3]
    
    
    if(target=="detailed"):
        #MIB-1 groups
        bins = [0, 5, 15, 30, 60, 100]  # Define the bins
        labels = ['1-5', '6-15', '15-30', '30-60','60-100']  # Define the group labels
        df['mib_group'] = pd.cut(df['MIB1_pnb_procent_preop_max'], bins=bins, labels=labels)

    elif(target in ["Helen_006"]):
        #MIB-1 groups RCT comparision
        bins = [0,30, 100]  # Define the bins
        labels = ['under 30%', 'over 30%']  # Define the group labels
        df['mib_group_30'] = pd.cut(df['MIB1_pnb_procent_preop_max'], bins=bins, labels=labels)
    
    #Tumor grade
    if(target  in ["Helen_006", "Neosphere"]):
        df.loc[df["Gradus_preop"].isna(), "Gradus_preop"] = "N/A"
        df.loc[df["Gradus_preop"].isin([1.0, 2.0]), "Gradus_preop"] = "1 or 2"
        df.loc[df["Gradus_preop"].isin([3.0]), "Gradus_preop"] = "3"

    #Multifocality
    df.loc[df["Multifocal_preop"].isna(), "Multifocal_preop"] = "N/A"

    #Node status
    df.loc[df["Node_cytology_preop"].isna(), "Node_cytology_preop"] = "N/A"

    #Tumor size mm
    df["Tumor_size_preop"] = df["Size_MRI_preop"].combine_first(df["Size_US_preop"])

    #Tumor size classes (we are not able to recognize T4)
    bins = [0, 20, 50, df["Tumor_size_preop"].max()]  # Define the bins
    labels = ['T1', 'T2', 'T3']  # Define the group labels
    df['Tumor_size_class_preop'] = pd.cut(df['Tumor_size_preop'], bins=bins, labels=labels).cat.add_categories("N/A")
    df.loc[df["Tumor_size_class_preop"].isna(), "Tumor_size_class_preop"] = "N/A"

    if(target=="detailed"):
        #Size, multifocality interconnection
        #df["Size_multifocal_status_inter"] = df["Tumor_size_class_preop"].astype(str) + ", " + df["Multifocal_preop"]
        #df.loc[df["Size_multifocal_status_inter"].str.contains("N/A"), "Size_multifocal_status_inter"] = "N/A"
        conditions = [
            ((df['Tumor_size_class_preop'].isin(['T1', 'T2'])) & (df["Multifocal_preop"]=="T")),
            ((df['Tumor_size_class_preop'].isin(['T3'])) & (df["Multifocal_preop"]=="T")),
            ((df['Tumor_size_class_preop'].isin(['T1', 'T2'])) & (df["Multifocal_preop"]=="F")),
            ((df['Tumor_size_class_preop'].isin(['T3'])) & (df["Multifocal_preop"]=="F"))
        
            ]
        choices = ['T1-T2 multifocal', 'T3 multifocal', 'T1-T2 unifocal', 'T3 unifocal', 'other']
        
        # Apply the conditions to create the new column
        df['Size_multifocal_status_inter'] = pd.Series(index=df.index)
        df.loc[conditions[0], 'Size_multifocal_status_inter'] = choices[0]
        df.loc[conditions[1], 'Size_multifocal_status_inter'] = choices[1]
        df.loc[conditions[2], 'Size_multifocal_status_inter'] = choices[2]
        df.loc[conditions[3], 'Size_multifocal_status_inter'] = choices[3]
        df.loc[((~conditions[0])&(~conditions[1])&(~conditions[2])&(~conditions[3])), 'Size_multifocal_status_inter'] = choices[4]
        
        

    elif(target in ["Helen_006","Neosphere"]):
        conditions = [
            (df['Tumor_size_class_preop'].isin(['T1', 'T2'])),
            (df['Tumor_size_class_preop'].isin(['T3']))
        
            ]
        choices = ['T1-T2', 'T3', 'other']
        
        # Apply the conditions to create the new column
        df['Tumor_size_status'] = pd.Series(index=df.index)
        df.loc[conditions[0], 'Tumor_size_status'] = choices[0]
        df.loc[conditions[1], 'Tumor_size_status'] = choices[1]
        df.loc[((~conditions[0])&(~conditions[1])), 'Tumor_size_status'] = choices[2]

    elif(target in ["Kristine"]):
        conditions = [
            (df['Tumor_size_class_preop'].isin(['T1'])),
            (df['Tumor_size_class_preop'].isin(['T2','T3']))
        
            ]
        choices = ['T1', 'T2-T3', 'N/A']
        
        # Apply the conditions to create the new column
        df['Tumor_size_status'] = pd.Series(index=df.index)
        df.loc[conditions[0], 'Tumor_size_status'] = choices[0]
        df.loc[conditions[1], 'Tumor_size_status'] = choices[1]
        df.loc[((~conditions[0])&(~conditions[1])), 'Tumor_size_status'] = choices[2]
        

    #ECOG
    df.loc[df["Ecog_preop"].isna(), "Ecog_preop"] = "N/A"

    if(target in ["Helen_006"]):
        df.loc[(df["Ecog_preop"].isin([1.0, 0.0])), "Ecog_preop"] = "0 or 1"
        df.loc[(df["Ecog_preop"].isin([2.0])), "Ecog_preop"] = "2"

    #IHC
    df.loc[df["HER2_IHC_pnb_grade_preop_max"].isin(["1+", "neg."]), "HER2_IHC_pnb_grade_preop_max"] = "1+ or neg."
    df.loc[df["HER2_IHC_pnb_grade_preop_max"].isna(), "HER2_IHC_pnb_grade_preop_max"] = "N/A"
    
    #ISH
    df.loc[df["CISH_pnb_preop"].isna(), "CISH_pnb_preop"] = "N/A"
    df.loc[df["CISH_pnb_preop"].str.contains("positiivinen"), "CISH_pnb_preop"] = "Positive"
    df.loc[(df["CISH_pnb_preop"]=="IHC c-erbB2 negatiivinen") | (df["CISH_pnb_preop"]=="c-erbB2 negatiivinen"), "CISH_pnb_preop"] = "Negative"

    if(target=="Helen_006"):
        
        conditions = [
            (df['HER2_IHC_pnb_grade_preop_max']==3),
            ((df['HER2_IHC_pnb_grade_preop_max']==2) & (df["CISH_pnb_preop"]=="Positive"))        
            ]

        choices = ['IHC3+', 'IHC2+ and ISH+', 'other']
        
        # Apply the conditions to create the new column
        df['IHC_ISH_inter'] = pd.Series(index=df.index)
        df.loc[conditions[0], 'IHC_ISH_inter'] = choices[0]
        df.loc[conditions[1], 'IHC_ISH_inter'] = choices[1]
        df.loc[((~conditions[0])&(~conditions[1])), 'IHC_ISH_inter'] = choices[2]

    #Bilateral
    df.loc[df["Bilateral"].isna(), "Bilateral"] = "False"
    df.loc[df["Bilateral"]==True, "Bilateral"] = "True"

    #StageII-III
    conditions = [
            (((df['Tumor_size_class_preop'].isin(['T1'])) & (df["Node_cytology_preop"]=="pos")) | ((df['Tumor_size_class_preop'].isin(['T2','T3'])))),       
            ]
    choices = ['II-III', 'other']
    
    # Apply the conditions to create the new column
    df['Stage_preop'] = pd.Series(index=df.index)
    df.loc[conditions[0], 'Stage_preop'] = choices[0]
    df.loc[((~conditions[0])), 'Stage_preop'] = choices[1]

    #biopsy nodal positive
    conditions = [
            (((df['Imusolmukenayte'].isin(['PNB', 'ONB'])) & (df["Node_cytology_preop"]=="pos"))),       
            ]
    choices = ['True', 'False']
    
    # Apply the conditions to create the new column
    df['biopsy_nodal_positive'] = pd.Series(index=df.index)
    df.loc[conditions[0], 'biopsy_nodal_positive'] = choices[0]
    df.loc[((~conditions[0])), 'biopsy_nodal_positive'] = choices[1]
    
    return df

    

# ##########
# Diagnostics

# Function to bootstrap confidence intervals
def bootstrap_ci(data, num_samples=1000, ci=95):
    boot_means = []
    for _ in range(num_samples):
        sample = np.random.choice(data, size=len(data), replace=True)
        boot_means.append(np.mean(sample))
    lower_bound = np.percentile(boot_means, (100 - ci) / 2)
    upper_bound = np.percentile(boot_means, 100 - (100 - ci) / 2)
    return lower_bound, upper_bound

def binned_residual_plot(df, x, categorical=False, y="residual", obs="y", n_bins=40, confidence=0.95):
    
    binned_residuals = df.copy()
        
    if not categorical:
        binned_residuals.sort_values(x)
        binned_residuals['bin'] = pd.qcut(binned_residuals[x], n_bins, labels=False, duplicates='drop')  # Equal-sized bins
                
        # Group by "bin" and calculate the mean and confidence intervals
        results = []
        for bin_label, group in binned_residuals.groupby('bin'):
            bin_mean = group[x].mean()
            residual_mean = group[y].mean()
            ci_lower, ci_upper = bootstrap_ci(group[y])
            results.append({
                'bin': bin_label,
                'bin_mean': bin_mean,
                'residual_mean':residual_mean,
                'yerr': ci_upper - ci_lower
            })
        
        # Convert results to a DataFrame
        results_df = pd.DataFrame(results).sort_values("bin_mean")
            
        # Scatter plot with confidence intervals
        plt.figure(figsize=(8, 6))
    
        plt.errorbar(
            results_df['bin_mean'], 
            results_df['residual_mean'], 
            yerr=results_df['yerr'], 
            fmt='o', color='blue', ecolor='gray', elinewidth=1.5, capsize=4, label="Binned Averages with CI"
                  )
    else:
        # Scatter plot with confidence intervals
        plt.figure(figsize=(8, 6))
        sns.barplot(x=x, y=y, data=binned_residuals, errorbar=("ci", 95))

        
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.8, label="Residual = 0")
    plt.title("Binned Scatterplot of Residual vs. " + x + " with Confidence Intervals")
    plt.xlabel("Average " + x + " (Binned)")
    plt.ylabel("Average Residual (Binned)")
    
    if categorical:
        plt.xticks(rotation=45)
    
    plt.legend()
    plt.grid()
    plt.show()




# #Real coefficients and posteriors

def _canon(s: str) -> str:
    """Canonicalize a name: lowercase, non-alnum -> '_', collapse/restrip underscores."""
    s = s.strip().lower()
    s = re.sub(r"[^0-9a-z]+", "_", s)
    s = re.sub(r"_+", "_", s)
    return s.strip("_")

def plot_betas_from_tuple_truth(idata, truth_dict, var_prefix="beta_", ncols=3):
    """
    Plot posterior betas for factor levels using a truth dict keyed by (var, lvl) tuples.
    
    Parameters
    ----------
    idata : arviz.InferenceData
        Posterior from PyMC.
    truth_dict : dict
        {(var, lvl): true_value, ...}
        Example: {("Age group med","Under 50"): 0.03, ("histology","Carcinoma ductale"): 0.47}
    var_prefix : str
        Prefix of the beta variables in idata (default "beta_").
    ncols : int
        Columns in the subplot grid.
    """
    # 1) Index all beta_* vectors in idata by a canonical factor name
    beta_vars = [v for v in idata.posterior.data_vars if v.startswith(var_prefix)]
    if not beta_vars:
        raise ValueError(f"No variables starting with '{var_prefix}' in idata.posterior")

    index = {}  # canon_factor -> {var, dim, levels, canon_levels}
    for v in beta_vars:
        da = idata.posterior[v]
        dims = [d for d in da.dims if d not in ("chain", "draw")]
        if len(dims) != 1:
            # Skip non-1D betas; adapt here if you have interactions etc.
            continue
        dim = dims[0]
        levels = list(map(str, da.coords[dim].values))
        canon_factor = _canon(v[len(var_prefix):])  # strip prefix and canonicalize
        index[canon_factor] = {
            "var": v,
            "dim": dim,
            "levels": levels,
            "canon_levels": [_canon(lvl) for lvl in levels],
        }

    # 2) Build the list of plots by matching (var, lvl) tuples
    to_plot = []  # (var_name, dim, lvl_text, truth_val)
    missing = []
    for (var, lvl), truth_val in truth_dict.items():
        canon_var = _canon(str(var))
        if canon_var not in index:
            missing.append((var, lvl, "factor-not-found"))
            continue
        rec = index[canon_var]
        canon_lvl = _canon(str(lvl))
        try:
            j = rec["canon_levels"].index(canon_lvl)
        except ValueError:
            missing.append((var, lvl, "level-not-found"))
            continue
        lvl_text = rec["levels"][j]  # original pretty label from idata coords
        to_plot.append((rec["var"], rec["dim"], lvl_text, float(truth_val)))

    if not to_plot:
        msg = "No matches between truth_dict (tuple keys) and idata betas."
        if missing:
            msg += f"\nFirst few misses: {missing[:5]}"
        raise ValueError(msg)

    # 3) Plot in a grid
    n = len(to_plot)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3.2*nrows))
    axes = np.array(axes).reshape(-1)  # safe even if nrows==1 or ncols==1

    for ax, (var_name, dim, lvl_text, truth_val) in zip(axes, to_plot):
        az.plot_posterior(
            idata,
            var_names=[var_name],
            coords={dim: [lvl_text]},
            point_estimate="mean",  # or "median"
            hdi_prob=0.95,
            ax=ax,
        )
        ax.axvline(truth_val, linestyle="--", linewidth=2, color="red")
        ax.set_title(f"{var_name}[{lvl_text}]\ntrue = {truth_val:+.3f}", fontsize=10)

    # Remove any unused axes
    for j in range(len(to_plot), len(axes)):
        fig.delaxes(axes[j])

    fig.suptitle("Posterior betas with true values (red dashed)", y=1.02, fontsize=13)
    fig.tight_layout()

    # Optional: report unmatched items for convenience
    if missing:
        print("Unmatched truth entries (showing up to 10):")
        for m in missing[:10]:
            print("  -", m)

    return fig

# ---------- Example usage ----------
# truth_tuple = {
#    ("Age group med", "Under 50"):  0.03,
#    ("Age group med", "Over 50"):  -0.03,
#    ("histology",     "Carcinoma ductale"):  0.47,
#    ("histology",     "carcinoma lobulare"): -0.15,
#    ("histology",     "other"):              -0.32,
#    # ...
# }
# fig = plot_betas_from_tuple_truth(idata, truth_tuple)
# plt.show()


#Plotting Target marginal posteriors
def plot_posterior_margins_vs_target_and_original(
    idata,
    L_blocks,
    N_target_marg_blocks,
    T_target,
    q_star,
):
    """
    Compare posterior-implied margins to:
      1) target margins (red dashed line)
      2) original sample-based margins (blue dashed line), scaled to T_target.

    Parameters
    ----------
    idata : arviz.InferenceData
        Posterior from PyMC, must contain "p_comp".
    L_blocks : list of 2D numpy arrays
        Each L_b has shape (K_b, C) mapping joint cells to margin categories.
    N_target_marg_blocks : list of 1D numpy arrays
        Each N_b has shape (K_b,), target marginal counts for block b.
    T_target : float or int
        Total population count used in Nhat = T_target * p_comp.
    q_star : 1D numpy array
        Original sample joint proportions over C cells (sum ~ 1).
        Used to compute original implied margins: (L_b @ q_star) * T_target.
    """

    # 1) Extract posterior samples for p_comp and convert to implied Nhat per cell
    p_post = idata.posterior["p_comp"].stack(sample=("chain", "draw")).values  # (C, S)
    Nhat_post = T_target * p_post  # (C, S)
    C, S = Nhat_post.shape

    # Basic sanity checks
    assert len(L_blocks) == len(N_target_marg_blocks), "Mismatch in number of blocks"
    assert q_star.shape[0] == C, "q_star length must match number of cells C"

    # 2) Loop over margin blocks
    for b, (L_b, N_b) in enumerate(zip(L_blocks, N_target_marg_blocks)):
        # L_b: (K_b, C), N_b: (K_b,)
        K_b = N_b.shape[0]

        # Posterior-implied margins: (K_b, S)
        Nhat_b_post = L_b @ Nhat_post

        # Original sample-based margins scaled to T_target: (K_b,)
        N_orig_b = (L_b @ q_star) * T_target

        print(f"\n=== Margin block {b} ===")
        print("Target marginals:   ", N_b)
        print("Original marginals: ", N_orig_b.round(1))

        # 3) Plot per category
        for k in range(K_b):
            samples = Nhat_b_post[k, :]      # S samples
            target_val = N_b[k]
            orig_val   = N_orig_b[k]

            fig, ax = plt.subplots(figsize=(8, 5))
            az.plot_posterior(
                samples,
                hdi_prob=0.95,
                point_estimate="mean",
                ax=ax,
            )
            # target (red)
            ax.axvline(target_val, color="red", linestyle="--", linewidth=2, label="Target margin")
            # original (blue)
            ax.axvline(orig_val, color="blue", linestyle="--", linewidth=2, label="Original (q⋅T_target)")

            ax.set_title(f"Block {b}, category {k}\n"
                         f"target={target_val:.1f}, original={orig_val:.1f}",
                         fontsize=10)
            ax.set_xlabel("Implied margin count Nhat_b[k]")
            ax.legend()
            plt.tight_layout()
            plt.show()

# ##########
# Raking methods

def build_loading_matrix_from_data(data, raking_vars, margin_formulas, marg_dict):
    """
    Build a loading matrix L and a poststratification contingency table
    directly from row-level data (e.g. individuals).

    Parameters
    ----------
    data : pd.DataFrame
        Raw row-level dataset with one row per individual and
        categorical raking variables as columns.

    raking_vars : list of str
        Names of the variables used for raking (e.g. ["age", "gender", "region"]).

    margin_formulas : list of lists
        Each element defines one margin as a list of variable names.
        Example: [["age"], ["gender"], ["region"], ["age", "gender"]]

    Returns
    -------
    L : np.ndarray
        Loading matrix of shape (D, J) mapping cell counts to marginal totals.
        D = total number of unique margin categories.
        J = number of poststratification cells.

    poststrat_table : pd.DataFrame
        DataFrame with all unique combinations of raking_vars and
        a 'Freq' column (cell counts).

    row_labels : list of str
        Names of margin categories (rows of L).

    col_labels : list of str
        Names of poststratification cells (columns of L).
    """
    # ------------------------------------------------------------------
    # 1️⃣ Build the poststratification contingency table
    # ------------------------------------------------------------------

    data_cp = data.copy()

    #let's remove the rows that do not exist in the target population
    for key, value in marg_dict.items():
        var, level = key.split(":", 1)

        if value==0:
            data_cp = data_cp[data_cp[var]!=level]
        
    poststrat_table = (
        data_cp[raking_vars]
        .astype("category")
        .apply(lambda x: x.cat.remove_unused_categories())
        .groupby(raking_vars)
        .size()
        .reset_index(name="Freq")
        .sort_values(raking_vars)
        .reset_index(drop=True)
    )

    # Label each cell combination
    col_labels = (
        poststrat_table[raking_vars].astype(str).agg(':'.join, axis=1).tolist()
    )
    J = len(col_labels)

    # ------------------------------------------------------------------
    # 2️⃣ Build loading matrix L
    # ------------------------------------------------------------------
    L_rows = []
    row_labels = []
    L_blocks = []

    for margin_vars in margin_formulas:
        # Group by each margin (e.g. ["age"], ["gender"], ["age", "gender"])
        grouped = poststrat_table.groupby(margin_vars).groups

        #One block
        L_b = []

        for margin_values, indices in grouped.items():
            row = np.zeros(J, dtype=int)
            row[list(indices)] = 1  # mark which cells belong to this margin

            # Human-readable margin label
            if isinstance(margin_values, tuple):
                label_parts = [f"{v}{val}" for v, val in zip(margin_vars, margin_values)]
                label = ":".join(label_parts)
            else:
                label = f"{margin_vars[0]}:{margin_values}"

            L_rows.append(row)
            row_labels.append(label)
            L_b.append(row)
            
        block = np.vstack(L_b)
        #Block-wise loading matrix
        L_blocks.append(block)

    #The whole loading matrix
    L = np.vstack(L_rows)

    return L, poststrat_table, row_labels, col_labels, L_blocks

#The following funtions are used to insert 0 values to the zero marginals that would not otherwise exist
def missing_label_indices(row_labels, row_labels_A):
    """
    Return the indices from row_labels whose values do not appear in row_labels_A.
    """
    set_A = set(row_labels_A)
    
    return [i for i, lbl in enumerate(row_labels) if lbl not in set_A]

def missing_label_indices(row_labels, row_labels_A):
    """
    Return the indices from row_labels whose values do not appear in row_labels_A.
    """
    set_A = set(row_labels_A)
    return [i for i, lbl in enumerate(row_labels) if lbl not in set_A]

def insert_missing_in_blocks(N_target_marg_blocks, missing_indices):
    """
    Insert zeros into an array-of-arrays structure (list of 1D arrays)
    at the positions given by missing_indices (long-format indexing).

    Parameters
    ----------
    N_target_marg_blocks : list of np.array
        Original block-wise margin arrays.
    missing_indices : list of int
        Long-format indices where new zeros should be inserted.

    Returns
    -------
    updated_blocks : list of np.array
        New blocks with zeros inserted at correct positions.
    """

    # Sort missing indices to handle in ascending order
    missing_indices = sorted(missing_indices)

    updated_blocks = []
    offset = 0
    mi_ptr = 0  # pointer into missing_indices

    for block in N_target_marg_blocks:
        block_len = len(block)
        new_block = block.tolist()  # work in Python list (easy insert)

        # Where does the next missing index fall relative to this block?
        while (
            mi_ptr < len(missing_indices)
            and missing_indices[mi_ptr] < offset + block_len
        ):
            # Local index inside this block:
            local_idx = missing_indices[mi_ptr] - offset

            # Insert zero BEFORE the element currently at local_idx
            new_block.insert(local_idx, 0)

            # Adjust block length and move to next missing index
            block_len += 1
            mi_ptr += 1

        updated_blocks.append(np.array(new_block))
        offset += block_len

    # If missing indices go beyond last block, append a new block
    while mi_ptr < len(missing_indices):
        # All remaining missing entries belong to a new block
        updated_blocks.append(np.array([0]))
        mi_ptr += 1

    return updated_blocks

# This function is needed to transfer the RCT marginal dict in to a block-list format

def marg_dict_to_blocks(marg_dict):
    result = []
    grouped = OrderedDict()
    
    for key, value in marg_dict.items():
        var, level = key.split(":", 1)
        if var not in grouped:
            grouped[var] = []
        grouped[var].append(value)
    
    result = list(grouped.values())

    return result

##################################
#MAIC analysis helper functions
def marg_dict_to_specs(
    dict_marg,
    *,
    drop="last",          # "last", "first", or explicit dict {var: level}
    sort_levels=True
):
    """
    Convert a dict of the form {'var:level': count, ...}
    into MAIC/entropy-balancing specs.

    Parameters
    ----------
    dict_marg : dict
        Keys like 'variable:level', values are counts.
    drop : {"last", "first", dict}
        Which level to drop (reference category) per variable.
    sort_levels : bool
        Whether to sort levels alphabetically for reproducibility.

    Returns
    -------
    specs : dict
        Compatible with build_marginal_frequency_constraints().
    """

    # 1) Parse keys into variable -> {level: count}
    by_var = defaultdict(dict)

    for k, v in dict_marg.items():
        if ":" not in k:
            raise ValueError(f"Key '{k}' must be of the form 'variable:level'")
        var, level = k.split(":", 1)
        by_var[var][level] = float(v)

    # 2) Build specs
    specs = {}

    for var, level_counts in by_var.items():
        levels = list(level_counts.keys())
        if sort_levels:
            levels = sorted(levels)

        counts = [level_counts[lv] for lv in levels]
        total = sum(counts)

        if total <= 0:
            raise ValueError(f"Total count for variable '{var}' must be positive")

        # convert counts -> proportions
        target_props = {lv: c / total for lv, c in zip(levels, counts)}

        # choose dropped level
        if isinstance(drop, dict):
            drop_level = drop[var]
        elif drop == "first":
            drop_level = levels[0]
        elif drop == "last":
            drop_level = levels[-1]
        else:
            raise ValueError("drop must be 'first', 'last', or a dict")

        specs[var] = {
            "type": "categorical",
            "levels": levels,
            "target": target_props,
            "drop": drop_level,
        }

    return specs

def remove_zero_marginals_MAIC(data_org, target_N, marg_dict):
    data_new = data_org.copy()
    marg_dict_new = {}

    #let's remove the rows that do not exist in the target population
    for key, value in marg_dict.items():
        var, level = key.split(":", 1)

        if value==0:
            data_new = data_new[data_new[var]!=level]
        elif value==target_N:
            continue
        else:
          marg_dict_new[key] = (value/target_N)

    return data_new, marg_dict_new
