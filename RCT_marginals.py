# RCT Marginals
##############################
# # Tryphaena group B 

#Tryphaena group B N
N_Tryp = 75

#Marginals
age_under_49_Tryp = 38
age_over_49_Tryp = 37

ecog_0_Tryp = 66
ecog_1_Tryp = 9
ecog_na_Tryp = 0

grade_1_Tryp = 2
grade_2_Tryp = 34
grade_3_Tryp = 26
grade_NA_Tryp = 13

ER_or_PR_pos_Tryp = 35
ER_and_PR_neg_Tryp = 40

her2_IHC_0_or_1plus_Tryp = 0
her2_IHC_2plus_Tryp = 1
her2_IHC_3plus_Tryp = 74

her2_ISH_pos_Tryp = 69
her2_ISH_neg_Tryp = 1
her2_ISH_NA_Tryp = 5

tumor_size_T1_T2_Tryp = 38
tumor_size_T3_Tryp = 37

##############################
# # Helen-006 group A
# Nab-paclitaxel, trastuzumab, and pertuzumab group

N_helenA = 332
N_pCR_helenA = 220

#Marginals
age_under_50_helenA = 166
age_over_50_helenA = 166

nodal_positive_helenA = 241
nodal_negative_helenA = 91

histology_ductal_helenA = 315
histology_lobular_helenA = 3
histology_other_helenA = 14

grade_1_helenA = 1
grade_2_helenA = 141
grade_3_helenA = 190
grade_NA_helenA = 0

ER_or_PR_pos_helenA = 196
ER_and_PR_neg_helenA = 136

her2_IHC_3plus_helenA = 233
her2_IHC_2plus_ISHpos_helenA = 99


tumor_size_T1_T2_helenA = 277
tumor_size_T3_helenA = 55

mib_under_30_helenA = 73
mib_over_30_helenA = 259


marginals_helenA = {'Age_group_med:Under 50':age_under_50_helenA,
 'Age_group_med:Over 50':age_over_50_helenA,
 'Gradus_preop:1.0':grade_1_helenA,
 'Gradus_preop:2.0':grade_2_helenA,
 'Gradus_preop:3.0':grade_3_helenA,
 'Gradus_preop:N/A':grade_NA_helenA,
 'histology:Carcinoma ductale':histology_ductal_helenA,
 'histology:Carcinoma lobulare':histology_lobular_helenA,
 'histology:other':histology_other_helenA,
 'mib_group_30:under 30%':mib_under_30_helenA,
 'mib_group_30:over 30%':mib_over_30_helenA,
 'Node_cytology_preop:N/A':0,
 'Node_cytology_preop:neg':nodal_negative_helenA,
 'Node_cytology_preop:pos':nodal_positive_helenA,
 'Tumor_size_status:T1-T2':tumor_size_T1_T2_helenA,
 'Tumor_size_status:T3':tumor_size_T3_helenA,
 'Tumor_size_status:other':0,
 'hormone_status:Negative':ER_and_PR_neg_helenA,
 'hormone_status:Positive':ER_or_PR_pos_helenA,
 'IHC_ISH_inter:IHC2+ and ISH+':her2_IHC_2plus_ISHpos_helenA,
 'IHC_ISH_inter:IHC3+':her2_IHC_3plus_helenA,
 'IHC_ISH_inter:other':0}

## Helen-006 group B
#NDocetaxel, carboplatin, trastuzumab, and pertuzumab group

N_helenB = 337
N_pCR_helenB = 194

#Marginals
age_under_50_helenB = 168
age_between_50_70 = 169
age_over_70_helenB = 0

nodal_positive_helenB = 247
nodal_negative_helenB = 90

histology_ductal_helenB = 317
histology_lobular_helenB = 4
histology_other_helenB = 16

grade_1_2_helenB = 143
grade_3_helenB = 194
grade_NA_helenB = 0

ER_or_PR_pos_helenB = 193
ER_and_PR_neg_helenB = 144

her2_IHC_3plus_helenB = 233
her2_IHC_2plus_ISHpos_helenB = 104


tumor_size_T1_T2_helenB = 281
tumor_size_T3_helenB = 56

mib_under_30_helenB = 75
mib_over_30_helenB = 262

bilateral_true_helenB = 0
bilateral_false_helenB = N_helenB

ECOG_0_1_helenB = N_helenB
ECOG_2_helenB = 0
ECOG_NA_helenB = 0

stage_II_III_helenB = N_helenB
stage_other_helenB = 0

marginals_helenB = {'Age_group_med:Under 50':age_under_50_helenB,
                    'Age_group_med:50 - 70':age_between_50_70,
                    'Age_group_med:Over 70':age_over_70_helenB,
                    'Gradus_preop:1 or 2':grade_1_2_helenB,
                    'Gradus_preop:3':grade_3_helenB,
                    'Gradus_preop:N/A':grade_NA_helenB,
                    'histology:Carcinoma ductale':histology_ductal_helenB,
                    'histology:Carcinoma lobulare':histology_lobular_helenB,
                    'histology:other':histology_other_helenB,
                    'mib_group_30:under 30%':mib_under_30_helenB,
                    'mib_group_30:over 30%':mib_over_30_helenB,
                    'Node_cytology_preop:N/A':0,
                    'Node_cytology_preop:neg':nodal_negative_helenB,
                    'Node_cytology_preop:pos':nodal_positive_helenB,
                    'Tumor_size_status:T1-T2':tumor_size_T1_T2_helenB,
                    'Tumor_size_status:T3':tumor_size_T3_helenB,
                    'Tumor_size_status:other':0,
                    #'Ecog_preop:0 or 1':ECOG_0_1_helenB,
                    #'Ecog_preop:2':ECOG_2_helenB,
                    #'Ecog_preop:N/A':ECOG_NA_helenB,
                    'hormone_status:Negative':ER_and_PR_neg_helenB,
                    'hormone_status:Positive':ER_or_PR_pos_helenB,
                    'IHC_ISH_inter:IHC2+ and ISH+':her2_IHC_2plus_ISHpos_helenB,
                    'IHC_ISH_inter:IHC3+':her2_IHC_3plus_helenB,
                    'IHC_ISH_inter:other':0,
                    'Bilateral:True':bilateral_true_helenB,
                    'Bilateral:False':bilateral_false_helenB,
                    'Stage_preop:II-III':stage_II_III_helenB,
                    'Stage_preop:other':stage_other_helenB
                   }

discont_NAT_helenB = 23
finalized_NAT_helenB = N_helenB - 23

marginals_helenB_discont = {'Age_group_med:Under 50':age_under_50_helenB,
                    'Age_group_med:50 - 70':age_between_50_70,
                    'Age_group_med:Over 70':age_over_70_helenB,
                    'Gradus_preop:1 or 2':grade_1_2_helenB,
                    'Gradus_preop:3':grade_3_helenB,
                    'Gradus_preop:N/A':grade_NA_helenB,
                    'histology:Carcinoma ductale':histology_ductal_helenB,
                    'histology:Carcinoma lobulare':histology_lobular_helenB,
                    'histology:other':histology_other_helenB,
                    'mib_group_30:under 30%':mib_under_30_helenB,
                    'mib_group_30:over 30%':mib_over_30_helenB,
                    'Node_cytology_preop:N/A':0,
                    'Node_cytology_preop:neg':nodal_negative_helenB,
                    'Node_cytology_preop:pos':nodal_positive_helenB,
                    'Tumor_size_status:T1-T2':tumor_size_T1_T2_helenB,
                    'Tumor_size_status:T3':tumor_size_T3_helenB,
                    'Tumor_size_status:other':0,
                    #'Ecog_preop:0 or 1':ECOG_0_1_helenB,
                    #'Ecog_preop:2':ECOG_2_helenB,
                    #'Ecog_preop:N/A':ECOG_NA_helenB,
                    'hormone_status:Negative':ER_and_PR_neg_helenB,
                    'hormone_status:Positive':ER_or_PR_pos_helenB,
                    'IHC_ISH_inter:IHC2+ and ISH+':her2_IHC_2plus_ISHpos_helenB,
                    'IHC_ISH_inter:IHC3+':her2_IHC_3plus_helenB,
                    'IHC_ISH_inter:other':0,
                    'Bilateral:True':bilateral_true_helenB,
                    'Bilateral:False':bilateral_false_helenB,
                    'Stage_preop:II-III':stage_II_III_helenB,
                    'Stage_preop:other':stage_other_helenB,
                    'NAT_discontinuation:True':discont_NAT_helenB,
                    'NAT_discontinuation:False':finalized_NAT_helenB
                   }

##############################
# # Neosphere group B
# Docetaxel, trastuzumab, and pertuzumab group

N_neosphereB = 107
N_pCR_neosphereB = 49

#Marginals
age_under_50_neosphereB = 54
age_50_77_neosphereB = 53
age_over_77_neosphere = 0

nodal_positive_neosphereB = 76
nodal_negative_neosphereB = 31

ER_or_PR_pos_neosphereB = 50
ER_and_PR_neg_neosphereB = 57


tumor_size_T1_neosphereB = 0
tumor_size_T2_neosphereB = 54
tumor_size_T3_neosphereB = 53
tumor_size_other_neosphereB = 0

bilateral_true_neosphereB = 0
bilateral_false_neosphereB = N_neosphereB

marginals_neosphereB = {'Age_group_med:Under 50':age_under_50_neosphereB,
 'Age_group_med:50 - 77':age_50_77_neosphereB,
 'Age_group_med:Over 77':age_over_77_neosphere,
 'Node_cytology_preop:N/A':0,
 'Node_cytology_preop:neg':nodal_negative_neosphereB,
 'Node_cytology_preop:pos':nodal_positive_neosphereB,
 'Tumor_size_class_preop:T1':tumor_size_T1_neosphereB,
 'Tumor_size_class_preop:T2':tumor_size_T2_neosphereB,
 'Tumor_size_class_preop:T3':tumor_size_T3_neosphereB,
 'Tumor_size_class_preop:N/A':0,
 'hormone_status:Negative':ER_and_PR_neg_neosphereB,
 'hormone_status:Positive':ER_or_PR_pos_neosphereB,
 'Bilateral:True':bilateral_true_neosphereB,
 'Bilateral:False':bilateral_false_neosphereB
}

discont_NAT_neosphereB = 5
finalized_NAT_neosphereB = N_neosphereB - 5

marginals_neosphereB_discont = {'Age_group_med:Under 50':age_under_50_neosphereB,
 'Age_group_med:50 - 77':age_50_77_neosphereB,
 'Age_group_med:Over 77':age_over_77_neosphere,
 'Node_cytology_preop:N/A':0,
 'Node_cytology_preop:neg':nodal_negative_neosphereB,
 'Node_cytology_preop:pos':nodal_positive_neosphereB,
 'Tumor_size_class_preop:T1':tumor_size_T1_neosphereB,
 'Tumor_size_class_preop:T2':tumor_size_T2_neosphereB,
 'Tumor_size_class_preop:T3':tumor_size_T3_neosphereB,
 'Tumor_size_class_preop:N/A':0,
 'hormone_status:Negative':ER_and_PR_neg_neosphereB,
 'hormone_status:Positive':ER_or_PR_pos_neosphereB,
 'Bilateral:True':bilateral_true_neosphereB,
 'Bilateral:False':bilateral_false_neosphereB,
 'NAT_discontinuation:True':discont_NAT_neosphereB,
 'NAT_discontinuation:False':finalized_NAT_neosphereB
}
##############################
# # KRISTINE group B
# Docetaxel, carboplatin, trastuzumab, and pertuzumab group

N_kristineB = 221
N_pCR_kristineB = 123

#Marginals
age_under_49_kristineB = 110
age_over_49_kristineB = 111

biopsy_nodal_positive_kristineB = 0
no_biopsy_nodal_positive_kristineB = N_kristineB

ER_or_PR_pos_kristineB= 138
ER_and_PR_neg_kristineB = 83


tumor_size_T1_kristineB = 0
tumor_size_T2_T3_kristineB = N_kristineB
tumor_size_other_kristineB = 0

bilateral_true_kristineB = 0
bilateral_false_kristineB = N_kristineB

marginals_kristineB = {'Age_group_med:Under 49':age_under_49_kristineB,
 'Age_group_med:Over 49':age_over_49_kristineB,
 'biopsy_nodal_positive:False':no_biopsy_nodal_positive_kristineB,
 'biopsy_nodal_positive:True':biopsy_nodal_positive_kristineB,
 'Tumor_size_status:T1':tumor_size_T1_kristineB,
 'Tumor_size_status:T2-T3':tumor_size_T2_T3_kristineB,
 'Tumor_size_status:N/A':0,
 'hormone_status:Negative':ER_and_PR_neg_kristineB,
 'hormone_status:Positive':ER_or_PR_pos_kristineB,
 'Bilateral:True':bilateral_true_kristineB,
 'Bilateral:False':N_kristineB
}

discont_NAT_kristineB = 18
finalized_NAT_kristineB = N_kristineB - 18

marginals_kristineB_discont = {'Age_group_med:Under 49':age_under_49_kristineB,
 'Age_group_med:Over 49':age_over_49_kristineB,
 'biopsy_nodal_positive:False':no_biopsy_nodal_positive_kristineB,
 'biopsy_nodal_positive:True':biopsy_nodal_positive_kristineB,
 'Tumor_size_status:T1':tumor_size_T1_kristineB,
 'Tumor_size_status:T2-T3':tumor_size_T2_T3_kristineB,
 'Tumor_size_status:N/A':0,
 'hormone_status:Negative':ER_and_PR_neg_kristineB,
 'hormone_status:Positive':ER_or_PR_pos_kristineB,
 'Bilateral:True':bilateral_true_kristineB,
 'Bilateral:False':N_kristineB,
 'NAT_discontinuation:True':discont_NAT_kristineB,
 'NAT_discontinuation:False':finalized_NAT_kristineB
}
