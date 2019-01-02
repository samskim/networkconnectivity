#list of summary statistics (traits to be analyzed)
declare -a LIST=("UKB_460K.body_BMIz" "UKB_460K.cov_EDU_YEARS" "UKB_460K.lung_FVCzSMOKE" "UKB_460K.cov_SMOKING_STATUS" 
     "UKB_460K.mental_NEUROTICISM" "UKB_460K.blood_WHITE_COUNT" "PASS_Years_of_Education2" "UKB_460K.bp_SYSTOLICadjMEDz" 
      "UKB_460K.body_HEIGHTz" "UKB_460K.other_MORNINGPERSON" "UKB_460K.body_WHRadjBMIz" "UKB_460K.lung_FEV1FVCzSMOKE" 
       "UKB_460K.repro_MENARCHE_AGE" "UKB_460K.blood_RED_COUNT" "UKB_460K.blood_PLATELET_COUNT" "UKB_460K.bmd_HEEL_TSCOREz" 
        "UKB_460K.blood_EOSINOPHIL_COUNT" "PASS_Schizophrenia" "UKB_460K.blood_RBC_DISTRIB_WIDTH" "PASS_Height1" "PASS_BMI1" 
	 "UKB_460K.disease_T2D" "PASS_AgeFirstBirth" "UKB_460K.disease_RESPIRATORY_ENT" "UKB_460K.body_BALDING1" "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP" 
	  "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED" "UKB_460K.disease_HI_CHOL_SELF_REP" "UKB_460K.repro_MENOPAUSE_AGE" "PASS_HDL" "UKB_460K.pigment_SUNBURN" 
	   "PASS_NumberChildrenEverBorn" "PASS_Anorexia" "PASS_LDL" "PASS_Crohns_Disease" "PASS_DS" "PASS_Ever_Smoked" "UKB_460K.pigment_HAIR" 
	    "PASS_Rheumatoid_Arthritis" "PASS_Type_2_Diabetes" "PASS_Autism" "UKB_460K.pigment_TANNING" "PASS_Ulcerative_Colitis" 
	     "UKB_460K.disease_DERMATOLOGY" "PASS_Coronary_Artery_Disease" "UKB_460K.disease_AID_SURE" "UKB_460K.pigment_SKIN");
	     
for TRAIT in "${LIST[@]}"
do
	python /random/directory/PartitionHeritability_baselineLD.py "\""$TRAIT"\"".sumstats "$1"
done

