merge_rawfeed <- function(path){ 

# sort INRAE FEED TABLE

FeedTables_I <- read_excel(paste(path, "feedTables.xlsx", sep = "/"), sheet = "INRAE_feedTable")
Feed_Name_I <- read_excel(paste(path, "category_names.xlsx", sep = "/"), sheet = "INRAE")
Feed_Table_I <- merge(FeedTables_I, Feed_Name_I, by = "FEED_ID_INRAE")

Feed_Table_I= subset(Feed_Table_I, 
           select = c ( 
             Feed_name,
             DM_av,	CP_av,	CF_av,	Cfat_av,	Ash_av,	Insol_ash_av,	NDF_av,	ADF_av,	Lignin_av, Starch_av,	Starch_enz_av,	Sugars_av,	GE_MJ_av,
             DE_growPig_av,	ME_growPig_av,	NE_growPig_av,	DE_adultPig_av,	ME_adultPig_av,	NE_adultPig_av,	ED_growPig_av,	ED_adultPig_av,	OMD_growPig_av,	OMD_adultPig_av,	
             ND_growPig_av,	ND_adultPig_av,	ND_ilealPig_av,	Fatdig_pig_av,	Pdig_pig_noPhyt_av,	Pdig_pig_wPhyt_av
             ))
Feed_Table_I <- subset(Feed_Table_I, Feed_Table_I$Feed_name != 'NA')

Feed_Table_I <- Feed_Table_I %>% 
  group_by(Feed_name) %>%
  summarise(across(, mean, na.rm = TRUE))

# sort NORFOR FEED TABLE

FeedTables_N <- read_excel(paste(path, "feedTables.xlsx", sep = "/"), sheet = "Norfor_feedTable")
Feed_Name_N <- read_excel(paste(path, "category_names.xlsx", sep = "/"), sheet = "Norfor")
Feed_Table_N <- merge(FeedTables_N, Feed_Name_N, by = "FormattedFeedStuffID")

Feed_Table_N= subset(Feed_Table_N, 
                     select = c ( 
                       Feed_name,
                       DM_N, Ash_N,	OM_N, I_CP_N, Cfat_N, NDF_N, iNDF_N, Starch_N,	iStarch_N, Cfiber_N,	Sugar_N, ADF_N,	Lignin_N
                     ))
Feed_Table_N <- subset(Feed_Table_N, Feed_Table_N$Feed_name != 'NA')

Feed_Table_N <- Feed_Table_N %>% 
  group_by(Feed_name) %>%
  summarise(across(, mean, na.rm = TRUE))

# SEGES (DANISH) FEED TABLE

FeedTables_S <- read_excel(paste(path, "feedTables.xlsx", sep = "/"), sheet = "SEGES_feedTable")
Feed_Name_S <- read_excel(paste(path, "category_names.xlsx", sep = "/"), sheet = "SEGES_Gris")
Feed_Table_S <- merge(FeedTables_S, Feed_Name_S, by = "Fodermid_delkode")

Feed_Table_S= subset(Feed_Table_S, 
                     select = c ( 
                       Feed_name,
                       DM_S, Ash_S, OM_S, CP_S,	Cfat_S, Starch_S,	Sugar_S, 	Cfiber_S,	Sol_DF_S,	Insol_DF_S,
                       P_dig_0fyt_S,	P_dig_60fyt_S,	P_dig_100fyt_S,	P_dig_150fyt_S,	P_dig_200fyt_S,	P_dig_250fyt_S,	P_dig_300fyt_S,	P_dig_350fyt_S,	P_dig_400fyt_S
                       ))
Feed_Table_S <- subset(Feed_Table_S, Feed_Table_S$Feed_name != 'NA')

Feed_Table_S <- Feed_Table_S %>% 
  group_by(Feed_name) %>%
  summarise(across(, mean, na.rm = TRUE))

# MERGE feed tables

Feed_Table_SI <- merge(Feed_Table_S, Feed_Table_I, by = "Feed_name")
Feed_Table_SNI <- merge(Feed_Table_SI, Feed_Table_N, by = "Feed_name" )

# export merged table
write.xlsx(Feed_Table_SNI, paste(path, "SNI.xlsx", sep = "/"))

digestabilities(path = path)
}






