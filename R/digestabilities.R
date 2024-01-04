digestabilities <- function(path){

SNI <- read_excel(paste(path, "SNI.xlsx", sep = "/"), sheet = "Sheet 1")

# 	Calculations:  	                                                                                  	 				                                 #
# CHEMICAL COMPOSITION______________________________________________________________________________________________ 

attach(SNI)

# DM:       Dry matter content in feed ingredients 
# Unit:     g/kg feed 
#           All data originates from SEGES feed table
SNI$DM      <- DM_S*10                                                                                            

# Ash:      Ash content in feed ingredients 
# Unit:     g/kg DM 
#           Data from SEGES feed table is prioritized. If unavailable, then the corresponding value from the Norfor feed table is used
SNI$Ash	    <- ifelse(is.na(SNI$Ash_S), Ash_N, Ash_S*10)                                                          

# OM:       Organic matter content in feed ingredients 
# Unit:     g/kg DM 
#           Data from SEGES feed table is prioritized. If unavailable, then the corresponding value content from the Norfor feed table is used
SNI$OM	    <- ifelse(is.na(SNI$OM_S), OM_N, OM_S)

# CP:       Crude protein content in feed ingredients 
# Unit:     g/kg DM 
#           All data originates from SEGES feed table
SNI$CP      <- CP_S*10

# CF:       Crude fat content in feed ingredients 
# Unit:     g/kg DM 
#           Data from SEGES feed table is prioritized. If unavailable, then the corresponding value from the Norfor feed table is used
SNI$CF    <- ifelse(is.na(SNI$Cfat_S), Cfat_N, Cfat_S*10)

# Starch:   Starch content in feed ingredients 
# Unit:     g/kg DM 
#           Data from SEGES feed table is prioritized. If unavailable, then the corresponding value from the Norfor feed table is used
SNI$Starch	<- ifelse(is.na(SNI$Starch_S), Starch_N, Starch_S)

# Sugar:    Sugar content in feed ingredients 
# Unit:     g/kg DM 
#           Data from SEGES feed table is prioritized. If unavailable, then the corresponding value from the Norfor feed table is used
SNI$Sugar   <- ifelse(is.na(SNI$Sugar_S), Sugar_N, Sugar_S) 

attach(SNI)

# StSugar:  Starch + Sugar content in feed ingredients 
# Unit:     g/kg DM 
SNI$StSugar <- Starch + Sugar

# ResFib:   Calculated Fiber content 
# Unit:     g/kg DM 
#           Calculated as the Dry matter (i.e., 1000 g DM) minus CP, CF, Starch, and Sugar. ResFib can not take negative values. Thus, calculated negative values are considered to be zero
SNI$ResFib  <- ifelse(1000 - Ash - CP - CF - Starch - Sugar < 0, 0, 1000 - Ash - CP - CF - Starch - Sugar)

# NDF:      Neutral detergent fiber content 
# Unit:     g/kg DM 
#           All data originates from Norfor feed table
SNI$NDF     <- NDF_N

# DIGESTABILITIES________________________________________________________________________________________

attach(SNI)

# OMdigGrow:          Organic matter digestibility of feed ingredients by GROWING pigs (Apparent total tract digestibilities)
# Unit:               g digested OM/g OM 
#____________         All data originates from INRAE feed table
SNI$OMdigGrow         <-  OMD_growPig_av/100

#____________
# OMdigAdult:         Organic matter digestibility of feed ingredients by ADULT pigs (Apparent total tract digestibilities)
# Unit:               g digested OM/g OM 
#____________         All data originates from INRAE feed table
SNI$OMdigAdult        <-  OMD_adultPig_av/100

#____________
# CFdig:              Crude CF digestibility of feed ingredients by all pigs (GROWING AND ADULTS) (Apparent total tract digestibilities)
# Unit:               g digested CF/g CF 
#____________         All data originates from INRAE feed table
SNI$CFdig             <-  Fatdig_pig_av/100
  
#____________
# CPdigGrow:          Crude protein digestibility of feed ingredients by GROWING pigs (Apparent total tract digestibilities)
# Unit:               g digested CP/g CP 
#____________         All data originates from INRAE feed table
SNI$CPdigGrow         <-  ND_growPig_av/100

#____________
# CPdigAdult:         Crude protein digestibility of feed ingredients by ADULT pigs (Apparent total tract digestibilities)
# Unit:               g digested CP/g CP 
#____________         All data originates from INRAE feed table
SNI$CPdigAdult        <-  ND_adultPig_av/100

# INDIGESTABILITIES____________________________________________________________________________________________________

attach(SNI)

# iCP:                Indigestible Crude protein content of feed ingredients (indigestible CP after 16 to 24 hour incubation in the Cow rumen following passage through the remaining part of the digestive tract) 
# Unit:               g/kg DM 
#____________         All data originates from Norfor feed table
SNI$iCP               <-  (I_CP_N/1000)*CP

# iNDF:               Indigestible neutral detergent fiber content of feed ingredients (indigestible NDF after 288 hour incubation in the Cow rumen) 
# Unit:               g/kg DM 
#____________         All data originates from Norfor feed table
SNI$iNDF              <-  (NDF_N/1000)*NDF

# iStSugar:           Indigestible starch content of feed ingredients (indigestible CP after 16 to 24 hour incubation in the Cow rumen following passage through the remaining part of the digestive tract)
# Unit:               g/kg DM 
#____________         All data originates from Norfor feed table
SNI$iStSugar          <-  (iStarch_N/1000)*StSugar


# DIGESTED/EXCRETED ORGANIC MATTER_______________________________________________________________________________________________________ 

attach(SNI)

# Dig_OMGrow:         Digested organic matter of feed ingredients by GROWING pigs 
# Unit:               g/kg DM 
#____________         Calculated from OM content and OM digestibility
SNI$Dig_OMGrow        <-  OM*OMdigGrow

# Ex_OMGrow:          Excreted organic matter of feed ingredients by GROWING pigs  
# Unit:               g/kg DM 
#____________         Calculated as the difference betwen OM content and digested OM
SNI$Ex_OMGrow         <-  OM-(OM*OMdigGrow)

# Dig_OMAdult:        Digested organic matter of feed ingredients by ADULT pigs
# Unit:               g/kg DM 
#____________         Calculated from OM content and OM digestibility
SNI$Dig_OMAdult       <-  OM*OMdigAdult

# Ex_OMGrow:          Excreted organic matter of feed ingredients by ADULT pigs  
# Unit:               g/kg DM 
#____________         Calculated as the difference between OM content and digested OM
SNI$Ex_OMAdult        <-  OM-(OM*OMdigAdult)

# DIGESTED/EXCRETED CRUDE CF__________________________________________________________________________________________________________________ 

attach(SNI)

# Dig_CF:             Digested crude CF of feed ingredients by ALL pigs
# Unit:               g/kg DM 
#____________         Calculated from CF content and CF digestibility
SNI$Dig_CF          <-  CF*CFdig


# Ex_CF:              Excreted Crude CF of feed ingredients by ALL pigs
# Unit:               g/kg DM 
#____________         Calculated from CF content and CF digestibility
SNI$Ex_CF            <-  CF-(CF*CFdig)


# DIGESTED/EXCRETED CRUDE PROTEIN_______________________________________________________________________________________________________________ 

attach(SNI)

# Dig_CPGrow:         Digested crude protein of feed ingredients by GROWING pigs
# Unit:               g/kg DM 
#____________         Calculated from CP content and CP digestibility
SNI$DigCPGrow         <-  CP*CPdigGrow

#____________
# Ex_CPGrow:          Excreted Crude protein of feed ingredients by GROWING pigs
# Unit:               g/kg DM 
#____________         Calculated from CP content and CP digestibility
SNI$Ex_CPGrow         <-  CP-(CP*CPdigGrow)

#____________
# iCP_Grow_fast:      Excreted Crude protein that is rapidly converted/metabolize during storage of feed ingredients by GROWING pigs
# Unit:               g/kg DM 
#____________         Assumed to equal difference between excreted CP and iCP content. The value i set to equal zero if iCP i greater than Ex_CPGrow
SNI$Ex_CPGrow_fast    <- ifelse(CP-(CP*CPdigGrow) < iCP, 0, (CP-(CP*CPdigGrow))-iCP)

#____________
# iCP_Grow_slow:      Excreted Crude protein that is slowly converted/metabolize during storage of feed ingredients by GROWING pigs
# Unit:               g/kg DM 
#____________         Assumed to equal iCP content originating from the Norfor feed table. The value i set to equal the total amount of excreted CP (Ex_CPGrow) if iCP i greater than Ex_CPGrow
SNI$Ex_CPGrow_slow     <- ifelse(CP-(CP*CPdigGrow) < iCP, CP-(CP*CPdigGrow), iCP)

# Dig_CPAdult:        Digested crude protein of feed ingredients by ADULT pigs
# Unit:               g/kg DM 
#____________         Calculated from CP content and CP digestibility
SNI$DigCPAdult        <-  CP*CPdigAdult

# Ex_CPAdult:         Excreted Crude protein of feed ingredients by ADULT pigs
# Unit:               g/kg DM 
#____________         Calculated from CP content and CP digestibility
SNI$Ex_CPAdult        <-  CP-(CP*CPdigAdult)

# iCP_Adult_fast:     Excreted Crude protein that is rapidly converted/metabolize during storage of feed ingredients by ADULT pigs
# Unit:               g/kg DM 
#____________         Assumed to equal difference between excreted CP and iCP content. The value i set to equal zero if iCP i greater than Ex_CPAdult
SNI$Ex_CPAdult_fast   <- ifelse(CP-(CP*CPdigAdult) < iCP, 0, (CP-(CP*CPdigAdult))-iCP)

# iCP_Adult_slow:     Excreted Crude protein that is slowly converted/metabolize during storage of feed ingredients by ADULT pigs
# Unit:               g/kg DM 
#____________         Assumed to equal iCP content originating from the Norfor feed table. The value i set to equal the total amount of excreted CP (Ex_CPGrow) if iCP i greater than Ex_CPGrow
SNI$Ex_CPAdult_slow    <- ifelse(CP-(CP*CPdigAdult) < iCP, CP-(CP*CPdigAdult), iCP)


# DIGESTED/EXCRETED STARCH AND SUGAR__________________________________________________________________________________________________________ 

attach(SNI)

# DigStSugar:         Digested Starch and sugar of feed ingredients by ALL pigs
# Unit:               g/kg DM 
#____________         Calculated from Starch and Sugar content and the indigestible amount of starch
SNI$DigStSugar        <- StSugar-iStSugar

# Ex_StSugar:         Excreted Starch and sugar of feed ingredients by ALL pigs
# Unit:               g/kg DM 
#____________         Assumed to equal iCP content originating from the Norfor feed table.
SNI$Ex_StSugar         <- iStSugar


# DIGESTED/EXCRETED RESFIB (Calculated fiber content)_____________________________________________________________________________________________________ 

attach(SNI)

# Dig_ResFibGrow:         Digested Calculated fiber (ResFib) of feed ingredients by GROWING pigs
# Unit:                   g/kg DM 
#                         Calculated as the difference beween Dig_OMGrow and digested CF, CP, Starch, and Sugar.  
#                         Dig_ResFibGrow can not take negative values. In that case, negative values are assumed to equal zero
#____________             Dig_ResFibGrow can not be grater than the content of ResFib in the feed. In that case, calculated values assumed to equal ResFib.      
SNI$Dig_ResFibGrow        <-  ifelse( Dig_OMGrow - DigCPGrow - Dig_CF - DigStSugar < 0, 0, 
                              ifelse( Dig_OMGrow - DigCPGrow - Dig_CF - DigStSugar >    ResFib, ResFib, Dig_OMGrow - DigCPGrow - Dig_CF - DigStSugar) )

# Dig_ResFibAdult:        Digested Calculated fiber (ResFib) of feed ingredients by GROWING pigs
# Unit:                   g/kg DM 
#                         Calculated as the difference beween Dig_OMAdult and digested CF, CP, Starch, and Sugar.  
#                         Dig_ResFibAdult can not take negative values. In that case, negative values are assumed to equal zero
#____________             Dig_ResFibAdult can not be grater than the content of ResFib in the feed. In that case, calculated values assumed to equal ResFib.
SNI$Dig_ResFibAdult       <-  ifelse( Dig_OMAdult - DigCPAdult - Dig_CF - DigStSugar < 0, 0, 
                              ifelse( Dig_OMAdult - DigCPAdult - Dig_CF - DigStSugar >  ResFib, ResFib, Dig_OMAdult - DigCPAdult - Dig_CF - DigStSugar))

attach(SNI)

# Dig_ResFib:             Excreted calculated fiber (ResFib) of feed ingredients by GROWING or Adult pigs
# Unit:                   g/kg DM 
#____________             Calculated as the difference beween residual fraction i feed and the estimated digested part of that fraction
SNI$Ex_ResFibGrow         <-  ResFib - Dig_ResFibGrow
SNI$Ex_ResFibAdult        <-  ResFib - Dig_ResFibAdult

# Ex_deg_ResFib:          Excreted calculated fiber (ResFib) of feed ingredients that that can be converted during storage 
# Unit:                   g/kg DM 
#                         Calculated as the excreted fraction minis the indigestible neutral detergent fiber content
#____________             Ex_deg_ResFibGrow can not take negative values. In that case, negative values are assumed to equal zero
SNI$Ex_deg_ResFibGrow     <-  ifelse( ResFib - Dig_ResFibGrow  - iNDF < 0 , 0, ResFib - Dig_ResFibGrow - iNDF )
SNI$Ex_deg_ResFibAdult    <-  ifelse( ResFib - Dig_ResFibAdult - iNDF < 0 , 0, ResFib - Dig_ResFibAdult - iNDF )

attach(SNI)

# Ex_nondeg_ResFib:       Excreted calculated fiber (ResFib) of feed ingredients with no or limited potential fro conversion/degradation during storage
# Unit:                   g/kg DM 
#                         Assumed to equal the content of indigestible  neutral detergent fiber (iNDF)
#____________             Ex_nondeg_ResFibGrow can not take values than the difference between the excreted fraction minus the excreted fraction that can be converted/degraded during storage.                         
SNI$Ex_nondeg_ResFibGrow  <-  ifelse( iNDF > Ex_ResFibGrow  - Ex_deg_ResFibGrow  , Ex_ResFibGrow  - Ex_deg_ResFibGrow,  iNDF)
SNI$Ex_nondeg_ResFibAdult <-  ifelse( iNDF > Ex_ResFibAdult - Ex_deg_ResFibAdult , Ex_ResFibAdult - Ex_deg_ResFibAdult, iNDF)

# fibdig:                 Calculated digestibility of the calculataed fiber fraction by GROWING and ADULT pigs
# Unit:                   g digestedfiber/g fiber 
#                         The proportion of digested amount of fiber relative to the total Calculated fiber content.
#____________             The 0.000001 was added to the denominator to allow calculation of ingredients with a fiber content of zero.
SNI$fibdigGrow             <-  Dig_ResFibGrow/(ResFib + 0.000001)
SNI$fibdigAdult            <-  Dig_ResFibAdult/(ResFib + 0.000001)

attach(SNI)

##################################################################################################################################################
# 	Variable names  	                                                                                 	 				                                 #
##################################################################################################################################################
#_________________________________________________________________________________________________________________________________________________
# COMPOSITION 
#_________________________________________________________________________________________________________________________________________________
# DM                      Dry matter                                                                                        g/kg feed
# Ash                     Ash                                                                                               g/kg DM
# OM                      Organic matter                                                                                    g/kg DM
# CP                      Crude protein                                                                                     g/kg DM    
# CF                      Crude CF                                                                                         g/kg DM
# Starch                  Starch                                                                                            g/kg DM
# Sugar                   Sugar                                                                                             g/kg DM  
# StSugar                 Starh + Sugar                                                                                     g/kg DM  
# ResFib                  Calculated total fiber / fiber calculated as a residual component                                 g/kg DM
# NDF                     Neutral detergent fiber                                                                           g/kg DM
#_________________________________________________________________________________________________________________________________________________
# DIGESTIBILITY 
#_________________________________________________________________________________________________________________________________________________
# OMdigGrow               Organic matter digestibility by growing pigs                                                      g digested OM/g OM
# OMdigAdult              Organic matter digestibility by adult pigs                                                        g digested OM/g OM
# CFdig                 Crude CF digestibility by all pigs                                                               g digested CF/g CF
# CPdigGrow               Crude protein digestibility by growing pigs                                                       g digested CP/g CP
# CPdigAdult              Crude protein digestibility by adult pigs                                                         g digested CP/g CP
#_________________________________________________________________________________________________________________________________________________
# INDIGESTIBILITY 
#_________________________________________________________________________________________________________________________________________________
# iCP                     Indigestible Crude protein (parameter from Norfor feed table)                                     g/kg DM
# iNDF                    Indigestible neutral detergent fiber (parameter from Norfor feed table)                           g/kg DM
# iStarch_N               Indigestible starch (parameter from Norfor feed table)                                            g/kg DM  
#_________________________________________________________________________________________________________________________________________________
# DIGESTED/EXCRETED 
#_________________________________________________________________________________________________________________________________________________
# DigStSugar              Digested starch and sugar by all pigs                                                             g/kg DM
# ExStSugar               Excreted starch and sugar by all pigs                                                             g/kg DM
# Dig_CF                  Digested crude CF by all pigs                                                                    g/kg DM  
# Ex_CF                   Excreted crude fat by all pigs                                                                    g/kg DM
# Dig_OMGrow              Digested organic matter by growing pigs                                                           g/kg DM
# Ex_OMGrow               Excreted organic matter by growing pigs                                                           g/kg DM   
# DigCPGrow               Digested Crude protein by growing pigs                                                            g/kg DM
# Ex_CPGrow               Excreted crude protein by growing pigs                                                            g/kg DM  
# Ex_CPGrow_fast          Excreted crude protein that is rapidly converted during storage                                   g/kg DM  
# Ex_CPGrow_slow          Excreted crude protein that is slowly converted during storage                                    g/kg DM
# Dig_ResFibGrow          Digested Calculated total fiber by growing pigs                                                   g/kg DM
# Ex_ResFibGrow           Excreted calculated total fiber by growing pigs                                                   g/kg DM
# Ex_deg_ResFibGrow       Excreted calculated total fiber that can be converted during storage                              g/kg DM
# Ex_nondeg_ResFibGrow    Excreted calculated total fiber with no or limited potential fro conversion during storage        g/kg DM

# Dig_OMAdult             Digested organic matter by adult pigs                                                             g/kg DM
# Ex_OMAdult              Excreted organic matter by adult pigs                                                             g/kg DM
# DigCPAdult              Digested Crude protein by adult pigs                                                              g/kg DM
# Ex_CPAdult              Excreted crude protein by adult pigs                                                              g/kg DM
# Ex_CPAdult_fast         Excreted crude protein that is rapidly converted during storage                                   g/kg DM  
# Ex_CPAdult_slow         Excreted crude protein that is slowly converted during storage                                    g/kg DM
# Dig_ResFibAdult         Digested Calculated total fiber by adult pigs                                                     g/kg DM
# Ex_ResFibAdult          Excreted calculated total fiber by adult pigs                                                     g/kg DM
# Ex_deg_ResFibAdult      Excreted calculated total fiber that can be converted during storage                              g/kg DM
# Ex_nondeg_ResFibAdult   Excreted calculated total fiber with no or limited potential fro conversion during storage        g/kg DM
# fibdigGrow              Calculated digestibility of the calculated fiber fraction (residual component) by adult pigs      g digested fiber/g fiber
# fibdigAdult             Calculated digestibility of the calculated fiber fraction (residual component) by adult pigs      g digested fiber/g fiber  

#################################################################################################################################################
# 	Table creation and export of data                                                                                  	 				                #
#################################################################################################################################################
#_________________________________________________________________________________________________________________________________________________
#   Selection of variables for table
#_________________________________________________________________________________________________________________________________________________
Table_SNI = subset(SNI, 
                    select = c ( 
         Feed_name,
         DM,            Ash,           OM,          CP,          CF,              Starch,           Sugar,           StSugar,         ResFib,               NDF,
         OMdigGrow,     OMdigAdult,    CFdig,       CPdigGrow,   CPdigAdult,
         iCP,           iNDF,          iStarch_N,
         DigStSugar,    Ex_StSugar,
         Dig_CF,        Ex_CF,
         Dig_OMGrow,    Ex_OMGrow,     DigCPGrow,   Ex_CPGrow,   Ex_CPGrow_fast,    Ex_CPGrow_slow,   Dig_ResFibGrow,  Ex_ResFibGrow,   Ex_deg_ResFibGrow,    Ex_nondeg_ResFibGrow,
         Dig_OMAdult,   Ex_OMAdult,    DigCPAdult,  Ex_CPAdult,  Ex_CPAdult_fast,   Ex_CPAdult_slow,  Dig_ResFibAdult, Ex_ResFibAdult,  Ex_deg_ResFibAdult,   Ex_nondeg_ResFibAdult,                       
         fibdigGrow,    fibdigAdult,   Ex_Ash,      Ex_P,        Ex_K,              Ex_TAN 
         ))

detach(SNI)

# Export data____________________________________________________________________________________________________________________________
write.xlsx(Table_SNI, paste(path, "Ex_feces.xlsx", sep = "/"))
}













