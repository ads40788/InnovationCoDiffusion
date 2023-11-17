require(sna)
require(igraph)
require(dils)
require(intergraph)
library(sna) 
library(asnipe)
library(igraph)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(stargazer)

####FUNCTIONS AND DATA

###########################
###PART 0.0 - LOAD DATA####
###########################

#setwd("C:/ProgramData/MySQL/MySQL Server 8.0/Data")
setwd("C:/Users/adsir/Dropbox/GitHub")
NetworkFrames_All<- read.csv('IntersectionFramesAll.csv',header=TRUE)
NetworkFrames_All$Lang2Perc <- NetworkFrames_All$L2Y1/NetworkFrames_All$BaseSample

####PART 0.1 = Define Functions

#function that turns dataframe into two weighted adjacency matricies
adjmatrify <- function(netframe,col_vals){
  netframe = droplevels(netframe)
  main_nodes = union(levels(netframe$Lang1),levels(netframe$Lang2))
  nodes = levels(netframe$Lang1)
  node_c = length(nodes)
  adj_matrix = matrix(0:0, nrow = node_c, ncol = node_c)
  for(x in 1:node_c){
    for(y in 1:node_c){
      Lang1 = nodes[x] 
      Lang2 = nodes[y]
      TargetValue = netframe[which((netframe['Lang1']==Lang1)& (netframe['Lang2']==Lang2)),][1,][[col_vals]]
      adj_matrix[x,y] = TargetValue
    }
  }
  rownames(adj_matrix) <- nodes
  return(adj_matrix)
}

adjmatrify_asym <- function(netframe,col_vals){
  netframe =droplevels(netframe)
  main_nodes_1 = levels(netframe$Lang1)
  main_nodes_2 = levels(netframe$Lang2)
  #nodes_1 = levels(netframe$Lang1)
  node_1_c = length(main_nodes_1)
  node_2_c = length(main_nodes_2)
  adj_matrix = matrix(0:0, nrow = node_1_c, ncol = node_2_c)
  for(x in 1:node_1_c){
    for(y in 1:node_2_c){
      Lang1 = main_nodes_1[x] 
      Lang2 = main_nodes_2[y]
      TargetValue = netframe[which((netframe['Lang1']==Lang1)& (netframe['Lang2']==Lang2)),][1,][[col_vals]]
      adj_matrix[x,y]=TargetValue
    }
  }
  rownames(adj_matrix) <- main_nodes_1
  colnames(adj_matrix) <- main_nodes_2
  return(adj_matrix)
}

#This calcualtes a measure of cosine similarity directly from Obs/Exp data.
cosineify <- function(netframe_cos,target_value,newcolumnname,asym=FALSE){
  #Get adjacency matrix from target value
  if(asym==FALSE){
    netframe_adj <- adjmatrify(netframe_cos,target_value)
  }
  else{
    netframe_adj <- adjmatrify_asym(netframe_cos,target_value)
  }
  #Replace diagonals with 0s
  netframe_adj[is.na(netframe_adj)==TRUE]<-0 
  #print(netframe_adj)
  #Count rows and create new matricies to store results of calculations
  rows = nrow(netframe_adj)
  #new_matrix = matrix(0:0,nrow = rows,ncol=rows)
  long_matrix =matrix(0:0,0,3)
  #iterate through cells and create matrix
  for(x in (1:rows)){
    for(y in (1:rows)){
      if(x!=y){
        #print(c(x,y))
        row1name = rownames(netframe_adj)[x]
        row2name = rownames(netframe_adj)[y]
        #print(c(row1name,row2name))
        row1vec <-netframe_adj[x,]
        row2vec <-netframe_adj[y,]
        row1norm <- 1/sqrt(sum(row1vec^2))
        row2norm <- 1/sqrt(sum(row2vec^2))
        equiv = (row1vec%*%row2vec)*row1norm*row2norm
        #print(equiv)
        #new_matrix[x,y] = equiv
        #rownames(new_matrix) = rownames(netframe_adj)
        #colnames(new_matrix) = rownames(netframe_adj)
        equivrow <-c(row1name,row2name,equiv)
        long_matrix =rbind(long_matrix,equivrow)
      }
    }
  }
  #Rename columns in long matrix and turn into dataframe
  colnames(long_matrix)<- c('Lang1','Lang2',newcolumnname)
  long_matrix <-as.data.frame(long_matrix)
  #Merge newly calculated cosine similarities into original matrix
  #netframe_new <-merge(netframe,long_matrix,by=c('Lang1','Lang2'))
  #netframe_new$newcolumnname <- as.numeric(as.character(netframe_new$newcolumnname))
  return(long_matrix)
}

###draws on cosinify returns final function
create_frame <-function(TheFrame,target_year,iv_cutoff,dv_cutoff){
  netframe_year = TheFrame[which(TheFrame$Year==target_year),]
  
  netframe_cutoff_low = netframe_year[((netframe_year$Lang1Perc>iv_cutoff) & (netframe_year$Lang2Perc>iv_cutoff)),]
  netframe_cutoff_high = netframe_year[((netframe_year$Lang1Perc>dv_cutoff) & (netframe_year$Lang2Perc>dv_cutoff)),]
  netframe_high_to_low = netframe_year[((netframe_year$Lang1Perc>dv_cutoff) & (netframe_year$Lang2Perc>iv_cutoff)& (netframe_year$Lang2Perc<dv_cutoff)),]
  
  matrix_cutoff_low <- cosineify(netframe_cutoff_low,'ObsExp','Cosine_LowCutoff')
  matrix_cutoff_high <- cosineify(netframe_cutoff_high,'ObsExp','Cosine_HighCutoff')
  matrix_high_to_low <- cosineify(netframe_high_to_low,'ObsExp','Cosine_HighToLow',asym=TRUE)
  
  netframe_year_final <- merge(netframe_year,matrix_cutoff_low,by=c('Lang1','Lang2'))
  netframe_year_final <- merge(netframe_year_final,matrix_cutoff_high,by=c('Lang1','Lang2'))
  netframe_year_final <- merge(netframe_year_final,matrix_high_to_low,by=c('Lang1','Lang2'))
  
  netframe_year_final$Cosine_LowCutoff <- as.numeric(as.character(netframe_year_final$Cosine_LowCutoff))
  netframe_year_final$Cosine_HighCutoff <- as.numeric(as.character(netframe_year_final$Cosine_HighCutoff))
  netframe_year_final$Cosine_HighToLow <- as.numeric(as.character(netframe_year_final$Cosine_HighToLow))
  return(netframe_year_final)
}

create_frame_simp<-function(TheFrame,target_year,cutoff){
  netframe_year = TheFrame[which(TheFrame$Year==target_year),]
  netframe_cutoff = netframe_year[((netframe_year$Lang1Perc>cutoff) & (netframe_year$Lang2Perc>cutoff)),]
  matrix_cutoff = cosineify(netframe_cutoff,'ObsExp','Cosine_')
  netframe_year_final = merge(netframe_cutoff,matrix_cutoff,by=c('Lang1','Lang2'))
}

filterFrame <- function(netframe,listoflanguages){
  netframe_filtered <- filter(netframe, netframe$Lang1 %in% listoflanguages)
  netframe_filtered <- filter(netframe_filtered, netframe_filtered$Lang2 %in% listoflanguages)
  return(netframe_filtered)
}

selectVariables <- function(data, element_column, weight_column, chosen_languages) {
  # Extract elements and weights from the data frame
  elements <- data[[element_column]]
  weights <- data[[weight_column]]
  # Randomly select indices based on weights
  if (length(elements)>chosen_languages){
    indices <- sample(length(elements), size = chosen_languages, replace = FALSE, prob = weights)
    selected_variables <- elements[indices]
  }
  else{
    selected_variables = elements
  }
  # Return selected variables
  return(selected_variables)
}

#####DATAFRAME FORM OUTPUT
QAP_Results_July_2023_Full <- data.frame(year = integer(),    # Create empty data frame
                                         cutoff_percentage = character(),
                                         ###output Piggybacking Model
                                         pbr_constant_coef = numeric(),pbr_cosine_coef = numeric(),pbr_obsexp_coef = numeric(),
                                         pbr_constant_tval = numeric(),pbr_cosine_tval = numeric(),pbr_obsexp_tval = numeric(),
                                         pbr_constant_pval1 = numeric(),pbr_cosine_pval1 = numeric(),pbr_obsexp_pval1 = numeric(),
                                         pbr_constant_pval2 = numeric(),pbr_cosine_pval2 = numeric(),pbr_obsexp_pval2 = numeric(),
                                         ###Output Cannibalism Model
                                         cbr_constant_coef = numeric(),cbr_cosine_coef = numeric(),cbr_obsexp_coef = numeric(),
                                         cbr_constant_tval = numeric(),cbr_cosine_tval = numeric(),cbr_obsexp_tval = numeric(),
                                         cbr_constant_pval1 = numeric(),cbr_cosine_pval1 = numeric(),cbr_obsexp_pval1 = numeric(),
                                         cbr_constant_pval2 = numeric(),cbr_cosine_pval2 = numeric(),cbr_obsexp_pval2 = numeric(),
                                         list_of_langauges = list(),
                                         stringsAsFactors = FALSE)


############################################
###Run Models Without Randomization ########
############################################


#OutputFrame = QAP_Results_July_2023_Full
### 1 - Create Frame Basic
#for(A_CUTOFF_PERCENTAGE in c(0.01)){
for(A_CUTOFF_PERCENTAGE in c(0.005,0.01,0.02)){
  #for(LANGS_PER_NETWORK in c(20)){
  #for(A_YEAR in 2010:2017){
  for(A_YEAR in 2010:2017){
    Netframe_Year<-create_frame_simp(NetworkFrames_All,A_YEAR,A_CUTOFF_PERCENTAGE)
    LangaugesWeights = unique(Netframe_Year[,c('Lang1','L1Y1')])
    variable_subset <-selectVariables(LangaugesWeights,'Lang1','L1Y1',10000)
    Netframe_Year_Filtered <- filterFrame(Netframe_Year,variable_subset)
    #variable_subset= unique(Netframe_Year[,c('Lang1','L1Y1')])
    #Netframe_Year_Filtered <- filterFrame(Netframe_Year,variable_subset)
    ###1 Create Network Variables
    NetworkFrames_CBR <-adjmatrify(Netframe_Year_Filtered, 'LogCBR')
    NetworkFrames_PBR <-adjmatrify(Netframe_Year_Filtered, 'LogPBR')
    NetworkFrames_ObsExp <-adjmatrify(Netframe_Year_Filtered, 'ObsExp')
    Netframe_Year_Filtered$'Cosine_'<-as.numeric(as.character(Netframe_Year_Filtered$'Cosine_'))
    NetworkFrames_Cos <-adjmatrify(Netframe_Year_Filtered, 'Cosine_')
    #### 2 Fit Models
    NetworkFrames_Predictors<-list(NetworkFrames_Cos,NetworkFrames_ObsExp)
    CBRModel <- netlm(NetworkFrames_CBR,NetworkFrames_Predictors, mode="graph", nullhyp="qap", test.statistic="t-value")
    PBRModel <- netlm(NetworkFrames_PBR,NetworkFrames_Predictors, mode="graph", nullhyp="qap", test.statistic="t-value")
    ####Print Output
    OutputRow<-c(A_YEAR,A_CUTOFF_PERCENTAGE,
                 PBRModel[1]$coefficients[1],PBRModel[1]$coefficients[2],PBRModel[1]$coefficients[3],
                 PBRModel[8]$tstat[1],PBRModel[8]$tstat[2],PBRModel[8]$tstat[3],
                 PBRModel[10]$pleeq[1],PBRModel[10]$pleeq[2],PBRModel[10]$pleeq[3],
                 PBRModel[12]$pgreqabs[1],PBRModel[12]$pgreqabs[2],PBRModel[12]$pgreqabs[3],
                 CBRModel[1]$coefficients[1],CBRModel[1]$coefficients[2],CBRModel[1]$coefficients[3],
                 CBRModel[8]$tstat[1],CBRModel[8]$tstat[2],CBRModel[8]$tstat[3],
                 CBRModel[10]$pleeq[1],CBRModel[10]$pleeq[2],CBRModel[10]$pleeq[3],
                 CBRModel[12]$pgreqabs[1],CBRModel[12]$pgreqabs[2],CBRModel[12]$pgreqabs[3]
    )
    OutputLength = nrow(QAP_Results_July_2023_Full)
    QAP_Results_July_2023_Full[OutputLength+1,]<-OutputRow
    print(c(A_CUTOFF_PERCENTAGE,A_YEAR))
  }
}

write.csv(QAP_Results_July_2023_Full,'prelim_data_july_2023.csv')


##########################################
##### Plot Models ########################
##########################################

FullResults <- QAP_Results_July_2023_Full[which(QAP_Results_July_2023_Full$cutoff_percentage >0.001),]

CBR_Coef_Cos <- FullResults[, c("year", "cbr_cosine_coef","cutoff_percentage")]
CBR_Coef_Cos['variable_type']='cosine (equivalence)'
colnames(CBR_Coef_Cos)[2]<-'Coefficient'

CBR_Coef_Obs <- FullResults[, c("year", "cbr_obsexp_coef","cutoff_percentage")]
CBR_Coef_Obs['variable_type']='observed/expected (cohesion)'
colnames(CBR_Coef_Obs)[2]<-'Coefficient'

CBR_Pval_Cos <- FullResults[, c("year", "cbr_cosine_pval1","cutoff_percentage")]
CBR_Pval_Cos['variable_type']='cosine (equivalence)'
colnames(CBR_Pval_Cos)[2]<-'PValue'

CBR_Pval_Obs <- FullResults[, c("year", "cbr_obsexp_pval1","cutoff_percentage")]
CBR_Pval_Obs['variable_type']='observed/expected (cohesion)'
colnames(CBR_Pval_Obs)[2]<-'PValue'

CBR_Coef <- rbind(CBR_Coef_Cos,CBR_Coef_Obs)
CBR_Pval <-rbind(CBR_Pval_Cos,CBR_Pval_Obs)

PBR_Coef_Cos <- FullResults[, c("year", "pbr_cosine_coef","cutoff_percentage")]
PBR_Coef_Cos['variable_type']='cosine (equivalence)'
colnames(PBR_Coef_Cos)[2]<-'Coefficient'

PBR_Coef_Obs <- FullResults[, c("year", "pbr_obsexp_coef","cutoff_percentage")]
PBR_Coef_Obs['variable_type']='observed/expected (cohesion)'
colnames(PBR_Coef_Obs)[2]<-'Coefficient'

PBR_Pval_Cos <- FullResults[, c("year", "pbr_cosine_pval1","cutoff_percentage")]
PBR_Pval_Cos['variable_type']='cosine (equivalence)'
colnames(PBR_Pval_Cos)[2]<-'PValue'

PBR_Pval_Obs <- FullResults[, c("year", "pbr_obsexp_pval1","cutoff_percentage")]
PBR_Pval_Obs['variable_type']='observed/expected (cohesion)'
colnames(PBR_Pval_Obs)[2]<-'PValue'

PBR_Coef <- rbind(PBR_Coef_Cos,PBR_Coef_Obs)
PBR_Pval <-rbind(PBR_Pval_Cos,PBR_Pval_Obs)

CBR_Coef$cutoff_adj = 0
CBR_Coef$cutoff_adj[CBR_Coef$cutoff_percentage==0.001]=-0.4
CBR_Coef$cutoff_adj[CBR_Coef$cutoff_percentage==0.005]=-0.3
CBR_Coef$cutoff_adj[CBR_Coef$cutoff_percentage==0.02]=0.3
CBR_Coef$year_adj = CBR_Coef$year+CBR_Coef$cutoff_adj

CBR_Pval$cutoff_adj = 0
CBR_Pval$cutoff_adj[CBR_Pval$cutoff_percentage==0.001]=-0.4
CBR_Pval$cutoff_adj[CBR_Pval$cutoff_percentage==0.005]=-0.3
CBR_Pval$cutoff_adj[CBR_Pval$cutoff_percentage==0.02]=0.3
CBR_Pval$year_adj = CBR_Pval$year+CBR_Pval$cutoff_adj

PBR_Coef$cutoff_adj = 0
PBR_Coef$cutoff_adj[PBR_Coef$cutoff_percentage==0.001]=-0.4
PBR_Coef$cutoff_adj[PBR_Coef$cutoff_percentage==0.005]=-0.3
PBR_Coef$cutoff_adj[PBR_Coef$cutoff_percentage==0.02]=0.3
PBR_Coef$year_adj = PBR_Coef$year+PBR_Coef$cutoff_adj

PBR_Pval$cutoff_adj = 0
PBR_Pval$cutoff_adj[PBR_Pval$cutoff_percentage==0.001]=-0.4
PBR_Pval$cutoff_adj[PBR_Pval$cutoff_percentage==0.005]=-0.3
PBR_Pval$cutoff_adj[PBR_Pval$cutoff_percentage==0.02]=0.3
PBR_Pval$year_adj = PBR_Pval$year+PBR_Pval$cutoff_adj


CBR_Coef_Plot <- ggplot(CBR_Coef, aes(x = year_adj, y = Coefficient, color = variable_type,shape=cutoff_percentage))+ 
  geom_point()+ geom_hline(yintercept=0, linetype="solid", color = "black") + 
  geom_vline(xintercept=c(2010.5,2011.5,2012.5,2013.5,2014.5,2015.5,2016.5), linetype="solid", color = "grey") + 
  scale_x_continuous(breaks=seq(2010,2017,1))+
  scale_fill_discrete() + labs(y = "Coefficient Size",x = "",color='Variable',shape='Cutoff Proportion') + 
  theme_bw() + theme(legend.position = "none")+ ggtitle("Competitive Diffusion")+ theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
#CBR_Coef_Plot

CBR_Pval_Plot <- ggplot(CBR_Pval, aes(x = year_adj, y = PValue, color = variable_type,shape=cutoff_percentage))+ 
  geom_point()+ geom_hline(yintercept=0.05, linetype="dashed", color = "blue")+
  geom_hline(yintercept=0.95, linetype="dotted", color = "red")+
  geom_vline(xintercept=c(2010.5,2011.5,2012.5,2013.5,2014.5,2015.5,2016.5), linetype="solid", color = "grey") +
  scale_x_continuous(breaks=seq(2010,2017,1))+
  scale_color_discrete(labels=c('Cosine (Equivalence)','Observed/Expected (Cohesion)')) +theme_bw()+ labs(y = "P-Value",x = "Year",color='Variable',shape='Cutoff Proportion')+
  theme(legend.position = "bottom")+ theme(plot.title = element_text(hjust = 0.5))+ guides(shape = "none")+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())

#CBR_Pval_Plot

PBR_Coef_Plot <- ggplot(PBR_Coef, aes(x = year_adj, y = Coefficient, color = variable_type,shape=cutoff_percentage))+ 
  geom_point()+ geom_hline(yintercept=0, linetype="solid", color = "black")+
  geom_vline(xintercept=c(2010.5,2011.5,2012.5,2013.5,2014.5,2015.5,2016.5), linetype="solid", color = "grey") +
  scale_x_continuous(breaks=seq(2010,2017,1))+
  scale_fill_discrete() +theme_bw()+ 
  labs(y = "Coefficient Size",x = "",color='Variable',shape='Cutoff Proportion')+theme(legend.position = "none")+
  ggtitle("Complementary Diffusion")+ theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())

#PBR_Coef_Plot

PBR_Pval_Plot <- ggplot(PBR_Pval, aes(x = year_adj, y = PValue, color = variable_type,shape=cutoff_percentage))+ 
  geom_point()+ geom_hline(yintercept=0.95, linetype="dashed", color = "blue")+
  geom_hline(yintercept=0.025, linetype="dotted", color = "red")+
  geom_hline(yintercept=0.975, linetype="dotted", color = "red")+
  geom_vline(xintercept=c(2010.5,2011.5,2012.5,2013.5,2014.5,2015.5,2016.5), linetype="solid", color = "grey") +
  scale_x_continuous(breaks=seq(2010,2017,1))+
  scale_fill_discrete() +theme_bw()+ labs(y = "P-Value",x = "Year",color='Variable',shape='Cutoff Proportion')+
  theme(legend.position = "bottom")+ theme(plot.title = element_text(hjust = 0.5))+ guides(color = "none")+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())

#PBR_Pval_Plot

grid.arrange(PBR_Coef_Plot,CBR_Coef_Plot,PBR_Pval_Plot,CBR_Pval_Plot,ncol=2)



############################
######Create Tables ########
############################

FullResults_Tabular <-FullResults
FullResults_Tabular$year_plus_1 = FullResults_Tabular$year+1

####Complementary Diffusion Models
#FullResults_Tabular_PBR[c('cutoff_percentage','year','year_plus_1','pbr_constant_coef','pbr_obsexp_coef','pbr_obsexp_tval','pbr_obsexp_pval1','pbr_cosine_coef','pbr_cosine_tval','pbr_cosine_pval1')]
#FullResults_Tabular_CBR[c('cutoff_percentage','year','year_plus_1','pbr_constant_coef','pbr_obsexp_coef','pbr_obsexp_tval','pbr_obsexp_pval1','pbr_cosine_coef','pbr_cosine_tval','pbr_cosine_pval1')]


FullResults_Tabular_PBR = as.data.frame(cbind(
  as.character(FullResults_Tabular$cutoff_percentage),
  as.character(FullResults_Tabular$year),
  as.character(FullResults_Tabular$year_plus_1),
  round(as.numeric(as.character(FullResults_Tabular$pbr_constant_coef)),digits=3),
  round(as.numeric(as.character(FullResults_Tabular$pbr_obsexp_coef)),digits=3),
  round(as.numeric(as.character(FullResults_Tabular$pbr_obsexp_tval)),digits=3),
  round(as.numeric(as.character(FullResults_Tabular$pbr_obsexp_pval1)),digits=3),
  round(as.numeric(as.character(FullResults_Tabular$pbr_cosine_coef)),digits=3),
  round(as.numeric(as.character(FullResults_Tabular$pbr_cosine_tval)),digits=3),
  round(as.numeric(as.character(FullResults_Tabular$pbr_cosine_pval1)),digits=3)))
names(FullResults_Tabular_PBR)<-c('Cutoff','Predictor Year','Outcome Year','Constant Coef.',
                                  'Cohes. Coef.','T-Stat', 'Sim. < Obs.',
                                  'Equiv. Coef.', 'T-Stat.','Sim. < Obs.')

xtable(FullResults_Tabular_PBR,file='')


FullResults_Tabular_CBR = as.data.frame(cbind(
  as.character(FullResults_Tabular$cutoff_percentage),
  as.character(FullResults_Tabular$year),
  as.character(FullResults_Tabular$year_plus_1),
  round(as.numeric(as.character(FullResults_Tabular$cbr_constant_coef)),digits=3),
  round(as.numeric(as.character(FullResults_Tabular$cbr_obsexp_coef)),digits=3),
  round(as.numeric(as.character(FullResults_Tabular$cbr_obsexp_tval)),digits=3),
  round(as.numeric(as.character(FullResults_Tabular$cbr_obsexp_pval1)),digits=3),
  round(as.numeric(as.character(FullResults_Tabular$cbr_cosine_coef)),digits=3),
  round(as.numeric(as.character(FullResults_Tabular$cbr_cosine_tval)),digits=3),
  round(as.numeric(as.character(FullResults_Tabular$cbr_cosine_pval1)),digits=3)))
names(FullResults_Tabular_CBR)<-c('Cutoff','Predictor Year','Outcome Year','Constant Coef.',
                                  'Cohes. Coef.','T-Stat', 'Sim. < Obs.',
                                  'Equiv. Coef.', 'T-Stat.','Sim. < Obs.')

xtable(FullResults_Tabular_CBR,file='')



####Table of Languages By Cutoff:
CutoffTable <- data.frame(Year =integer(),
                          Cutoff = numeric(),
                          Count= integer(),
                          stringsAsFactors = FALSE
)
CutoffPercentages = c(0.005,0.01,0.02)
for(A_YEAR in 2010:2017){
  for(A_CUTOFF_PERCENTAGE in CutoffPercentages){
    Netframe_Year<-create_frame_simp(NetworkFrames_All,A_YEAR,A_CUTOFF_PERCENTAGE)
    LangaugeCount = nrow(unique(Netframe_Year[,c('Lang1','L1Y1')]))
    CutoffRow=c(A_YEAR,A_CUTOFF_PERCENTAGE,LangaugeCount)
    OutputLength = nrow(CutoffTable)
    CutoffTable[OutputLength+1,]<-CutoffRow
  }
}
reshape(CutoffTable, idvar = "Year", timevar = "Cutoff", direction = "wide")



