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

### LOAD DATA ###

#setwd("C:/ProgramData/MySQL/MySQL Server 8.0/Data")
setwd("C:/Users/adsir/Dropbox/GitHub")
NetworkFrames_All<- read.csv('IntersectionFramesAll.csv',header=TRUE)
NetworkFrames_All$Lang2Perc <- NetworkFrames_All$L2Y1/NetworkFrames_All$BaseSample


### DEFINE FUNCTIONS ###

#First, a function to scrub down dataframe to relevant year and languages
create_frame_simp<-function(TheFrame,target_year,cutoff){
  netframe_year = TheFrame[which(TheFrame$Year==target_year),]
  netframe_cutoff = netframe_year[((netframe_year$Lang1Perc>cutoff) & (netframe_year$Lang2Perc>cutoff)),]
  matrix_cutoff = cosineify(netframe_cutoff,'ObsExp','Cosine_')
  netframe_year_final = merge(netframe_cutoff,matrix_cutoff,by=c('Lang1','Lang2'))
  netframe_year_final$Cosine_ <- as.numeric(as.character(netframe_year_final$Cosine_))
  return(netframe_year_final)
}

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
      adj_matrix[x,y]=TargetValue
    }
  }
  rownames(adj_matrix) <- nodes
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

#Visualization functions for PBR/CBR
piggyback_plot <- function(Year,Cutoff,high_cutoff,low_cutoff){
  TargetFrame <-eval(as.name(paste("NetworkFrame",toString(Year),toString(Cutoff),sep='_')))
  ####Default Color Parameters
  color_high = 'purple'
  color_medium = 'grey'
  color_low = 'green'
  shape_high = 3
  shape_medium = 2
  shape_low = 1
  ####Assign Percentile-Based Categories To PBR Values
  high_cutoff_pig = quantile(TargetFrame$LogPBR,high_cutoff)
  low_cutoff_pig = quantile(TargetFrame$LogPBR,low_cutoff)
  TargetFrame$PiggybackingCategory <- cut(TargetFrame$LogPBR,
                                          breaks=c(min(TargetFrame$LogPBR)-.01,low_cutoff_pig,high_cutoff_pig,max(TargetFrame$LogPBR)+.01),
                                          labels=c('Low PBR', 'Medium PBR', 'High PBR'))
  #####Calculate Means For Density Plot Lines
  mean_piggyback_df <- TargetFrame %>%
    group_by(PiggybackingCategory) %>%
    summarise(meanCos=mean(Cosine_),meanOE = mean(ObsExp))
  ####Density Plots
  ####Cosine Similarity vs. PBR
  cos_pig_density <- ggplot(TargetFrame, aes(x=Cosine_, fill=PiggybackingCategory)) + 
    geom_density(alpha=.5) + 
    geom_vline(data = mean_piggyback_df, aes(xintercept = meanCos,color = PiggybackingCategory), size=1.5)+
    scale_fill_manual(values = c(color_high,color_medium,color_low)) + 
    scale_colour_manual(values = c(color_high,color_medium,color_low)) + 
    theme_bw()+theme(legend.position = "empty")+xlab(paste(toString(Year),'Cosine Similarity (Equivalence)'))
  ####Observed/Expected Ratio vs. PBR
  oe_pig_density <- ggplot(TargetFrame, aes(y=ObsExp, fill=PiggybackingCategory)) + 
    geom_density(alpha=.5) + 
    geom_hline(data = mean_piggyback_df, aes(yintercept = meanOE,color = PiggybackingCategory), size=1.5)+
    scale_fill_manual(values = c(color_high,color_medium,color_low)) + 
    scale_colour_manual(values = c(color_high,color_medium,color_low)) + 
    theme_bw()+theme(legend.position = "empty")+ylab(paste(toString(Year),'Observed/Expected (Cohesion)'))
  #Scatter Plot
  pig_scatter <-ggplot(TargetFrame)+ 
    geom_point(aes(x=Cosine_, y=ObsExp,color = PiggybackingCategory,shape=PiggybackingCategory),size =2,alpha=0.75)+
    scale_color_manual(values = c(color_high,color_medium,color_low)) + 
    scale_shape_manual(values = c(shape_high,shape_medium,shape_low)) +
    theme_bw()+theme(legend.position = "empty")+xlab(paste(toString(Year),'Cosine Similarity (Equivalence)'))+ylab(paste(toString(Year),'Observed/Expected (Cohesion)'))
  empty <- ggplot(TargetFrame)+geom_point(aes(colour=PiggybackingCategory,shape=PiggybackingCategory),x=NA,y=NA,size=4 )+
    scale_color_manual(name=paste(toString(Year),'to',toString(Year+1),':\nPiggybacking Diffusion Rate'),values = c(color_high,color_medium,color_low)) +
    scale_shape_manual(name=paste(toString(Year),'to',toString(Year+1),':\nPiggybacking Diffusion Rate'),values = c(shape_high,shape_medium,shape_low)) +
    theme(axis.ticks=element_blank(), 
          panel.background=element_blank(), 
          axis.text.x=element_blank(), axis.text.y=element_blank(),           
          axis.title.x=element_blank(), axis.title.y=element_blank())+
    labs(color= paste(toString(Year),'to',toString(Year+1),':\nPiggybacking Diffusion Rate'))
  grid.arrange(cos_pig_density, empty, pig_scatter, oe_pig_density, ncol=2, nrow=2, widths=c(5, 2), heights=c(2, 5))
}

cannibal_plot<-function(Year,Cutoff,high_cutoff,low_cutoff){
  #high_cutoff and low_cutoff correspond  to categorization of rates.
  TargetFrame <-eval(as.name(paste("NetworkFrame",toString(Year),toString(Cutoff),sep='_')))
  ####Default Color Parameters
  color_high = 'purple'
  color_medium = 'grey'
  color_low = 'green'
  shape_high = 3
  shape_medium = 2
  shape_low = 1  
  ####Assign Percentile-Based Categories To CBR Values
  high_cutoff_can = quantile(TargetFrame$LogCBR,high_cutoff)
  low_cutoff_can = quantile(TargetFrame$LogCBR,low_cutoff)
  TargetFrame$CannibalismCategory <- cut(TargetFrame$LogCBR,
                                         breaks=c(min(TargetFrame$LogCBR)-.01,
                                                  low_cutoff_can,high_cutoff_can, max(TargetFrame$LogCBR+.01)),
                                         labels=c('Low CBR', 'Medium CBR', 'High CBR'))
  #####Calculate Means For Density Plot Lines
  mean_cannibal_df <- TargetFrame %>%
    group_by(CannibalismCategory) %>%
    summarise(meanCos=mean(Cosine_),meanOE = mean(ObsExp))
  ####Density Plots
  cos_can_density <- ggplot(TargetFrame, aes(x=Cosine_, fill=CannibalismCategory)) + 
    geom_density(alpha=.5) + 
    geom_vline(data = mean_cannibal_df, aes(xintercept = meanCos,color = CannibalismCategory), size=1.5)+
    scale_fill_manual(values = c(color_high,color_medium,color_low)) + 
    scale_colour_manual(values = c(color_high,color_medium,color_low)) + 
    theme_bw()+theme(legend.position = "empty")+xlab(paste(toString(Year),'Cosine Similarity (Equivalence)'))
  
  oe_can_density <- ggplot(TargetFrame, aes(y=ObsExp, fill=CannibalismCategory)) + 
    geom_density(alpha=.5) + 
    geom_hline(data = mean_cannibal_df, aes(yintercept = meanOE,color = CannibalismCategory), size=1.5)+
    scale_fill_manual(values = c(color_high,color_medium,color_low)) + 
    scale_colour_manual(values = c(color_high,color_medium,color_low)) + 
    theme_bw()+theme(legend.position = "empty")+ylab(paste(toString(Year),'Observed/Expected (Cohesion)'))
  
  ####Scatter Plot
  can_scatter <-ggplot(TargetFrame)+ 
    geom_point(aes(x=Cosine_, y=ObsExp,color = CannibalismCategory,shape=CannibalismCategory),size =2,alpha=0.75)+
    scale_color_manual(values = c(color_high,color_medium,color_low)) + 
    scale_shape_manual(values = c(shape_high,shape_medium,shape_low)) +
    theme_bw()+theme(legend.position = "empty") +xlab(paste(toString(Year),'Cosine Similarity (Equivalence)'))+ylab(paste(toString(Year),'Observed/Expected (Cohesion)'))
  
  ####Legend
  
  empty <- ggplot(TargetFrame)+geom_point(aes(colour=CannibalismCategory,shape=CannibalismCategory),x=NA,y=NA,size=4 )+
    scale_color_manual(name=paste(toString(Year),'to',toString(Year+1),':\nCannibalistic Diffusion Rate'),values = c(color_high,color_medium,color_low)) +
    scale_shape_manual(name=paste(toString(Year),'to',toString(Year+1),':\nCannibalistic Diffusion Rate'),values = c(shape_high,shape_medium,shape_low)) +
    theme(axis.ticks=element_blank(), 
          panel.background=element_blank(), 
          axis.text.x=element_blank(), axis.text.y=element_blank(),           
          axis.title.x=element_blank(), axis.title.y=element_blank())+
    labs(color= paste(toString(Year),'to',toString(Year+1),':\nCannibalistic Diffusion Rate'))
  grid.arrange(cos_can_density, empty, can_scatter, oe_can_density, ncol=2, nrow=2, widths=c(5, 2), heights=c(2, 5))
}


### DEFINE LANGAUGES CUTOFF IN EQUIVALENCE CALCULATIONS ###

NETWORK_PERCENTAGE_CUTOFF = 0.005

### ITERATE THROUGH AND CREATE FRAMES YEAR BY YEAR FOR THE GIVEN CUTOFF ###

for (TARGET_YEAR in 2010:2018){
  print(TARGET_YEAR)
  calculated_frame = create_frame_simp(NetworkFrames_All,TARGET_YEAR,NETWORK_PERCENTAGE_CUTOFF)
  assign(paste("NetworkFrame",toString(TARGET_YEAR),toString(NETWORK_PERCENTAGE_CUTOFF),sep='_'),calculated_frame)
}

### IMPLEMENT VISUALIZATIONS ###

piggyback_plot(2014,0.005,0.75,0.25)
cannibal_plot(2014,0.005,0.75,0.25)