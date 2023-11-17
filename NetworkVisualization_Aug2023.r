#####NETWORK VISUALIZATION CODE######

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
library(ggraph)

##########################################################
########## 0. LOAD DATA (Same as Main Analysis File)######
##########################################################

setwd("C:/Users/adsir/Dropbox/GitHub")
NetworkFrames_All<- read.csv('IntersectionFramesAll.csv',header=TRUE)
NetworkFrames_All$Lang2Perc <- NetworkFrames_All$L2Y1/NetworkFrames_All$BaseSample

##############################################################
########## 1. FUNCTIONS (Copied from RemodelingJune2023.r)####
##############################################################

create_frame_simp<-function(TheFrame,target_year,cutoff){
  netframe_year = TheFrame[which(TheFrame$Year==target_year),]
  netframe_cutoff = netframe_year[((netframe_year$Lang1Perc>cutoff) & (netframe_year$Lang2Perc>cutoff)),]
  matrix_cutoff = cosineify(netframe_cutoff,'ObsExp','Cosine_')
  netframe_year_final = merge(netframe_cutoff,matrix_cutoff,by=c('Lang1','Lang2'))
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

#####################################################
########## 2. CREATE NETWORK DATAFRAME ##############
#####################################################

####Select the year and minimum percentage for inclusion.
VIZ_YEAR = 2014
VIZ_CUTOFF = 0.02
Netframe_Year_VIZ<-create_frame_simp(NetworkFrames_All,VIZ_YEAR,VIZ_CUTOFF)

LangaugesWeights_VIZ = unique(Netframe_Year_VIZ[,c('Lang1','L1Y1')])
LangaugesWeights_VIZ$LogVolume = log(LangaugesWeights_VIZ$L1Y1)

variable_subset_VIZ <-selectVariables(LangaugesWeights_VIZ,'Lang1','L1Y1',10000)
Netframe_Year_Filtered_VIZ <- filterFrame(Netframe_Year_VIZ,variable_subset_VIZ)

######################################################
#########  3. CREATE NETWORK NODES AND EDGES #########
######################################################

#####CREATE VECTOR FOR EDGE CATEGORIZATION
high_cutoff_net = 0.75
low_cutoff_net = 0.25

###CREATE COLUMNS WITH NEW NEWTORK VISUALIZATION DETAIILS

###Creating cutoff values for high and low for variables
######Dependent Variables
high_cutoff_ObsExp_net = quantile(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'ObsExp')),high_cutoff_net)
low_cutoff_ObsExp_net = quantile(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'ObsExp')),low_cutoff_net)
high_cutoff_Cosine_net = quantile(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'Cosine')),high_cutoff_net)
low_cutoff_Cosine_net = quantile(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'Cosine')),low_cutoff_net)

######Independent Variables
high_cutoff_PBR_net = quantile(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'LogPBR')),high_cutoff_net)
low_cutoff_PBR_net = quantile(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'LogPBR')),low_cutoff_net)
high_cutoff_CBR_net = quantile(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'LogCBR')),high_cutoff_net)
low_cutoff_CBR_net = quantile(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'LogCBR')),low_cutoff_net)


####OBS EXP CATEGORIES
Netframe_Year_Filtered_VIZ$LogPBR_Color<- cut(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'LogPBR')),
                                                 breaks=c(min(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'LogPBR')))-.01,
                                                          low_cutoff_PBR_net,high_cutoff_PBR_net,
                                                          max(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'LogPBR')))+.01),
                                                 labels=c('(1) High PBR','(2) Medium PBR','(3) Low PBR'))

Netframe_Year_Filtered_VIZ$LogCBR_Color<- cut(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'LogCBR')),
                                              breaks=c(min(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'LogCBR')))-.01,
                                                       low_cutoff_CBR_net,high_cutoff_CBR_net,
                                                       max(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'LogCBR')))+.01),
                                              labels=c('(1) High CBR','(2) Medium CBR','(3) Low CBR'))

Netframe_Year_Filtered_VIZ$ObsExp_Color<- cut(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'ObsExp')),
                                                 breaks=c(min(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'ObsExp')))-.01,
                                                          low_cutoff_ObsExp_net,high_cutoff_ObsExp_net,
                                                          max(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'ObsExp')))+.01),
                                                 labels=c('(1) High Cohesion','(2) Medium Cohesion','(3) Low Cohesion'))

Netframe_Year_Filtered_VIZ$Cosine_Color<- cut(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'Cosine')),
                                                 breaks=c(min(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'Cosine')))-.01,
                                                          low_cutoff_Cosine_net,high_cutoff_Cosine_net,
                                                          max(as.numeric(as.character(Netframe_Year_Filtered_VIZ$'Cosine')))+.01),
                                                 labels=c('(1) High Equivalence','(2) Medium Equivalence','(3) Low Equivalence'))

###############################################################
#### 4. CREATE GRAPH AND ASSIGN NODE AND EDGE ATTRIBUTES ######
###############################################################

###Graph 1 Has Underlying Structure of ObsExp
Netframe_Year_Filtered_VIZ_OE <- Netframe_Year_Filtered_VIZ
Netframe_Year_Filtered_VIZ_OE['weight'] <- Netframe_Year_Filtered_VIZ_OE['ObsExp']
viz_edges_OE = Netframe_Year_Filtered_VIZ_OE[c('Lang1','Lang2','weight')]

graph1 <- graph_from_data_frame(viz_edges_OE, directed = TRUE, vertices = NULL)
E(graph1)$weight <- as.numeric(as.character(Netframe_Year_Filtered_VIZ_OE$'ObsExp'))
E(graph1)$weight_color_CBR <- as.numeric(as.character(Netframe_Year_Filtered_VIZ_OE$'LogCBR'))
E(graph1)$weight_color_PBR <- as.numeric(as.character(Netframe_Year_Filtered_VIZ_OE$'LogPBR'))
V(graph1)$node_size = sqrt(as.numeric(as.character(LangaugesWeights_VIZ$L1Y1)))/3

E(graph1)$weight_color_PBR_category <- as.character(Netframe_Year_Filtered_VIZ_OE$LogPBR_Color)
E(graph1)$weight_color_CBR_category <- as.character(Netframe_Year_Filtered_VIZ_OE$LogCBR_Color)
E(graph1)$weight_color_Cosine_category <- as.character(Netframe_Year_Filtered_VIZ_OE$Cosine_Color)
E(graph1)$weight_color_ObsExp_category <- as.character(Netframe_Year_Filtered_VIZ_OE$ObsExp_Color)


###Graph 2 Has Underlying Structure of Cosine
Netframe_Year_Filtered_VIZ_Cos <- Netframe_Year_Filtered_VIZ
Netframe_Year_Filtered_VIZ_Cos['weight'] <- Netframe_Year_Filtered_VIZ_Cos['Cosine']
viz_edges_Cos = Netframe_Year_Filtered_VIZ_Cos[c('Lang1','Lang2','weight')]

graph2 <- graph_from_data_frame(viz_edges_Cos, directed = TRUE, vertices = NULL)
E(graph2)$weight <- as.numeric(as.character(Netframe_Year_Filtered_VIZ_Cos$'Cosine'))
E(graph2)$weight_color_CBR <- as.numeric(as.character(Netframe_Year_Filtered_VIZ_Cos$'LogCBR'))
E(graph2)$weight_color_PBR <- as.numeric(as.character(Netframe_Year_Filtered_VIZ_Cos$'LogPBR'))
V(graph2)$node_size = sqrt(as.numeric(as.character(LangaugesWeights_VIZ$L1Y1)))/3
E(graph2)$cat_color_CBR <- as.numeric(as.character(Netframe_Year_Filtered_VIZ_Cos$'LogCBR'))

E(graph2)$weight_color_PBR_category <- as.character(Netframe_Year_Filtered_VIZ_Cos$LogPBR_Color)
E(graph2)$weight_color_CBR_category <- as.character(Netframe_Year_Filtered_VIZ_Cos$LogCBR_Color)
E(graph2)$weight_color_Cosine_category <- as.character(Netframe_Year_Filtered_VIZ_Cos$Cosine_Color)
E(graph2)$weight_color_ObsExp_category <- as.character(Netframe_Year_Filtered_VIZ_Cos$ObsExp_Color)


######################################################
########   4. BUILD AND PLOT NETWORKS ################
######################################################

####adjust_layout takes a network layout and stops nodes from crowding on top of each other ####

adjust_layout<-function(a_graph,iterations=500,min_dist=0.15,repel_dist=0.25){
  set.seed(1)
  coordinates <- cbind(a_graph$data$x,a_graph$data$y)
  node_count <- nrow(coordinates)
  for (move in 1:iterations){
    selected_nodes <- sample(1:node_count,2,replace=FALSE)
    x1=coordinates[selected_nodes[[1]],1]
    y1=coordinates[selected_nodes[[1]],2]
    x2=coordinates[selected_nodes[[2]],1]
    y2=coordinates[selected_nodes[[2]],2]
    #print(c(x1,y1,x2,y2))
    xdist = abs(x1-x2)
    ydist = abs(y1-y2)
    distance = sqrt(xdist**2+ydist**2)
    #print(distance)
    mid_x= (x1+x2)/2
    mid_y =(y1+y2)/2
    #print(distance,min_dist)
    distance_to_move = repel_dist - distance
    if (distance < min_dist){
      if(x1<x2){
        new_x1 = mid_x - 0.5*xdist*((repel_dist+.001)/(distance+.001))
        new_x2 = mid_x + 0.5*xdist*((repel_dist+.001)/(distance+.001))
      }
      else{
        new_x1 = mid_x + 0.5*xdist*((repel_dist+.001)/(distance+.001))
        new_x2 = mid_x - 0.5*xdist*((repel_dist+.001)/(distance+.001))
      }
      if(y1<y2){
        new_y1 = mid_y - 0.5*ydist*((repel_dist+.001)/(distance+.001))
        new_y2 = mid_y + 0.5*ydist*((repel_dist+.001)/(distance+.001))
      }
      else{
        new_y1 = mid_y + 0.5*ydist*((repel_dist+.001)/(distance+.001))
        new_y2 = mid_y - 0.5*ydist*((repel_dist+.001)/(distance+.001))
      }
    #else{}
    coordinates[selected_nodes[[1]],1]=new_x1
    coordinates[selected_nodes[[1]],2]=new_y1
    coordinates[selected_nodes[[2]],1]=new_x2
    coordinates[selected_nodes[[2]],2]=new_y2
    }
  }
    new_graph = a_graph
    new_graph$data$x = coordinates[,1]
    new_graph$data$y = coordinates[,2]
    return(new_graph)
}

###ObsExp Networks

#### VIZ SUBSTEP A1 - Generate Graph with 1st Layout and Color by PBR

g<-ggraph(graph1,type="spring") +
  geom_edge_fan(aes(width = abs(weight),edge_color = weight_color_PBR_category, edge_alpha=weight_color_PBR_category)) +
  #scale_edge_color_viridis()+
  scale_edge_color_manual(values=c("darkgreen","grey90","purple"))+
  scale_edge_alpha_manual(values=c(0.65,0.25,0.65))+
  geom_node_point(aes(size = node_size), color='cyan') +
  #geom_node_text(aes(label = name), color = "black", border="black",nudge_y=-.015,repel=FALSE,size = 9) +
  geom_node_label(aes(label = name), color = "black", border="black",nudge_y=-.015,nudge_x=0.015,alpha=.8,label.padding=0.25,repel=TRUE,size = 6) +
  theme_graph()

#### VIZ SUBSTEP A2 - Adjust Coordinates to force between for Distance and Save Coordinates

#adjust_layout(g)
#new_x = g$data$x
#new_y = g$data$y

#### VIZ SUBSTEP A3 - Create Visualization 1 (PBR) and Clean Up Legends

obsexp_PBR_net = adjust_layout(g,2500,0.15,0.2)+ guides(size = "none",edge_width="none",edge_alpha="none")+theme(legend.position = 'bottom')+labs(edge_colour='Piggybacking Diffusion Ratio')
obsexp_PBR_net = obsexp_PBR_net+ guides(size = "none",edge_width="none",edge_alpha="none")
obsexp_PBR_net

#### VIZ SUBSTEP A4 - Create Graph with 1st Layout and Color by CBR, Cos, ObsExp

g2<-ggraph(graph1,layout="eigen",type="adjacency") +
  geom_edge_fan(aes(width = abs(weight),edge_color = weight_color_CBR_category, edge_alpha=weight_color_CBR_category)) +
  #scale_edge_color_viridis()+
  scale_edge_color_manual(values=c("darkgreen","grey90","purple"))+
  scale_edge_alpha_manual(values=c(0.65,0.25,0.65))+
  geom_node_point(aes(size = node_size), color='cyan') +
  geom_node_label(aes(label = name), color = "black", border="black",nudge_y=-.015,nudge_x=0.015,alpha=.8,label.padding=0.25,repel=TRUE,size = 6) +
  theme_graph()

g3<-ggraph(graph1,layout="eigen",type="adjacency") +
  geom_edge_link(aes(width = abs(weight),edge_color = weight_color_Cosine_category, edge_alpha=weight_color_Cosine_category)) +
  #scale_edge_color_viridis()+
  scale_edge_color_manual(values=c("darkgreen","grey90","purple"))+
  scale_edge_alpha_manual(values=c(0.65,0.25,0.65))+
  geom_node_point(aes(size = node_size), color='cyan') +
  geom_node_label(aes(label = name), color = "black", border="black",nudge_y=-.015,nudge_x=0.015,alpha=.8,label.padding=0.25,repel=TRUE,size = 6) +
  theme_graph()

g4<-ggraph(graph1,layout="eigen",type="adjacency") +
  geom_edge_link(aes(width = abs(weight),edge_color = weight_color_ObsExp_category, edge_alpha=weight_color_ObsExp_category)) +
  #scale_edge_color_viridis()+
  scale_edge_color_manual(values=c("darkgreen","grey90","purple"))+
  scale_edge_alpha_manual(values=c(0.65,0.25,0.65))+
  geom_node_point(aes(size = node_size), color='cyan') +
  geom_node_label(aes(label = name), color = "black", border="black",nudge_y=-.015,nudge_x=0.015,alpha=.8,label.padding=0.25,repel=TRUE,size = 6) +
  theme_graph()

#### VIZ SUBSTEP A5 - Update Coordinates to Match PBR Network
g2$data$x = obsexp_PBR_net$data$x
g2$data$y = obsexp_PBR_net$data$y

g3$data$x = obsexp_PBR_net$data$x
g3$data$y = obsexp_PBR_net$data$y

g4$data$x = obsexp_PBR_net$data$x
g4$data$y = obsexp_PBR_net$data$y

#### VIZ SUBSTEP A6 - Create Remaining ObsExp Visualizations and Clean Up Legends

obsexp_CBR_net = g2+ guides(size = "none",edge_width="none",edge_alpha="none",edge="none")+theme(legend.position = 'bottom')+labs(edge_colour='Cannibalistic Diffusion Ratio')
#obsexp_CBR_net

obsexp_COS_net = g3+ guides(size = "none",edge_width="none",edge_alpha="none",edge="none")+theme(legend.position = 'bottom')+labs(edge_colour='Functional Equivalence')
#obsexp_COS_net

obsexp_OE_net = g4+ guides(size = "none",edge_width="none",edge_alpha="none",edge="none")+theme(legend.position = 'bottom')+labs(edge_colour='Functional Cohesion')
#obsexp_OE_net




##############COSINE NETWORKS################

#### VIZ STEP B1 - Generate Graph with 1st Layout and Color by PBR

g_cos<-ggraph(graph2,weights=weight) +
  geom_edge_fan(aes(width = abs(weight),edge_color = weight_color_PBR_category, edge_alpha=weight_color_PBR_category)) +
  scale_edge_color_manual(values=c("darkgreen","grey90","purple"))+
  scale_edge_alpha_manual(values=c(0.65,0.25,0.65))+
  geom_node_point(aes(size = node_size), color='cyan') +
  geom_node_label(aes(label = name), color = "black", border="black",nudge_y=-.015,nudge_x=0.015,alpha=.8,label.padding=0.25,repel=TRUE,size = 6) +
  theme_graph()

#### VIZ SUBSTEP B3 - Create Visualization 1 (PBR) and Clean Up Legends

cosine_PBR_net = adjust_layout(g_cos,2500,250,400)
cosine_PBR_net = cosine_PBR_net+ guides(size = "none",edge_width="none",edge_alpha="none")+theme(legend.position = 'bottom')+labs(edge_colour='Piggybacking Diffusion Ratio')
cosine_PBR_net


#### VIZ SUBSTEP B4 - Create Graph with 1st Layout and Color by CBR, Cos, ObsExp

g2_cos<-ggraph(graph2,layout="eigen",type="adjacency") +
  geom_edge_fan(aes(width = abs(weight),edge_color = weight_color_CBR_category, edge_alpha=weight_color_CBR_category)) +
  #scale_edge_color_viridis()+
  scale_edge_color_manual(values=c("darkgreen","grey90","purple"))+
  scale_edge_alpha_manual(values=c(0.75,0.25,0.75))+
  geom_node_point(aes(size = node_size), color='cyan') +
  geom_node_label(aes(label = name), color = "black", border="black",nudge_y=-.015,nudge_x=0.015,alpha=.8,label.padding=0.25,repel=TRUE,size = 6) +
  theme_graph()

g3_cos<-ggraph(graph2,layout="eigen",type="adjacency") +
  geom_edge_fan(aes(width = abs(weight),edge_color = weight_color_Cosine_category, edge_alpha=weight_color_Cosine_category)) +
  #scale_edge_color_viridis()+
  scale_edge_color_manual(values=c("darkgreen","grey90","purple"))+
  scale_edge_alpha_manual(values=c(0.75,0.25,0.75))+
  geom_node_point(aes(size = node_size), color='cyan') +
  geom_node_label(aes(label = name), color = "black", border="black",nudge_y=-.015,nudge_x=0.015,alpha=.8,label.padding=0.25,repel=TRUE,size = 6) +
  theme_graph()

g4_cos<-ggraph(graph2,layout="eigen",type="adjacency") +
  geom_edge_link(aes(width = abs(weight),edge_color = weight_color_ObsExp_category, edge_alpha=weight_color_ObsExp_category)) +
  #scale_edge_color_viridis()+
  scale_edge_color_manual(values=c("darkgreen","grey90","purple"))+
  scale_edge_alpha_manual(values=c(0.75,0.25,0.75))+
  geom_node_point(aes(size = node_size), color='cyan') +
  geom_node_label(aes(label = name), color = "black", border="black",nudge_y=-.015,nudge_x=0.015,alpha=.8,label.padding=0.25,repel=TRUE,size = 6) +
  theme_graph()


#### VIZ SUBSTEP B5 - Update Coordinates to Match PBR Network

g2_cos$data$x = cosine_PBR_net$data$x
g2_cos$data$y = cosine_PBR_net$data$y

g3_cos$data$x = cosine_PBR_net$data$x
g3_cos$data$y = cosine_PBR_net$data$y

g4_cos$data$x = cosine_PBR_net$data$x
g4_cos$data$y = cosine_PBR_net$data$y

#### VIZ SUBSTEP B6 - Create Remaining ObsExp Visualizations and Clean Up Legends

cosine_CBR_net = g2_cos + guides(size = "none",edge_width="none",edge_alpha="none",edge="none")+theme(legend.position = 'bottom')+labs(edge_colour='Cannibalistic Diffusion Ratio')
#cosine_CBR_net

cosine_COS_net = g3_cos + guides(size = "none",edge_width="none",edge_alpha="none",edge="none")+theme(legend.position = 'bottom')+labs(edge_colour='Functional Equivalence')
#cosine_COS_net

cosine_OE_net = g4_cos + guides(size = "none",edge_width="none",edge_alpha="none",edge="none")+theme(legend.position = 'bottom')+labs(edge_colour='Functional Cohesion')
#cosine_OE_net




####COMPOSITE VISUALIZATIONS####
figure <- ggarrange(obsexp_PBR_net,obsexp_CBR_net,
                    labels = c("",""),
                    ncol = 1, nrow = 2)
figure

figure <- ggarrange(cosine_PBR_net,cosine_CBR_net,
                    labels = c("",""),
                    ncol = 1, nrow = 2)
figure
