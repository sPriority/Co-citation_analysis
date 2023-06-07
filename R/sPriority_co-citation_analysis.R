### Synthesis Paper: Co-citation Analysis
### Authors: Benjamin M Delory & Elle M Barnes
### Last Updated: June 2, 2023

library(curl)
library(readr)
library(tidyverse)
library(readxl)
library(bibliometrix)
library(ggpubr)
library(dplyr)
library(igraph)
library(brainGraph)
library(RCy3)
library(ggrepel)

##### Set path #####

setwd("C:/Users/Delory/Documents/GitHub/Co-citation_analysis")

##### Generate Co-citation network #####

M<-convert2df(file=c("Data/Results_WoS_1-500.txt",
                      "Data/Results_WoS_501-902.txt"),
              format = "plaintext")

NetMatrix <- biblioNetwork(M, 
                           analysis = "co-citation", 
                           network = "references",
                           shortlabel = FALSE)


set.seed(123) #For reproducibility

net=networkPlot(NetMatrix, 
                n = 1000, 
                Title = "", 
                type = "fruchterman",
                label=FALSE,
                labelsize=1,
                label.cex=FALSE,
                label.color=FALSE,
                label.n=0,
                halo=FALSE,
                cluster="louvain",
                community.repulsion = 0.1,
                size=1,
                size.cex=FALSE,
                remove.multiple=TRUE,
                remove.isolates=TRUE,
                edges.min=2,
                verbose=FALSE)

#Create dataset for plotting results

data_plot<-cbind(net$cluster_res, net$layout)

#Add color information

data_plot$color<-igraph::V(net$graph)$color
                              
colnames(data_plot)[6:7]<-c("x", "y")

#Plot results using cluster as a grouping factor 
    # cluster <- ggplot(data_plot, aes(x=x, y=y, col=factor(cluster)))+
    #     theme_bw()+
    #     geom_point()+
    #   ggtitle("Colored by cluster")

#Plot results using color as a grouping factor 
    # color.fac <- ggplot(data_plot, aes(x=x, y=y, col=factor(color)))+
    #   theme_bw()+
    #   geom_point()+
    #   ggtitle("Color by V$color")
    
    # ggarrange(cluster, color.fac,
    #          ncol=2, nrow=1)

#Export high resolution figure

pdf("co-citation-network_1.pdf",
    width=15/cm(1),
    height=15/cm(1))

igraph::plot.igraph(net$graph)
legend("topleft", 
       legend=unique(data_plot$cluster), 
       col=unique(igraph::V(net$graph)$color),
       pch=16,
       bty="n")

dev.off()

#Export network to Cytoscape

cytoscapePing()

#Export network to Cytoscape
createNetworkFromIgraph(net$graph,"co-citation")

#Export table with clustered references

write.csv(data_plot,
          "co-citation-network_2.csv",
          row.names = FALSE)

##### Calculate network attributes #####

data.clust <- net$cluster_res

partco <- part_coeff(net$graph, memb = data.clust$cluster) # participation coefficient
within.moddeg <- within_module_deg_z_score(net$graph, memb = data.clust$cluster) # within module degree   

rn.brainGraph <- as.data.frame(do.call(cbind, list(participation = partco, connectivity = within.moddeg)))
name <- rownames(rn.brainGraph)
rownames(rn.brainGraph) <- NULL
data <- cbind(name,rn.brainGraph)

# assign network roles based on Guimera & Amaral (2005)

data$role <- with(data, ifelse(connectivity > 2.5, 'Hub',
                               ifelse(connectivity < 2.5 & participation < 0.05, 'Ultra-peripheral', 
                                      ifelse(connectivity < 2.5 & participation < 0.62, 'Peripheral',
                                             ifelse(connectivity < 2.5 & participation < 0.8, 'Connector', 'Kinless')))))

write.csv(data, "cocitation.brainGraph.csv") # export results


categories <- c("Connector", "Hub", "Peripheral", "Ultra-peripheral")
color <- c("#a7bf9c", "#33454e", "#81a9ad", "#e2e2d2")

# plot each citation by connectivity x participation & label hubs

ggplot(data, aes(x=participation, y=connectivity, color=role)) +
  geom_point(size=2) +
  scale_color_manual(breaks=categories, values=color)+
  ylab("within-module degree") +
  xlab("among-module connectivity \n (participation coefficient)") +
  theme_classic() +
  geom_hline(yintercept = 2.5, linetype = "dashed") +
  geom_vline(xintercept = 0.62, linetype = "dashed")+
  geom_label_repel(aes(label=ifelse(connectivity>2.5,as.character(name),'')),
                   size = 1.7,
                   segment.color = 'grey50')

# export image 

pdf("cocitation.roles_1.pdf",
    width=11,
    height=8)

ggplot(data, aes(x=participation, y=connectivity, color=role)) +
  geom_point(size=2) +
  scale_color_manual(breaks=categories, values=color)+
  ylab("within-module degree") +
  xlab("among-module connectivity \n (participation coefficient)") +
  theme_classic() +
  geom_hline(yintercept = 2.5, linetype = "dashed") +
  geom_vline(xintercept = 0.62, linetype = "dashed")+
  geom_label_repel(aes(label=ifelse(connectivity>2.5,as.character(name),'')),
                   size = 1.7,
                   segment.color = 'grey50')

dev.off()


# identify connections within vs across modules

el <- data.frame(as_edgelist(net$graph)) # extract list of edges (pairs of vertices)
data.clust <- net$cluster_res

el <- el %>% 
  rename("vertex" = "X1")
el.join1 <- inner_join(el, data.clust, by = "vertex") # merge cluster info with vertex name for v1 (aka source)
el.join1 <- el.join1 %>% 
  rename("V1" = "vertex",
         "V1.clust" = "cluster")

el.join1 <- el.join1 %>% 
  rename("vertex" = "X2")
el.join2 <- inner_join(el.join1, data.clust, by = "vertex") # merge cluster info with vertex name for v2 (aka target)
el.join.final <- el.join2 %>% 
  rename("V2" = "vertex",
         "V2.clust" = "cluster")

# identify if co-citation is between papers within same cluster or across clusters
el.join.final$within.v.across <- with(el.join.final, ifelse(V1.clust == V2.clust, 'Within', 'Across'))

# write.csv(el.join.final, "cocitation.edges.csv") # export results


##### Stacked bar charts of within vs. across cluster co-citations #####

el.join.final$No.Edge <- c(1) # add value of 1 for each row = 1 edge

connect.clust <- el.join.final %>%
  mutate(V1.clust = recode(V1.clust, '1' = "Theory, etc.", '2' = "Animals", '3' =  "Plants", '4' = "Evolution", '5' = "Parasitology", '6' = "Polar" ))

connect.clust <- connect.clust %>%
  mutate(V2.clust = recode(V2.clust, '1' = "Theory, etc.", '2' = "Animals", '3' =  "Plants", '4' = "Evolution", '5' = "Parasitology", '6' = "Polar" ))

### Across vs Within by cluster

connect.clust.sum <- connect.clust %>%
  group_by(V1.clust, within.v.across) %>%
  summarise(sum = sum(No.Edge))

connect.clust %>%
  group_by(V1.clust) %>%
  summarise(sum = sum(No.Edge)) # to get total number of co-citations by cluster for x-axis annotation


# stacked bar plot
ggplot(connect.clust.sum, aes(x = factor(V1.clust, levels = c("Theory, etc.", "Animals", "Plants", "Evolution", "Parasitology", "Polar")),
                                     y = sum, fill=within.v.across)) + 
  scale_fill_manual(values=c("#c7e9b4", "#2c7fb8")) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  labs(y = "Proportion of Co-citations")+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position="top",
        legend.justification="right")+
  scale_x_discrete(labels=c("Theory\n (15,849)", "Animals\n (3,018)", "Plants\n (10,061)", "Evolution\n (2,356)", "Parasitol.\n (399)", "Polar\n (3,086)"))

### Edge identity by cluster

connect.clust.sum2 <- connect.clust %>%
  group_by(V1.clust, V2.clust) %>%
  summarise(sum = sum(No.Edge))

print(tbl_df(connect.clust.sum2), n=25) # show full df

connect.clust.sum2.across <- connect.clust.sum2[-c(1,7,13,18,25),] # remove rows that show total within cluster, except Polar

connect.clust.sum2.across %>%
  group_by(V1.clust) %>%
  summarise(sum = sum(sum)) # to get total number of across-cluster co-citations by cluster for x-axis annotation

# stacked bar plot
ggplot(connect.clust.sum2.across, aes(x = factor(V1.clust, levels = c("Theory, etc.", "Animals", "Plants", "Evolution", "Parasitology", "Polar")),
                          y = sum, fill=factor(V2.clust, levels = c("Theory, etc.", "Animals", "Plants", "Evolution", "Parasitology", "Polar")))) + 
  scale_fill_manual(values=c("#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#0c2c84")) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  labs(y = "Proportion of Co-citations")+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position="top",
        legend.justification="right") +
  scale_x_discrete(labels=c("Theory\n (5,705)", "Animals\n (1,309)", "Plants\n (1,415)", "Evolution\n (693)", "Parasitol.\n (24)", "Polar\n (0)"))
