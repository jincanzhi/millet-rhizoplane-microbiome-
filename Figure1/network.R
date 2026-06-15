pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(ggplot2)
library(ggClusterNet)
library(igraph)
library(microeco)

dir.create('network')

set.seed(113)

MAG_table <- read.csv('./data/MAG_abundance_sum.csv', row.names = 1,check.names = F)
metadata<-read.csv(file="./data/metadata.csv",row.names = 1,header=T)


##RP_network
RP_table <- MAG_table[1:17]

MAG_dataset_all <- microtable$new(sample_table = metadata, otu_table = RP_table)
MAG_dataset_all$tidy_dataset()
MAG_dataset_all

MAG_network_All <- trans_network$new(dataset = MAG_dataset_all, cor_method = 'spearman',use_WGCNA_pearson_spearman = T, filter_thres = 0, nThreads = 18)
write.csv(MAG_network_All$res_cor_p$cor, './network/RP_WGCNA_spearman_cor_raw.csv')
write.csv(MAG_network_All$res_cor_p$p, './network/RP_WGCNA_spearman_p_raw.csv')

MAG_network_All$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
MAG_network_All$cal_module(method = 'cluster_fast_greedy')
MAG_network_All$save_network(filepath = './network/RP_WGCNA_spearman_network.gexf')

MAG_network_All$cal_network_attr()
write.csv(MAG_network_All$res_network_attr, './network/RP_network_property_All.csv')
MAG_network_All$get_node_table()
write.csv(MAG_network_All$res_node_table, './network/RP_network_node_property_All.csv', row.names = F)
MAG_network_All$get_edge_table()
write.csv(MAG_network_All$res_edge_table, './network/RP_network_edge_property_All.csv', row.names = F)


MAG_network_cor_all <- read.csv('./network/RP_WGCNA_spearman_cor_raw.csv', row.names = 1)
MAG_network_node_all <- read.csv('./network/RP_network_node_property_All.csv')
MAG_network_edge_all <- read.csv('./network/RP_network_edge_property_All.csv')

network_igraph <- graph_from_data_frame(MAG_network_edge_all, directed = F, vertices = MAG_network_node_all)
network_net_property <- net_properties(network_igraph)

network_node_property <- node_properties(network_igraph)

network_zipi <- ZiPiPlot(igraph = network_igraph, method = 'cluster_fast_greedy')
ggsave('./network/RP_network_zipi.pdf', network_zipi[[1]], width = 8, height = 6)
network_zipi[[1]]

network_hub_info <- data.frame(network_node_property,
                               network_zipi[[2]][row.names(network_node_property),])
write.csv(network_hub_info, './network/RP_network_hub_info.csv')



##RS_network
RS_table <- MAG_table[18:34]

MAG_dataset_all <- microtable$new(sample_table = metadata, otu_table = RS_table)
MAG_dataset_all$tidy_dataset()
MAG_dataset_all

MAG_network_All <- trans_network$new(dataset = MAG_dataset_all, cor_method = 'spearman',use_WGCNA_pearson_spearman = T, filter_thres = 0, nThreads = 18)
write.csv(MAG_network_All$res_cor_p$cor, './network/RS_WGCNA_spearman_cor_raw.csv')
write.csv(MAG_network_All$res_cor_p$p, './network/RS_WGCNA_spearman_p_raw.csv')

MAG_network_All$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
MAG_network_All$cal_module(method = 'cluster_fast_greedy')
MAG_network_All$save_network(filepath = './network/RS_WGCNA_spearman_network.gexf')

MAG_network_All$cal_network_attr()
write.csv(MAG_network_All$res_network_attr, './network/RS_network_property_All.csv')
MAG_network_All$get_node_table()
write.csv(MAG_network_All$res_node_table, './network/RS_network_node_property_All.csv', row.names = F)
MAG_network_All$get_edge_table()
write.csv(MAG_network_All$res_edge_table, './network/RS_network_edge_property_All.csv', row.names = F)


MAG_network_cor_all <- read.csv('./network/RS_WGCNA_spearman_cor_raw.csv', row.names = 1)
MAG_network_node_all <- read.csv('./network/RS_network_node_property_All.csv')
MAG_network_edge_all <- read.csv('./network/RS_network_edge_property_All.csv')

network_igraph <- graph_from_data_frame(MAG_network_edge_all, directed = F, vertices = MAG_network_node_all)
network_net_property <- net_properties(network_igraph)

network_node_property <- node_properties(network_igraph)

network_zipi <- ZiPiPlot(igraph = network_igraph, method = 'cluster_fast_greedy')
ggsave('./network/RS_network_zipi.pdf', network_zipi[[1]], width = 8, height = 6)
network_zipi[[1]]

network_hub_info <- data.frame(network_node_property,
                               network_zipi[[2]][row.names(network_node_property),])
write.csv(network_hub_info, './network/RS_network_hub_info.csv')




##All_network
MAG_dataset_all <- microtable$new(sample_table = metadata, otu_table = MAG_table)
MAG_dataset_all$tidy_dataset()
MAG_dataset_all

MAG_network_All <- trans_network$new(dataset = MAG_dataset_all, cor_method = 'spearman',use_WGCNA_pearson_spearman = T, filter_thres = 0, nThreads = 18)
write.csv(MAG_network_All$res_cor_p$cor, './network/Whole_WGCNA_spearman_cor_raw.csv')
write.csv(MAG_network_All$res_cor_p$p, './network/Whole_WGCNA_spearman_p_raw.csv')

MAG_network_All$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
MAG_network_All$cal_module(method = 'cluster_fast_greedy')
MAG_network_All$save_network(filepath = './network/Whole_WGCNA_spearman_network.gexf')

MAG_network_All$cal_network_attr()
write.csv(MAG_network_All$res_network_attr, './network/Whole_network_property_All.csv')
MAG_network_All$get_node_table()
write.csv(MAG_network_All$res_node_table, './network/Whole_network_node_property_All.csv', row.names = F)
MAG_network_All$get_edge_table()
write.csv(MAG_network_All$res_edge_table, './network/Whole_network_edge_property_All.csv', row.names = F)


MAG_network_cor_all <- read.csv('./network/Whole_WGCNA_spearman_cor_raw.csv', row.names = 1)
MAG_network_node_all <- read.csv('./network/Whole_network_node_property_All.csv')
MAG_network_edge_all <- read.csv('./network/Whole_network_edge_property_All.csv')

network_igraph <- graph_from_data_frame(MAG_network_edge_all, directed = F, vertices = MAG_network_node_all)
network_net_property <- net_properties(network_igraph)

network_node_property <- node_properties(network_igraph)

network_zipi <- ZiPiPlot(igraph = network_igraph, method = 'cluster_fast_greedy')
ggsave('./network/Whole_network_zipi.pdf', network_zipi[[1]], width = 8, height = 6)
network_zipi[[1]]

network_hub_info <- data.frame(network_node_property,
                               network_zipi[[2]][row.names(network_node_property),])
write.csv(network_hub_info, './network/Whole_network_hub_info.csv')

