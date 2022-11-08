
library(dplyr)
library(ggplot2)
library(corrplot)

GSE <- read.table("GSE19860_series_matrix 2.txt", header = T, fill = T, skip = 61, sep = "\t")
GSE <- GSE[-54676,]

GPL <- read.table("GPL570-55999.txt", fill = T, quote = '\\', header = T, skip = 16, sep='\t')

data <- GSE
data$gene_name <- GPL$Gene.Symbol[match(GPL$ID,GSE$ID_REF)] # объединенный датасет из двух файлов

length(unique(data$gene_name)) #количество уникальных названий генов

# Анализ экспрессии одного гена (APC)
apc_data <- data[data$gene_name == 'APC', 2:41] #датафрейм, содержащий строки с геном 'APC'
apc_mean_df<-as.numeric(sapply(apc_data, FUN=mean)) # для каждого пациента ищем среднее значение уровня экспрессии гена 'APC'

# Function hist гистограмма

hist(apc_mean_df, breaks = 30, col = 'purple', xlab = 'Gene expression level', main = 'Histogram of expression level of gene "APC"')

# Function qplot cxcl1
CXCL1_data <- data[data$gene_name == 'CXCL1', 2:41] #датафрейм, содержащий строки с геном 'APC'
qplot(apc_mean_df, CXCL1_numeric, geom = c("point", "line"), main = 'Dependence APC expression level from CXCL1 gene expresion', xlim = 'APC', xlab = 'APC gene expression level', ylab = 'CXCL1 gene expression level')
#xlab ylab названия осей

# Function ggplot2 + geom histogram

apc_df <- as.data.frame(apc_mean_df)
colnames(apc_df) <- 'Gene_expression_value'

ggplot(apc_df, aes(x = Gene_expression_value)) + geom_histogram(binwidth = 3, color="yellow", fill="pink") + labs(title="Gene expression level of 'APC' histogram plot", x ="Gene expression level of APC", y = "Frequency")

ggplot(apc_df, aes(x = Gene_expression_value)) + geom_density(color="black", fill="pink") + labs(title="Gene expression level of 'APC' density plot",x="Gene expression level of APC", y = "Frequency")

# Зависимость экспрессии одного гена от другого

CXCL1_data <- data[data$gene_name == 'CXCL1', 2:41] # датафрейм, содержащий строки с геном 'APC'
CXCL1_numeric<-as.numeric(CXCL1_data) # для каждого пациента ищем среднее значение уровня экспрессии гена 'APC'
df_with_2_genes <- as.data.frame(cbind(CXCL1_numeric, apc_mean_df)) # создаем таблицу с уровнем экспрессии двух генов 'APC' и 'CXCL1'
colnames(df_with_2_genes) <- c('CXCL1', 'APC')

# график зависимости экспрессии одного гена от другого

ggplot(df_with_2_genes, aes(x = CXCL1, y= APC)) + geom_point(aes(colour= APC)) + ggtitle("My scatter plot") + labs(title="Dependence of gene expression level of 'CXCL1' from 'APC' scatter plot",x="CXCL1", y = "APC") 

#Таблица c  AAK1  ZMYM6 POLH
AAK1_data <- data[data$gene_name == 'AAK1', 2:41] # создаем датафрейм, содержащий строки с геном 'AAK1'
AAK1_mean_df<-as.data.frame(sapply(AAK1_data, FUN=mean)) # для каждого пациента ищем среднее значение уровня экспрессии гена 'APC'
AAK1_vector <- as.data.frame(rep('AAK1',40)) #повторить название 40 раз тк у нас 40 строчек
AAK1_df <- cbind(AAK1_vector,AAK1_mean_df)
AAK1_df$sample_names <- rownames(AAK1_df)
colnames(AAK1_df) <- c('gene_name','gene_expression_level', 'sample_name')

ZMYM6_data <- data[data$gene_name == 'ZMYM6', 2:41] # датафрейм, содержащий строки с геном 'ZMYM6'
ZMYM6_mean_df<-as.data.frame(sapply(ZMYM6_data, FUN=mean)) # для каждого пациента ищем среднее значение уровня экспрессии гена 'APC'
ZMYM6_vector <- as.data.frame(rep('ZMYM6',40))
ZMYM6_df <- cbind(ZMYM6_vector,ZMYM6_mean_df)
ZMYM6_df$sample_names <- rownames(ZMYM6_df)
colnames(ZMYM6_df) <- c('gene_name','gene_expression_level', 'sample_name')

POLH_data <- data[data$gene_name == 'POLH', 2:41] #датафрейм, содержащий строки с геном 'POLH'
POLH_mean_df<-as.data.frame(sapply(POLH_data, FUN=mean)) # для каждого пациента ищем среднее значение уровня экспрессии гена 'APC'
POLH_vector <- as.data.frame(rep('POLH',40))
POLH_df <- cbind(POLH_vector,POLH_mean_df)
POLH_df$sample_names <- rownames(POLH_df)
colnames(POLH_df) <- c('gene_name','gene_expression_level', 'sample_name')

data_3_genes <- rbind(AAK1_df,ZMYM6_df,POLH_df)

#boxplot 3 genes

boxplot(log10(gene_expression_level) ~ gene_name, data_3_genes, xlab = 'Gene name', ylab = 'Gene expression level', main="Boxplot of gene expresion level of different genes", col.main="purple", col="pink")

# ggplot 3 genes

ggplot(data_3_genes, aes(x= gene_expression_level, y=sample_name, color =gene_name )) + geom_point() + xlab("Gene expression level") + ylab("Sample name") + ggtitle('Gene expression level')

ggplot(data_3_genes, aes(x = gene_expression_level,  y = sample_name, color = gene_name)) + geom_boxplot() + ggtitle('Gene expression level') + xlab("Gene expression level") + ylab("Sample name")

# корреляцию Спирмена между образцами по первым 5000 генам 
data_unique <- data %>% distinct(gene_name, .keep_all=TRUE) # датасет с уникальными значениями переменной 'gene_name'

data_5000 <- data_unique[1:21,2:41]

cor_res <-  cor(data_5000)
corrplot(order= 'hclust', cor_res, method = 'circle')
corrplot( order= 'hclust',method = 'number')

cor_res1 <-  cor(data_5000[1:10,1:10]) #таблица с меньшим колличеством данных
corrplot(cor_res1, method = 'circle')

corrplot.mixed(cor_res1, order = 'alphabet') 

# extract metadata
md <- read.csv("GSE19860_series_matrix 2.txt", header = T, fill = T, sep = "\t", skip = 27)[1:32,]
table(unlist(md[10,-1]))
number_of_patients<-ncol(md)-1
numbers_for_FL_non_responders <- which(md[10,-1]=="treatment response: FL_Non_responder")
numbers_for_FL_BV_nonresponders <- which(md[10,-1]=="treatment response: FL_Non_responder, BV_Non_reponder")
numbers_for_FL_nonresponders_BV_responder <- which(md[10,-1]=="treatment response: FL_Non_responder, BV_Responder")
numbers_for_FL_responders <- which(md[10,-1]=="treatment response: FL_Responder")
numbers_for_FL_responders_BV_non_responder <- which(md[10,-1]=="treatment response: FL_Responder, BV_Non_reponder")
numbers_for_FL_responders_BV_responder <- which(md[10,-1]=="treatment response: FL_Responder, BV_Responder")

color_vector<-rep("",number_of_patients)
color_vector[numbers_for_FL_non_responders]<-"green"
color_vector[numbers_for_FL_BV_nonresponders]<-"red"
color_vector[numbers_for_FL_nonresponders_BV_responder]<-"pink"
color_vector[numbers_for_FL_responders]<-"purple"
color_vector[numbers_for_FL_responders_BV_non_responder ]<-"blue"
color_vector[numbers_for_FL_responders_BV_responder]<-"orange"

#3d PCA visualization
p<-prcomp(t(data_unique [1:20,2:41]))
t(data_unique [1:20,2:41])
pov<-p$sdev^2/sum(p$sdev^2)
s3d <- scatterplot3d(p$x[,1], p$x[,2], p$x[,3], angle=60, pch = 16, type="p",xlab=paste0("PC1 (",round(pov[1]*100,digits = 2),"%)"),ylab=paste0("PC2 (",round(pov[2]*100,digits = 2),"%)"),zlab=paste0("PC3 (",round(pov[3]*100,digits = 2),"%)"),cex.symbols = 0.5, cex.axis = 0.6, cex.lab = 0.8, color = color_vector)

length(color_vector)
dim(p$x)
#2D
plot(p$x[,1], p$x[,2], pch = 19, type="p",xlab=paste0("PC1 (",round(pov[1]*100,digits = 2),"%)"),ylab=paste0("PC2 (",round(pov[2]*100,digits = 2),"%)"),cex.symbols = 0.5, cex.axis = 0.6, cex.lab = 0.8, col = color_vector, mar=c(1,2,1,1),main="PCA visualization Genes")

#heatmaps
df2 <- scale(data_unique [1:20,2:41])
heatmap(df2, scale = "none")
col<- colorRampPalette(c("red", "white", "blue"))(256)
library("RColorBrewer")
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
# Use RColorBrewer color palette names
library("RColorBrewer")
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)

library("pheatmap")
pheatmap(df2, cutree_rows = 4)
library("d3heatmap")
d3heatmap(scale(data_unique), colors = "RdYlBu",
          k_row = 4, # Number of groups in rows
          k_col = 2 # Number of groups in columns
)
library(dendextend)
# order for rows
Rowv  <- df2 %>% scale %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k = 3) %>% set("branches_lwd", 1.2) %>%
  ladderize
# Order for columns: We must transpose the data
Colv  <- df2 %>% scale %>% t %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k = 2, value = c("orange", "blue")) %>%
  set("branches_lwd", 1.2) %>%
  ladderize
library("d3heatmap")
d3heatmap(scale(df2), colors = "RdBu",
          Rowv = Rowv, Colv = Colv)
if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
Heatmap(df2, 
        name = "genes", #title of legend
        column_title = "Expression", row_title = "Probes",
        row_names_gp = gpar(fontsize = 7) 
        )        
library(circlize)
mycols <- colorRamp2(breaks = c(-2, 0, 2), 
                     colors = c("green", "white", "red"))
Heatmap(df2, name = "expressions", col = mycols)
library("circlize")
library("RColorBrewer")
Heatmap(df2, name = "genes",
        col = colorRamp2(c(-2, 0, 2), brewer.pal(n=3, name="RdBu")))

# Heatmap 1
ht1 = Heatmap(df2, name = "ht1", km = 2,
              column_names_gp = gpar(fontsize = 9))
# Heatmap 2
ht2 = Heatmap(df2, name = "ht2", 
              col = circlize::colorRamp2(c(-2, 0, 2), c("yellow", "white", "black")),
              column_names_gp = gpar(fontsize = 9))

ht1 + ht2

# Compute distances and hierarchical clustering
dd <- dist(scale(t(GSE [1:2000, 2:41])), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

plot(hc)
plot(hc,  cex = 0.6)
plot(hc, hang = -1, cex = 0.6)

hcd <- as.dendrogram(hc)
# Default plot
plot(hcd, type = "rectangle", ylab = "Probes")
# Triangle plot
plot(hcd, type = "triangle", ylab = "Probest")
# Zoom in to the first dendrogram
plot(hcd, xlim = c(1, 20), ylim = c(1,8))
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "pink")
# Customized plot; remove labels
plot(hcd, ylab = "Height", nodePar = nodePar, leaflab = "none")
# Horizontal plot
plot(hcd,  xlab = "Height",
     nodePar = nodePar, horiz = TRUE)
plot(hcd,  xlab = "Probes", nodePar = nodePar, 
     edgePar = list(col = 2:3, lwd = 2:1))


DDR1 <- c(unlist(data_unique[1,-1]))
DDR1_RFC2 <- c(unlist(data_unique[2,-1]))

#Student | Wilcox test
t.test(DDR1,DDR1_RFC2, TRUE)
wilcox.test(DDR1,DDR1_RFC2)

pv <- NA [1:10]
pv_w <- NA[1:10]
for (i in 1:100){
  gene_expr_for_resp<-unlist(aggregated_expr[i,numbers_for_responders+1])
  gene_expr_for_nonresp<-unlist(aggregated_expr[i,numbers_for_nonresponders+1])}

  
  # if both vectors is not constants and not na
  
  if(!is.constant(gene_expr_for_resp) & !is.constant(gene_expr_for_nonresp)& sum(is.na(gene_expr_for_nonresp))==0 & sum(is.na(gene_expr_for_resp))==0 ){
    t1<-
      t.test(gene_expr_for_resp,gene_expr_for_nonresp)
    pv[i] <-t1$p.value
    t2<-
      wilcox.test(gene_expr_for_resp,gene_expr_for_nonresp)
    pv_w[i] <-t2$p.value
  }
}

install.packages("forecast")
library(forecast)
for (i in 1:100){
  gene_expr_for_resp<-unlist(GSE[i,numbers_for_FL_BV_nonresponders])
  
  gene_expr_for_nonresp<-unlist(GSE[i,numbers_for_FL_non_responders])
  
  gene_expr_for_nonresp2<-unlist(GSE[i,numbers_for_FL_nonresponders_BV_responder])
  
  gene_expr_for_nonresp3<-unlist(GSE[i,numbers_for_FL_responders])
  
  gene_expr_for_nonresp4<-unlist(GSE[i,numbers_for_FL_responders_BV_non_responder])
  
  gene_expr_for_nonresp5<-unlist(GSE[i, numbers_for_FL_responders_BV_responder])
  
  # if both vectors is not constants and not na
  
  #if(!is.constant(gene_expr_for_resp) & !is.constant(gene_expr_for_nonresp)& sum(is.na(gene_expr_for_nonresp))==0 & sum(is.na(gene_expr_for_resp))==0 ){
   # t1<-
    #  t.test(gene_expr_for_resp,gene_expr_for_nonresp)
  #  pv[i] <-t1$p.value
    #t2<-
     # wilcox.test(gene_expr_for_resp,gene_expr_for_nonresp)
  #  pv_w[i] <-t2$p.value
  }
}
pv
pv_w
# поправка на значения 
pv_adjusted <- p.adjust(pv, method="fdr")
# поправка на значения 
pw_adjusted <- p.adjust(pv_w, method="fdr")
pv_adjusted
pw_adjusted

sum(na.rm = TRUE, pv_adjusted < 0.05)
sum(na.rm = TRUE, pw_adjusted < 0.05)
