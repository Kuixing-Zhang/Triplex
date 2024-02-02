###信号倒着画

##该用ggplot2
# 加载包
library(tidyverse)
library(stringr)
library(reshape2)

setwd("D:\\南大生科院\\Dr.chen\\triplex\\photoBQQ\\补实验\\S1_END_seq")

##读取数据
data_path_pattern2 <- "DDX5_SMARCA5_PSIP1_S1ENDSeq_shared_CTAG_center_2K"
data_path_pattern3 <- "DNMT1_TOP2A_MATR3_S1ENDSeq_shared_CTAG_center_2K"
data_path_pattern1 <- "HNRNPK_PCBP1_PCBP2_S1ENDSeq_shared_CTAG_center_2k"
data_path_KM12 <- "KM12_FiveCellineHDNAmirror_shared_CT_center"

CT_pattern1 <- read_delim(data_path_pattern1,
                        delim = "\t",
                        skip = 1,
                        col_names = F) %>%
  na.omit()
CT_pattern2 <- read_delim(data_path_pattern2,
                        delim = "\t",
                        skip = 1,
                        col_names = F) %>%
  na.omit()
CT_pattern3 <- read_delim(data_path_pattern3,
                          delim = "\t",
                          skip = 1,
                          col_names = F) %>%
  na.omit()

CT_KM12 <- read_delim(data_path_KM12,
                      delim = "\t",
                      skip = 1,
                      col_names = F) %>%
  na.omit()


# 保留区间id和信号值列
## protein
data_agg_merge_CT_pattern1 <- CT_pattern1[,-c(1:3,5,6)] #去除第1,2,3,5,6列
data_agg_merge_CT_pattern1_mean <- colMeans(data_agg_merge_CT_pattern1[,-1]) %>%
  as.data.frame()
data_agg_merge_CT_pattern1_mean$bin_number <- c(1:400)
data_agg_merge_CT_pattern1_mean$type[1:400] <- "HNRNPK"
data_agg_merge_CT_pattern1_mean$type[401:800] <- "PCBP1"
data_agg_merge_CT_pattern1_mean$type[801:1200] <- "PCBP2"
#data_agg_merge_CT_strand_mean$type[1201:1600] <- "KM12_fwd"
#data_agg_merge_CT_strand_mean$type[1601:2000] <- "KM12_rev"

## KM12 bw
data_agg_merge_CT_KM12 <- CT_KM12[,-c(1:3,5,6)] #去除第1,2,3,5,6列
data_agg_merge_CT_KM12_mean <- colMeans(data_agg_merge_CT_KM12[,-1]) %>%
  as.data.frame()
data_agg_merge_CT_KM12_mean$bin_number <- c(1:400)
data_agg_merge_CT_KM12_mean$type[1:400] <- "Forward"
data_agg_merge_CT_KM12_mean$type[401:800] <- "Reverse"


# 添加列名
colnames(data_agg_merge_CT_pattern1_mean) <- c("signal", "bin_number","type")
colnames(data_agg_merge_CT_KM12_mean) <- c("signal", "bin_number","type")
# plot
data_plot_pattern1 <- data_agg_merge_CT_pattern1_mean
data_plot_pattern1$type <- factor(data_plot_pattern1$type, levels  = unique(data_plot_pattern1$type))

data_plot_KM12 <- data_agg_merge_CT_KM12_mean
data_plot_KM12$type <- factor(data_plot_KM12$type, levels  = unique(data_plot_KM12$type))


# 画图参数
y_label <- "signal (CPM)"
title_label <- ""
y_limit <- c(round(min(data_plot_KM12$signal)),
             round(max(data_plot_KM12$signal)) + 1)
y_axis <- seq(round(min(data_plot_KM12$signal)),
              round(max(data_plot_KM12$signal))+1,
              2)
x_label <- "Repeat Center (kb)"
legend_label <- ""
fill_value <- c("#AAB1FF","#004CFF","#29FFCE","#000000","#808080")
size <- 4
y <- data_plot_KM12$signal
x <- data_plot_KM12$bin_number

x_limit <- c(0,max(data_plot_KM12$bin_number))
x_axis <- seq(0,max(data_plot_KM12$bin_number),
              max(data_plot_KM12$bin_number)/4)
x_axis_pattern1 <- seq(0,max(data_plot_pattern1$bin_number),
                       max(data_plot_pattern1$bin_number)/4)

legend_position <- c(0.8,0.7)
legend_direction <- "vertical"

top.mar=0.6
right.mar=0.6
bottom.mar=0.6
left.mar=0.6

mytheme<-theme_classic()+
  theme(text=element_text(family = "sans",colour ="black",size = 20),  # 文本字体为sans，颜色为灰色30%，大小为12
        legend.text=element_text(colour ="black",size = 15),  # 图例文本颜色为灰色30%，大小为10
        legend.title=element_text(colour ="black",size = 15),  # 图例标题颜色为灰色30%，大小为10
        legend.key.size=unit(6,units = "mm"),  # 图例键的大小为6毫米
        legend.position="top", # 图例位置在上侧
        axis.line = element_line(size = 0.8,colour = "gray30"),  # 坐标轴线的大小为0.4，颜色为灰色30%
        axis.ticks = element_line(size = 0.8,colour = "gray30"),  # 坐标轴刻度线的大小为0.4，颜色为灰色30%
        axis.ticks.length = unit(1.8,units = "mm"),  # 坐标轴刻度线长度为1.5毫米
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),  # 绘图区域的边距为top.mar、right.mar、bottom.mar、left.mar
                         units="inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 

# 画图
# 计算第一个y轴的范围
library(ggplot2)
library(lubridate)
p1 <- ggplot(data = data_plot_pattern1, aes(x = data_plot_pattern1$bin_number, colors=type)) +
  geom_line(aes(y = ifelse(type %in% c("HNRNPK", "PCBP1", "PCBP2"), signal, NA))) +
  scale_x_continuous(breaks = x_axis_pattern1,
                     labels = c("-2","-1","Center","1","2")) +
  #scale_color_manual(values= c("#67c2a3" ,"#29abe2", "#e889bd")) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA, color = 'black'), 
        axis.text.y = element_text(color = 'blue'), axis.ticks.y = element_line(color = 'blue'), 
        axis.title.y = element_text(color = 'blue')) +
  labs(y = 'MS protein signal (CPM)')
p1
p2 <- ggplot(data = data_plot_KM12, aes(x = x,colors=type)) +
  geom_line(aes(y = ifelse(type %in% c("Forward", "Reverse"), signal, NA))) +
  scale_x_continuous(breaks = x_axis,
                     labels = c("-2","-1","Center","1","2")) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA, color = 'black'), 
        axis.text.y = element_text(color = 'blue'), axis.ticks.y = element_line(color = 'blue'), 
        axis.title.y = element_text(color = 'blue')) +
  labs(y = 'HDNA signal (CPM)')
p2
y2_plot(p1, p2)

ggsave("Pattern1_shared_S1ENDSeq.pdf")
ggsave("Pattern2_shared_S1ENDSeq.pdf")
ggsave("Pattern3_shared_S1ENDSeq.pdf")

library(gtable)
library(grid)

# 定义函数，用于组合 ggplot2 绘图结果构建双坐标轴，参考自：
# https://stackoverflow.com/questions/36754891/ggplot2-adding-secondary-y-axis-on-top-of-a-plot

y2_plot <- function(p1, p2) {
  p1 <- ggplotGrob(p1)
  p2 <- ggplotGrob(p2)
  
  # Get the location of the plot panel in p1.
  # These are used later when transformed elements of p2 are put back into p1
  pp <- c(subset(p1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  p1 <- gtable_add_grob(p1, p2$grobs[[which(p2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  # Get the y axis title from p2
  index <- which(p2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- p2$grobs[[index]]                # Extract that grob
  ylab <- hinvert_title_grob(ylab)         # Swap margins and fix justifications
  
  # Put the transformed label on the right side of p1
  p1 <- gtable_add_cols(p1, p2$widths[p2$layout[index, ]$l], pp$r)
  p1 <- gtable_add_grob(p1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from p2 (axis line, tick marks, and tick mark labels)
  index <- which(p2$layout$name == 'axis-l')  # Which grob
  yaxis <- p2$grobs[[index]]                  # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of p1
  p1 <- gtable_add_cols(p1, p2$widths[p2$layout[index, ]$l], pp$r)
  p1 <- gtable_add_grob(p1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(p1)
}

