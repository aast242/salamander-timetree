# Library
library(ggplot2)

# Dummy data

smatrix = read.table("supermatrix_long.tsv",
                     header = F, sep = "\t", col.names=c("gene", "sp", "pres",
                                                         "V4"))

# Heatmap 
hmap = ggplot(smatrix, aes(x=sp, y=gene, fill= pres)) + 
  geom_tile() +
  scale_fill_gradient(low=rgb(247/255, 251/255, 255/255),
                      high=rgb(8/255, 48/255, 107/255)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave("supermatrix.png", plot = hmap, device="png", limitsize = FALSE,
       width=24, height=16, dpi=600, bg="transparent")
