taxaplot = function(N,grouping_info,radist,ylabel="Proportions") {
  
  library(RColorBrewer)
  getPalette = colorRampPalette(brewer.pal(8, "Set3"))
    
  x=radist[,order(colSums(radist),decreasing=T)]
  taxa_list<-colnames(x)[1:N]
  #remove "No hits" and add it to others
  taxa_list<-taxa_list[!grepl("No hits",taxa_list)]
  N<-length(taxa_list)
  
  #Generate a new table with everything added to Others
  new_x<-data.frame(x[,colnames(x) %in% taxa_list],Others=rowSums(x[,!colnames(x) %in% taxa_list]))
  
  
  #You can change the Type=grouping_info[,1] should you desire any other grouping of panels
  df<-NULL
  
  for (i in 1:dim(new_x)[2]){
    tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),
                    Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),
                    Value=new_x[,i],Type=grouping_info[,1])
    
    if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
  }
  
  print(unique(df$Taxa))
  df$Taxa = factor(df$Taxa, levels = unique(df$Taxa))
  
  #colours = c(getPalette(N), "black")


  if(N>29){
    colours <- c(getPalette(N-29), rainbow(5),"#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                 "#0075DC", "#F0A3FF", "#993F00","#4C005C","#2BCE48","#FFCC99",
                 "#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380",
                 "#FFA8BB","#426600","#5EF1F2","#00998F",
                 "black","grey");
  }
  else if(N>24){
    colours <- c(getPalette(N-24), "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                 "#0075DC", "#F0A3FF", "#993F00","#4C005C","#2BCE48","#FFCC99",
                 "#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380",
                 "#FFA8BB","#426600","#5EF1F2","#00998F",
                 "black","grey");
                 
  }
  else if(N>16){
      colours <- c(getPalette(N-8), "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                   "#0075DC", "grey");
  }
  else {
    colours = c(getPalette(N),"grey")
  }
  
  
  
  library(ggplot2)
  library(grid)
  p<-ggplot(df, aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
  p<-p+scale_fill_manual(values=colours[1:(N+1)])
  p<-p+theme_bw()+ylab(ylabel)
  p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines"))
  p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))  
  print(p) 
}
