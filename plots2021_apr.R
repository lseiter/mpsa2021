library(igraph)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggraph)
library(graphlayouts)
library(tidyverse)

t<- as.data.table(read.csv("dataSets/tmp6.csv", stringsAsFactors = FALSE))

t<-t[lengthShock==3]
t<-t[sizeShock>=-0.1 & sizeShock <=0.1]

dim(t)
setnames(t, "winner", "candidate")  #oldname, newname
setnames(t, "sizeShock", "size_shock")  #oldname, newname
setnames(t, "lengthShock", "length_shock")  #oldname, newname

t$decay = as.factor(t$decay)
t$true_believer = as.factor(t$true_believer)
#t$seed = as.factor(t$seed)
t$seed <- factor(t$seed, levels=c( "0.18", "0.2", "0.22"))
t$candidate <- factor(t$candidate)
t$size_shock <- as.factor(t$size_shock)
t$length_shock <- as.factor(t$length_shock)

#nodecay<-t[decay==0.0] 
#nodecay$seed <- factor(nodecay$seed, levels=c( "0.18", "0.2", "0.22"))
#trueB<-t[true_believer==1]
#nodecaytrueB<-nodecay[true_believer==1]

winnerCount <- t %>%
  group_by(candidate, true_believer, seed, size_shock, length_shock) %>%
  summarize(N = n()) %>%
  mutate(freq = N / 50 *100 )%>%
  ungroup() %>%
  complete( candidate, true_believer, seed, size_shock, length_shock,
            fill = list(N = 0, freq = 0))


p<-ggplot(winnerCount,aes(seed, freq, fill=candidate))+
  geom_bar(stat="identity", width=0.5, position='dodge') + 
  geom_text(aes(label=freq),size=4,position=position_dodge(width=0.5),vjust=-0.25) +
  xlab("Candidate 1 Seed Proportion") +ylab("Win")+ 
  facet_grid(size_shock  ~ true_believer + length_shock, labeller=label_both)+
  theme( panel.spacing = unit(1, "lines"))+ 
  theme(text = element_text(size=12))+
  ggtitle("BA Graph 1000 Nodes, Candidate 2 and 3 Seed Proportion 0.15")
p  + ylim(c(0,120))






#first plot in paper, no decay, all true believers

winnerCount <- nodecaytrueB %>%
  group_by(candidate, seed) %>%
  summarize(N = n()) %>%
  mutate(freq = N / 100 )%>%
  ungroup() %>%
  complete( candidate, seed,
           fill = list(N = 0, freq = 0))

p<- ggplot(winnerCount,aes(seed, freq, fill=candidate))+
  geom_bar(stat="identity" ,width=0.5)  + 
  xlab("Candidate 1 Seed Proportion") +ylab("Win Proportion") +
  theme_bw() +
  theme(text = element_text(size=20))+
  ggtitle(t$configuration[1] 
  )
p + scale_fill_grey() 

#second plot, tb
winnerCount <- nodecay %>%
  group_by(candidate, true_believer, seed) %>%
  summarize(N = n()) %>%
  mutate(freq = N / 100 )%>%
  ungroup() %>%
  complete( candidate, true_believer, seed,
            fill = list(N = 0, freq = 0))

p<-ggplot(winnerCount,aes(seed,freq, fill=candidate))+
  geom_bar(stat="identity", width=0.5) + 
  xlab("Candidate 1 Seed Proportion") +ylab("Win Proportion")+
  facet_grid(. ~ true_believer ,labeller=label_both)+
  theme_bw()+
  theme( panel.spacing = unit(1, "lines")) +
  theme(text = element_text(size=20))+ 
  ggtitle(t$configuration[1])
p + scale_fill_grey() 


#third plot decay

winnerCount <- t %>%
  group_by(candidate, true_believer, seed, decay) %>%
  summarize(N = n()) %>%
  mutate(freq = N / 100 )%>%
  ungroup() %>%
  complete( candidate, true_believer, seed, decay,
            fill = list(N = 0, freq = 0))

p<-ggplot(winnerCount,aes(seed, freq, fill=candidate))+
  geom_bar(stat="identity", width=0.5) + 
  xlab("Candidate 1 Seed Proportion") +ylab("Win Proportion")+ 
  facet_grid(decay ~ true_believer, labeller=label_both)+
  theme_bw()+
  theme( panel.spacing = unit(1, "lines"))+ 
  theme(text = element_text(size=20))+
  ggtitle(t$configuration[1])
p   + scale_fill_grey()



#fourth compare candidate1 wins
candidate1Wins<-t[candidate=='1']
winnerCount<- candidate1Wins[,.N,by=.( true_believer, seed, decay, configuration)]
winnerCount$N = winnerCount$N/100

p<-ggplot(winnerCount,aes(x=seed, y=N, group=configuration))+
  geom_line(aes(color=configuration), size=1)+scale_x_discrete(expand = c(0.05, 0))+
 # scale_linetype_manual(values=c("dotted","longdash","twodash", "solid"))+
  scale_color_manual(values=c('#FF0000','#00FF00','#0000FF','#FFFF00','#FF00FF', '#000000'))+
  geom_point() +
  xlab("Candidate 1 Seed Proportion") + 
  facet_grid(decay ~ true_believer, labeller=label_both)+
  theme_bw()+
  theme(text = element_text(size=18), legend.position="bottom",  
        axis.title.y = element_blank(),
        panel.spacing = unit(1, "lines"))
p


#fifth no decay
candidate1Wins<-nodecay[candidate=='1']
winnerCount<- candidate1Wins[,.N,by=.( true_believer, seed,  configuration)]
winnerCount$N = winnerCount$N/100

p<-ggplot(winnerCount,aes(x=seed, y=N, group=configuration))+
  geom_line(aes(color=configuration), size=1)+geom_point()+
 # scale_linetype_manual(values=c("dotted","longdash","twodash", "solid"))+
  scale_color_manual(values=c('#FF0000','#00FF00','#0000FF','#FFFF00','#FF00FF', '#000000'))+
  xlab("Candidate 1 Seed Proportion") +ylab("Win Proportion")+ 
  facet_grid(. ~ true_believer, labeller=label_both)+theme_bw()+
  theme(text = element_text(size=16), legend.position="bottom")+
  ggtitle("Proportion of Candidate 1 Wins")
p




g<- sample_pa( n=1000, directed=FALSE)
g<- sample_pa( n=1000, m=4,directed=FALSE)
mean(degree(g))

ggraph(g,layout="stress")+
  geom_edge_link(width=0.2,colour="grey")+
  geom_node_point(col="black",size=0.3)+
  theme_graph()


#compare neighbors

winnerCount <- nodecaytrueB %>%
  group_by(candidate, seed, configuration) %>%
  summarize(N = n()) %>%
  mutate(freq = N / 100 )%>%
  ungroup() %>%
  complete( candidate, seed, configuration,
            fill = list(N = 0, freq = 0))

p<- ggplot(winnerCount,aes(seed, freq, fill=candidate))+
  geom_bar(stat="identity" ,width=0.5)  + 
  xlab("Candidate 1 Seed Proportion") +ylab("Win Proportion") +
  facet_grid(. ~ configuration, labeller=label_both) +
  theme_bw() +
  theme(text = element_text(size=16))+
  ggtitle(t$configuration[1] 
  )
p + scale_fill_grey() 


