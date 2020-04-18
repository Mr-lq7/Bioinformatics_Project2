setwd("C:/Users/linqingquan/Desktop/GSE_121787");
Sys.setlocale('LC_ALL','Chinese');

GPL_table = read.table('GPL21185-21174.txt',sep="\t",comment.char="#",
                       stringsAsFactors=F,header=T,fill=TRUE,quote="");

GSE_matrix <- read.table('GSE121787_series_matrix.txt',sep="\t",comment.char="!",
                         stringsAsFactors=F,header=T,fill=TRUE);
ID_Sybmol = GPL_table[,c(1,6)];
colnames(ID_Sybmol)[2]="Symbol";

Exp=merge(ID_Sybmol,GSE_matrix,by.x="ID",by.y="ID_REF",all=T);
Exp=Exp[,-1];
View(Exp)

Exp=Exp[Exp$Symbol != "",];
Exp=na.omit(Exp);

Exp1=data.frame(Exp[-grep("/",Exp$"Symbol"),]);
meanfun <- function(x) {
  x1 <- data.frame(unique(x[,1]));
  colnames(x1) <- c("Symbol");
  for (i in 2:ncol(x)) {
    x2 <- data.frame(tapply(x[,i],x[,1],mean));
    x2[,2] <- rownames(x2);
    colnames(x2) <- c(colnames(x)[i],"Symbol");
    x1 <- merge(x1,x2,by.x="Symbol",by.y="Symbol",all=FALSE);
  }
  return(x1);
}


Exp2 <- meanfun(Exp1);

par(cex=0.7);
n.sample=ncol(Exp2[,-1]);
if(n.sample>40)par(cex=0.5);
cols <- rainbow(n.sample*1.2);
boxplot(Exp2[,-1],col=cols,main="expression value",las=2)

write.table(Exp2,"Exp_Original.txt",row.names=F,quote=F,sep="\t")

row.names(Exp2)=Exp2[,1];
Exp2=log(Exp2[,-1]);

par(cex=0.7);
n.sample=ncol(Exp2);
if(n.sample>40)par(cex=0.5);
cols <- rainbow(n.sample*1.2);
boxplot(Exp2,col=cols,main="expression value",las=2);

Symbol=row.names(Exp2);
Exp_test=cbind(Symbol,Exp2);
write.table(Exp_test,"Exp.txt",row.names=F,quote=F,sep="\t")


top_diff <- read.table('top_diff.txt',sep="\t",comment.char="!", stringsAsFactors=F,header=T,fill=TRUE);
f2=merge(ID_Sybmol,top_diff,by.x="ID",by.y="ID",all=T)
f3=merge(f2,Exp_test,by.x = "Symbol",by.y = "Symbol",all=F);
f4=na.omit(f3)
f5 <- data.frame(unique(f4[,1]))
colnames(f5) <- c("ss")
f6=merge(f4,f5,by.x="Symbol",by.y = "ss",all=F)
f6 <- f4[!duplicated(f4$Symbol),]
write.table(f6,"merge_diff_gene.txt",row.names=F,quote=F,sep="\t")
            

