#220228
##1.create tx2gene*.csv
in PC conda base R:
```
library("biomaRt")
setwd("H:/work/main/bioinfo")

mart1 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host = 'https://asia.ensembl.org')
mart2 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl",host = 'https://asia.ensembl.org')
df1 <- biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id","ensembl_transcript_id","refseq_mrna","refseq_ncrna"),mart = mart1)
df2 <- biomaRt::getBM(attributes = c("external_gene_name","ensembl_gene_id","ensembl_transcript_id","refseq_mrna","refseq_ncrna"),mart = mart2)
refseq1<-ifelse(df1$refseq_mrna!="",df1$refseq_mrna,df1$refseq_ncrna)
refseq2<-ifelse(df2$refseq_mrna!="",df2$refseq_mrna,df2$refseq_ncrna)

t2g<-data.frame(transcript=df1$ensembl_transcript_id,gene=df1$hgnc_symbol)
t2g<-t2g[t2g$gene!="",]
write.csv(t2g,file="./resource/tx2gene/tx2gene_ens2hgnc_hsapiens.csv",row.names=F,quote=F)
t2g<-data.frame(transcript=df2$ensembl_transcript_id,gene=df2$external_gene_name)
t2g<-t2g[t2g$gene!="",]
write.csv(t2g,file="./resource/tx2gene/tx2gene_ens2symbol_mmusculus.csv",row.names=F,quote=F)

t2g<-data.frame(transcript=refseq1,gene=df1$hgnc_symbol)
t2g<-t2g[t2g$gene!="" & t2g$transcript!="",]
write.csv(t2g,file="./resource/tx2gene/tx2gene_ref2hgnc_hsapiens.csv",row.names=F,quote=F)
t2g<-data.frame(transcript=refseq2,gene=df2$external_gene_name)
t2g<-t2g[t2g$gene!="" & t2g$transcript!="",]
write.csv(t2g,file="./resource/tx2gene/tx2gene_ref2symbol_mmusculus.csv",row.names=F,quote=F)

t2g<-data.frame(transcript=df1$ensembl_transcript_id,gene=df1$ensembl_gene_id)
write.csv(t2g,file="./resource/tx2gene/tx2gene_ens2ens_hsapiens.csv",row.names=F,quote=F)
t2g<-data.frame(transcript=df2$ensembl_transcript_id,gene=df2$ensembl_gene_id)
write.csv(t2g,file="./resource/tx2gene/tx2gene_ens2ens_mmusculus.csv",row.names=F,quote=F)

inte<-data.frame(tx=df1$ensembl_transcript_id,ens=df1$ensembl_gene_id,symbol=df1$hgnc_symbol)
gene<-ifelse(inte$symbol!="",inte$symbol,inte$ens)
t2g<-data.frame(transcript=df1$ensembl_transcript_id,gene=gene)
write.csv(t2g,file="./resource/tx2gene/tx2gene_ens2mix_hsapiens.csv",row.names=F,quote=F)
inte<-data.frame(tx=df2$ensembl_transcript_id,ens=df2$ensembl_gene_id,symbol=df2$external_gene_name)
gene<-ifelse(inte$symbol!="",inte$symbol,inte$ens)
t2g<-data.frame(transcript=df2$ensembl_transcript_id,gene=gene)
write.csv(t2g,file="./resource/tx2gene/tx2gene_ens2mix_mmusculus.csv",row.names=F,quote=F)
```