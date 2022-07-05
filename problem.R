###############################################
###############################################
#Reading files

fil<- read.csv(file = "77_cancer_proteomes_CPTAC_itraq.csv")
fil_2<-read.csv(file = "clinical_data_breast_cancer.csv")
fil_3<-read.csv(file = "PAM50_proteins.csv")
######part 1
# we select column REFSEQ freom file PAM50 to make condition on it with the coloumns
#in the transposed clinical data in correlation and t_test instead of making join

vec_ref_seq_pamfile<-c()
vec_ref_seq_pamfile<-fil_3$RefSeqProteinID 
#1.1
#extracting column complete tcga from the clinical file to make some changes on it
#we start from column 4 because the 3 coloumns doesn't have the genes names
#we make change on genes names on cancer file to be suitable with genes in clinical data

listPatient2<-fil_2$Complete.TCGA.ID

column_names<-as.list(colnames(fil))
column_names<-column_names[4:83]
column_names
#subset of col names in file cancer to compare  it by rows name in col Complete TCGA
#ID in clinical file 
listPatient<-c()
for ( i in column_names){
  part1<-substr(i, nchar(i)-3, nchar(i))
  part3<-substr(i,4, nchar(i)-7)
  part2<-substr(i,1,2)
  listPatient<-c(listPatient,paste0(part1,"-",part2,"-",part3))
}
#replace patient columns with correct format in file proteome
for(i in 1:length(listPatient)){
names(fil)[i+3]<-listPatient[i]
print(names(fil)[i+3])
}
#1.2
fil=na.omit(fil)
tansposedProteome<-as.data.frame(t(fil[4:83]))
list_refseq=fil$RefSeq_accession_number
list_refseq

rownames_trans=rownames(tansposedProteome)
rownames(tansposedProteome)=NULL  #remove rownames as its numbers
tansposedProteome=cbind(rownames_trans,tansposedProteome)
tansposedProteome
# making colum names of transposedproteome with genes names and
#the first coloum with complete.TCGA.ID
proteome_colums=c("Complete.TCGA.ID",list_refseq)
colnames(tansposedProteome)=NULL
colnames(tansposedProteome)=proteome_colums

library(dplyr)
#select colum complete and herefinal from file cancer that are interested in it 
clinicalcolumns<-as.data.frame(cbind(fil_2$Complete.TCGA.ID,fil_2$HER2.Final.Status))
colnames(clinicalcolumns)=NULL
colnames(clinicalcolumns)=c("Complete.TCGA.ID","HER2.Final.Status")

joined_result<-inner_join(clinicalcolumns,tansposedProteome,by="Complete.TCGA.ID")
#########################################################################
#####part2
# convert to numeric
# convert every value in colum her to not having error when make t_test and corr
#2.1
vec_class<-c()
for (i in joined_result$HER2.Final.Status){
  if(i=='Negative'){
    vec_class=c(vec_class,-1)
  }else if(i=='Positive'){
    vec_class=c(vec_class,1)
  }else{
    vec_class=c(vec_class,0)
  }
}
#bn3'yrha fe ele joined_result nfsha
joined_result$HER2.Final.Status=vec_class
vec_correlation<-c()
col.names.joined.result=colnames(joined_result)#extract colnames from joinedresult

#print(joined_result[768]  %in% col.names.joined.result ) 
#print(col.names.joined.result %in% vec_ref_seq_pamfile)

# perform cor with every colums in joinedresult and her coloum but
#before it cheak if the coloum in the coloum refseq in pamfile instead of
#making JOIN with this file 
#length joined 6914
for (j in 3:length(joined_result)){
  if( col.names.joined.result[j] %in% vec_ref_seq_pamfile){
    result_cor<-cor(x=joined_result[j],y=joined_result$HER2.Final.Status,method ='pearson')
    vec_correlation<-c(vec_correlation,result_cor)
  }
}
#2.2
# order them and function order return indexes not values 
result_oreder<-order(x=vec_correlation,decreasing = TRUE)
result_oreder
#2.3
threshold=0.2
filtered_vec_correlation<-c()
filtered_vec_columns<-c()
# comparing it with threshold to filliter them and save the coloums names from proteme coloums the coloums of transpose cancer
for(i in result_oreder){
  if(abs(vec_correlation[i])>threshold){
    filtered_vec_correlation=c(filtered_vec_correlation,vec_correlation[i])
  filtered_vec_columns=c(filtered_vec_columns,proteome_colums[i])  }
}
#################################################
####part3

##3.1
#spilit that according herfinal coloum that contain + = -
res_split<-split(joined_result,joined_result$HER2.Final.Status)
neg_split<-res_split$`-1`
pos_split<-res_split$`1`
#combing the - + to perform t_test
#res_combin<-rbind(neg_split,pos_split)
vec_p_value<-c()
#t_test
##3.2
#perform t_test with every colums in rescombin but before
#it cheak if the coloum in the coloum refseq in pamfile instead of making JOIN with this file 

for(i in 3:length(pos_split)){
  if( col.names.joined.result[i] %in% vec_ref_seq_pamfile){
    #x<-t.test(res_combin[,i]~res_combin$HER2.Final.Status)
    x<-t.test(x=abs(neg_split),y = pos_split)
    vec_p_value<-c(vec_p_value,x$p.value)
  }
}
vec_p_value
##3.3
res_ordered<-order(x=vec_p_value,decreasing = TRUE)
filtered_t_test_columns<-c()
for(i in res_ordered){
  if(vec_p_value[i]<0.05){
    filtered_t_test_columns=c(filtered_t_test_columns,proteome_colums[i])
    
  }
}
##4.4
## all features in t-test abears in correlation features and there numbers are 8
for(i in filtered_vec_columns){
  for(j in filtered_t_test_columns ){
    if(i==j){
      print(i)
    }
  }
}
