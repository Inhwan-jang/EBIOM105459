rm(list=ls())

#Parameters setup
work.dir = "~/source code/Data_Fig3c_d"  # route of directory

setwd(work.dir)
library(ggplot2)
library(survival)




#example of input file
cox.file = c("Liquid_NK_total.txt","Liquid_NK_top30.txt","Liquid_NK_bottom70.txt")  # top20-bottom80 | top30-bottom70  | top40-bottom60  | top50-bottom50 
type.list = c("Liquid")  
out.file = "Liquid_70_30"

cox.group = c("Total","NKrich","NKpoor") #NK, CD4T, or CD8T



  fore_rmeta.labels.list = list()
  
  for(n in 1:length(cox.file)){  
    
    cox.data= read.table(cox.file[n], header = T, sep ="\t",stringsAsFactors = F)
    
    
    sample.ID = cox.data[,1]
    cox.OS = cox.data[,c(2,3)]
    
    cox.type = cox.data[,4]
    cox.expr = cox.data[,-c(1,2,3,4)]
    cox.input.gene = cox.expr
    
    cox.input = cbind(cox.OS,cox.input.gene)
    
    cancer.ID = as.character(unique(cox.type))
    
    coefficients=list()
    fore_rmeta = list()
    
    
  
     

     

   result.cox = coxph(Surv(OS.month,OS) ~  NCR3LG1 + MICA + MICB + ULBP1 + ULBP2 + ULBP3 + BAG6 + 
                            CADM1 + CD70 + HSPG2 + ICAM1 + PCNA + 
                            PDGFD +  RAET1E + RAET1G + RAET1L  + TNFSF9 + VIM,data = cox.input)

       
   result.cox.summary = summary(result.cox)
   coefficients = cbind(result.cox$coefficients,result.cox.summary$coefficient[,5])
   fore_rmeta = cbind(result.cox.summary$conf.int[,1],result.cox.summary$conf.int[,3],result.cox.summary$conf.int[,4])

    
    # ggforest(result.cox)
    coefficients.output =as.data.frame(coefficients)
    coefficients.output = cbind(rownames(coefficients.output),coefficients.output)
    coefficients.output[,2] = as.numeric(coefficients.output[,2])
    coefficients.output[,3] = as.numeric(coefficients.output[,3])
    
    
    
    ######forest plot by cancer type #####
    

    fore_rmeta.frame = as.data.frame(fore_rmeta)
    colnames(fore_rmeta.frame) = c("HR","lower","upper")
    

    Pval = c(coefficients.output[,3])
    lower = round(fore_rmeta.frame[,2],2)
    upper = round(fore_rmeta.frame[,3],2)
    Type = c(rownames(fore_rmeta.frame))
    HR = c(round(fore_rmeta.frame[,1],2))
    group = rep(cox.group[n],length(HR))
    
    fore_rmeta.labels = as.data.frame(cbind(Type,HR,lower,upper,group,Pval))
    
    
    fore_rmeta.labels.list[[n]] = fore_rmeta.labels 
  }

  fore_rmeta.labels.fin = NULL
  for(j in 1:length(cox.file)){
    
    fore_rmeta.labels.fin = rbind(fore_rmeta.labels.fin,fore_rmeta.labels.list[[j]])
    
  }
  
write.table(fore_rmeta.labels.fin,file = paste(out.file,type.list,"multi_cox.txt"),row.names=F,sep = "\t", quote=F, col.names = T)

