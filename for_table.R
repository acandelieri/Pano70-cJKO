rm(list=ls()); graphics.off(); cat("\014")
library(transport)

df <- NULL

files <- list.files("JKO_results",full.names=T)

for( f in files ) {
  
  res <- readRDS(f)
  
  ot <- transport(pp(res$X0),pp(res$Xks[[1]]),p=2,"primaldual")
  T.cost <- wasserstein(pp(res$X0),pp(res$Xks[[1]]),p=2,tplan=ot)
  for( k in 2:length(res$Xks) ) {
    ot <- transport(pp(res$Xks[[k-1]]),pp(res$Xks[[k]]),p=2,"primaldual")
    T.cost <- T.cost + wasserstein(pp(res$Xks[[k-1]]),pp(res$Xks[[k]]),p=2,tplan=ot)
  }
    
  ot <- transport(pp(res$X0),pp(res$Xks[[length(res$Xks)]]),p=2,"primaldual") 
  FW2 <- wasserstein(pp(res$X0),pp(res$Xks[[length(res$Xks)]]),p=2,tplan=ot)
  T.quality <- round((T.cost - FW2)/FW2,6)
  
  tmp <- unlist(strsplit(unlist(strsplit(f,"/"))[2],"_"))
  df <- rbind( df, data.frame(case=tmp[2],
                              algo=tmp[1],
                              obj.type=tmp[3],
                              h_or_eps=res$h,
                              F.final=round(FW2,6),
                              iterations=length(res$Xks),
                              runtime=sum(res$runTimes),
                              T.quality=round(T.quality,6),
                              stringsAsFactors=F) )
}


files <- list.files("cJKO_results",full.names=T)
for( f in files ) { 
  
  res <- readRDS(f)
  
  ot <- transport(pp(res$X0),pp(res$Xks[[1]]),p=2,"primaldual")
  T.cost <- wasserstein(pp(res$X0),pp(res$Xks[[1]]),p=2,tplan=ot)
  for( k in 2:length(res$Xks) ) {
    ot <- transport(pp(res$Xks[[k-1]]),pp(res$Xks[[k]]),p=2,"primaldual")
    T.cost <- T.cost + wasserstein(pp(res$Xks[[k-1]]),pp(res$Xks[[k]]),p=2,tplan=ot)
  }
  
  ot <- transport(pp(res$X0),pp(res$Xks[[length(res$Xks)]]),p=2,"primaldual") 
  T.quality <- T.cost - wasserstein(pp(res$X0),pp(res$Xks[[length(res$Xks)]]),p=2,tplan=ot)
  
  
  tmp <- unlist(strsplit(unlist(strsplit(f,"/"))[2],"_"))
  df <- rbind( df, data.frame(case=tmp[2],
                              algo=tmp[1],
                              obj.type=tmp[3],
                              h_or_eps=res$eps,
                              F.final=round(res$objs[length(res$objs)],6),
                              iterations=length(res$Xks),
                              runtime=sum(res$runTimes),
                              T.quality=round(T.quality,6),
                              stringsAsFactors=F) )
}

print(df[order(df$case),])
