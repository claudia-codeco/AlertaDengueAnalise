# Gerador de Alerta - modelo claudia v1 
#========================================

#--------------------------
#Alerta por APS em 4 niveis
#--------------------------

#**Verde (atividade baixa)** 
#- se temperatura < 22 graus por 3 semanas 
#- se atividade de tweet for normal (nao aumentada)
#- ausencia de transmissao sustentada
#- se incidencia < 100:100.000

#**Amarelo (Alerta)**
#- se temperatura > 22C por mais de 3 semanas
#- se atividade de tweet aumentar

#**Laranja (Transmissao sustentada)**
#- se numero reprodutivo >1, por 3 semanas

#**Vermelho (atividade alta)**
#- se incidencia > 100:100.000

source("fun/Rt.r")
source("fun/data2SE.r")
source("fun/corrigecasos.r")


# Abre os dados mais recentes
#dadosAPS <- paste("../",dadosAPS,sep="")
#d<-read.csv(dadosAPS)


# Funcao getData pega a tabela mais atual 
# comecando na semana epidemiologica ini (default 201001)
# e corrige os dados de notificacao por causa do atraso (modelo lognormal): atualmente unica opcao 
getData <- function(dadosAPS,correction="lognormal",ini=201001){
  d<-read.csv(dadosAPS)
  dd<-subset(d,SE>=ini)
  listaAPS<-unique(dd$APS)
  if (correction=="lognormal") for(i in 1:10) dd$casosm[dd$APS==listaAPS[i]] <-corrigecasos(dd$casos[dd$APS==listaAPS[i]])
  tot<-aggregate(dd$casosm,by=list(dd$SE),sum)
  names(tot)<-c("SE","casosestimados")
  message("correcao da incidencia pelo atraso de notificacao")
  print(tail(tot))
  dd 
}
d<-getData(dadosAPS)


geraAlerta <- function(d,fun="original",pop="tabelas/populacao2010porAPS_RJ.csv",tempcrit=22){
  # Aplica o modelo original de alerta (desenvolvido pela claudia em 2014)
  
  listaAPS<-unique(d$APS)
  # data.frame para colocar os resultados
  SE<-d$SE[d$APS==listaAPS[1]]
  d2 <- expand.grid(SE=SE,APS=listaAPS)
  d2<-merge(d2,d[,c("SE","APS","data","tweets","estacao","casos","casosm","tmin")],by=c("SE","APS")) # precisa especificar colunas?

  # agregar dados de populacao (cuidado, desorganiza tudo!) 
  pop<-read.csv(file=pop)
  d2<-merge(d2,pop)
  d2<-d2[order(d2$APS,d2$SE),]

  #"Alerta de Temperatura minima semanal > 22 graus por 3 semanas"
  #Temperatura > Tcrit
  detcli <- function(temp,tempcrit=tempcrit,lag=3){
    t1<-as.numeric(temp>=tempcrit)
    le <- length(t1)
    ac <- t1[lag:le]
    for(i in 1:(lag-1)) ac <- ac+t1[(lag-i):(le-i)]
    c(rep(NA,(lag-1)),ac)
  }

  d2$alertaCli <- NA
  for(i in 1:10) d2$alertaCli[which(d2$APS==listaAPS[i])] <-detcli(d2$tmin[d2$APS==listaAPS[i]])

  # Alerta de Tweet usa Rt beta
  # Rtw = crescimento significativo de tweet na ultima semana (repete os valores para todas as APS)
  d2$Rtw <-NA  # Rt do tweet
  d2$ptw1 <-NA # Prob(Rt tweet >1)
  d2$Rtwlr <-NA # lim inf IC Rt tweet
  d2$Rtwur <-NA  # lim sup IC Rt tweet

  for(i in 1:10) d2[d2$APS==listaAPS[i],c("Rtw","ptw1","Rtwlr","Rtwur")] <- Rt.beta(d2$tweets[d2$APS==listaAPS[i]])

  # Rt > 1 se pr(Rt >1) > pcrit por lag weeks  # serve tanto para Rt de tweet como para casos
  Rtgreat1 <- function(p1,pcrit=0.8,lag=0){
    t1<-as.numeric(p1>pcrit)
    le <- length(t1)
    if (lag==0) ac<-t1[1:le]
    if (lag>0) {
          ini <- lag # ini <- (lag+1)
          ac <- t1[ini:le]
          for(i in 1:(ini-1)) ac <- ac+t1[(ini-i):(le-i)]
          ac<-c(rep(NA,(ini-1)),ac)
          }
    ac
  }
  # Calcula se Rt(tweet) > 1  
  for(i in 1:10) d2$twgreat1[d2$APS==listaAPS[i]] <-Rtgreat1(d2$ptw1[d2$APS==listaAPS[i]],pcrit=0.85,lag=0)
  # alertaRtweet = acumulado de Rt>1 por 3 vezes
  for(i in 1:10) d2$alertaRtweet[d2$APS==listaAPS[i]] <-Rtgreat1(d2$ptw1[d2$APS==listaAPS[i]],pcrit=0.85,lag=3)

  # calculo do Rt de dengue
  d2$Rt <-NA
  d2$pRt1 <-NA
  d2$Rtlr <-NA
  d2$Rtur <-NA

  for(i in 1:10) d2[d2$APS==listaAPS[i],c("Rt","pRt1","Rtlr","Rtur")] <- Rt.beta((d2$casosm[d2$APS==listaAPS[i]]))

  #  imputa os P(Rt) de casos faltantes com p(Rtw) dos tweets
  d2$pRti <- d2$pRt1
  d2$pRti[is.na(d2$pRt1)]<-d2$ptw1[is.na(d2$pRt1)]

  # Rt(dengue) > 1 se Pr>0.85 
  for(i in 1:10) d2$Rtgreat1[d2$APS==listaAPS[i]] <-Rtgreat1(d2$pRti[d2$APS==listaAPS[i]],pcrit=0.85,lag=0)
  # alertaRt = acumulado de Rt>1 por 3 vezes
  for(i in 1:10) d2$alertaRt[d2$APS==listaAPS[i]] <-Rtgreat1(d2$pRti[d2$APS==listaAPS[i]],pcrit=0.85,lag=3)

  # Casos:
  # onde houver dados faltantes, completar calculando o valor esperado considerando a tendencia de variacao das ultimas 3 semanas.
  # Isso é, calcula a média de variacao (Rt) das ultimas 3 semanas e multiplica pelo valor observado de casos. 

  fillCasos <- function(casos, Rt){
     casos_est <- casos
     n <- which(is.na(casos_est))
     if(length(n)>0) casos_est[n]<-casos[(n-2)]*mean(Rt[(n-2)],na.rm=TRUE)
     casos_est
  }
  d2$casos_est<-NA
  for(i in 1:10) d2$casos_est[d2$APS==listaAPS[i]]<-fillCasos((d2$casos[d2$APS==listaAPS[i]]),
                                                            d2$Rt[d2$APS==listaAPS[i]])
  d2$inc <- d2$casos_est/d2$Pop2010*100000
  d2$alertaCasos <- as.numeric(d2$inc>100)

  #Resultado
  # - alertaClima = 1, se temperatura > 22C por mais de 3 semanas
  # - alertaTweet = 1, se Tweet com tendencia de aumento 
  # - alertaTransmissao = 1, se casos com tendencia de aumento
  # - alertaCasos = 1, se Incidencia > 100 por 100 mil

  le = length(d2$tmin[d2$APS==listaAPS[1]])

  #funcao que define cor
  def.cor<-function(d2v){
    # d2v = dados de uma ap
    # 1 = verde, 2=amarelo, 3 =laranja, 4 = vermelho
    les = dim(d2v)[1]
    d2v$cor <-NA
    d2v$cor[intersect(6:les,which(d2v$alertaCli<3 & d2v$alertaRtweet<3 & d2v$alertaRt<3 & d2v$alertaCasos==0))]<-1
    d2v$cor[intersect(6:les,which(d2v$alertaCli>=3 | d2v$alertaRtweet>=3))]<-2
    d2v$cor[intersect(6:les,which(d2v$alertaCli>=3| d2v$alertaRtweet>=3)+1)]<-2 # inercia para desligar
    d2v$cor[intersect(6:les,which(d2v$alertaCli>=3| d2v$alertaRtweet>=3)+2)]<-2 # inercia para desligar
    d2v$cor[intersect(6:les,which(d2v$alertaCli>=3| d2v$alertaRtweet>=3)+3)]<-2 # inercia para desligar
    d2v$cor[intersect(6:les,which(d2v$alertaRt>=3))]<-3
    d2v$cor[intersect(6:les,which(d2v$alertaRt>=3)+1)]<-3 # inercia para desligar
    d2v$cor[intersect(6:les,which(d2v$alertaRt>=3)+2)]<-3 # inercia para desligar
    d2v$cor[intersect(6:les,which(d2v$alertaRt>=3)+3)]<-3 # inercia para desligar
    d2v$cor[intersect(6:les,which(d2v$alertaCasos==1))]<-4
    d2v$cor[intersect(6:les,which(d2v$alertaCasos==1)+1)]<-4  # inercia para desligar
    d2v$cor[intersect(6:les,which(d2v$alertaCasos==1)+2)]<-4   # inercia para desligar
    d2v$cor[intersect(6:les,which(d2v$alertaCasos==1)+3)]<-4
    d2v
  }

  d2$cor<-NA
  for(i in 1:10) d2[d2$APS==listaAPS[i],]<-def.cor(d2[d2$APS==listaAPS[i],])

  # codificacao para o relatorio
  d2$alertaTweet <- ifelse(d2$alertaRtweet>=3,1,0)
  d2$alertaClima <- ifelse(d2$alertaCli>=3,1,0)
  d2$alertaTransmissao <- ifelse(d2$alertaRt>=3,1,0)
  d2$nivel<-"nulo"
  d2$nivel[d2$cor==1]<-"verde"
  d2$nivel[d2$cor==2]<-"amarelo"
  d2$nivel[d2$cor==3]<-"laranja"
  d2$nivel[d2$cor==4]<-"vermelho"
  d2
}

# Aplica funcao geraAlerta aos dados
d2<-geraAlerta(d)

# Funcoes de log e salvamento de dados a alterar

options(digits = 2)
for(i in 1:10) {
  message(paste ("Estado do Alerta de ",listaAPS[i]))
  print(tail(d2[d2$APS==listaAPS[i],c(1,2,6,7,8,10,11,12,16,17,18,23,25,26,31)],n=4))
}
        
outputfile = paste("alerta/alertaAPS_",max(d2$SE),".csv",sep="")
message("alerta salvo em ",outputfile)

write.table(d2,file=outputfile,row.names=FALSE,sep=",")

nalerta <- outputfile
knit(input="geraAlerta/relatorio.rmd",quiet=TRUE,envir=new.env())

message("relatorio.md gerado. Abri-lo no RStudio e rodar PREVIEW PDF")
message("fim")