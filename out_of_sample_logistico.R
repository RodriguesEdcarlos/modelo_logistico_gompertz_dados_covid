#out of sample

options(scipen = 999999)
#MUDANDO O DIRETÓRIO PARA PASTA ONDE SE ENCONTRAM OS DADOS DA COVID NO BRASIL
setwd('insira o diretório dos arquivos aqui')
install.packages("lubridate")

library(lubridate)
# importando a base de dados:
dados1 = read.csv2('HIST_PAINEL_COVIDBR_2020_Parte1_26set2022.csv',sep = ';', 
                   encoding = 'UTF-8',header = TRUE)

dados2 = read.csv2('HIST_PAINEL_COVIDBR_2020_Parte2_26set2022.csv',sep = ';', 
                   encoding = 'UTF-8',header = TRUE)

dados3 = read.csv2('HIST_PAINEL_COVIDBR_2021_Parte1_26set2022.csv',sep = ';', 
                   encoding = 'UTF-8',header = TRUE)

dados4 = read.csv2('HIST_PAINEL_COVIDBR_2021_Parte2_26set2022.csv',sep = ';', 
                   encoding = 'UTF-8',header = TRUE)

dados5 = read.csv2('HIST_PAINEL_COVIDBR_2022_Parte1_26set2022.csv',sep = ';', 
                   encoding = 'UTF-8',header = TRUE)

dados6 = read.csv2('HIST_PAINEL_COVIDBR_2022_Parte2_26set2022.csv',sep = ';', 
                   encoding = 'UTF-8',header = TRUE)

# empilhar dados

dados = rbind(dados1,dados2,dados3,dados4,dados5,dados6)

head(dados)
tail(dados)
dim(dados)
#dados$data = ymd(dados$data)
## tratamento da base até chegar na cidade de Guarulhos
dados_guarulhos2 = dados[which(dados$estado=='SP' 
                              & dados$municipio=='Guarulhos' 
                              & ymd(dados$data) >= ymd("2020-03-28") 
                              & ymd(dados$data) <= ymd("2022-06-15") 
                              & dados$obitosAcumulado != 0),]

y = dados_guarulhos2$obitosAcumulado

T_ = length(y)
T_
t_ = seq(1:T_)

a.c = max(y)*1.48
a.c

b.c = -log(a.c/min(y)-1)
b.c

g.c = 0

sigma.c = sd(y)
sigma.c
theta.c = c(a.c,b.c,g.c,sigma.c)

func_ll = function(theta){
  
  ##defino os coeficientes
  a <- theta[1] # alpha
  b <- theta[2] # beta
  g <- theta[3] # gama
  sigma <- theta[4]
  
  ## cira-se o vetor de residuos
  
  u = y - a/(1+exp(-b-g*t_))
  
  
  ## escrevo minha função de logartimo da verossimilhança (ln_L),
  ## com uma distribuição normal
  ln_L = -(-0.5*T_*log(2*pi)-T_*log(abs(sigma))-0.5*sum(u^2)/(sigma^2))
  return(ln_L)
  
  ## (-) na frente para maximizá-la 
}

lnL_max = nlm(func_ll,theta.c, steptol = 1e-15,
              gradtol = 1e-10,print.level = 1,iterlim = 1000000,
              ndigit = 12, fscale = 1000)

#alpha hat
a.h = lnL_max$estimate[1]
a.h
#beta hat
b.h = lnL_max$estimate[2]
b.h
# gama hat
g.h = lnL_max$estimate[3]
g.h
# sigma hat
sigma.h = lnL_max$estimate[4]

y.h = a.h/(1+exp(-b.h-g.h*t_))

## gráfico


plot(y~t_,type="p",main="ÓBITOS Acumulados de COVID-19 na Cidade de Guarulhos - In Sample",
     ylab = 'Óbitos', xlab ='Dias' )
lines(y.h~t_,col = "blue")

library("numDeriv")
# p-valor para teste de significância de 95%
hess=hessian(func_ll,lnL_max$estimate)

var_covar_ll = solve(hess)

ep.ll = (diag(var_covar_ll)^0.5)[1:(length(diag(var_covar_ll))-1)]

theta.hat = lnL_max$estimate[1:3]
theta.hat
z = theta.hat/ep.ll

pvalor = 2*(1-pnorm(abs(z)))



resultado.ll = cbind(theta.hat,ep.ll,z,pvalor)
colnames(resultado.ll)=c("Coef","Erro Padrão","Z","P-valor")
resultado.ll

#métricas de erro
MSE = sum((y-y.h)^2)/length(t_)

RMSE = MSE^0.5

MAE = sum(abs(y-y.h))/length(t_)

MAPE = sum(abs(y-y.h)/y)/length(t_)


metricas.ll = cbind(MSE,RMSE,MAE,MAPE)  
metricas.ll

## Treinando o modelo logístico
teste = dados[which(dados$estado=='SP' 
                               & dados$municipio=='Guarulhos' 
                               & ymd(dados$data) >= ymd("2022-06-16") 
                               & ymd(dados$data) <= ymd("2022-09-23") 
                               & dados$obitosAcumulado != 0),]

y_test <- teste$obitosAcumulado
T_test= length(y_test)

t_test =length(t_)+seq(1:100)
t_test


y.h.teste = a.h/(1+exp(-b.h-g.h*t_test))
y.h.teste

par(mfrow = c(1, 1))
plot(y~t_,type="p",main="ÓBITOS Acumulados de COVID-19 na Cidade de Guarulhos - In Sample Logístico",
     ylab = 'Óbitos', xlab ='Dias' )
lines(y.h~t_,col = "blue")

plot(y_test~t_test,type="p",main="ÓBITOS Acumulados de COVID-19 na Cidade de Guarulhos - Out of Sample Modelo Logístico",
     ylab = 'Óbitos', xlab ='Dias' )
lines(y.h.teste~t_test,col = "blue")



MSE.teste = sum((y_test-y.h.teste)^2)/length(t_test)
MSE.teste

RMSE.teste = MSE.teste^0.5
RMSE.teste

MAE.teste = sum(abs(y_test-y.h.teste))/length(t_test)
MAE.teste

MAPE.teste = sum(abs(y_test-y.h.teste)/y_test)/length(t_)
MAPE.teste


metricas.teste = cbind(MSE.teste,RMSE.teste,MAE.teste, MAPE.teste)
metricas.teste
