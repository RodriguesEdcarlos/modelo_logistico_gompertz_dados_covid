### Edcarlos Rodrigues Carneiro Silva
### Cidade: Guarulhos
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
dados_guarulhos = dados[which(dados$estado=='SP' 
                              & dados$municipio=='Guarulhos' 
                              & ymd(dados$data) >= ymd("2020-03-28") 
                              & ymd(dados$data) <= ymd("2022-09-23") 
                              & dados$obitosAcumulado != 0),]



dados_guarulhos$data <- ymd(dados_guarulhos$data)

plot(x = dados_guarulhos$data, y = dados_guarulhos$obitosAcumulado)

head(dados_guarulhos)
tail(dados_guarulhos)
dim(dados_guarulhos)
write.table(dados_guarulhos, file = 'E:/dados_garulhos.csv', sep = ';')
# implementação do modelo logistico:
y = dados_guarulhos$obitosAcumulado

T_ = length(y)
T_
t_ = seq(1:T_)

a.c = max(y)
a.c

b.c = -log(a.c/min(y)-1)
b.c

g.c = 0

sigma.c = sd(y)
sigma.c
theta.c = c(a.c,b.c,g.c,sigma.c)

# agora defino funções que serão de meu interesse

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

## defino chutes iniciais

a.c = max(y)*2.5
b.c = -log(a.c/min(y)-1)
b.c
g.c = 0
sigma.c = sd(y)
theta.c = c(a.c,b.c,g.c,sigma.c)

## iterlim = a quantidade de tentativas
lnL_max = nlm(func_ll,theta.c, steptol = 1e-12,
              gradtol = 1e-7,print.level = 1,iterlim = 1000000,
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


plot(y~t_,type="p",main="ÓBITOS Acumulados de COVID-19 na Cidade de Guarulhos - Logístico",
     ylab = 'Óbitos', xlab ='Dias' )
lines(y.h~t_,col = "blue")


library("numDeriv")
# calculo dos p-valores para teste de significância de 95%
hess=hessian(func_ll,lnL_max$estimate)
hess

var_covar_ll = solve(hess)
var_covar_ll

ep.ll = (diag(var_covar_ll)^0.5)[1:(length(diag(var_covar_ll))-1)]
ep.ll



theta.hat = lnL_max$estimate[1:3]
theta.hat
z = theta.hat/ep.ll
z

pvalor = 2*(1-pnorm(abs(z)))
pvalor


resultado.ll = cbind(theta.hat,ep.ll,z,pvalor)
colnames(resultado.ll)=c("Coef","Erro Padrão","Z","P-valor")
resultado.ll


# p-valor: hipotese nula é falsa

## beta 1
# a probabilidade de estar errado ao rejeitar a hipotese nula  
#de que o verdadeiro beta 1 é igual a zero, por máxima verossimilhança é de 98,88%

## fazer a mesma interpretação para o restante dos parâmetros

# <= 0.05 p-valor é estatisticamente significativo

## gráfico


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


plot(y~t_,type="p",main="ÓBITOS Acumulados de COVID-19 na Cidade de Guarulhos - Logístico",
     ylab = 'Dias', xlab ='óbitos' )
lines(y.h~t_,col = "blue")


#1.11
# logistico métricas de avaliação de erro
MSE = sum((y-y.h)^2)/length(t_)
MSE

RMSE = MSE^0.5
RMSE

MAE = sum(abs(y-y.h))/length(t_)
MAE

MAPE = sum(abs(y-y.h)/y)/length(t_)
MAPE


metricas.ll = cbind(MSE,RMSE,MAE,MAPE)  
metricas.ll

# Gerando previsões de óbitos acumulados para os próximos 1, 2 e 30 dias 
vetor = c(1:100)
previsoes= c(1:100)*0
for (i in vetor) {
  previsoes[i] = a.h/(1+exp(-b.h-g.h*(length(t_)+vetor[i])))
}
previsoes = c(previsoes[1], previsoes[2], previsoes[30])
previsoes=matrix(previsoes,nrow = 1, ncol = 3)
colnames(previsoes)=c('Previsão para 1 dia', 'Previsaão para 2 dias', 'Previsão para 30 dias')
previsoes


