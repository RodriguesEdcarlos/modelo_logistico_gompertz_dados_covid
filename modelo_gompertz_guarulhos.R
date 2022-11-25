### Edcarlos Rodrigues Carneiro Silva
### Cidade: Guarulhos

# implementação do modelo de Gompertz:

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

y = dados_guarulhos$obitosAcumulado

tail(y)
T_ = length(y)

T_

t_ = seq(1:T_)


func_LL.G =  function(theta){
  
  ##defino os coeficientes
  a <- theta[1] # alpha
  b <- theta[2] # beta
  g <- theta[3] # gama
  sigma <- theta[4]
  
  ## cira-se o vetor de residuos
  
  u = y- a*exp(-exp(-b-g*t_))
  
  
  ## escrevo minha função de logartimo da verossimilhança (ln_L),
  ## com uma distribuição normal
  ln_L = -(-0.5*T_*log(2*pi)-T_*log(abs(sigma))-0.5*sum(u^2)/(sigma^2))
  return(ln_L)
  
  ## (-) na frente para maximizá-la  no nlm
  
}

a.c.g.61 = max(y)*1.61

b.c.g.61 = -log(log(a.c.g.61)-log(min(y)))

g.c.g.61 = 0

sigma.c.g.61 = (var(y)*(T_-1)/T_)^0.5


theta.c.g.61 = c(a.c.g.61,b.c.g.61,g.c.g.61,sigma.c.g.61)


lnL_max.g.61 = nlm(func_LL.G,theta.c.g.61, 
                steptol = 1e-20,
                gradtol = 1e-10,print.level = 1,iterlim = 1000000,
                ndigit = 12, fscale = 1000)


theta.c.g.61

-lnL_max.g.61$minimum

#alpha hat
a.h.g.61 = lnL_max.g.61$estimate[1]
#beta hat
b.h.g.61 = lnL_max.g.61$estimate[2]

# gama hat
g.h.g.61 = lnL_max.g.61$estimate[3]

# sigma hat
sigma.h.g.61 = lnL_max.g.61$estimate[4]

y.h.g.61 = a.h.g.61*exp(-exp(-b.h.g.61-g.h.g.61*t_))

## gráfico
plot(y~t_,type="s",main="ÓBITOS Acumulados de COVID-19 na Cidade de Guarulhos - Gompertz",
     ylab = 'Óbitos', xlab ='Dias')
lines(y.h.g.61~t_,col = "red")


# Gerando p-valores para teste de signifiância de 95%
library(numDeriv)
hess.g=hessian(func_LL.G,lnL_max.g.61$estimate)
hess.g

var_covar_ll.g = solve(hess.g)
var_covar_ll.g

ep.g = (diag(var_covar_ll.g)^0.5)[1:(length(diag(var_covar_ll.g))-1)]
ep.g



theta.hat.g = lnL_max.g.61$estimate[1:3]
theta.hat.g

z.g = theta.hat.g/ep.g
z.g

pvalor.g = 2*(1-pnorm(abs(z.g)))
pvalor.g


resultado.g = cbind(theta.hat.g,ep.g,z.g,pvalor.g)

colnames(resultado.g)=c('Coef', 'Erro Padrao', 'Razao Z', 'P-valor')
resultado.g


# métricas de avaliação de erro modelo de gompertz
MSE.g = sum((y-y.h.g.61)^2)/length(t_)

RMSE.g = MSE.g^0.5


MAE.g = sum(abs(y-y.h.g.61))/length(t_)


MAPE.g = sum(abs(y-y.h.g.61)/y)/length(t_)
MAPE.g

metricas.g = cbind(MSE.g,RMSE.g,MAE.g, MAPE.g)
metricas.g


# Gerando previsões de óbitos acumulados para os próximos 1, 2 e 30 dias 
vetor = c(1:100)
previsoes.g= c(1:100)*0
for (i in vetor) {
  previsoes.g[i] = a.h.g.61*exp(-exp(-b.h.g.61-g.h.g.61*(length(t_)+i)))
}
previsoes.g = c(previsoes.g[1], previsoes.g[2], previsoes.g[30])
previsoes.g=matrix(previsoes.g,nrow = 1, ncol = 3)
colnames(previsoes.g)=c('Previsão para 1 dia', 'Previsaão para 2 dias', 'Previsão para 30 dias')
previsoes.g
