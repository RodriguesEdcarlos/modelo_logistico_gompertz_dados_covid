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

## tratamento da base até chegar na cidade de Guarulhos
dados_guarulhos2 = dados[which(dados$estado=='SP' 
                               & dados$municipio=='Guarulhos' 
                               & ymd(dados$data) >= ymd("2020-03-28") 
                               & ymd(dados$data) <= ymd("2022-06-15") 
                               & dados$obitosAcumulado != 0),]


# implementação do modelo de gompertz:
y = dados_guarulhos2$obitosAcumulado
T_ = length(y)


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

a.c.g.5 = max(y)*1.5

b.c.g.5 = -log(log(a.c.g.5)-log(min(y)))

g.c.g.5 = 0

sigma.c.g.5 = (var(y)*(T_-1)/T_)^0.5


theta.c.g.5 = c(a.c.g.5,b.c.g.5,g.c.g.5,sigma.c.g.5)



lnL_max.g.5 = nlm(func_LL.G,theta.c.g.5, 
                   steptol = 1e-20,
                   gradtol = 1e-10,print.level = 1,iterlim = 1000000,
                   ndigit = 12, fscale = 1000)



#alpha hat
a.h.g.5 = lnL_max.g.5$estimate[1]
#beta hat
b.h.g.5 = lnL_max.g.5$estimate[2]

# gama hat
g.h.g.5 = lnL_max.g.5$estimate[3]

# sigma hat
sigma.h.g.5 = lnL_max.g.5$estimate[4]

y.h.g.5 = a.h.g.5*exp(-exp(-b.h.g.5-g.h.g.5*t_))

## gráfico
plot(y~t_,type="s",main="ÓBITOS Acumulados de COVID-19 na Cidade de Guarulhos - In Sample Gompertz",
     ylab = 'Óbitos', xlab ='Dias')
lines(y.h.g.5~t_,col = "red")





library(numDeriv)
hess.g=hessian(func_LL.G,lnL_max.g.5$estimate)


var_covar_ll.g = solve(hess.g)

ep.g = (diag(var_covar_ll.g)^0.5)[1:(length(diag(var_covar_ll.g))-1)]



theta.hat.g = lnL_max.g.5$estimate[1:3]

z.g = theta.hat.g/ep.g

pvalor.g = 2*(1-pnorm(abs(z.g)))


resultado.g = cbind(theta.hat.g,ep.g,z.g,pvalor.g)

colnames(resultado.g)=c('Coef', 'Erro Padrao', 'Razao Z', 'P-valor')
resultado.g
# modelo de gompertz
MSE.g = sum((y-y.h.g.5)^2)/length(t_)

RMSE.g = MSE.g^0.5


MAE.g = sum(abs(y-y.h.g.5))/length(t_)


MAPE.g = sum(abs(y-y.h.g.5)/y)/length(t_)


metricas.g = cbind(MSE.g,RMSE.g,MAE.g, MAPE.g)
metricas.g

vetor = c(1:100)
previsoes.g= c(1:100)*0
for (i in vetor) {
  previsoes.g[i] = a.h.g.5*exp(-exp(-b.h.g.5-g.h.g.5*(length(t_)+i)))
}
previsoes.g = c(previsoes.g[1], previsoes.g[2], previsoes.g[30])
previsoes.g=matrix(previsoes.g,nrow = 1, ncol = 3)
colnames(previsoes.g)=c('Previsão para 1 dia', 'Previsaão para 2 dias', 'Previsão para 30 dias')
previsoes.g


# gerando dados para teste do modelo
teste = dados[which(dados$estado=='SP' 
                    & dados$municipio=='Guarulhos' 
                    & ymd(dados$data) >= ymd("2022-06-16") 
                    & ymd(dados$data) <= ymd("2022-09-23") 
                    & dados$obitosAcumulado != 0),]

y_test.g <- teste$obitosAcumulado

t_test.g =length(t_)+seq(1:100)


# gerando previsões para os próximos 100 dias e comparar com os dados reais observados
y.h.g.5.teste = a.h.g.5*exp(-exp(-b.h.g.5-g.h.g.5*t_test.g))

par(mfrow = c(1, 1))

# gráfico modelo treinado
plot(y~t_,type="p",main="ÓBITOS Acumulados de COVID-19 na Cidade de Guarulhos - In Sample Gompertz",
     ylab = 'Óbitos', xlab ='Dias' )
lines(y.h.g.5~t_,col = "red")


# gráfico com os dados de teste
plot(y_test.g~t_test.g,type="p",main="ÓBITOS Acumulados de COVID-19 na Cidade de Guarulhos - Out of Sample Gompertz",
     ylab = 'Óbitos', xlab ='Dias', ylim= c(5000,6000))
lines(y.h.g.5.teste~t_test.g,col = "red")



# métricas de erro
MSE.g.teste = sum((y_test.g-y.h.g.5.teste)^2)/length(y_test.g)


RMSE.g.teste = MSE.g.teste^0.5


MAE.g.teste = sum(abs(y_test.g-y.h.g.5.teste))/length(y_test.g)


MAPE.g.teste = sum(abs(y_test.g-y.h.g.5.teste)/y_test.g)/length(y_test.g)

metricas.g.teste = cbind(MSE.g.teste,RMSE.g.teste,MAE.g.teste, MAPE.g.teste)
metricas.g.teste
