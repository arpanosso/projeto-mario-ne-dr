
# função para plotar diferentes modelos -----------------------------------
plot_my_models <- function(modelo_1,modelo_2,modelo_3){
  sqr.f1<-round(attr(modelo_1, "SSErr"),4); c01<-round(modelo_1$psill[[1]],4); c0_c11<-round(sum(modelo_1$psill),4);a1<-round(modelo_1$range[[2]],2)
  sqr.f2<-round(attr(modelo_2, "SSErr"),4); c02<-round(modelo_2$psill[[1]],4); c0_c12<-round(sum(modelo_2$psill),4);a2<-round(3*modelo_2$range[[2]],2)
  sqr.f3<-round(attr(modelo_3, "SSErr"),4); c03<-round(modelo_3$psill[[1]],4); c0_c13<-round(sum(modelo_3$psill),4);a3<-round(modelo_3$range[[2]]*(3^.5),2)

  df_aux <- vari_exp |>
    mutate(
      gamma_m1 = ifelse(dist <= a1, c01 + (c0_c11-c01)*(3/2*(dist/a1)-1/2*(dist/a1)^3),c0_c11),
      gamma_m2 = c02 + (c0_c12-c02)*(1-exp(-3*(dist/a2))),
      gamma_m3 = c03 + (c0_c13-c03)*(1-exp(-(dist/a3)^2)),
      residuo_total = (gamma-mean(gamma))^2,
      residuo_mod_1 = (gamma - gamma_m1)^2,
      residuo_mod_2 = (gamma - gamma_m2)^2,
      residuo_mod_3 = (gamma - gamma_m3)^2
    ) |>
    summarise(
      r2_1=(sum(residuo_total) - sum(residuo_mod_1))/sum(residuo_total),
      r2_2=(sum(residuo_total) - sum(residuo_mod_2))/sum(residuo_total),
      r2_3=(sum(residuo_total) - sum(residuo_mod_3))/sum(residuo_total),
    )
  r21<-as.vector(round(df_aux[1],4))
  r22<-as.vector(round(df_aux[2],4))
  r23<-as.vector(round(df_aux[3],4))
  print(plot(vari_exp,
             model=modelo_1,
             col=1,pl=F,
             pch=16,
             cex=1.2,cex.main=7,
             ylab=list("Semivariância",cex=1.3),
             xlab=list("Distância de Separação h (m)",cex=1.3),
             main =paste("Esf(C0= ",c01,"; C0+C1= ",
                         c0_c11, "; a= ", a1,"; r2 = ",
                         r21,")",sep="")))
  modelo_3$range[[2]] <- a3
  print(plot(vari_exp,model=modelo_2, col=1,pl=F,pch=16,cex=1.2,cex.main=7,ylab=list("Semivariância",cex=1.3),xlab=list("Distância de Separação h (m)",cex=1.3),main =paste("Exp(C0= ",c02,"; C0+C1= ", c0_c12, "; a= ", a2,"; r2 = ", r22,")",sep="")))
  print(plot(vari_exp,model=modelo_3, col=1,pl=F,pch=16,cex=1.2,cex.main=7,ylab=list("Semivariância",cex=1.3),xlab=list("Distância de Separação h (m)",cex=1.3),main =paste("Gau(C0= ",c03,"; C0+C1= ", c0_c13, "; a= ", a3,"; r2 = ", r23,")",sep="")))
}
# função para validação cruzada -------------------------------------------
my_cross_validation <- function(df,formula,m1,m2,m3){
  conjunto_validacao <- df |>
    as_tibble() |>
    sample_n(300)
  coordinates(conjunto_validacao) = ~x + y
  modelos<-list(m1,m2,m3)
  for(j in 1:3){
    est<-0
    # vari<-as.character(form)[2]
    for(i in 1:nrow(conjunto_validacao)){
      valid <- krige(formula=formula, conjunto_validacao[-i,], conjunto_validacao, model=modelos[[j]])
      est[i]<-valid$var1.pred[i]
    }
    obs<-as.data.frame(conjunto_validacao)[,3]
    RMSE<-round((sum((obs-est)^2)/length(obs))^.5,3)
    mod<-lm(obs~est)
    b<-round(mod$coefficients[2],3)
    se<-round(summary(mod)$coefficients[4],3)
    r2<-round(summary(mod)$r.squared,3)
    a<-round(mod$coefficients[1],3)
    plot(est,obs,xlab="Estimado", ylab="Observado",pch=j,col="blue",
         main=paste("Modelo = ",modelos[[j]][2,1],"; Coef. Reg. = ", b, " (SE = ",se, ", r2 = ", r2,")\ny intersept = ",a,"RMSE = ",RMSE ))
    abline(lm(obs~est))
    abline(0,1,lty=3)
  }
}
