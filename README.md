
<!-- README.md is generated from README.Rmd. Please edit that file -->

# APRENDIZADO DE MÁQUINA PARA PREDIÇÃO DA VARIABILIDADE ESPAÇO-TEMPORAL DA PRODUÇÃO DE CANA-DE-AÇÚCAR

## Carregando pacotes

``` r
library(tidymodels)
library(tidyverse)
library(modeldata)
library(ggridges)
library(corrplot)
library(brulee)
library(keras3)
library(ggpubr)
library(vegan)
library(gstat)
library(ISLR)
library(vip)
library(sp)
source("R/my-functions.R")
theme_set(theme_bw())
```

## Lendo o banco de dados

``` r
data_set <- read_rds("data/sugarcane-soil.rds")  |> 
  rename_with(~ str_replace(., "_2$", ""), ends_with("_2")) |>  
    rename(TCH = tch_real,
           Ca = ca,
           `m%` = m,
           SB = sb,
           CTC = ctc,
           P = p_resina,
           `H+Al` = h_al,
           `V%` = v,
           K = k,
           Mg = mg,
           MO = mo,
           pH = ph_cacl2,
           S = s) 
glimpse(data_set)
```

## 1. Análise Estatística Exploratória

Tabela com a estatística descritiva por unidade (CAT e POT) para as
variáveis numéricas (n, min, q1, mediana, média, q3, max, dp, epm, cv,
cassimetria, ccurtose) para as variáveis `tch_real` até `s_2`…nos
diferentes tempos, pra caracterizarmos a variabilidade total do conjunto
de dados

criando a função para a análise exploratória dos dados

``` r
estat_desc <- function(x){
  set.seed(1235)
  xs <- sample(x,500)
  normal <- nortest::lillie.test(xs)
  log_norm <- nortest::lillie.test(log(xs+1))
  n <- length(x)
  m <- mean(x,na.rm = TRUE)
  md <- median(x)
  mini <- min(x,na.rm = TRUE)
  q1 <- quantile(x,.25)
  q3 <- quantile(x,.75)
  maxi <- max(x,na.rm = TRUE)
  dp <- sd(x,na.rm = TRUE)
  epm <- dp/sqrt(n)
  cv <- 100*dp/m
  ass <- agricolae::skewness(x)
  curt <- agricolae::kurtosis(x)
  c(n,mini,q1,md,m,q3,maxi,dp,epm,cv,ass,curt,normal$p.value,log_norm$p.value)
}
estat_names <- c("n","Min","Q1","Med","Média","Q3",
                 "Max","DP","EP","CV","Skn","Krt","norm","log_norm")
```

Tabela de estatística descritiva presente no arquivo do excel na pasta
`output`.

``` r
# data_set |> 
#   group_by(ano,unidade) |> 
#   reframe(
#     across(.cols = TCH:S,
#            .fns = estat_desc,
#            .names = "{.col}")
#   )|> 
#   add_column(estatistica = rep(estat_names,2*3)) |> 
#   relocate(estatistica) |> 
#   writexl::write_xlsx("output/estat-desc.xlsx")
```

Boxplot para comparação

``` r
var_names <- data_set |> select(TCH:S) |> names()
map(var_names, ~{
    x<-data_set |> pull(!!sym(.x))
  q3 <- quantile(x,0.99)
  data_set |> 
    filter(!!sym(.x) <= q3) |>  
  ggplot(aes(x=as_factor(ano),y=!!sym(.x),fill=unidade)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#009E73", "#D55E00")) +
  labs(x="Ano", y=.x)
})
```

``` r
map(var_names, ~{
  x<-data_set |> pull(!!sym(.x))
  q3 <- quantile(x,0.99)
  data_set |> 
    filter(!!sym(.x) <= q3) |> 
    ggplot(aes(x=!!sym(.x),fill=unidade)) +
    geom_histogram(color="black") +
    scale_fill_manual(values=c("#009E73", "#D55E00")) +
    labs(x="Ano", y=.x)+
    facet_wrap(~ano,ncol=2)
})
```

``` r
map(var_names, ~{
  x<-data_set |> pull(!!sym(.x))
  q3 <- quantile(x,0.99)
  data_set |> 
    filter(!!sym(.x) <= q3) |> 
    ggplot(aes(y=as_factor(ano))) +
    geom_density_ridges(rel_min_height = 0.03,
                        aes(x=!!sym(.x), fill=unidade),
                        alpha = .6, color = "black"
    ) +
    scale_fill_cyclical(values = c("#ff8080","#238B45"),
                        name = "classe", guide = "legend") +
    labs(y="Ano", x=.x)+
    theme_ridges()
})
```

## 2. Matriz de correlação linear (corrplot) fazer por ano (2016, 2017 e 2018) e por UNIDADE

Para CAT

``` r
map(2016:2018, ~{data_set |> filter(unidade == "CAT") |> 
  filter(ano == .x) |> 
  select(TCH:S) |> 
  cor() |> 
  corrplot::corrplot( method = "color",
         outline = T,,
         addgrid.col = "darkgray",cl.pos = "r", tl.col = "black",
         tl.cex = 1, cl.cex = 1, type = "upper", bg="azure2",
         diag = FALSE,
         addCoef.col = "black",
         cl.ratio = 0.2,
         cl.length = 5,
         number.cex = 0.8)}) 
```

Para POT

``` r
map(2016:2018, ~{data_set |> filter(unidade == "POT") |> 
  filter(ano == .x) |> 
  select(TCH:S) |> 
  cor() |> 
  corrplot::corrplot( method = "color",
         outline = T,,
         addgrid.col = "darkgray",cl.pos = "r", tl.col = "black",
         tl.cex = 1, cl.cex = 1, type = "upper", bg="azure2",
         diag = FALSE,
         addCoef.col = "black",
         cl.ratio = 0.2,
         cl.length = 5,
         number.cex = 0.8)}) 
```

### 2.1 Análise de agrupamento hierárquico, independente do ano

``` r
da <- data_set |> 
  group_by(x,y) |> 
  summarise(
    across(.cols = TCH:S,
           .fns = mean,
           .names = "{.col}"),
    .groups = "drop"
  ) |> select(-x,-y)

da_pad <- decostand(da, 
                  method = "standardize",
                  na.rm=TRUE)
da_pad_euc <- vegdist(da_pad, "euclidean") 
da_pad_euc_ward<-hclust(da_pad_euc, method="ward.D")
plot(da_pad_euc_ward, 
     ylab="Distancia Euclidiana",
     xlab="Acessos", hang=-1,
     col="blue", las=1,
     cex=.6,lwd=1.5);box()
grupo <- cutree(da_pad_euc_ward,4)
```

### 2.2 Análise de Componentes Principais

``` r
print("======== Análise de Componentes Principais ========== ")
pca <-  prcomp(da_pad,scale.=TRUE)
# Autovalores
eig<-pca$sdev^2
print("==== Autovalores ====")
print(round(eig,3))
print("==== % da variância explicada ====")
ve<-eig/sum(eig)
print(round(ve,4))
print("==== % da variância explicada acumulada ====")
print(round(cumsum(ve),4)*100)
print("==== Poder Discriminante ====")
mcor<-cor(da_pad,pca$x)
corrplot(mcor)
print("==== screeplot ====")
screeplot(pca)
abline(h=1)
```

``` r
pc1V<-cor(da_pad,pca$x)[,1]/sd(cor(da_pad,pca$x)[,1])
pc2V<-cor(da_pad,pca$x)[,2]/sd(cor(da_pad,pca$x)[,2])
pc3V<-cor(da_pad,pca$x)[,3]/sd(cor(da_pad,pca$x)[,3])
pc1c<-pca$x[,1]/sd(pca$x[,1])
pc2c<-pca$x[,2]/sd(pca$x[,2])
pc3c<-pca$x[,3]/sd(pca$x[,3])
nv<-ncol(da) # número de variáveis utilizadas na análise

# gráfico biplot
bip<-data.frame(pc1c,pc2c,pc3c,grupo)
texto <- data.frame(
  x = pc1V,
  y = pc2V,
  z = pc3V,
  label = names(da)
)

bi_plot <- bip |> 
  ggplot(aes(x=pc1c,y=pc2c,color = as_factor(grupo)))+
  geom_point() + 
  theme_minimal() +
  # scale_shape_manual(values=16:18)+
  # scale_color_manual(values=c("#009E73", "#999999","#D55E00")) +
  annotate(geom="text", x=pc1V, y=pc2V, label=names(pc1V),
              color="black")+
  geom_vline(aes(xintercept=0),
             color="black", size=1)+
  geom_hline(aes(yintercept=0),
             color="black", size=1)+
  annotate(geom="segment",
           x=rep(0,length(da)),
           xend=texto$x,
           y=rep(0,length(da)),
           yend=texto$y,color="black",lwd=.5)+
  geom_label(data=texto,aes(x=x,y=y,label=label),
             color="black",angle=0,fontface="bold",size=4,fill="white")+
  labs(x=paste("CP1 (",round(100*ve[1],2),"%)",sep=""),
       y=paste("CP2 (",round(100*ve[2],2),"%)",sep=""),
       color="",shape="")+
  theme(legend.position = "top")
bi_plot +
  coord_cartesian(
    xlim = c(-4,3),
    ylim = c(-3,8)
  )
```

``` r
# gráfico biplot
bi_plot_2 <- bip |> 
  ggplot(aes(x=pc1c,y=pc3c,color = as_factor(grupo)))+
  geom_point() + 
  theme_minimal() +
  # scale_shape_manual(values=16:18)+
  # scale_color_manual(values=c("#009E73", "#999999","#D55E00")) +
  annotate(geom="text", x=pc1V, y=pc3V, label=names(pc1V),
              color="black")+
  geom_vline(aes(xintercept=0),
             color="black", size=1)+
  geom_hline(aes(yintercept=0),
             color="black", size=1)+
  annotate(geom="segment",
           x=rep(0,length(da)),
           xend=texto$x,
           y=rep(0,length(da)),
           yend=texto$z,color="black",lwd=.5)+
  geom_label(data=texto,aes(x=x,y=z,label=label),
             color="black",angle=0,fontface="bold",size=4,fill="white")+
  labs(x=paste("CP1 (",round(100*ve[1],2),"%)",sep=""),
       y=paste("CP3 (",round(100*ve[3],2),"%)",sep=""),
       color="",shape="")+
  theme(legend.position = "top")
bi_plot_2 +
  coord_cartesian(
    xlim = c(-4,3),
    ylim = c(-3,8)
  )
```

``` r
print("==== Tabela da correlação dos atributos com cada PC ====")
    ck<-sum(pca$sdev^2>=0.98)
    tabelapca<-vector()
    for( l in 1:ck) tabelapca<-cbind(tabelapca,mcor[,l])
    colnames(tabelapca)<-paste(rep(c("PC"),ck),1:ck,sep="")
    pcat<-round(tabelapca,3)
    tabelapca<-tabelapca[order(abs(tabelapca[,1])),]
    print(tabelapca)
writexl::write_xlsx(data.frame(tabelapca) |>
                      add_column(nome = row.names(tabelapca)) |>
                      relocate(nome), "output/pca-table.xlsx")
```

``` r
 data_set |> 
  group_by(x,y) |> 
  summarise(
    across(.cols = TCH:S,
           .fns = mean,
           .names = "{.col}"),
    .groups = "drop"
  ) |> 
  add_column(grupo) |> 
  ggplot(aes(x=x,y=y,color=as_factor(grupo))) +
  geom_point() # +
  # scale_shape_manual(values=16:18)+
  # scale_color_manual(values=c("#009E73", "#999999","#D55E00")) 
```

``` r
data_set |> 
  group_by(x,y) |> 
  summarise(
    across(.cols = TCH:S,
           .fns = mean,
           .names = "{.col}"),
    .groups = "drop"
  )|> 
  add_column(grupo) |> 
  select(x,y,grupo) |> 
  left_join(
    data_set |> 
  group_by(unidade,x,y) |> 
  summarise(
    across(.cols = TCH:S,
           .fns = mean,
           .names = "{.col}"),
    .groups = "drop"),
    by = c("x","y")
  ) |> 
  mutate(
    unidade = as.factor(unidade),
    grupo = as.factor(grupo)) |> 
  group_by(unidade, grupo) |> 
  summarise(
    n = n()
  ) |> 
  ggplot(aes(x=unidade,y=n,fill=grupo)) +
  geom_col(color="black") +
  scale_fill_viridis_d(option = "magma")
```

``` r
data_set |> 
  group_by(unidade,x,y) |> 
  summarise(
    across(.cols = TCH:S,
           .fns = mean,
           .names = "{.col}"),
    .groups = "drop"
  ) |> 
  filter(unidade == "POT") |>
  ggplot(aes(x,y)) +
  geom_point()
```

## 3 Aprendizado de Máquina Estatístico

### Definindo a Base de treino e teste

``` r
data_set_ml <- data_set  |>  
  janitor::clean_names() |> 
  drop_na() |> 
  filter(tch >= 25)
tch_initial_split <- initial_split(data_set_ml, prop = 0.80)
tch_train <- training(tch_initial_split)
# tch_test <- testing(tch_initial_split)
# visdat::vis_miss(tch_test)
tch_train  %>% 
  ggplot(aes(x=tch, y=..density..))+
  geom_histogram(bins = 30, color="black",  fill="lightgray")+
  geom_density(alpha=.05,fill="red")+
  theme_bw() +
  labs(x="tch - treino", y = "Densidade")
```

``` r
tch_testing <- testing(tch_initial_split)
tch_testing  |> 
  ggplot(aes(x=tch, y=..density..))+
  geom_histogram(bins = 30, color="black",  fill="lightgray")+
  geom_density(alpha=.05,fill="blue")+
  theme_bw() +
  labs(x="tch - teste", y = "Densidade")
```

``` r
tch_train |> 
  select(tch:s) |> 
  cor() |> 
  corrplot::corrplot( method = "color",
         outline = T,,
         addgrid.col = "darkgray",cl.pos = "r", tl.col = "black",
         tl.cex = 1, cl.cex = 1, type = "upper", bg="azure2",
         diag = FALSE,
         addCoef.col = "black",
         cl.ratio = 0.2,
         cl.length = 5,
         number.cex = 0.8)
```

### Definindo a `receita` da análise

``` r
tch_recipe <- recipe(tch ~ ., 
                      data = tch_train  |>  
            select(unidade,ambiente,textura,
                   corte,tch:s,variedade) 
) |>  
  step_naomit() %>%  
  step_zv(all_predictors())  |> 
  step_impute_median(where(is.numeric)) |> # inputação da mediana nos   
  step_normalize(all_numeric_predictors()) |> # padronização dos dados  
  step_novel(all_nominal_predictors())  |># trata categorias novas 
  step_dummy(all_nominal_predictors()) # variáveis categóricas
bake(prep(tch_recipe), new_data = NULL)
```

### Definindo a reamostragem

``` r
tch_resamples <- vfold_cv(tch_train, v = 5) 
```

## REDE NEURAL ARTIFICIAL

### Definição do Modelo de RNA - MultiLayer Perceptron

``` r
tch_nn_model <- mlp() |>  # margin sempre para regressão
  set_mode("regression") |> 
  set_engine("nnet")
```

### Definindo os parâmetros de tunagem

``` r
tch_nn_model <- mlp(
  hidden_units = tune(),
  penalty = tune(),
  epochs = tune()
  ) |>  # margin sempre para regressão
  set_mode("regression") |> 
  set_engine("nnet")
```

### Workflow e tunagem

``` r
tch_nn_wf <- workflow()   |> 
  add_model(tch_nn_model) |> 
  add_recipe(tch_recipe)

grid_nn <- grid_regular(
  hidden_units(range = c(2, 30)), ## tentar até 250
  penalty(range = c(-7, 0), trans = scales::log10_trans()), ## no máximo 30
  epochs(range = c(3, 20)),
  levels = c(3, 3, 3)
)

tch_nn_tune_grid <- tune_grid(
  tch_nn_wf,
  resamples = tch_resamples,
  grid = grid_nn,
  metrics = metric_set(rmse)
)
autoplot(tch_nn_tune_grid)
```

### Coletando métricas

``` r
collect_metrics(tch_nn_tune_grid)
tch_nn_tune_grid |> 
  show_best(metric = "rmse", n = 6)
```

### Desempenho do modelo final

``` r
tch_nn_best_params <- select_best(tch_nn_tune_grid, metric = "rmse")
tch_nn_wf <- tch_nn_wf |> 
  finalize_workflow(tch_nn_best_params)
tch_nn_last_fit <- last_fit(tch_nn_wf, tch_initial_split)

## Criando os preditos
tch_test_preds <- bind_rows(
  collect_predictions(tch_nn_last_fit)  |> 
    mutate(modelo = "nn"))

tch_test <- testing(tch_initial_split)

tch_test_preds |> 
  ggplot(aes(x=.pred, y=tch)) +
  geom_point()+
  theme_bw() +
  geom_smooth(method = "lm") +
  stat_regline_equation(ggplot2::aes(
  label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~"))) +
  geom_abline (slope=1, linetype = "dashed", color="Red")
```

``` r
tch_nn_last_fit_model <- tch_nn_last_fit$.workflow[[1]]$fit$fit
vip(tch_nn_last_fit_model,
    aesthetics = list(color = "black", fill = "orange")) +
    theme(axis.text.y=element_text(size=rel(1.5)),
          axis.text.x=element_text(size=rel(1.5)),
          axis.title.x=element_text(size=rel(1.5))
          ) +
  theme_bw()
```

### Principais Métricas

``` r
da <- tch_test_preds |> 
  filter(tch > 0, .pred>0 )

my_r <- cor(da$tch,da$.pred)
my_r2 <- my_r*my_r
my_mse <- Metrics::mse(da$tch,da$.pred)
my_rmse <- Metrics::rmse(da$tch,
                         da$.pred)
my_mae <- Metrics::mae(da$tch,da$.pred)
my_mape <- Metrics::mape(da$tch,da$.pred)*100

vector_of_metrics <- c(r=my_r, R2=my_r2, MSE=my_mse, RMSE=my_rmse, MAE=my_mae, MAPE=my_mape)
print(data.frame(vector_of_metrics))
```

## 4 Análise geoestatítica

### 4.1 Criando os contornos

#### Contorno para POT

``` r
contorno_pot <- data_set |> 
  group_by(unidade,x,y) |> 
  summarise(
    n=n(),
    .groups = "drop"
  ) |> 
  filter(unidade == "POT") |> 
  mutate(
    x1=mean(x),
    y1=mean(y),
    ang_rad = round(atan2(y-y1,x-x1),1),
    distancia = sqrt((x - x1)^2 + (y - y1)^2),
  ) |> 
  group_by(ang_rad) |> 
  mutate(
    max_dist = max(distancia)
  ) |> 
  ungroup() |> 
  mutate(flag_dist =distancia == max_dist) |> 
  filter(flag_dist)  |> 
  arrange(ang_rad) |> 
  select(x,y) |>   
  as.data.frame()
contorno_pot <- rbind(contorno_pot, contorno_pot[1,]) 

ggplot(contorno_pot, aes(x = x, y = y)) +
  geom_polygon(fill = "lightblue", color = "black") +
  coord_equal() +  # Mantém proporção real dos eixos
  theme_minimal()
```

#### Malha refinada para estimativa em POT

``` r
x<-contorno_pot$x
y<-contorno_pot$y
dis <- 0.005 #Distância entre pontos
grid_pot <- expand.grid(X=seq(min(x),max(x),dis), Y=seq(min(y),max(y),dis)) |> 
    mutate(flag = def_pol(X,Y,contorno_pot)) |>  
  filter(flag) |> select(-flag)
gridded(grid_pot) = ~ X + Y
plot(grid_pot) 
points(x,y,col="red",pch=4)
```

#### Contorno para CAT

``` r
contorno_cat <- data_set |> 
  group_by(unidade,x,y) |> 
  summarise(
    n=n(),
    .groups = "drop"
  ) |> 
  filter(unidade == "CAT") |> 
  mutate(
    x1=mean(x),
    y1=mean(y),
    ang_rad = round(atan2(y-y1,x-x1),1),
    distancia = sqrt((x - x1)^2 + (y - y1)^2),
  ) |> 
  group_by(ang_rad) |> 
  mutate(
    max_dist = max(distancia)
  ) |> 
  ungroup() |> 
  mutate(flag_dist =distancia == max_dist) |> 
  filter(flag_dist)  |> 
  arrange(ang_rad) |> 
  select(x,y) |>   
  as.data.frame()
contorno_cat <- rbind(contorno_cat, contorno_cat[1,]) 

ggplot(contorno_cat, aes(x = x, y = y)) +
  geom_polygon(fill = "pink", color = "black") +
  coord_equal() +  # Mantém proporção real dos eixos
  theme_minimal()
```

#### Malha refinada para estimativa em CAT

``` r
x<-contorno_cat$x
y<-contorno_cat$y
dis <- 0.005 #Distância entre pontos
grid_cat <- expand.grid(X=seq(min(x),max(x),dis), Y=seq(min(y),max(y),dis)) |> 
    mutate(flag = def_pol(X,Y,contorno_cat)) |>  
  filter(flag) |> select(-flag)
gridded(grid_cat) = ~ X + Y
plot(grid_cat) 
points(x,y,col="red",pch=4)
```

#### Contorno para CAT e POT (contorno total)

``` r
contorno_total <- rbind(contorno_cat |> 
                    mutate(unidade="CAT"),
                  contorno_pot |> 
                    mutate(unidade="POT"))
contorno_total |> ggplot(aes(x = x, y = y, group = unidade, fill = unidade)) +
  geom_polygon(color = "black", alpha = 0.5) +
  coord_equal() +
  theme_minimal()+
  geom_point(data = data_set |> 
               filter(ano == 2018),
             aes(x,y,color=unidade))
```

### 4.2 Comparação - Variáveis Originais e Log transformadas

``` r
map(var_names, ~{
  # x<-data_set |> pull(!!sym(.x))
  # x<-log(x+1)
  # q3 <- quantile(x,0.99)
  new_var <- paste0("log(", .x,"+1)")
  data_set |> 
    mutate(!!new_var := log(!!sym(.x) +1)) |> 
    select(!!sym(.x), !!sym(new_var)) |> 
    pivot_longer(cols = c(!!sym(.x), !!sym(new_var)),
                 names_to = "variavel",
                 values_to = "valor") |> 
    ggplot(aes(valor,fill = variavel)) +
    geom_histogram(color="black") +
    # labs(x="Ano", y=paste0("log(",.x,"+1)")) +
      facet_wrap(~variavel, scale="free") +
    scale_fill_brewer(palette = "Dark2")
})
```

### 4.3 Análise geoestatística CAT

#### Separando o banco e transformar se necessário

``` r
ano_analise <- 2018
variavel <- "m%"
unidade_analise <- "CAT"
data_set_aux <- data_set |> 
  filter(
    ano == ano_analise,
    unidade == unidade_analise) |> 
  select(x,y,variavel) |> 
  rename(z = !!sym(variavel)) |> 
  group_by(x,y) |> 
  summarise(
    z = mean(z,na.rm = TRUE),
    .groups = "drop"
  ) |> mutate(z = log(z+1)) ## TRANSFORMAR SE NECESSÁRIO
glimpse(data_set_aux)
```

#### Criando a fórmula para a análise

``` r
coordinates(data_set_aux) = ~ x + y  
form <- z ~ 1 # fórmula da função variogram
```

#### Construção do Semivariograma Experimental

``` r
vari_exp <- variogram(form, data = data_set_aux,
                      cressie = FALSE,
                      cutoff = 0.50, # distância máxima do semivariograma
                      width = .01) # distância entre pontos
vari_exp  |>  
 ggplot(aes(x=dist, y=gamma)) +
 geom_point() +
 labs(x="lag (º)",
      y=expression(paste(gamma,"(h)")))
```

#### Ajustando os modelodos

``` r
patamar=0.1
alcance=0.02
epepita=0
modelo_1 <- fit.variogram(vari_exp,vgm(patamar,"Sph",alcance,epepita))
modelo_2 <- fit.variogram(vari_exp,vgm(patamar,"Exp",alcance,epepita))
modelo_3 <- fit.variogram(vari_exp,vgm(patamar,"Gau",alcance,epepita))
plot_my_models(modelo_1,modelo_2,modelo_3)
```

#### Validação Cruzada

``` r
# my_cross_validation(data_set_aux,form,modelo_1,modelo_2,modelo_3)
```

#### Seleção do melhor modelo

O melhor modelo é aquele que apresenta um coeficiente de regressão o
mais próximo de $1$, o intercepto o mais próximo de $0$, com o menor
valor de $RMSE$ possível.

``` r
modelo <- modelo_1 ## sempre modificar
```

#### Definido o melhor modelo, precisamos guardar os valores dos parâmetros.

Salvamos os parâmetros, e a figura o semivariograma + modelo ajustado

``` r
model <- modelo |> slice(2) |> pull(model)
rss <- round(attr(modelo, "SSErr"),4)
c0 <- round(modelo$psill[[1]],4)
c0_c1 <- round(sum(modelo$psill),4)
a <- ifelse(model == "Gau", round(modelo$range[[2]]*(3^.5),2),
            ifelse(model == "Exp",round(3*modelo$range[[2]],2),
            round(modelo$range[[2]],2)))
r2 <- vari_exp |> add_column( model = model, a=a, c0 = c0,
                                  c0_c1 = c0_c1) |>
    mutate(
      gamma_m = ifelse(model == "Sph",
        ifelse(dist <= a, c0 + (c0_c1 - c0) * (3/2 * (dist/a) - 1/2 * (dist/a)^3),c0_c1), ifelse(model == "Exp", c0 + (c0_c1-c0)*(1-exp(-3*(dist/a))),c0 + (c0_c1-c0)*(1-exp(-(dist/a)^2)))),
      residuo_total = (gamma-mean(gamma))^2,
      residuo_mod = (gamma - gamma_m)^2
    ) |>
    summarise(
      r2=(sum(residuo_total) - sum(residuo_mod))/sum(residuo_total)
    ) |> pull(r2)

tibble(
  ano_analise, unidade_analise, variavel, model, c0, c0_c1, a, rss, r2
) |> mutate(gde = c0/c0_c1, .after = "a") |>
  write_csv(paste0("output/best-fit/",ano_analise,"-",
                   unidade_analise,"-",str_remove(variavel,"%"),".csv"))

ls_csv <- list.files("output/best-fit/",full.names = TRUE,pattern = ".csv")
map_df(ls_csv, read_csv) |>
  writexl::write_xlsx("output/semivariogram-models.xlsx")
png(filename = paste0("output/semivariogram-img/semivar-",
                      ano_analise,"-",
                   unidade_analise,"-",str_remove(variavel,"%"),".png"),
    width = 800, height = 600)
plot(vari_exp,model=modelo,cex.lab=2, col=1,pl=F,pch=16,cex=2.2,ylab=list("Semivariância",cex=2.3),xlab=list("Distância de Separação h (m)",cex=2.3,cex.axis=4))
dev.off()
```

#### Krigragem ordinária (KO)

Estimar o atributo para locais não amostrados.

``` r
# ko_variavel <- krige(formula=form, data_set_aux, grid_cat, model=modelo,
#     block=c(0,0),
#     nsim=0,
#     na.action=na.pass,
#     debug.level=-1,
#     )
```

#### Salva o mapa e os valores estimados

``` r
# mapa <- as_tibble(ko_variavel) |> 
#   ggplot(aes(x=X, y=Y)) + 
#   geom_tile(aes(fill = var1.pred)) +
#   # scale_fill_gradient(low = "yellow", high = "blue") + 
#   scale_fill_viridis_c() +
#   coord_equal() + 
#   labs(fill=variavel,
#        x="Longitude",
#        y="Latitude")
# mapa
# ggsave(paste0("output/krigagem/krg-",ano_analise,"-",
#                    unidade_analise,"-",str_remove(variavel,"%"),".png"))
```

``` r
# df <- ko_variavel |> 
#   as_tibble() |> 
#   mutate(var1.var = sqrt(var1.var)) |> 
#   rename(
#     !!variavel := var1.pred,
#     !!paste0(variavel,"_sd") := var1.var,
#   )
# write_rds(df,paste0("output/krigagem/estimativa-",ano_analise,"-",
#                    unidade_analise,"-",str_remove(variavel,"%"),".rds"))
```

### 4.4 Análise geoestatística para POT

#### Separando o banco e transformar se necessário

``` r
unidade_analise <- "POT"
data_set_aux <- data_set |> 
  filter(
    ano == ano_analise,
    unidade == unidade_analise) |> 
  select(x,y,variavel) |> 
  rename(z = !!sym(variavel)) |> 
  group_by(x,y) |> 
  summarise(
    z = mean(z,na.rm = TRUE),
    .groups = "drop"
  ) |> mutate(z = log(z+1)) ## TRANSFORMAR SE NECESSÁRIO
glimpse(data_set_aux)
```

#### Criando a fórmula para a análise

``` r
coordinates(data_set_aux) = ~ x + y  
form <- z ~ 1 # fórmula da função variogram
```

#### Construção do Semivariograma Experimental

``` r
vari_exp <- variogram(form, data = data_set_aux,
                      cressie = FALSE,
                      cutoff = 0.05, # distância máxima do semivariograma
                      width = .004) # distância entre pontos
vari_exp  |>  
 ggplot(aes(x=dist, y=gamma)) +
 geom_point() +
 labs(x="lag (º)",
      y=expression(paste(gamma,"(h)")))
```

#### Ajustando os modelodos

``` r
patamar=600
alcance=0.02
epepita=0
modelo_1 <- fit.variogram(vari_exp,vgm(patamar,"Sph",alcance,epepita))
modelo_2 <- fit.variogram(vari_exp,vgm(patamar,"Exp",alcance,epepita))
modelo_3 <- fit.variogram(vari_exp,vgm(patamar,"Gau",alcance,epepita))
plot_my_models(modelo_1,modelo_2,modelo_3)
```

#### Validação Cruzada

``` r
# my_cross_validation(data_set_aux,form,modelo_1,modelo_2,modelo_3)
```

#### Seleção do melhor modelo

O melhor modelo é aquele que apresenta um coeficiente de regressão o
mais próximo de $1$, o intercepto o mais próximo de $0$, com o menor
valor de $RMSE$ possível.

``` r
modelo <- modelo_1 ## sempre modificar
```

#### Definido o melhor modelo, precisamos guardar os valores dos parâmetros.

Salvamos os parâmetros, e a figura o semivariograma + modelo ajustado

``` r
model <- modelo |> slice(2) |> pull(model)
rss <- round(attr(modelo, "SSErr"),4)
c0 <- round(modelo$psill[[1]],4)
c0_c1 <- round(sum(modelo$psill),4)
a <- ifelse(model == "Gau", round(modelo$range[[2]]*(3^.5),2),
            ifelse(model == "Exp",round(3*modelo$range[[2]],2),
            round(modelo$range[[2]],2)))
r2 <- vari_exp |> add_column( model = model, a=a, c0 = c0,
                                  c0_c1 = c0_c1) |>
    mutate(
      gamma_m = ifelse(model == "Sph",
        ifelse(dist <= a, c0 + (c0_c1 - c0) * (3/2 * (dist/a) - 1/2 * (dist/a)^3),c0_c1), ifelse(model == "Exp", c0 + (c0_c1-c0)*(1-exp(-3*(dist/a))),c0 + (c0_c1-c0)*(1-exp(-(dist/a)^2)))),
      residuo_total = (gamma-mean(gamma))^2,
      residuo_mod = (gamma - gamma_m)^2
    ) |>
    summarise(
      r2=(sum(residuo_total) - sum(residuo_mod))/sum(residuo_total)
    ) |> pull(r2)

tibble(
  ano_analise, unidade_analise, variavel, model, c0, c0_c1, a, rss, r2
) |> mutate(gde = c0/c0_c1, .after = "a") |>
  write_csv(paste0("output/best-fit/",ano_analise,"-",
                   unidade_analise,"-",str_remove(variavel,"%"),".csv"))

ls_csv <- list.files("output/best-fit/",full.names = TRUE,pattern = ".csv")
map_df(ls_csv, read_csv) |>
  writexl::write_xlsx("output/semivariogram-models.xlsx")
png(filename = paste0("output/semivariogram-img/semivar-",
                      ano_analise,"-",
                   unidade_analise,"-",str_remove(variavel,"%"),".png"),
    width = 800, height = 600)
plot(vari_exp,model=modelo,cex.lab=2, col=1,pl=F,pch=16,cex=2.2,ylab=list("Semivariância",cex=2.3),xlab=list("Distância de Separação h (m)",cex=2.3,cex.axis=4))
dev.off()
```

#### Krigragem ordinária (KO)

Estimar o atributo para locais não amostrados.

``` r
# ko_variavel <- krige(formula=form, data_set_aux, grid_pot, model=modelo,
#     block=c(0,0),
#     nsim=0,
#     na.action=na.pass,
#     debug.level=-1,
#     )
```

#### Salva o mapa e os valores estimados

``` r
# mapa <- as_tibble(ko_variavel) |>
#   ggplot(aes(x=X, y=Y)) +
#   geom_tile(aes(fill = var1.pred)) +
#   # scale_fill_gradient(low = "yellow", high = "blue") +
#   scale_fill_viridis_c() +
#   coord_equal() +
#   labs(fill=variavel,
#        x="Longitude",
#        y="Latitude")
# mapa
# ggsave(paste0("output/krigagem/krg-",ano_analise,"-",
#                    unidade_analise,"-",str_remove(variavel,"%"),".png"))
```

``` r
# df <- ko_variavel |>
#   as_tibble() |>
#   mutate(var1.var = sqrt(var1.var)) |>
#   rename(
#     !!variavel := var1.pred,
#     !!paste0(variavel,"_sd") := var1.var,
#   )
# write_rds(df,paste0("output/krigagem/estimativa-",ano_analise,"-",
#                    unidade_analise,"-",str_remove(variavel,"%"),".rds"))
```
