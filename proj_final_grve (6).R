# PROJET GRVE 
# STEPHAN Corentin, LECOEUR Bertrand, DROIN Charlotte

# ------------------------------------------------------------------------------
# ÉTAPE 1 : PRÉPARATION 
# ------------------------------------------------------------------------------
# On charge les "boîtes à outils" nécessaires :
# - readxl : Pour lire le fichier Excel fourni.
# - ismev & evd : Les références pour la Théorie des Valeurs Extrêmes (GEV, GPD).
# - copula : Pour modéliser la dépendance entre X et Y.
# - MASS : Pour des outils statistiques standards.
suppressPackageStartupMessages({
  if (!require("readxl")) install.packages("readxl")
  library(readxl)
  if (!require("ismev")) install.packages("ismev")
  library(ismev)
  if (!require("evd")) install.packages("evd")
  library(evd)
  if (!require("copula")) install.packages("copula")
  library(copula)
  library(MASS)
})


# ------------------------------------------------------------------------------
# ÉTAPE 2 : IMPORTATION ET NETTOYAGE DES DONNÉES
# ------------------------------------------------------------------------------
cat("\n 1 - CHARGEMENT DES DONNÉES\n")

# Lecture sécurisée le tryCatch empêche le script de planter si le fichier est absent.
data <- tryCatch(
  read_excel("Base_X_Y.xlsx"),
  error = function(e) stop("ERREUR: 'Base_X_Y.xlsx' introuvable.")
)

# On renomme les colonnes pour être sûr de les appeler correctement (Date, X, Y).
colnames(data) <- c("Date", "X", "Y") 

# Suppression des lignes vides pour éviter les erreurs de calcul.
data <- na.omit(data)

# Conversion des formats : Dates en format temporel, Montants en numérique.
data$Date <- as.Date(data$Date)
X <- as.numeric(data$X)
Y <- as.numeric(data$Y)

# On garde le nombre total d'observations (n) et on extrait les années pour la GEV.
n <- nrow(data)
annees <- format(data$Date, "%Y")


################################################################################
# PARTIE A : ANALYSE SÉRIE PAR SÉRIE
################################################################################
cat(" 2 - CALCULS EN COURS\n")

# ------------------------------------------------------------------------------
# ANALYSE DE LA SÉRIE X
# ------------------------------------------------------------------------------

# 1. Modèle Loi Normale
# On calcule la moyenne et l'écart-type classiques.
# On en déduit le quantile 99.5% théorique 
mu_X <- mean(X); sd_X <- sd(X)
q995_norm_X <- qnorm(0.995, mu_X, sd_X)
q995_emp_X <- quantile(X, 0.995)

# 2. Modèle GEV (Generalized Extreme Value) Méthode des Maxima par Blocs
# On extrait la perte maximale observée pour chaque année.
X_max_ann <- as.numeric(tapply(X, annees, max))

# On ajuste la loi GEV sur ces maxima annuels.
# L'option "BFGS" est une méthode d'optimisation numérique robuste.
suppressWarnings({ 
  fit_gev_X <- tryCatch(
    gev.fit(X_max_ann, show = FALSE), 
    error = function(e) gev.fit(X_max_ann, show = FALSE, method = "BFGS")
  ) 
})
# On récupère les paramètres (Position, Echelle, Forme) pour calculer le quantile.
loc_X <- fit_gev_X$mle[1]; sc_X <- fit_gev_X$mle[2]; sh_X <- fit_gev_X$mle[3]
q995_gev_X <- qgev(0.995, loc_X, sc_X, sh_X)


# Affichage
cat("Probabilité mensuelle équivalente =", p_month, "\n")
cat("Quantile GEV mensuel équivalent =", q995_gev_X_month, "\n")


# 3. Modèle GPD (Generalized Pareto Distribution) Méthode POT (Seuil)
tester_seuils_gpd <- function(x,
                              probs = seq(0.80, 0.98, by = 0.05),  # mets by=0.01 si tu veux tous les 1%
                              methode = c("ismev", "evd"),
                              show = FALSE) {
  methode <- match.arg(methode)
  n <- length(x)
  
  # seuils candidats
  u_grid <- as.numeric(quantile(x, probs = probs, type = 8, na.rm = TRUE))
  
  res <- data.frame(
    prob = probs,
    u = u_grid,
    nexc = NA_integer_,
    p_exc = NA_real_,
    beta = NA_real_,
    xi = NA_real_,
    ok = FALSE
  )
  
  for (k in seq_along(u_grid)) {
    u <- u_grid[k]
    
    fit <- tryCatch({
      if (methode == "ismev") {
        ismev::gpd.fit(x, threshold = u, show = show)
      } else {
        tmp <- evd::fpot(x, threshold = u, std.err = FALSE)
        list(mle = tmp$estimate, nexc = tmp$nat)
      }
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      res$nexc[k] <- fit$nexc
      res$p_exc[k] <- fit$nexc / n
      res$beta[k] <- fit$mle[1]
      res$xi[k] <- fit$mle[2]
      res$ok[k] <- TRUE
    }
  }
  
  return(res)
}

res_X <- tester_seuils_gpd(X, probs = seq(0.80, 0.98, by = 0.01), methode = "ismev")
# La zone de linéarité du Mean Excess Plot débute dès le quantile 90 %, ce qui 
# justifie le choix d’un seuil plus bas afin de réduire la variance de l’estimation.
u_X <- quantile(X, 0.90, type = 8)

# Bloc de sécurité technique :
# Parfois, l'optimisation mathématique "ismev" échoue (matrice singulière).
# Si cela arrive, on bascule automatiquement sur "evd::fpot" qui est plus tolérante.
fit_gpd_X <- tryCatch({
  ismev::gpd.fit(X, threshold = u_X, show = FALSE)
}, error = function(e) {
  tmp <- evd::fpot(X, threshold = u_X, std.err = FALSE)
  list(mle = tmp$estimate, nexc = tmp$nat)
})

# Extraction des paramètres GPD (Scale, Shape) et calcul du quantile 99.5%.
beta_X <- fit_gpd_X$mle[1]; xi_X <- fit_gpd_X$mle[2]
nb_exces_X <- fit_gpd_X$nexc; prob_exces_X <- nb_exces_X / n
q995_gpd_X <- u_X + (beta_X / xi_X) * ( ((1 - 0.995) / prob_exces_X)^(-xi_X) - 1 )

# La vraie perte de retour 200 ANS
# 200 ans = 2400 mois. La probabilité cible est 1 - 1/2400
prob_200_ans <- 1 - (1 / (200 * 12))
perte_200_ans_X <- u_X + (beta_X / xi_X) * ( ((1 - prob_200_ans) / prob_exces_X)^(-xi_X) - 1 )
# ------------------------------------------------------------------------------
# ANALYSE DE LA SÉRIE Y (Même logique que X)
# ------------------------------------------------------------------------------

# 1. Modèle Normal Y
mu_Y <- mean(Y); sd_Y <- sd(Y)
q995_norm_Y <- qnorm(0.995, mu_Y, sd_Y)
q995_emp_Y <- quantile(Y, 0.995)

# 2. Modèle GEV Y Méthode des Maxima par Blocs
# On extrait la perte maximale observée pour chaque année.
Y_max_ann <- as.numeric(tapply(Y, annees, max))
suppressWarnings({ 
  fit_gev_Y <- tryCatch(
    gev.fit(Y_max_ann, show = FALSE), 
    error = function(e) gev.fit(Y_max_ann, show = FALSE, method = "BFGS")
  ) 
})
# On récupère les paramètres (Position, Echelle, Forme) pour calculer le quantile.
loc_Y <- fit_gev_Y$mle[1]; sc_Y <- fit_gev_Y$mle[2]; sh_Y <- fit_gev_Y$mle[3]
q995_gev_Y <- qgev(0.995, loc_Y, sc_Y, sh_Y)

# 3. Modèle GPD Y 

res_Y <- tester_seuils_gpd(
  Y,
  probs = seq(0.80, 0.98, by = 0.01),
  methode = "ismev"
)

res_Y

u_Y <- quantile(Y, 0.90, type = 8)

# Bloc de sécurité technique :
# Parfois, l'optimisation mathématique "ismev" échoue (matrice singulière).
# Si cela arrive, on bascule automatiquement sur "evd::fpot" qui est plus tolérante.
fit_gpd_Y <- tryCatch({
  ismev::gpd.fit(Y, threshold = u_Y, show = FALSE)
}, error = function(e) {
  tmp <- evd::fpot(Y, threshold = u_Y, std.err = FALSE)
  list(mle = tmp$estimate, nexc = tmp$nat)
})

# Extraction des paramètres GPD (Scale, Shape) et calcul du quantile 99.5%.
beta_Y <- fit_gpd_Y$mle[1]
xi_Y   <- fit_gpd_Y$mle[2]

nb_exces_Y   <- fit_gpd_Y$nexc
prob_exces_Y <- nb_exces_Y / n

q995_gpd_Y <- u_Y + (beta_Y / xi_Y) *
  ( ((1 - 0.995) / prob_exces_Y)^(-xi_Y) - 1 )

# La vraie perte de retour 200 ans
# On réutilise la probabilité calculée plus haut (1 - 1/2400)
perte_200_ans_Y <- u_Y + (beta_Y / xi_Y) * ( ((1 - prob_200_ans) / prob_exces_Y)^(-xi_Y) - 1 )

# ============================================================================== 
# PARTIE B : DÉPENDANCE (COPULES) ET SIMULATION 
# ============================================================================== 

cat(" 3 - CALCULS COPULES & SIMULATION\n") 
# 1. Transformation en Rangs (Pseudo-observations) 
# On transforme les montants en probabilités uniformes [0,1] pour supprimer l'effet 
# de l'échelle et ne garder que la structure de corrélation pure. 
U <- rank(X) / (n + 1) 
V <- rank(Y) / (n + 1) 
data_cop <- cbind(U, V) 

# 2. Calibrage des Copules 
# On teste 4 familles célèbres pour voir laquelle colle le mieux aux données. 
cop_models <- list(Gumbel=gumbelCopula(dim=2), Frank=frankCopula(dim=2), Normal=normalCopula(dim=2), Clayton=claytonCopula(dim=2)) 
res_aic <- c(); res_bic <- c(); fits <- list() 
loglik_indep <- sum(log(U * V)) 
res_aic["Indep"] <- -2 * loglik_indep 
res_bic["Indep"] <- -2 * loglik_indep 

for(nom in names(cop_models)) { 
  suppressWarnings({ fits[[nom]] <- fitCopula(cop_models[[nom]], 
  data_cop, method = "ml") }) 
  
# On stocke les critères AIC/BIC (plus c'est bas, meilleur est le modèle). 
  res_aic[nom] <- AIC(fits[[nom]]) 
  res_bic[nom] <- BIC(fits[[nom]]) 
}
  
# Classement des copules du meilleur au moins bon. 
tab_res <- data.frame(Copule = names(res_aic), AIC = round(res_aic, 2), BIC = round(res_bic, 2)) 
tab_res <- tab_res[order(tab_res$AIC), ] 
  
# 3. Simulation de Monte-Carlo 
# On choisit Gumbel (souvent la meilleure pour les extrêmes) et on récupère son paramètre Theta. 
theta_gumbel <- coef(fits[["Gumbel"]]) 
cop_final <- gumbelCopula(theta_gumbel, dim = 2)
# On simule 120 000 couples (U,V) fictifs (10 000 années x 12 mois). 
nsim <- 10000 * 12 
UV_sim <- rCopula(nsim, cop_final) 
  
# 4. Fonction Inverse
# Cette fonction transforme les probas simulées [0,1] en montants réels. 
# - Si la proba est extrême (> seuil), on utilise la formule GPD inverse. 
# - Si la proba est normale (< seuil), on utilise l'historique empirique. 
inv_cdf <- function(p, dat, u, b, x_shape) { 
  res <- numeric(length(p))
  p_seuil <- mean(dat <= u) 
  idx <- p > p_seuil 
  res[idx] <- u + (b / x_shape) * ( ((1 - p[idx]) / (1 - p_seuil))^(-x_shape) - 1 ) 
  res[!idx] <- quantile(dat, p[!idx], na.rm = TRUE) 
  return(res) 
} 
  
# Application de la transformation inverse sur les données simulées. 
X_sim <- inv_cdf(UV_sim[,1], X, u_X, beta_X, xi_X)
Y_sim <- inv_cdf(UV_sim[,2], Y, u_Y, beta_Y, xi_Y)
  
# 5. Agrégation et Calcul Final 
# On additionne X+Y, on groupe par paquets de 12 (Années) et on prend le Max annuel.
n_extra <- (12 - length(X_sim) %% 12) %% 12
Z_ann <- matrix(c(X_sim + Y_sim, rep(NA, n_extra)), ncol = 12, byrow = TRUE) 
Z_max <- apply(Z_ann, 1, max, na.rm = TRUE) 
# Le quantile 99.5% de ces max annuels donne la perte de retour 200 ans. 
perte_200 <- quantile(Z_max, 0.995)



# ==============================================================================
# PARTIE GRAPHIQUE
# ==============================================================================
cat(" 4 - GENERATION DES GRAPHIQUES\n")

plot_normal_diag <- function(date_vec, serie, mu_hat, sd_hat, titre="Serie") {
  dev.new(title = paste("Normal diagnostics -", titre))
  par(mfrow=c(2,2))
  plot(date_vec, serie, type="l",
       main=paste(titre, "(mensuelle)"), xlab="Date", ylab=titre)
  hist(serie, breaks=30, freq=FALSE,
       main=paste("Histogramme", titre, "+ densite Normale"), xlab=titre)
  curve(dnorm(x, mean=mu_hat, sd=sd_hat), add=TRUE, lwd=2)
  qqnorm(serie, main="QQ-plot vs Normale")
  qqline(serie)
  plot(density(serie), main="Densite empirique vs Normale", xlab=titre)
  curve(dnorm(x, mean=mu_hat, sd=sd_hat), add=TRUE, lwd=2)
  legend("topright", legend=c("Empirique","Normale"), lwd=c(1,2), bty="n")
  par(mfrow=c(1,1))
}

plot_gev_diag_safe <- function(fit_gev, titre="") {
  dev.new(title = paste("GEV diagnostics -", titre))
  # ismev::gev.diag fait plusieurs graphes de diagnostic
  tryCatch({
    ismev::gev.diag(fit_gev)
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("GEV diag impossible:", e$message))
  })
}

# ------------------------------------------------------------------
# Loi normale : X et Y
# ------------------------------------------------------------------
plot_normal_diag(data$Date, X, mu_X, sd_X, titre="X")
plot_normal_diag(data$Date, Y, mu_Y, sd_Y, titre="Y")

# ------------------------------------------------------------------
# Loi GEV : maxima annuels (X et Y)
# ------------------------------------------------------------------
plot_gev_diag_safe(fit_gev_X, titre="X (maxima annuels)")
plot_gev_diag_safe(fit_gev_Y, titre="Y (maxima annuels)")

# ------------------------------------------------------------------
# Graphes de stabilité : xi(u), beta(u), nexc(u) - X (GPD)
# ------------------------------------------------------------------
dev.new()
plot(res_X$u, res_X$xi, type="b", pch=19,
     xlab="Seuil u", ylab="xi (shape)",
     main="Stabilité du paramètre xi (série X)")
grid()

dev.new()
plot(res_X$u, res_X$beta, type="b", pch=19,
     xlab="Seuil u", ylab="beta (scale)",
     main="Stabilité du paramètre beta (série X)")
grid()

dev.new()
plot(res_X$u, res_X$nexc, type="b", pch=19,
     xlab="Seuil u", ylab="Nombre d'excès",
     main="Nombre d'excès au-dessus de u (série X)")
grid()
# Le seuil 90 % est retenu car il correspond à une zone de stabilité des paramètres 
# GPD, tout en garantissant un nombre suffisant d’excès.

# ------------------------------------------------------------------
# Graphes de stabilité : xi(u), beta(u), nexc(u) - Y (GPD)
# ------------------------------------------------------------------
dev.new()
plot(res_Y$u, res_Y$xi, type="b", pch=19,
     xlab="Seuil u", ylab="xi (shape)",
     main="Stabilité du paramètre xi (série Y)")
grid()

dev.new()
plot(res_Y$u, res_Y$beta, type="b", pch=19,
     xlab="Seuil u", ylab="beta (scale)",
     main="Stabilité du paramètre beta (série Y)")
grid()

dev.new()
plot(res_Y$u, res_Y$nexc, type="b", pch=19,
     xlab="Seuil u", ylab="Nombre d'excès",
     main="Nombre d'excès au-dessus de u (série Y)")
grid()
# Le seuil 90 % est retenu car il correspond à une zone de stabilité des paramètres
# GPD, tout en garantissant un nombre suffisant d’excès.


# Mean Excess Plot X. Permet de verifier visuellement si le seuil u est bon
dev.new(title = "1. Mean Excess Plot - X")
mrlplot(X, main = "Mean Excess Plot (Serie X)")
abline(v = u_X, col = "red", lty = 2)
legend("topright", "Seuil choisi", col = "red", lty = 2, bty = "n")

# Mean Excess Plot Y. 
dev.new(title = "2. Mean Excess Plot - Y")
mrlplot(Y, main = "Mean Excess Plot (Serie Y)")
abline(v = u_Y, col = "red", lty = 2)
legend("topright", "Seuil choisi", col = "red", lty = 2, bty = "n")

# Fenetre 3 : Pseudo-Observations. Visualisation brute de la dependance des rangs.
dev.new(title = "3. Pseudo-Observations")
plot(U, V, pch = 19, cex = 0.4, col = "blue", main = "Pseudo-Observations (Rangs)", xlab = "Rang X", ylab = "Rang Y")
grid()

# Fenetre 4 : Comparaison des Copules.
# On compare le nuage de points reel (noir) avec des simulations des modeles (bleu).
# Cela aide a valider visuellement le choix de la copule (ex: forme en pointe pour Gumbel).
dev.new(title = "4. Comparaison Copules")
par(mfrow = c(2, 3))
plot(U, V, pch = 19, cex = 0.2, col = "black", main = "DONN?ES R?ELLES", xlab = "", ylab = "")
for(nom in names(fits)) {
  s_tmp <- rCopula(2000, fits[[nom]]@copula)
  plot(s_tmp, pch = 19, cex = 0.2, col = "blue", main = paste("Simul:", nom), xlab = "", ylab = "")
}
plot(matrix(runif(4000), ncol = 2), pch = 19, cex = 0.2, col = "gray", main = "Simul: Indep", xlab = "", ylab = "")

# Fenetre 5 : Resultat Final. Distribution des pertes maximales annuelles simul?es.
dev.new(title = "5. R?sultat Final")
par(mfrow = c(1, 1))
hist(Z_max, breaks = 50, col = "orange", border = "white", prob = TRUE, main = "Distribution Pertes Annuelles (X+Y)", xlab = "Montant")
lines(density(Z_max), col = "darkred", lwd = 2)
abline(v = perte_200, col = "red", lwd = 3, lty = 2)
legend("topright", legend = paste("VaR 200 ans\n", format(round(perte_200, 0), big.mark = " ")), text.col = "red", bty = "n")

################################################################################
# AFFICHAGE FINAL 
################################################################################
cat("\n\n################################################################")
cat("\n                 RÉSULTATS DE L'ANALYSE                       ")

cat("\n--- SÉRIE X ---\n")
cat("Moyenne (Mean) :", round(mu_X, 2), "\n")
cat("Ecart-type (Sd) :", round(sd_X, 2), "\n")
cat("\nNormale - Quantile 99.5% :", round(q995_norm_X, 2), "\n")
cat("\nGEV - Paramètres : Location =", round(loc_X, 2), "| Scale =", round(sc_X, 2), "| Shape =", round(sh_X, 4), "\n")
cat("GEV - Quantile 99.5% :", round(q995_gev_X, 2), "\n")
cat("GPD - Différents seuils testés (grille POT) :\n")
print(res_X)
cat("\nGPD - Paramètres : Seuil (u) =", round(u_X, 2), "| Scale (beta) =", round(beta_X, 2), "| Shape (xi) =", round(xi_X, 4), "\n")
cat("GPD - Quantile 99.5% :", round(q995_gpd_X, 2), "\n")

cat("\n--- SÉRIE Y ---\n")
cat("Moyenne (Mean) :", round(mu_Y, 2), "\n")
cat("Ecart-type (Sd) :", round(sd_Y, 2), "\n")
cat("\nNormale - Quantile 99.5% :", round(q995_norm_Y, 2), "\n")
cat("\nGEV - Paramètres : Location =", round(loc_Y, 2), "| Scale =", round(sc_Y, 2), "| Shape =", round(sh_Y, 4), "\n")
cat("GEV - Quantile 99.5% :", round(q995_gev_Y, 2), "\n")
cat("GPD - Différents seuils testés (grille POT) :\n")
print(res_Y)
cat("\nGPD - Paramètres : Seuil (u) =", round(u_Y, 2), "| Scale (beta) =", round(beta_Y, 2), "| Shape (xi) =", round(xi_Y, 4), "\n")
cat("GPD - Quantile 99.5% :", round(q995_gpd_Y, 2), "\n")

cat("\n---  COMPARAISON DES LOIS (QUANTILES 99.5%) ---\n")
res_lois <- data.frame(
  Modele = c("Empirique", "Normale", "GEV", "GPD"),
  Serie_X = format(round(c(q995_emp_X, q995_norm_X, q995_gev_X, q995_gpd_X), 2), nsmall = 2),
  Serie_Y = format(round(c(q995_emp_Y, q995_norm_Y, q995_gev_Y, q995_gpd_Y), 2), nsmall = 2)
)
print(res_lois)

cat("\n--- NIVEAUX DE RETOUR 200 ANS ---\n")
cat("La perte de retour 200 ans pour X est :", format(round(perte_200_ans_X, 2), big.mark=" "), "\n")
cat("La perte de retour 200 ans pour Y est :", format(round(perte_200_ans_Y, 2), big.mark=" "), "\n")
# ---------------------------------------------------------------------

cat("\n--- . PARAMÈTRES RETENUS (GPD) ---\n")
cat("X -> Seuil u:", round(u_X, 2), "| Scale (beta):", round(beta_X, 2), "| Shape (xi):", xi_X, "\n")
cat("Y -> Seuil u:", round(u_Y, 2), "| Scale (beta):", round(beta_Y, 2), "| Shape (xi):", xi_Y, "\n")

cat("\n--- . TABLEAU DE SÉLECTION DES COPULES (AIC/BIC) ---\n")
print(tab_res)
cat(">>> Copule retenue pour simulation : Gumbel (Theta =", round(theta_gumbel, 3), ")\n")

cat("\n----------------------------------------------------------------\n")
cat(">>> RÉSULTAT FINAL DU PROJET\n")
cat(">>> Perte Annuelle Combinée (Retour 200 ans) :", format(round(perte_200, 1), big.mark = " "), "\n")
cat("----------------------------------------------------------------\n")

