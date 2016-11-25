## according to Martinez & Reyes-Valdes 2014
Att0 <- apply(peaklist[, 20:79][1:5], 1, mean)
Att72 <- apply(peaklist[, 20:79][6:10], 1, mean)
Obt0 <- apply(peaklist[, 20:79][11:15], 1, mean)
Obt72 <- apply(peaklist[, 20:79][16:20], 1, mean)
Clev0 <- apply(peaklist[, 20:79][21:25], 1, mean)
Clev72 <- apply(peaklist[, 20:79][26:30], 1, mean)
Quad0 <- apply(peaklist[, 20:79][31:35], 1, mean)
Quad72 <- apply(peaklist[, 20:79][36:40], 1, mean)
X10270 <- apply(peaklist[, 20:79][41:45], 1, mean)
X102772 <- apply(peaklist[, 20:79][46:50], 1, mean)
X571260 <- apply(peaklist[, 20:79][51:55], 1, mean)
X5712672 <- apply(peaklist[, 20:79][56:59], 1, mean) ## do not use 60

## relative frequency pij 
relFreq <- function(speciesFeat) {speciesFeat / sum(speciesFeat)}
relFreqAtt0 <- relFreq(Att0)
relFreqAtt72 <- relFreq(Att72)
relFreqObt0 <- relFreq(Obt0)
relFreqObt72 <- relFreq(Obt72)
relFreqClev0 <- relFreq(Clev0)
relFreqClev72 <- relFreq(Clev72)
relFreqQuad0 <- relFreq(Quad0)
relFreqQuad72 <- relFreq(Quad72)
relFreqX10270 <- relFreq(X10270)
relFreqX102772 <- relFreq(X102772)
relFreqX571260 <- relFreq(X571260)
relFreqX5712672 <- relFreq(X5712672)

## H_j
H_j <- function(relFreq) {- sum(relFreq * log2(relFreq))}
(H_att0 <- H_j(relFreqAtt0))
(H_att72 <- H_j(relFreqAtt72))
(H_obt0 <- H_j(relFreqObt0))
(H_obt72 <- H_j(relFreqObt72))
(H_clev0 <- H_j(relFreqClev0))
(H_clev72 <- H_j(relFreqClev72))
(H_quad0 <- H_j(relFreqQuad0))
(H_quad72 <- H_j(relFreqQuad72))
(H_x10270 <- H_j(relFreqX10270))
(H_x102772 <- H_j(relFreqX102772))
(H_x571260 <- H_j(relFreqX571260))
(H_x5712672 <- H_j(relFreqX5712672))
H_js <- c(H_att0, H_att72, H_obt0, H_obt72, H_clev0, H_clev72, H_quad0, 
          H_quad72, H_x10270, H_x102772, H_x571260, H_x5712672)

H_js0 <- c(H_att0, H_obt0, H_clev0, H_quad0, H_x10270, H_x571260 )
H_js72 <- c(H_att72, H_obt72, H_clev72, H_quad72, H_x102772, H_x5712672)

## p_i, average frequency 
p_i <- function(mm) {1 / dim(mm)[2] * apply(mm, 1, sum)}
mm <- matrix(c(relFreqAtt0, relFreqAtt72, relFreqObt0, relFreqObt72, relFreqClev0, relFreqClev72,
               relFreqQuad0, relFreqQuad72, relFreqX10270, relFreqX102772, 
               relFreqX571260, relFreqX5712672), ncol = 12)
mm0 <- matrix(c(relFreqAtt0, relFreqObt0, relFreqClev0, 
                relFreqQuad0, relFreqX10270, 
                relFreqX571260), ncol = 6)
mm72 <- matrix(c(relFreqAtt72, relFreqObt72, relFreqClev72, 
                 relFreqQuad72, relFreqX102772, 
                 relFreqX5712672), ncol = 6)
prob_i <- p_i(mm)
prob_i0 <- p_i(mm0)
prob_i72 <- p_i(mm72)

## gene/metabolite specificity: S_i = 1/t (sum_j=1^t of p_ij/p_i * log2(p_ij/p_i))
S_i <- function(mm, prob_i) {
    1 / dim(mm)[2] * apply(apply(mm, 2, function(x) {x / prob_i * log2(x / prob_i)}), 1, sum)}
metaboliteSpecificity <- S_i(mm, prob_i)
metaboliteSpecificity0 <- S_i(mm0, prob_i0)
metaboliteSpecificity72 <- S_i(mm72, prob_i72)
## 'S_i will give a value of zero if the gene (metabolite) is transcribed (detected) at the same
## frequency in all tissues and a maximum value of log2(t) if the
## gene (metabolite) is exclusively expressed (detected) in a single tissue"

## tissue specialisation
## delta_j varies between zero (= all metabolites are completely unspecific)
## and log2(t) when all metabolites in the tissue are not present anywhere else
delta_j <- function(relFreq, S_i) {sum(relFreq * S_i)}
(delta_att0 <- delta_j(relFreqAtt0, metaboliteSpecificity))
(delta_att72 <- delta_j(relFreqAtt72, metaboliteSpecificity))
(delta_obt0 <- delta_j(relFreqObt0, metaboliteSpecificity))
(delta_obt72 <- delta_j(relFreqObt72, metaboliteSpecificity))
(delta_clev0 <- delta_j(relFreqClev0, metaboliteSpecificity))
(delta_clev72 <- delta_j(relFreqClev72, metaboliteSpecificity))
(delta_quad0 <- delta_j(relFreqQuad0, metaboliteSpecificity))
(delta_quad72 <- delta_j(relFreqQuad72, metaboliteSpecificity))
(delta_x10270 <- delta_j(relFreqX10270, metaboliteSpecificity))
(delta_x102772 <- delta_j(relFreqX102772, metaboliteSpecificity))
(delta_x571260 <- delta_j(relFreqX571260, metaboliteSpecificity))
(delta_x5712672 <- delta_j(relFreqX5712672, metaboliteSpecificity))
## only for 0h
delta_att0_0 <- delta_j(relFreqAtt0, metaboliteSpecificity0)
delta_obt0_0 <- delta_j(relFreqObt0, metaboliteSpecificity0)
delta_clev0_0 <- delta_j(relFreqClev0, metaboliteSpecificity0)
delta_quad0_0 <- delta_j(relFreqQuad0, metaboliteSpecificity0)
delta_x10270_0 <- delta_j(relFreqX10270, metaboliteSpecificity0)
delta_x571260_0 <- delta_j(relFreqX571260, metaboliteSpecificity0)
## only for 72h
delta_att0_72 <- delta_j(relFreqAtt72, metaboliteSpecificity72)
delta_obt0_72 <- delta_j(relFreqObt72, metaboliteSpecificity72)
delta_clev0_72 <- delta_j(relFreqClev72, metaboliteSpecificity72)
delta_quad0_72 <- delta_j(relFreqQuad72, metaboliteSpecificity72)
delta_x10270_72 <- delta_j(relFreqX102772, metaboliteSpecificity72)
delta_x571260_72 <- delta_j(relFreqX5712672, metaboliteSpecificity72)

delta_js <- c(delta_att0, delta_att72, delta_obt0, delta_obt72, delta_clev0, delta_clev72,
              delta_quad0, delta_quad72, delta_x10270, delta_x102772, delta_x571260, delta_x5712672)
delta_js0 <- c(delta_att0, delta_obt0, delta_clev0, delta_quad0, delta_x10270, delta_x571260)
delta_js72 <- c(delta_att72, delta_obt72, delta_clev72, delta_quad72, delta_x102772, delta_x5712672)
plot(H_js, delta_js, col = c("blue", "blue", "red", "red", "yellow", "yellow", "green", "green", "darkgrey", "darkgrey", "lightgrey", "lightgrey"), xlab = "diversity", ylab = "specificity")
plot(H_js72, delta_js72, col = c("blue", "red", "yellow", "green", "darkgrey", "lightgrey"),
     xlab = "diversity", ylab = "specificity", main = "72h")
plot(H_js0, delta_js0, col = c("blue", "red", "yellow", "green", "darkgrey", "lightgrey"), 
     xlab = "diversity", ylab = "specificity", main = "0h")

## look at specialised metabolites for N. clevelandii and N. quadrivalvis (seperately) and see 
## how they differ with respect to N. attenuata and N. obtusifolia after induction (72h)

## N. clevelandii
peaklistClev72 <- peaklist[, c(20:79)[c(6:10, 16:20, 26:30)]]

## calculate Hj
H_jsClev72 <- c(H_att72, H_obt72, H_clev72)

mmClev72 <- matrix(c(relFreqAtt72, relFreqObt72, relFreqClev72), ncol = 3)
prob_iClev72 <- p_i(mmClev72)
metaboliteSpecificityClev72 <- S_i(mmClev72, prob_iClev72)

specialisedMetabolitesClev72 <- peaklist[which(metaboliteSpecificityClev72 > 0.6),c(1,4)]
##cbind(peaklist[, c(1,4,8:13)], cbind(Att0, Att72, Obt0, Obt72, Clev0, Clev72))[as.numeric(rownames(specialisedMetabolitesClev72)),]

## tissue specialisation
delta_j(relFreqAtt72, metaboliteSpecificityClev72)
delta_j(relFreqObt72, metaboliteSpecificityClev72)
delta_j(relFreqClev72, metaboliteSpecificityClev72)

## N. quadrivalvis
H_jsQuad72 <- c(H_att72, H_obt72, H_quad72)
mmQuad72 <- matrix(c(relFreqAtt72, relFreqObt72, relFreqQuad72), ncol = 3)
prob_iQuad72 <- p_i(mmQuad72)
metaboliteSpecificityQuad72 <- S_i(mmQuad72, prob_iQuad72)

specialisedMetabolitesQuad72 <- peaklist[which(metaboliteSpecificityQuad72 > 0.4),c(1,4)]
cbind(peaklist[, c(1,4,8:11, 14:15)], cbind(Att0, Att72, Obt0, Obt72, Quad0, Quad72))[as.numeric(rownames(specialisedMetabolitesQuad72)),]

## specificity
delta_j(relFreqAtt72, metaboliteSpecificityQuad72)
delta_j(relFreqObt72, metaboliteSpecificityQuad72)
delta_j(relFreqClev72, metaboliteSpecificityQuad72)

## N. x obt 10/27
H_jsX102772 <- c(H_att72, H_obt72, H_x102772)
mmX102772 <- matrix(c(relFreqAtt72, relFreqObt72, relFreqX102772), ncol = 3)
prob_iX102772 <- p_i(mmX102772)
metaboliteSpecificityX102772 <- S_i(mmX102772, prob_iX102772)

specialisedMetabolitesX102772 <- peaklist[which(metaboliteSpecificityX102772 > 0.6),c(1,4)]
cbind(peaklist[, c(1,4,8:11, 16:17)], cbind(Att0, Att72, Obt0, Obt72, X10270, X102772))[as.numeric(rownames(specialisedMetabolitesX102772)),]

## tissue specialisation
delta_j(relFreqAtt72, metaboliteSpecificityX102772)
delta_j(relFreqObt72, metaboliteSpecificityX102772)
delta_j(relFreqClev72, metaboliteSpecificityX102772)

## N. x obt 57/126
H_jsX5712672 <- c(H_att72, H_obt72, H_x5712672)
mmX5712672 <- matrix(c(relFreqAtt72, relFreqObt72, relFreqX5712672), ncol = 3)
prob_iX5712672 <- p_i(mmX5712672)
metaboliteSpecificityX5712672 <- S_i(mmX5712672, prob_iX5712672)

specialisedMetabolitesX5712672 <- peaklist[which(metaboliteSpecificityX5712672 > 0.6),c(1,4)]
cbind(peaklist[, c(1,4,8:11, 16:17)], cbind(Att0, Att72, Obt0, Obt72, X571260, X5712672))[as.numeric(rownames(specialisedMetabolitesX5712672)),]
specialisedMetX5712672_matrix <- as.matrix(peaklist[as.numeric(rownames(specialisedMetabolitesX5712672)),  c("mz", "rt", colnames(peaklist)[cols[c(6:10, 16:20, 56:60)]])])
colnames(specialisedMetX5712672_matrix) <- c("mz", "rt", rep("N. attenuata", 5), rep("N. obtusifolia", 5), rep("N. x obtusiata 57126", 5))
rownames(specialisedMetX5712672_matrix) <- paste(round(specialisedMetX5712672_matrix[, "mz"], 3), round(specialisedMetX5712672_matrix[, "rt"], 2), sep = "/")
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("blue", "white","red"))(n = 299)
col_breaks <-  c(seq(1,1.05,length=100),  # for blue
                 seq(1.051,3.66,length=100), # for white
                 seq(3.661,900,length=100)) 
heatmap.2(specialisedMetX5712672_matrix[, -c(1:2)], col = my_palette, breaks = col_breaks,
          Colv="NA", trace = "none", dendrogram = "row")
## specialised metabolites in N. clev 72 compared to progenitors: 
peaklist[3990, c("mz", "rt", colnames(peaklist)[cols[c(6:10, 16:20, 46:50)]])] ## also present in N. attenuata and N. obtusifolia but much higer in Nxobt1027
peaklist[4003, c("mz", "rt", colnames(peaklist)[cols[c(6:10, 16:20, 46:50)]])] ## also present in N. attenuata and N. obtusifolia but much higher in Nxobt1027
peaklist[4019, c("mz", "rt", colnames(peaklist)[cols[c(6:10, 16:20, 46:50)]])] 
peaklist[4020, c("mz", "rt", colnames(peaklist)[cols[c(6:10, 16:20, 46:50)]])]
peaklist[4208, c("mz", "rt", colnames(peaklist)[cols[c(6:10, 16:20, 46:50)]])]
peaklist[9359, c("mz", "rt", colnames(peaklist)[cols[c(6:10, 16:20, 46:50)]])]
peaklist[9460, c("mz", "rt", colnames(peaklist)[cols[c(6:10, 16:20, 46:50)]])]

## tissue specialisation
delta_j(relFreqAtt72, metaboliteSpecificityX102772)
delta_j(relFreqObt72, metaboliteSpecificityX102772)
delta_j(relFreqClev72, metaboliteSpecificityX102772)