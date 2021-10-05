# Libraries
library(tidyverse)

# Make a vector of all amino acids
AA <- c( "G", "A", "V", "L", "I", "S", "T", "C", "M", "F", "Y", "W", "P", "N", "Q", "D", "E", "K", "H", "R")

# Total amount of possible combinations?
N <- 20^4 + 20 ^ 3 + 20 ^ 2 + 20 ^ 1
N

###############################################################################
##### 1 AMINO ACIDS
# All possible combinations
Combinations_1 <- crossing(AA1 = AA)

# Giant for-loop shenanigans
Combinations_weight_1 <- data.frame("AA1" = character(), "Mass" = character(), "AtomicComp" = character())
for(rows in 1:nrow(Combinations_1)){
  Weight <- 0
  C <- 0
  H <- 0
  N <- 0
  O <- 0
  S <- 0
  
  polypeptide <- Combinations_1[rows, 1]
  
  # Second looptydoo
  for(AminoAcid in polypeptide){
    if(AminoAcid == "G"){
      Weight <- Weight + 75.07
      C <- C + 2
      H <- H + 5
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "A"){
      Weight <- Weight + 89.09
      C <- C + 3
      H <- H + 7
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "V"){
      Weight <- Weight + 117.15
      C <- C + 5
      H <- H + 11
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "L"){
      Weight <- Weight + 131.18
      C <- C + 6
      H <- H + 13
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "I"){
      Weight <- Weight + 131.18
      C <- C + 6
      H <- H + 13
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "S"){
      Weight <- Weight + 105.09
      C <- C + 3
      H <- H + 7
      N <- N + 1
      O <- O + 3
    } else if(AminoAcid == "T"){
      Weight <- Weight + 119.12
      C <- C + 4
      H <- H + 9
      N <- N + 1
      O <- O + 3
    } else if(AminoAcid == "C"){
      Weight <- Weight + 121.15
      C <- C + 3
      H <- H + 7
      N <- N + 1
      O <- O + 2
      S <- S + 1
    } else if(AminoAcid == "M"){
      Weight <- Weight + 149.21
      C <- C + 5
      H <- H + 11
      N <- N + 1
      O <- O + 2
      S <- S + 1
    } else if(AminoAcid == "F"){
      Weight <- Weight + 165.19
      C <- C + 9
      H <- H + 11
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "Y"){
      Weight <- Weight + 181.19
      C <- C + 9
      H <- H + 11
      N <- N + 1
      O <- O + 3
    } else if(AminoAcid == "W"){
      Weight <- Weight + 204.23
      C <- C + 11
      H <- H + 12
      N <- N + 2
      O <- O + 2
    } else if(AminoAcid == "P"){
      Weight <- Weight + 115.13
      C <- C + 5
      H <- H + 9
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "N"){
      Weight <- Weight + 132.12
      C <- C + 4
      H <- H + 8
      N <- N + 2
      O <- O + 3
    } else if(AminoAcid == "Q"){
      Weight <- Weight + 146.15
      C <- C + 5
      H <- H + 10
      N <- N + 2
      O <- O + 3
    } else if(AminoAcid == "D"){
      Weight <- Weight + 133.10
      C <- C + 4
      H <- H + 7
      N <- N + 1
      O <- O + 4
    } else if(AminoAcid == "E"){
      Weight <- Weight + 147.13
      C <- C + 5
      H <- H + 9
      N <- N + 1
      O <- O + 4
    } else if(AminoAcid == "K"){
      Weight <- Weight + 146.19
      C <- C + 6
      H <- H + 14
      N <- N + 2
      O <- O + 2
    } else if(AminoAcid == "H"){
      Weight <- Weight + 155.16
      C <- C + 6
      H <- H + 9
      N <- N + 3
      O <- O + 2
    } else if(AminoAcid == "R"){
      Weight <- Weight + 174.20
      C <- C + 6
      H <- H + 14
      N <- N + 4
      O <- O + 2
    } 
  }
  
  # Add to new DF
  AtomComp <- paste0("C", C, "H", H, "N", N, "O", O, "S", S)
  Temp <- data.frame(polypeptide, Weight, AtomComp)
  names(Temp) <- c("AA1", "Mass","AtomicComp")
  Combinations_weight_1 <- rbind(Combinations_weight_1, Temp)
  rm(Temp, Weight, C, H, N, O, S)
}

Combinations_weight_1 %>% unite(AA1, col = "Peptide", sep = "") -> Combinations_final_1

# Making the error interval
Combinations_final_1 %>% mutate("LowerError" = Mass - 0.005, "UpperError" = Mass + 0.005) -> Combinations_final_1

rm(Combinations_1)
rm(Combinations_weight_1)



###############################################################################
##### 2 AMINO ACIDS
# All possible combinations
Combinations_2 <- crossing(AA1 = AA, AA2 = AA)


# Giant for-loop shenanigans
Combinations_weight_2 <- data.frame("AA1" = character(),"AA2" = character(), "Mass" = character(), "AtomicComp" = character())
for(rows in 1:nrow(Combinations_2)){
  Weight <- 0
  C <- 0
  H <- 0
  N <- 0
  O <- 0
  S <- 0
  
  polypeptide <- Combinations_2[rows, 1:2]
  
  # Second looptydoo
  for(AminoAcid in polypeptide){
    if(AminoAcid == "G"){
      Weight <- Weight + 75.07
      C <- C + 2
      H <- H + 5
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "A"){
      Weight <- Weight + 89.09
      C <- C + 3
      H <- H + 7
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "V"){
      Weight <- Weight + 117.15
      C <- C + 5
      H <- H + 11
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "L"){
      Weight <- Weight + 131.18
      C <- C + 6
      H <- H + 13
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "I"){
      Weight <- Weight + 131.18
      C <- C + 6
      H <- H + 13
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "S"){
      Weight <- Weight + 105.09
      C <- C + 3
      H <- H + 7
      N <- N + 1
      O <- O + 3
    } else if(AminoAcid == "T"){
      Weight <- Weight + 119.12
      C <- C + 4
      H <- H + 9
      N <- N + 1
      O <- O + 3
    } else if(AminoAcid == "C"){
      Weight <- Weight + 121.15
      C <- C + 3
      H <- H + 7
      N <- N + 1
      O <- O + 2
      S <- S + 1
    } else if(AminoAcid == "M"){
      Weight <- Weight + 149.21
      C <- C + 5
      H <- H + 11
      N <- N + 1
      O <- O + 2
      S <- S + 1
    } else if(AminoAcid == "F"){
      Weight <- Weight + 165.19
      C <- C + 9
      H <- H + 11
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "Y"){
      Weight <- Weight + 181.19
      C <- C + 9
      H <- H + 11
      N <- N + 1
      O <- O + 3
    } else if(AminoAcid == "W"){
      Weight <- Weight + 204.23
      C <- C + 11
      H <- H + 12
      N <- N + 2
      O <- O + 2
    } else if(AminoAcid == "P"){
      Weight <- Weight + 115.13
      C <- C + 5
      H <- H + 9
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "N"){
      Weight <- Weight + 132.12
      C <- C + 4
      H <- H + 8
      N <- N + 2
      O <- O + 3
    } else if(AminoAcid == "Q"){
      Weight <- Weight + 146.15
      C <- C + 5
      H <- H + 10
      N <- N + 2
      O <- O + 3
    } else if(AminoAcid == "D"){
      Weight <- Weight + 133.10
      C <- C + 4
      H <- H + 7
      N <- N + 1
      O <- O + 4
    } else if(AminoAcid == "E"){
      Weight <- Weight + 147.13
      C <- C + 5
      H <- H + 9
      N <- N + 1
      O <- O + 4
    } else if(AminoAcid == "K"){
      Weight <- Weight + 146.19
      C <- C + 6
      H <- H + 14
      N <- N + 2
      O <- O + 2
    } else if(AminoAcid == "H"){
      Weight <- Weight + 155.16
      C <- C + 6
      H <- H + 9
      N <- N + 3
      O <- O + 2
    } else if(AminoAcid == "R"){
      Weight <- Weight + 174.20
      C <- C + 6
      H <- H + 14
      N <- N + 4
      O <- O + 2
    } 
  }
  # Adjust for water being split off
  Weight <- Weight - 1 * 18.00988
  H <- H - 2
  O <- O - 1
  
  # Add to new DF
  AtomComp <- paste0("C", C, "H", H, "N", N, "O", O, "S", S)
  Temp <- data.frame(polypeptide, Weight, AtomComp)
  names(Temp) <- c("AA1", "AA2", "Mass","AtomicComp")
  Combinations_weight_2 <- rbind(Combinations_weight_2, Temp)
  rm(Weight, Temp, C, H, N, O, S)
}

Combinations_weight_2 %>% unite(AA1:AA2, col = "Peptide", sep = "") -> Combinations_final_2

# Making the error interval
Combinations_final_2 %>% mutate("LowerError" = Mass - 0.00005, "UpperError" = Mass + 0.00005) -> Combinations_final_2

rm(Combinations_2)
rm(Combinations_weight_2)



###############################################################################
##### 3 AMINO ACIDS
# All possible combinations
Combinations_3 <- crossing(AA1 = AA, AA2 = AA, AA3 = AA)


# Giant for-loop shenanigans
Combinations_weight_3 <- data.frame("AA1" = character(),"AA2" = character(),"AA3" = character(), "Mass" = character(), "AtomicComp" = character())
for(rows in 1:nrow(Combinations_3)){
  Weight <- 0
  C <- 0
  H <- 0
  N <- 0
  O <- 0
  S <- 0
  
  polypeptide <- Combinations_3[rows, 1:3]
  
  # Second looptydoo
  for(AminoAcid in polypeptide){
    if(AminoAcid == "G"){
      Weight <- Weight + 75.07
      C <- C + 2
      H <- H + 5
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "A"){
      Weight <- Weight + 89.09
      C <- C + 3
      H <- H + 7
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "V"){
      Weight <- Weight + 117.15
      C <- C + 5
      H <- H + 11
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "L"){
      Weight <- Weight + 131.18
      C <- C + 6
      H <- H + 13
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "I"){
      Weight <- Weight + 131.18
      C <- C + 6
      H <- H + 13
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "S"){
      Weight <- Weight + 105.09
      C <- C + 3
      H <- H + 7
      N <- N + 1
      O <- O + 3
    } else if(AminoAcid == "T"){
      Weight <- Weight + 119.12
      C <- C + 4
      H <- H + 9
      N <- N + 1
      O <- O + 3
    } else if(AminoAcid == "C"){
      Weight <- Weight + 121.15
      C <- C + 3
      H <- H + 7
      N <- N + 1
      O <- O + 2
      S <- S + 1
    } else if(AminoAcid == "M"){
      Weight <- Weight + 149.21
      C <- C + 5
      H <- H + 11
      N <- N + 1
      O <- O + 2
      S <- S + 1
    } else if(AminoAcid == "F"){
      Weight <- Weight + 165.19
      C <- C + 9
      H <- H + 11
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "Y"){
      Weight <- Weight + 181.19
      C <- C + 9
      H <- H + 11
      N <- N + 1
      O <- O + 3
    } else if(AminoAcid == "W"){
      Weight <- Weight + 204.23
      C <- C + 11
      H <- H + 12
      N <- N + 2
      O <- O + 2
    } else if(AminoAcid == "P"){
      Weight <- Weight + 115.13
      C <- C + 5
      H <- H + 9
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "N"){
      Weight <- Weight + 132.12
      C <- C + 4
      H <- H + 8
      N <- N + 2
      O <- O + 3
    } else if(AminoAcid == "Q"){
      Weight <- Weight + 146.15
      C <- C + 5
      H <- H + 10
      N <- N + 2
      O <- O + 3
    } else if(AminoAcid == "D"){
      Weight <- Weight + 133.10
      C <- C + 4
      H <- H + 7
      N <- N + 1
      O <- O + 4
    } else if(AminoAcid == "E"){
      Weight <- Weight + 147.13
      C <- C + 5
      H <- H + 9
      N <- N + 1
      O <- O + 4
    } else if(AminoAcid == "K"){
      Weight <- Weight + 146.19
      C <- C + 6
      H <- H + 14
      N <- N + 2
      O <- O + 2
    } else if(AminoAcid == "H"){
      Weight <- Weight + 155.16
      C <- C + 6
      H <- H + 9
      N <- N + 3
      O <- O + 2
    } else if(AminoAcid == "R"){
      Weight <- Weight + 174.20
      C <- C + 6
      H <- H + 14
      N <- N + 4
      O <- O + 2
    } 
  }
  # Adjust for water being split off
  Weight <- Weight - 2 * 18.00988
  H <- H - 4
  O <- O - 2
  
  # Add to new DF
  AtomComp <- paste0("C", C, "H", H, "N", N, "O", O, "S", S)
  Temp <- data.frame(polypeptide, Weight, AtomComp)
  names(Temp) <- c("AA1", "AA2", "AA3","Mass","AtomicComp")
  Combinations_weight_3 <- rbind(Combinations_weight_3, Temp)
  rm(Weight, Temp, C, H, N, O, S)
}

Combinations_weight_3 %>% unite(AA1:AA3, col = "Peptide", sep = "") -> Combinations_final_3

# Making the error interval
Combinations_final_3 %>% mutate("LowerError" = Mass - 0.00005, "UpperError" = Mass + 0.00005) -> Combinations_final_3

rm(Combinations_3)
rm(Combinations_weight_3)

############################################################################3
##### 4 AMINO ACIDS
# All possible combinations
Combinations_4 <- crossing(AA1 = AA, AA2 = AA, AA3 = AA, AA4 = AA)

# Giant for-loop shenanigans
Combinations_weight_4 <- data.frame("AA1" = character(),"AA2" = character(),"AA3" = character(),"AA4" = character(), "Mass" = character(), "AtomicComp" = character())
for(rows in 1:nrow(Combinations_4)){
  Weight <- 0
  C <- 0
  H <- 0
  N <- 0
  O <- 0
  S <- 0
  
  polypeptide <- Combinations_4[rows, 1:4]
  
  # Second looptydoo
  for(AminoAcid in polypeptide){
    if(AminoAcid == "G"){
      Weight <- Weight + 75.07
      C <- C + 2
      H <- H + 5
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "A"){
      Weight <- Weight + 89.09
      C <- C + 3
      H <- H + 7
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "V"){
      Weight <- Weight + 117.15
      C <- C + 5
      H <- H + 11
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "L"){
      Weight <- Weight + 131.18
      C <- C + 6
      H <- H + 13
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "I"){
      Weight <- Weight + 131.18
      C <- C + 6
      H <- H + 13
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "S"){
      Weight <- Weight + 105.09
      C <- C + 3
      H <- H + 7
      N <- N + 1
      O <- O + 3
    } else if(AminoAcid == "T"){
      Weight <- Weight + 119.12
      C <- C + 4
      H <- H + 9
      N <- N + 1
      O <- O + 3
    } else if(AminoAcid == "C"){
      Weight <- Weight + 121.15
      C <- C + 3
      H <- H + 7
      N <- N + 1
      O <- O + 2
      S <- S + 1
    } else if(AminoAcid == "M"){
      Weight <- Weight + 149.21
      C <- C + 5
      H <- H + 11
      N <- N + 1
      O <- O + 2
      S <- S + 1
    } else if(AminoAcid == "F"){
      Weight <- Weight + 165.19
      C <- C + 9
      H <- H + 11
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "Y"){
      Weight <- Weight + 181.19
      C <- C + 9
      H <- H + 11
      N <- N + 1
      O <- O + 3
    } else if(AminoAcid == "W"){
      Weight <- Weight + 204.23
      C <- C + 11
      H <- H + 12
      N <- N + 2
      O <- O + 2
    } else if(AminoAcid == "P"){
      Weight <- Weight + 115.13
      C <- C + 5
      H <- H + 9
      N <- N + 1
      O <- O + 2
    } else if(AminoAcid == "N"){
      Weight <- Weight + 132.12
      C <- C + 4
      H <- H + 8
      N <- N + 2
      O <- O + 3
    } else if(AminoAcid == "Q"){
      Weight <- Weight + 146.15
      C <- C + 5
      H <- H + 10
      N <- N + 2
      O <- O + 3
    } else if(AminoAcid == "D"){
      Weight <- Weight + 133.10
      C <- C + 4
      H <- H + 7
      N <- N + 1
      O <- O + 4
    } else if(AminoAcid == "E"){
      Weight <- Weight + 147.13
      C <- C + 5
      H <- H + 9
      N <- N + 1
      O <- O + 4
    } else if(AminoAcid == "K"){
      Weight <- Weight + 146.19
      C <- C + 6
      H <- H + 14
      N <- N + 2
      O <- O + 2
    } else if(AminoAcid == "H"){
      Weight <- Weight + 155.16
      C <- C + 6
      H <- H + 9
      N <- N + 3
      O <- O + 2
    } else if(AminoAcid == "R"){
      Weight <- Weight + 174.20
      C <- C + 6
      H <- H + 14
      N <- N + 4
      O <- O + 2
    } 
  }
  # Adjust for water being split off
  Weight <- Weight - 3 * 18.00988
  H <- H - 6
  O <- O - 3
  
  # Add to new DF
  AtomComp <- paste0("C", C, "H", H, "N", N, "O", O, "S", S)
  Temp <- data.frame(polypeptide, Weight, AtomComp)
  names(Temp) <- c("AA1", "AA2", "AA3", "AA4", "Mass","AtomicComp")
  Combinations_weight_4 <- rbind(Combinations_weight_4, Temp)
  rm(Weight, Temp, C, H, N, O, S)
}

Combinations_weight_4 %>% unite(AA1:AA4, col = "Peptide", sep = "") -> Combinations_final_4

# Making the error interval
Combinations_final_4 %>% mutate("LowerError" = Mass - 0.00005, "UpperError" = Mass + 0.00005) -> Combinations_final_4
rm(AminoAcid, rows, Combinations_weight_4, Combinations_4, polypeptide)




# Merging all to one.
Combinations <- data.frame("Peptide" = character(),"Mass" = double(),"LowerError" = double(),"UpperError" = double(),"AtomicComp" = character())
Combinations <- bind_rows(Combinations, Combinations_final_1)
Combinations <- bind_rows(Combinations, Combinations_final_2)
Combinations <- bind_rows(Combinations, Combinations_final_3)
Combinations <- bind_rows(Combinations, Combinations_final_4)

Combinations %>% gather(key = "Weight", value = "Mass", LowerError, UpperError, Mass) -> test

########### Data analyses
# Distribution of mass
windows()
Combinations %>% ggplot(aes(x = Mass)) + geom_histogram(fill = "dodgerblue", alpha = 0.5) + theme_minimal() + labs(title = "Distribution of mass of polypeptides", subtitle = "Up to a length of 4 amino acids") + scale_x_continuous(limits = c(50, 800))

# Plot error intervals??


# Count the occurences of masses.
Combinations %>% select(Mass) %>% count(factor(Mass), sort = TRUE) %>% slice(1:5)


