# Calculates the probability that B is better than A and expected loss for A and B
# Also returns samples used for calculations
{
  # Arguments:
  #   Alpha and beta parameters of beta distributions for A and B
  #   k and theta parameters of gamma distributions for A and B
  #   alpha = number of successes + 1
  #   beta = number of failures + 1
  #   k = payer count + 1
  #   theta = 1/(1 + revenue)
  # Assumptions:
  #   payer ~ Bernoulli(lambda)
  #   lambda ~ Beta(alpha, beta)
  #   revenue ~ Exp(omega)
  #   omega ~ Gamma(k, theta)
  #   ARPU = E(payer)*E(revenue) = lambda*1/omega
  #   all varianles are calculated for A and B group
  # Return value:
  #   probBbeatsA = probability that the B version has higher ARPU than the A version
  #   expLossA = expected loss on condition that we select A but B is better
  #   expLossB = expected loss on condition that we select B but A is better
  #   expLossC = expected loss on condition that we select C but A is better
  # Example of use:
  #   If version A converted 180 out of 10000 users with revenue 5000 euros and
  #     version B converted 200 out of 11000 users with revenue 6000 euros then use
  #     bayes_arpu(
  #       180 + 1, 10000 - 180 + 1, 180 + 1, 1/(1 + 5000),
  #       200 + 1, 11000 - 200 + 1, 200 + 1, 1/(1 + 6000)
  #     )
  # Implemented based on VWO white paper
  # https://cdn2.hubspot.net/hubfs/310840/VWO_SmartStats_technical_whitepaper.pdf
}
bayes_arpu <- function(
  alphaA, betaA,
  kA, thetaA,
  alphaB, betaB,
  kB, thetaB,
  alphaC, betaC,
  kC, thetaC,
  alphaD, betaD,
  kD, thetaD,
  MSamples
) {
  if(alphaA <= 0 || betaA <= 0 || alphaB <= 0 || betaB <= 0 || alphaC <= 0 || betaC <= 0 || alphaD <= 0 || betaD <= 0 ||
     kA <= 0 || thetaA <= 0 || kB <= 0 || thetaB <= 0 || kC <= 0 || thetaC <= 0 || kD <= 0 || thetaD <= 0) {
    probBbeatsA <- 0
    probCbeatsA <- 0
    probDbeatsA <- 0
    probCbeatsB <- 0
    probDbeatsB <- 0
    probDbeatsC <- 0
    expLossA <- 0
    expLossB <- 0
    expLossC <- 0
    expLossD <- 0
    lambdaA <- 0
    lambdaB <- 0
    lambdaC <- 0
    lambdaD <- 0
    omegaA <- 0
    omegaB <- 0
    omegaC <- 0
    omegaD <- 0
  } else {
    lambdaA <- rbeta(MSamples, alphaA, betaA)
    lambdaB <- rbeta(MSamples, alphaB, betaB)
    lambdaC <- rbeta(MSamples, alphaC, betaC)
    lambdaD <- rbeta(MSamples, alphaD, betaD)
    omegaA <- rgamma(MSamples, shape = kA, scale = thetaA)
    omegaB <- rgamma(MSamples, shape = kB, scale = thetaB)
    omegaC <- rgamma(MSamples, shape = kC, scale = thetaC)
    omegaD <- rgamma(MSamples, shape = kD, scale = thetaD)
    
    convProbBbeatsA <- sum(lambdaB > lambdaA)/MSamples
    diffTemp <- lambdaB - lambdaA
    convExpLossA_AB <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossB_AB <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbAbeatsB <- sum(lambdaA > lambdaB)/MSamples
    diffTemp <- lambdaA - lambdaB
    convExpLossA_AB2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossB_AB2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbCbeatsA <- sum(lambdaC > lambdaA)/MSamples
    diffTemp <- lambdaC - lambdaA
    convExpLossA_AC <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossC_AC <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbAbeatsC <- sum(lambdaA > lambdaC)/MSamples
    diffTemp <- lambdaA - lambdaC
    convExpLossA_AC2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossC_AC2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbDbeatsA <- sum(lambdaD > lambdaA)/MSamples
    diffTemp <- lambdaD - lambdaA
    convExpLossA_AD <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossD_AD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbAbeatsD <- sum(lambdaA > lambdaD)/MSamples
    diffTemp <- lambdaA - lambdaD
    convExpLossA_AD2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossD_AD2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbCbeatsB <- sum(lambdaC > lambdaB)/MSamples
    diffTemp <- lambdaC - lambdaB
    convExpLossB_BC <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossC_BC <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbBbeatsC <- sum(lambdaB > lambdaC)/MSamples
    diffTemp <- lambdaB - lambdaC
    convExpLossB_BC2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossC_BC2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbDbeatsB <- sum(lambdaD > lambdaB)/MSamples
    diffTemp <- lambdaD - lambdaB
    convExpLossB_BD <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossD_BD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbBbeatsD <- sum(lambdaB > lambdaD)/MSamples
    diffTemp <- lambdaB - lambdaD
    convExpLossB_BD2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossD_BD2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbDbeatsC <- sum(lambdaD > lambdaC)/MSamples
    diffTemp <- lambdaD - lambdaC
    convExpLossC_CD <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossD_CD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbCbeatsD <- sum(lambdaC > lambdaD)/MSamples
    diffTemp <- lambdaC - lambdaD
    convExpLossC_CD2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossD_CD2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbBbeatsA <- sum(1/omegaB > 1/omegaA)/MSamples
    diffTemp <- 1/omegaB - 1/omegaA
    revExpLossA_AB <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossB_AB <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbAbeatsB <- sum(1/omegaA > 1/omegaB)/MSamples
    diffTemp <- 1/omegaA - 1/omegaB
    revExpLossA_AB2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossB_AB2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbCbeatsA <- sum(1/omegaC > 1/omegaA)/MSamples
    diffTemp <- 1/omegaC - 1/omegaA
    revExpLossA_AC <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossC_AC <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbAbeatsC <- sum(1/omegaA > 1/omegaC)/MSamples
    diffTemp <- 1/omegaA - 1/omegaC
    revExpLossA_AC2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossC_AC2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbDbeatsA <- sum(1/omegaD > 1/omegaA)/MSamples
    diffTemp <- 1/omegaD - 1/omegaA
    revExpLossA_AD <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossD_AD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbAbeatsD <- sum(1/omegaA > 1/omegaD)/MSamples
    diffTemp <- 1/omegaA - 1/omegaD
    revExpLossA_AD2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossD_AD2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbCbeatsB <- sum(1/omegaC > 1/omegaB)/MSamples
    diffTemp <- 1/omegaC - 1/omegaB
    revExpLossB_BC <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossC_BC <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbBbeatsC <- sum(1/omegaB > 1/omegaC)/MSamples
    diffTemp <- 1/omegaB - 1/omegaC
    revExpLossB_BC2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossC_BC2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbDbeatsB <- sum(1/omegaD > 1/omegaB)/MSamples
    diffTemp <- 1/omegaD - 1/omegaB
    revExpLossB_BD <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossD_BD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbBbeatsD <- sum(1/omegaB > 1/omegaD)/MSamples
    diffTemp <- 1/omegaB - 1/omegaD
    revExpLossB_BD2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossD_BD2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbDbeatsC <- sum(1/omegaD > 1/omegaC)/MSamples
    diffTemp <- 1/omegaD - 1/omegaC
    revExpLossC_CD <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossD_CD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbCbeatsD <- sum(1/omegaC > 1/omegaD)/MSamples
    diffTemp <- 1/omegaC - 1/omegaD
    revExpLossC_CD2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossD_CD2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbBbeatsA <- sum(lambdaB/omegaB > lambdaA/omegaA)/MSamples
    diffTemp <- lambdaB/omegaB - lambdaA/omegaA
    arpuExpLossA_AB <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossB_AB <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbAbeatsB <- sum(lambdaA/omegaA > lambdaB/omegaB)/MSamples
    diffTemp <- lambdaA/omegaA - lambdaB/omegaB
    arpuExpLossA_AB2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossB_AB2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbCbeatsA <- sum(lambdaC/omegaC > lambdaA/omegaA)/MSamples
    diffTemp <- lambdaC/omegaC - lambdaA/omegaA
    arpuExpLossA_AC <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossC_AC <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbAbeatsC <- sum(lambdaA/omegaA > lambdaC/omegaC)/MSamples
    diffTemp <- lambdaA/omegaA - lambdaC/omegaC
    arpuExpLossA_AC2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossC_AC2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbDbeatsA <- sum(lambdaD/omegaD > lambdaA/omegaA)/MSamples
    diffTemp <- lambdaD/omegaD - lambdaA/omegaA
    arpuExpLossA_AD <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossD_AD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbAbeatsD <- sum(lambdaA/omegaA > lambdaD/omegaD)/MSamples
    diffTemp <- lambdaA/omegaA - lambdaD/omegaD
    arpuExpLossA_AD2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossD_AD2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbCbeatsB <- sum(lambdaC/omegaC > lambdaB/omegaB)/MSamples
    diffTemp <- lambdaC/omegaC - lambdaB/omegaB
    arpuExpLossB_BC <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossC_BC <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbBbeatsC <- sum(lambdaB/omegaB > lambdaC/omegaC)/MSamples
    diffTemp <- lambdaB/omegaB - lambdaC/omegaC
    arpuExpLossB_BC2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossC_BC2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbDbeatsB <- sum(lambdaD/omegaD > lambdaB/omegaB)/MSamples
    diffTemp <- lambdaD/omegaD - lambdaB/omegaB
    arpuExpLossB_BD <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossD_BD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbBbeatsD <- sum(lambdaB/omegaB > lambdaD/omegaD)/MSamples
    diffTemp <- lambdaB/omegaB - lambdaD/omegaD
    arpuExpLossB_BD2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossD_BD2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbDbeatsC <- sum(lambdaD/omegaD > lambdaC/omegaC)/MSamples
    diffTemp <- lambdaD/omegaD - lambdaC/omegaC
    arpuExpLossC_CD <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossD_CD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbCbeatsD <- sum(lambdaC/omegaC > lambdaD/omegaD)/MSamples
    diffTemp <- lambdaC/omegaC - lambdaD/omegaD
    arpuExpLossC_CD2 <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossD_CD2 <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
  }
  list(
    convProbBbeatsA = convProbBbeatsA,
    convProbAbeatsB = convProbAbeatsB,
    convProbCbeatsA = convProbCbeatsA,
    convProbAbeatsC = convProbAbeatsC,
    convProbDbeatsA = convProbDbeatsA,
    convProbAbeatsD = convProbAbeatsD,
    convProbCbeatsB = convProbCbeatsB,
    convProbBbeatsC = convProbBbeatsC,
    convProbDbeatsB = convProbDbeatsB,
    convProbBbeatsD = convProbBbeatsD,
    convProbDbeatsC = convProbDbeatsC,
    convProbCbeatsD = convProbCbeatsD,
    convExpLossA_AB = convExpLossA_AB,
    convExpLossB_AB = convExpLossB_AB,
    convExpLossA_AC = convExpLossA_AC,
    convExpLossC_AC = convExpLossC_AC,
    convExpLossA_AD = convExpLossA_AD,
    convExpLossD_AD = convExpLossD_AD,
    convExpLossB_BC = convExpLossB_BC,
    convExpLossC_BC = convExpLossC_BC,
    convExpLossB_BD = convExpLossB_BD,
    convExpLossD_BD = convExpLossD_BD,
    convExpLossC_CD = convExpLossC_CD,
    convExpLossD_CD = convExpLossD_CD,
    convExpLossA_AB2 = convExpLossA_AB2,
    convExpLossB_AB2 = convExpLossB_AB2,
    convExpLossA_AC2 = convExpLossA_AC2,
    convExpLossC_AC2 = convExpLossC_AC2,
    convExpLossA_AD2 = convExpLossA_AD2,
    convExpLossD_AD2 = convExpLossD_AD2,
    convExpLossB_BC2 = convExpLossB_BC2,
    convExpLossC_BC2 = convExpLossC_BC2,
    convExpLossB_BD2 = convExpLossB_BD2,
    convExpLossD_BD2 = convExpLossD_BD2,
    convExpLossC_CD2 = convExpLossC_CD2,
    convExpLossD_CD2 = convExpLossD_CD2,
    
    revProbBbeatsA = revProbBbeatsA,
    revProbAbeatsB = revProbAbeatsB,
    revProbCbeatsA = revProbCbeatsA,
    revProbAbeatsC = revProbAbeatsC,
    revProbDbeatsA = revProbDbeatsA,
    revProbAbeatsD = revProbAbeatsD,
    revProbCbeatsB = revProbCbeatsB,
    revProbBbeatsC = revProbBbeatsC,
    revProbDbeatsB = revProbDbeatsB,
    revProbBbeatsD = revProbBbeatsD,
    revProbDbeatsC = revProbDbeatsC,
    revProbCbeatsD = revProbCbeatsD,
    revExpLossA_AB = revExpLossA_AB,
    revExpLossB_AB = revExpLossB_AB,
    revExpLossA_AC = revExpLossA_AC,
    revExpLossC_AC = revExpLossC_AC,
    revExpLossA_AD = revExpLossA_AD,
    revExpLossD_AD = revExpLossD_AD,
    revExpLossB_BC = revExpLossB_BC,
    revExpLossC_BC = revExpLossC_BC,
    revExpLossB_BD = revExpLossB_BD,
    revExpLossD_BD = revExpLossD_BD,
    revExpLossC_CD = revExpLossC_CD,
    revExpLossD_CD = revExpLossD_CD,
    revExpLossA_AB2 = revExpLossA_AB2,
    revExpLossB_AB2 = revExpLossB_AB2,
    revExpLossA_AC2 = revExpLossA_AC2,
    revExpLossC_AC2 = revExpLossC_AC2,
    revExpLossA_AD2 = revExpLossA_AD2,
    revExpLossD_AD2 = revExpLossD_AD2,
    revExpLossB_BC2 = revExpLossB_BC2,
    revExpLossC_BC2 = revExpLossC_BC2,
    revExpLossB_BD2 = revExpLossB_BD2,
    revExpLossD_BD2 = revExpLossD_BD2,
    revExpLossC_CD2 = revExpLossC_CD2,
    revExpLossD_CD2 = revExpLossD_CD2,
    
    arpuProbBbeatsA = arpuProbBbeatsA,
    arpuProbAbeatsB = arpuProbAbeatsB,
    arpuProbCbeatsA = arpuProbCbeatsA,
    arpuProbAbeatsC = arpuProbAbeatsC,
    arpuProbDbeatsA = arpuProbDbeatsA,
    arpuProbAbeatsD = arpuProbAbeatsD,
    arpuProbCbeatsB = arpuProbCbeatsB,
    arpuProbBbeatsC = arpuProbBbeatsC,
    
    arpuProbDbeatsB = arpuProbDbeatsB,
    arpuProbBbeatsD = arpuProbBbeatsD,
    arpuProbDbeatsC = arpuProbDbeatsC,
    arpuProbCbeatsD = arpuProbCbeatsD,
    
    
    arpuExpLossA_AB = arpuExpLossA_AB,
    arpuExpLossB_AB = arpuExpLossB_AB,
    arpuExpLossA_AC = arpuExpLossA_AC,
    arpuExpLossC_AC = arpuExpLossC_AC,
    arpuExpLossA_AD = arpuExpLossA_AD,
    arpuExpLossD_AD = arpuExpLossD_AD,
    arpuExpLossB_BC = arpuExpLossB_BC,
    arpuExpLossC_BC = arpuExpLossC_BC,
    arpuExpLossB_BD = arpuExpLossB_BD,
    arpuExpLossD_BD = arpuExpLossD_BD,
    arpuExpLossC_CD = arpuExpLossC_CD,
    arpuExpLossD_CD = arpuExpLossD_CD,
    arpuExpLossA_AB2 = arpuExpLossA_AB2,
    arpuExpLossB_AB2 = arpuExpLossB_AB2,
    arpuExpLossA_AC2 = arpuExpLossA_AC2,
    arpuExpLossC_AC2 = arpuExpLossC_AC2,
    arpuExpLossA_AD2 = arpuExpLossA_AD2,
    arpuExpLossD_AD2 = arpuExpLossD_AD2,
    arpuExpLossB_BC2 = arpuExpLossB_BC2,
    arpuExpLossC_BC2 = arpuExpLossC_BC2,
    arpuExpLossB_BD2 = arpuExpLossB_BD2,
    arpuExpLossD_BD2 = arpuExpLossD_BD2,
    arpuExpLossC_CD2 = arpuExpLossC_CD2,
    arpuExpLossD_CD2 = arpuExpLossD_CD2,
    sampleLambdaA = lambdaA,
    sampleLambdaB = lambdaB,
    sampleLambdaC = lambdaC,
    sampleLambdaD = lambdaD,
    sampleOmegaA = omegaA,
    sampleOmegaB = omegaB,
    sampleOmegaC = omegaC,
    sampleOmegaD = omegaD
  )
}

# Calculates the HDI interval from a sample of representative values
# estimated as shortest credible interval
{
  # Arguments:
  #   sampleVec is a vector of representative values from a probability
  #     distribution.
  #   credMass is a scalar between 0 and 1, indicating the mass within
  #     the credible interval that is to be estimated.
  # Return value:
  #   Highest density interval (HDI) limits in a vector
}
hdi_of_sample <- function(sampleVec, credMass = 0.95) {
  sortedPts <- sort(sampleVec)
  sortedPtsLength <- length(sortedPts)
  if(sortedPtsLength >= 3) {
    ciIdxInc <- min(ceiling(credMass*sortedPtsLength), sortedPtsLength - 1)
    nCIs <- sortedPtsLength - ciIdxInc
    ciWidth <- rep(0, nCIs)
    for (i in 1:nCIs) {
      ciWidth[i] <- sortedPts[i + ciIdxInc] - sortedPts[i]
    }
    HDImin <- sortedPts[which.min(ciWidth)]
    HDImax <- sortedPts[which.min(ciWidth) + ciIdxInc]
    HDIlim <- c(HDImin, HDImax)
  } else {
    HDIlim <- c(min(sortedPts), max(sortedPts))
  }
  return(HDIlim)
}

# Calculates the probability that B is better than A
{
  # Arguments:
  #   Alpha and beta parameters of beta distributions for A and B
  #   alpha = number of successes + 1
  #   beta = number of failures + 1
  # Assumptions:
  #   p_A ~ Beta(alpha_A, beta_A)
  #   p_B ~ Beta(alpha_B, beta_B)
  # Return value:
  #   Pr(p_B > p_A) = probability that the B version is better than the A version
  # Example of use:
  #   If version A converted 180 out of 1000 users and
  #     version B converted 200 out of 950 users type
  #     prob_B_beats_A(180 + 1, 1000 - 180 + 1, 200 + 1, 950 - 200 + 1)
  # Implemented based on Evan Miller`s blog post
  # http://www.evanmiller.org/bayesian-ab-testing.html#binary_ab
}
prob_B_beats_A <- function(alphaA, betaA, alphaB, betaB) {
  if(alphaA <= 0 || betaA <= 0 || alphaB <= 0 || betaB <= 0)
    result <- 0
  else {
    result <- 1
    for (i in 0:(alphaA - 1)) {
      result <- result - 
        exp(
          lbeta(alphaB + i, betaB + betaA) -
            log(betaA + i) -
            lbeta(1 + i, betaA) -
            lbeta(alphaB, betaB)
        )
    }
  }
  result
}

prob_A_beats_B <- function(alphaB, betaB, alphaA, betaA) {
  if(alphaB <= 0 || betaB <= 0 || alphaA <= 0 || betaA <= 0)
    result <- 0
  else {
    result <- 1
    for (i in 0:(alphaB - 1)) {
      result <- result - 
        exp(
          lbeta(alphaA + i, betaA + betaB) -
            log(betaB + i) -
            lbeta(1 + i, betaB) -
            lbeta(alphaA, betaA)
        )
    }
  }
  result
}

prob_C_beats_A <- function(alphaA, betaA, alphaC, betaC) {
  if(alphaA <= 0 || betaA <= 0 || alphaC <= 0 || betaC <= 0)
    result <- 0
  else {
    result <- 1
    for (i in 0:(alphaA - 1)) {
      result <- result - 
        exp(
          lbeta(alphaC + i, betaC + betaA) -
            log(betaA + i) -
            lbeta(1 + i, betaA) -
            lbeta(alphaC, betaC)
        )
    }
  }
  result
}

prob_A_beats_C <- function(alphaC, betaC, alphaA, betaA) {
  if(alphaC <= 0 || betaC <= 0 || alphaA <= 0 || betaA <= 0)
    result <- 0
  else {
    result <- 1
    for (i in 0:(alphaC - 1)) {
      result <- result - 
        exp(
          lbeta(alphaA + i, betaA + betaC) -
            log(betaC + i) -
            lbeta(1 + i, betaC) -
            lbeta(alphaA, betaA)
        )
    }
  }
  result
}

prob_D_beats_A <- function(alphaA, betaA, alphaD, betaD) {
  if(alphaA <= 0 || betaA <= 0 || alphaD <= 0 || betaD <= 0)
    result <- 0
  else {
    result <- 1
    for (i in 0:(alphaA - 1)) {
      result <- result - 
        exp(
          lbeta(alphaD + i, betaD + betaA) -
            log(betaA + i) -
            lbeta(1 + i, betaA) -
            lbeta(alphaD, betaD)
        )
    }
  }
  result
}

prob_A_beats_D <- function(alphaD, betaD, alphaA, betaA) {
  if(alphaD <= 0 || betaD <= 0 || alphaA <= 0 || betaA <= 0)
    result <- 0
  else {
    result <- 1
    for (i in 0:(alphaD - 1)) {
      result <- result - 
        exp(
          lbeta(alphaA + i, betaA + betaD) -
            log(betaD + i) -
            lbeta(1 + i, betaD) -
            lbeta(alphaA, betaA)
        )
    }
  }
  result
}

prob_C_beats_B <- function(alphaB, betaB, alphaC, betaC) {
  if(alphaB <= 0 || betaB <= 0 || alphaC <= 0 || betaC <= 0)
    result <- 0
  else {
    result <- 1
    for (i in 0:(alphaB - 1)) {
      result <- result - 
        exp(
          lbeta(alphaC + i, betaC + betaB) -
            log(betaB + i) -
            lbeta(1 + i, betaB) -
            lbeta(alphaC, betaC)
        )
    }
  }
  result
}

prob_B_beats_C <- function(alphaC, betaC, alphaB, betaB) {
  if(alphaC <= 0 || betaC <= 0 || alphaB <= 0 || betaB <= 0)
    result <- 0
  else {
    result <- 1
    for (i in 0:(alphaC - 1)) {
      result <- result - 
        exp(
          lbeta(alphaB + i, betaB + betaC) -
            log(betaC + i) -
            lbeta(1 + i, betaC) -
            lbeta(alphaB, betaB)
        )
    }
  }
  result
}

prob_D_beats_B <- function(alphaB, betaB, alphaD, betaD) {
  if(alphaB <= 0 || betaB <= 0 || alphaD <= 0 || betaD <= 0)
    result <- 0
  else {
    result <- 1
    for (i in 0:(alphaB - 1)) {
      result <- result - 
        exp(
          lbeta(alphaD + i, betaD + betaB) -
            log(betaB + i) -
            lbeta(1 + i, betaB) -
            lbeta(alphaD, betaD)
        )
    }
  }
  result
}

prob_B_beats_D <- function(alphaD, betaD, alphaB, betaB) {
  if(alphaD <= 0 || betaD <= 0 || alphaB <= 0 || betaB <= 0)
    result <- 0
  else {
    result <- 1
    for (i in 0:(alphaD - 1)) {
      result <- result - 
        exp(
          lbeta(alphaB + i, betaB + betaD) -
            log(betaD + i) -
            lbeta(1 + i, betaD) -
            lbeta(alphaB, betaB)
        )
    }
  }
  result
}

prob_D_beats_C <- function(alphaC, betaC, alphaD, betaD) {
  if(alphaC <= 0 || betaC <= 0 || alphaD <= 0 || betaD <= 0)
    result <- 0
  else {
    result <- 1
    for (i in 0:(alphaC - 1)) {
      result <- result - 
        exp(
          lbeta(alphaD + i, betaD + betaC) -
            log(betaC + i) -
            lbeta(1 + i, betaC) -
            lbeta(alphaD, betaD)
        )
    }
  }
  result
}

prob_C_beats_D <- function(alphaD, betaD, alphaC, betaC) {
  if(alphaD <= 0 || betaD <= 0 || alphaC <= 0 || betaC <= 0)
    result <- 0
  else {
    result <- 1
    for (i in 0:(alphaD - 1)) {
      result <- result - 
        exp(
          lbeta(alphaC + i, betaC + betaD) -
            log(betaD + i) -
            lbeta(1 + i, betaD) -
            lbeta(alphaC, betaC)
        )
    }
  }
  result
}
