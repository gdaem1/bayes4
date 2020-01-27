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
    
    convProbCbeatsA <- sum(lambdaC > lambdaA)/MSamples
    diffTemp <- lambdaC - lambdaA
    convExpLossA_AC <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossC_AC <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbDbeatsA <- sum(lambdaD > lambdaA)/MSamples
    diffTemp <- lambdaD - lambdaA
    convExpLossA_AD <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossD_AD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbCbeatsB <- sum(lambdaC > lambdaB)/MSamples
    diffTemp <- lambdaC - lambdaB
    convExpLossB_BC <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossC_BC <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbDbeatsB <- sum(lambdaD > lambdaB)/MSamples
    diffTemp <- lambdaD - lambdaB
    convExpLossB_BD <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossD_BD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    convProbDbeatsC <- sum(lambdaD > lambdaC)/MSamples
    diffTemp <- lambdaD - lambdaC
    convExpLossC_CD <- sum(diffTemp*(diffTemp > 0))/MSamples
    convExpLossD_CD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbBbeatsA <- sum(1/omegaB > 1/omegaA)/MSamples
    diffTemp <- 1/omegaB - 1/omegaA
    revExpLossA_AB <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossB_AB <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbCbeatsA <- sum(1/omegaC > 1/omegaA)/MSamples
    diffTemp <- 1/omegaC - 1/omegaA
    revExpLossA_AC <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossC_AC <- sum(-diffTemp*(-diffTemp > 0))/MSamples
  
    revProbDbeatsA <- sum(1/omegaD > 1/omegaA)/MSamples
    diffTemp <- 1/omegaD - 1/omegaA
    revExpLossA_AD <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossD_AD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
  
    revProbCbeatsB <- sum(1/omegaC > 1/omegaB)/MSamples
    diffTemp <- 1/omegaC - 1/omegaB
    revExpLossB_BC <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossC_BC <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbDbeatsB <- sum(1/omegaD > 1/omegaB)/MSamples
    diffTemp <- 1/omegaD - 1/omegaB
    revExpLossB_BD <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossD_BD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    revProbDbeatsC <- sum(1/omegaD > 1/omegaC)/MSamples
    diffTemp <- 1/omegaD - 1/omegaC
    revExpLossC_CD <- sum(diffTemp*(diffTemp > 0))/MSamples
    revExpLossD_CD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbBbeatsA <- sum(lambdaB/omegaB > lambdaA/omegaA)/MSamples
    diffTemp <- lambdaB/omegaB - lambdaA/omegaA
    arpuExpLossA_AB <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossB_AB <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbCbeatsA <- sum(lambdaC/omegaC > lambdaA/omegaA)/MSamples
    diffTemp <- lambdaC/omegaC - lambdaA/omegaA
    arpuExpLossA_AC <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossC_AC <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbDbeatsA <- sum(lambdaD/omegaD > lambdaA/omegaA)/MSamples
    diffTemp <- lambdaD/omegaD - lambdaA/omegaA
    arpuExpLossA_AD <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossD_AD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbCbeatsB <- sum(lambdaC/omegaC > lambdaB/omegaB)/MSamples
    diffTemp <- lambdaC/omegaC - lambdaB/omegaB
    arpuExpLossB_BC <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossC_BC <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbDbeatsB <- sum(lambdaD/omegaD > lambdaB/omegaB)/MSamples
    diffTemp <- lambdaD/omegaD - lambdaB/omegaB
    arpuExpLossB_BD <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossD_BD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
    arpuProbDbeatsC <- sum(lambdaD/omegaD > lambdaC/omegaC)/MSamples
    diffTemp <- lambdaD/omegaD - lambdaC/omegaC
    arpuExpLossC_CD <- sum(diffTemp*(diffTemp > 0))/MSamples
    arpuExpLossD_CD <- sum(-diffTemp*(-diffTemp > 0))/MSamples
    
  }
  list(
    convProbBbeatsA = convProbBbeatsA,
    convProbCbeatsA = convProbCbeatsA,
    convProbDbeatsA = convProbDbeatsA,
    convProbCbeatsB = convProbCbeatsB,
    convProbDbeatsB = convProbDbeatsB,
    convProbDbeatsC = convProbDbeatsC,
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
    revProbBbeatsA = revProbBbeatsA,
    revProbCbeatsA = revProbCbeatsA,
    revProbDbeatsA = revProbDbeatsA,
    revProbCbeatsB = revProbCbeatsB,
    revProbDbeatsB = revProbDbeatsB,
    revProbDbeatsC = revProbDbeatsC,
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
    arpuProbBbeatsA = arpuProbBbeatsA,
    arpuProbCbeatsA = arpuProbCbeatsA,
    arpuProbDbeatsA = arpuProbDbeatsA,
    arpuProbCbeatsB = arpuProbCbeatsB,
    arpuProbDbeatsB = arpuProbDbeatsB,
    arpuProbDbeatsC = arpuProbDbeatsC,
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
  #   Highest density iterval (HDI) limits in a vector
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
  if(alphaC <= 0 || beta <= 0 || alphaD <= 0 || betaD <= 0)
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

# Creates plotly chart of approximate density based on a sample data
plot_sample_density <- function(
  sampleDataA, hdiA = c(-0.1, 0.1),
  sampleDataB, hdiB = c(-0.1, 0.1),
  sampleDataC, hdiC = c(-0.1, 0.1),
  printPlot = TRUE, ...
) {
  if(printPlot) {
    densA <- density(sampleDataA, bw = 'nrd')
    plotDataA <- data.frame(x = densA$x, y = densA$y)
    xA <- plotDataA$x
    yA <- plotDataA$y
    xAreaA <- plotDataA[plotDataA$x >= hdiA[1] & plotDataA$x <= hdiA[2], ]$x
    yAreaA <- plotDataA[plotDataA$x >= hdiA[1] & plotDataA$x <= hdiA[2], ]$y
    densB <- density(sampleDataB, bw = 'nrd')
    plotDataB <- data.frame(x = densB$x, y = densB$y)
    xB <- plotDataB$x
    yB <- plotDataB$y
    xAreaB <- plotDataB[plotDataB$x >= hdiB[1] & plotDataB$x <= hdiB[2], ]$x
    yAreaB <- plotDataB[plotDataB$x >= hdiB[1] & plotDataB$x <= hdiB[2], ]$y
    densC <- density(sampleDataC, bw = 'nrd')
    plotDataC <- data.frame(x = densC$x, y = densC$y)
    xC <- plotDataC$x
    yC <- plotDataC$y
    xAreaC <- plotDataC[plotDataC$x >= hdiC[1] & plotDataC$x <= hdiC[2], ]$x
    yAreaC <- plotDataC[plotDataC$x >= hdiC[1] & plotDataC$x <= hdiC[2], ]$y
    plot_ly(
      x = xA,
      y = yA,
      type = 'scatter',
      mode = 'lines',
      name = 'ARPU B - A',
      line = list(color = '#008B00'),
      hoverinfo = 'none',
      # text = ~paste(sprintf('%.3g%%', x*100)),
      showlegend = FALSE
    ) %>% add_trace(
      x = xAreaA,
      y = yAreaA,
      type = 'scatter',
      mode = 'lines',
      fill = 'tozeroy',
      fillcolor = 'rgba(0, 139, 0, 0.2)',
      line = list(color = '#008B00'),
      hoverinfo = 'none',
      # text = ~paste(sprintf('%.3g%%', xAreaA*100)),
      showlegend = FALSE
    ) %>% layout(
      x = xB,
      y = yB,
      type = 'scatter',
      mode = 'lines',
      name = 'ARPU C - A',
      line = list(color = '##ff2f00'),
      hoverinfo = 'none',
      # text = ~paste(sprintf('%.3g%%', x*100)),
      showlegend = FALSE
    ) %>% add_trace(
      x = xAreaB,
      y = yAreaB,
      type = 'scatter',
      mode = 'lines',
      fill = 'tozeroy',
      fillcolor = 'rgba(255, 47, 0, 0.2)',
      line = list(color = '#ff2f00'),
      hoverinfo = 'none',
      # text = ~paste(sprintf('%.3g%%', xAreaB*100)),
      showlegend = FALSE
    ) %>% layout(
      x = xC,
      y = yC,
      type = 'scatter',
      mode = 'lines',
      name = 'ARPU C - B',
      line = list(color = '#4493c9'),
      hoverinfo = 'none',
      # text = ~paste(sprintf('%.3g%%', x*100)),
      showlegend = FALSE
    ) %>% add_trace(
      x = xAreaC,
      y = yAreaC,
      type = 'scatter',
      mode = 'lines',
      fill = 'tozeroy',
      fillcolor = 'rgba(68, 147, 201, 0.2)',
      line = list(color = '#4493c9'),
      hoverinfo = 'none',
      # text = ~paste(sprintf('%.3g%%', xArea*100)),
      showlegend = FALSE
    ) %>% layout(
      title = 'Approximate distribution of difference of ARPU A, B, and C',
      title = list(fontsize = 14),
      xaxis = list(showgrid = FALSE, fixedrange = TRUE),
      yaxis = list(showline = FALSE, showticklabels = FALSE, showgrid = FALSE, fixedrange = TRUE),
      legend = list(x = 0.9, y = 0.9)
    ) %>% config(displayModeBar = FALSE)
    
  } else {
    plot_ly(type = 'scatter', mode = 'markers')
  }
}

# Creates plotly chart of approximate densities based on a sample data
plot_sample_densities <- function(
  sampleDataA, hdiA = c(-0.1, 0.1),
  sampleDataB, hdiB = c(-0.1, 0.1),
  sampleDataC, hdiC = c(-0.1, 0.1),
  printPlot = TRUE, ...
) {
  if(printPlot) {
    densA <- density(sampleDataA, bw = 'nrd')
    plotDataA <- data.frame(x = densA$x, y = densA$y)
    xA <- plotDataA$x
    yA <- plotDataA$y
    xAreaA <- plotDataA[plotDataA$x >= hdiA[1] & plotDataA$x <= hdiA[2], ]$x
    yAreaA <- plotDataA[plotDataA$x >= hdiA[1] & plotDataA$x <= hdiA[2], ]$y
    densB <- density(sampleDataB, bw = 'nrd')
    plotDataB <- data.frame(x = densB$x, y = densB$y)
    xB <- plotDataB$x
    yB <- plotDataB$y
    xAreaB <- plotDataB[plotDataB$x >= hdiB[1] & plotDataB$x <= hdiB[2], ]$x
    yAreaB <- plotDataB[plotDataB$x >= hdiB[1] & plotDataB$x <= hdiB[2], ]$y
    densC <- density(sampleDataC, bw = 'nrd')
    plotDataC <- data.frame(x = densC$x, y = densC$y)
    xC <- plotDataC$x
    yC <- plotDataC$y
    xAreaC <- plotDataC[plotDataC$x >= hdiC[1] & plotDataC$x <= hdiC[2], ]$x
    yAreaC <- plotDataC[plotDataC$x >= hdiC[1] & plotDataC$x <= hdiC[2], ]$y
    plot_ly(
      x = xA,
      y = yA,
      type = 'scatter',
      mode = 'lines',
      name = 'Sample A',
      line = list(color = '#008B00'),
      hoverinfo = 'none',
      showlegend = TRUE
    ) %>% add_trace(
      x = xAreaA,
      y = yAreaA,
      type = 'scatter',
      mode = 'lines',
      fill = 'tozeroy',
      fillcolor = 'rgba(0, 139, 0, 0.2)',
      line = list(color = '#008B00'),
      hoverinfo = 'none',
      showlegend = FALSE
    ) %>% add_trace(
      x = xB,
      y = yB,
      type = 'scatter',
      mode = 'lines',
      name = 'Sample B',
      line = list(color = '#ff2f00'),
      hoverinfo = 'none',
      showlegend = TRUE
    ) %>% add_trace(
      x = xAreaB,
      y = yAreaB,
      type = 'scatter',
      mode = 'lines',
      name = 'Sample B',
      fill = 'tozeroy',
      fillcolor = 'rgba(255, 47, 0, 0.2)',
      line = list(color = '#ff2f00'),
      hoverinfo = 'none',
      showlegend = FALSE
    ) %>% layout(
      x = xC,
      y = yC,
      type = 'scatter',
      mode = 'lines',
      name = 'Sample C',
      line = list(color = '#4493c9'),
      hoverinfo = 'none',
      showlegend = FALSE
    ) %>% add_trace(
      x = xAreaC,
      y = yAreaC,
      type = 'scatter',
      mode = 'lines',
      name = 'Sample C',
      fill = 'tozeroy',
      fillcolor = 'rgba(68, 147, 201, 0.2)',
      line = list(color = '#4493c9'),
      hoverinfo = 'none',
      showlegend = TRUE
    ) %>% layout(
      title = 'Distribution of ARPU A, B, and C',
      title = list(fontsize = 14),
      xaxis = list(showgrid = FALSE, fixedrange = TRUE),
      yaxis = list(showline = FALSE, showticklabels = FALSE, showgrid = FALSE, fixedrange = TRUE),
      legend = list(x = 0.9, y = 0.9),
      showlegend = TRUE
    ) %>% config(displayModeBar = FALSE)
  } else {
    plot_ly(type = 'scatter', mode = 'markers')
  }
}
