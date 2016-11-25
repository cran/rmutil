\name{rmutil}
\title{Utilities for Repeated Measurements Library}
\alias{rmutil}
\description{
\code{\link[rmutil]{\%^\%}}			  Power of a Matrix

\code{\link[rmutil]{covariates}}		Extract Covariate Matrices from a Data Object

\code{\link[rmutil]{covind}}			Nesting Indicator for Observations within Individuals in a Data Object

\code{\link[rmutil]{dbetabinom}}		Density of Beta Binomial Distribution

\code{\link[rmutil]{dboxcox}}			Density of Box-Cox Distribution

\code{\link[rmutil]{dburr}}			Density of Burr Distribution

\code{\link[rmutil]{ddoublebinom}}		Density of Double Binomial Distribution

\code{\link[rmutil]{ddoublepois}}		Density of Double Poisson Distribution

\code{\link[rmutil]{delta}}			Extract Units of Measurement Vector from a Data Object

\code{\link[rmutil]{dftorep}}			Transform a Dataframe to a repeated Object

\code{\link[rmutil]{dgammacount}}		Density of Gamma Count Distribution

\code{\link[rmutil]{dgextval}}			Density of Generalized Extreme Value Distribution

\code{\link[rmutil]{dggamma}}			Density of Generalized Gamma Distribution

\code{\link[rmutil]{dginvgauss}}		Density of Generalized Inverse Gaussian Distribution

\code{\link[rmutil]{dglogis}}			Density of Generalized Logistic Distribution

\code{\link[rmutil]{dgweibull}}			Density of Generalized Weibull Distribution

\code{\link[rmutil]{dhjorth}}			Density of Hjorth Distribution

\code{\link[rmutil]{dinvgauss}}			Density of Inverse Gaussian Distribution

\code{\link[rmutil]{dlaplace}}			Density of Laplace Distribution

\code{\link[rmutil]{dlevy}}			Density of Levy Distribution

\code{\link[rmutil]{dmultbinom}}		Density of Multiplicative Binomial Distribution

\code{\link[rmutil]{dmultpois}}			Density of Multiplicative Poisson Distribution

\code{\link[rmutil]{dpareto}}			Density of Pareto Distribution

\code{\link[rmutil]{dpowexp}}			Density of Power Exponential Distribution

\code{\link[rmutil]{dpvfpois}}			Density of Power Variance Function Poisson Distribution

\code{\link[rmutil]{dsimplex}}			Density of Simplex Distribution

\code{\link[rmutil]{dskewlaplace}}		Density of Skew Laplace Distribution

\code{\link[rmutil]{finterp}}			Formula Interpreter

\code{\link[rmutil]{fmobj}}			Object Finder in Formulae

\code{\link[rmutil]{fnenvir}}			Check Covariates and Parameters of a Function

\code{\link[rmutil]{formula}}			Extract Formula Used to Create Time-constant Covariate Matrix in a Data Object

\code{\link[rmutil]{gauss.hermite}}		Calculate Gauss-Hermite Quadrature Points

\code{\link[rmutil]{gettvc}}			Create Time-varying Covariates

\code{\link[rmutil]{int}}			Vectorized One-dimensional Numerical Integration

\code{\link[rmutil]{int2}}			Vectorized Two-dimensional Numerical Integration

\code{\link[rmutil]{iprofile}}			Produce Individual Time Profiles for Plotting

\code{\link[rmutil]{lin.diff.eqn}}		Solution of Autonomous Linear Differential Equations

\code{\link[rmutil]{lvna}}			Create a Repeated Object Leaving NAs

\code{\link[rmutil]{mexp}}			Matrix Exponentiation

\code{\link[rmutil]{mprofile}}			Produce Marginal Time Profiles for Plotting

\code{\link[rmutil]{names}}			Extract Names of Covariates from a Data Object

\code{\link[rmutil]{nesting}}			Extract Nesting Indicators from a Data Object

\code{\link[rmutil]{nobs}}			Extract Number of Observations per Individual from a Data Object

\code{\link[rmutil]{pbetabinom}}		Distribution Function of Beta Binomial Distribution

\code{\link[rmutil]{pboxcox}}			Distribution Function of Box-Cox Distribution

\code{\link[rmutil]{pburr}}			Distribution Function of Burr Distribution

\code{\link[rmutil]{pdoublebinom}}		Distribution Function of Double Binomial Distribution

\code{\link[rmutil]{pdoublepois}}		Distribution Function of Double Poisson Distribution

\code{\link[rmutil]{pgammacount}}		Distribution Function of Gamma Count Distribution

\code{\link[rmutil]{pgextval}}			Distribution Function of Generalized Extreme Value Distribution

\code{\link[rmutil]{pggamma}}			Distribution Function of Generalized Gamma Distribution

\code{\link[rmutil]{pginvgauss}}		Distribution Function of Generalized Inverse Gaussian Distribution

\code{\link[rmutil]{pglogis}}			Distribution Function of Generalized Logistic Distribution

\code{\link[rmutil]{pgweibull}}			Distribution Function of Generalized Weibull Distribution

\code{\link[rmutil]{phjorth}}			Distribution Function of Hjorth Distribution

\code{\link[rmutil]{pinvgauss}}			Distribution Function of Inverse Gaussian Distribution

\code{\link[rmutil]{pkpd}}			Pharmacokinetic Model Functions

\code{\link[rmutil]{plaplace}}			Distribution Function of Laplace Distribution

\code{\link[rmutil]{plevy}}			Distribution Function of Levy Distribution

\code{\link[rmutil]{plot.residuals}}		Plot Residuals for Carma

\code{\link[rmutil]{pmultbinom}}		Distribution Function of Multiplicative Binomial Distribution

\code{\link[rmutil]{pmultpois}}			Distribution Function of Multiplicative Poisson Distribution

\code{\link[rmutil]{ppareto}}			Distribution Function of Pareto Distribution

\code{\link[rmutil]{ppowexp}}			Distribution Function of Power Exponential Distribution

\code{\link[rmutil]{ppvfpois}}			Distribution Function of Power Variance Function Poisson Distribution

\code{\link[rmutil]{psimplex}}			Distribution Function of Simplex Distribution

\code{\link[rmutil]{pskewlaplace}}		Distribution Function of Skew Laplace Distribution

\code{\link[rmutil]{qbetabinom}}		Quantiles of Beta Binomial Distribution

\code{\link[rmutil]{qboxcox}}			Quantiles of Box-Cox Distribution

\code{\link[rmutil]{qburr}}			Quantiles of Burr Distribution

\code{\link[rmutil]{qdoublebinom}}		Quantiles of Double Binomial Distribution

\code{\link[rmutil]{qdoublepois}}		Quantiles of Double Poisson Distribution

\code{\link[rmutil]{qgammacount}}		Quantiles of Gamma Count Distribution

\code{\link[rmutil]{qgextval}}			Quantiles of Generalized Extreme Value Distribution

\code{\link[rmutil]{qggamma}}			Quantiles of Generalized Gamma Distribution

\code{\link[rmutil]{qginvgauss}}		Quantiles of Generalized Inverse Gaussian Distribution

\code{\link[rmutil]{qglogis}}			Quantiles of Generalized Logistic Distribution

\code{\link[rmutil]{qgweibull}}			Quantiles of Generalized Weibull Distribution

\code{\link[rmutil]{qhjorth}}			Quantiles of Hjorth Distribution

\code{\link[rmutil]{qinvgauss}}			Quantiles of Inverse Gaussian Distribution

\code{\link[rmutil]{qlaplace}}			Quantiles of Laplace Distribution

\code{\link[rmutil]{qlevy}}			Quantiles of Levy Distribution

\code{\link[rmutil]{qmultbinom}}		Quantiles of Multiplicative Binomial Distribution

\code{\link[rmutil]{qmultpois}}			Quantiles of Multiplicative Poisson Distribution

\code{\link[rmutil]{qpareto}}			Quantiles of Pareto Distribution

\code{\link[rmutil]{qpowexp}}			Quantiles of Power Exponential Distribution

\code{\link[rmutil]{qpvfpois}}			Quantiles of Power Variance Function Poisson Distribution

\code{\link[rmutil]{qsimplex}}			Quantiles of Simplex Distribution

\code{\link[rmutil]{qskewlaplace}}		Quantiles of Skew Laplace Distribution

\code{\link[rmutil]{rbetabinom}}		Random Number Generation for Beta Binomial Distribution

\code{\link[rmutil]{rboxcox}}			Random Number Generation for Box-Cox Distribution

\code{\link[rmutil]{rburr}}			Random Number Generation for Burr Distribution

\code{\link[rmutil]{rdoublebinom}}		Random Number Generation for Double Binomial Distribution

\code{\link[rmutil]{rdoublepois}}		Random Number Generation for Double Poisson Distribution

\code{\link[rmutil]{read.list}}			Read a List of Matrices of Unbalanced Repeated Measurements from a File

\code{\link[rmutil]{read.rep}}			Read a Rectangular Data Set from a File to Create a repeated Object

\code{\link[rmutil]{read.surv}}			Read a List of Vectors of Event Histories from a File

\code{\link[rmutil]{response}}			Extract Response Vector from a Data Object

\code{\link[rmutil]{restovec}}			Create a Response Object

\code{\link[rmutil]{rgammacount}}		Random Number Generation for Gamma Count Distribution

\code{\link[rmutil]{rgextval}}			Random Number Generation for Generalized Extreme Value Distribution

\code{\link[rmutil]{rggamma}}			Random Number Generation for Generalized Gamma Distribution

\code{\link[rmutil]{rginvgauss}}		Random Number Generation for Generalized Inverse Gaussian Distribution

\code{\link[rmutil]{rglogis}}			Random Number Generation for Generalized Logistic Distribution

\code{\link[rmutil]{rgweibull}}			Random Number Generation for Generalized Weibull Distribution

\code{\link[rmutil]{rhjorth}}			Random Number Generation for Hjorth Distribution

\code{\link[rmutil]{rinvgauss}}			Random Number Generation for Inverse Gaussian Distribution

\code{\link[rmutil]{rlaplace}}			Random Number Generation for Laplace Distribution

\code{\link[rmutil]{rlevy}}			Random Number Generation for Levy Distribution

\code{\link[rmutil]{rmna}}			Create a Repeated Object

\code{\link[rmutil]{rmultbinom}}		Random Number Generation for Multiplicative Binomial Distribution

\code{\link[rmutil]{rmultpois}}			Random Number Generation for Multiplicative Poisson Distribution

\code{\link[rmutil]{rpareto}}			Random Number Generation for Pareto Distribution

\code{\link[rmutil]{rpowexp}}			Random Number Generation for Power Exponential Distribution

\code{\link[rmutil]{rpvfpois}}			Random Number Generation for Power Variance Function Poisson Distribution

\code{\link[rmutil]{rsimplex}}			Random Number Generation for Simplex Distribution

\code{\link[rmutil]{rskewlaplace}}		Random Number Generation for Skew Laplace Distribution

\code{\link[rmutil]{runge.kutta}}		Runge-Kutta Method for Solving Differential Equations

\code{\link[rmutil]{tcctomat}}			Create a Time-constant Covariate (tccov) Object

\code{\link[rmutil]{times}}			Extract Times Vector from a Data Object

\code{\link[rmutil]{transform}}			Transform Variables in a Data Object

\code{\link[rmutil]{tvctomat}}			Create a Time-varying Covariate (tvcov) Object

\code{\link[rmutil]{wr}}			Find the Response Vector and Design Matrix for a Model Formula
}
\keyword{documentation}
