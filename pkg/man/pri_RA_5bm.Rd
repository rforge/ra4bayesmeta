\name{pri_RA_5bm}
\alias{pri_RA_5bm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior reference analysis based on a data frame using 5 benchmarks
%%  ~~function to do ... ~~
}
\description{
Computes a table of Hellinger distance estimates between the
actual heterogeneity prior(s) specified and
four benchmark heterogeneity priors
proposed in the Supplementary Material of Ott et al. (2021).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pri_RA_5bm(df, tau.prior=list(function(x) dhalfnormal(x, scale=1)),
           m_J=NA, M_J=NA, upper.J=3, digits.J=2,
           m_inf=NA, M_inf=NA, rlmc0=0.0001, rlmc1=0.9999,
           mu.mean=0, mu.sd=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df}{
    data frame with one column "y" containing the (transformed) effect estimates for the individual studies and
    one column "sigma" containing the standard errors of these estimates. 
%%     ~~Describe \code{x} here~~
}
\item{tau.prior}{
    list of prior specifications, which are \emph{either} functions returning the probability densities of the actual priors of interest for the heterogeneity parameter tau \emph{or} character strings specifying priors implemented in the \code{bayesmeta} function. See the documentation of the argument \code{tau.prior} of the \code{bayesmeta} function for details.
}
%\item{tau.prior}{
%    list of functions returning the probability densities of the actual priors of interest for the heterogeneity parameter tau.
%}
\item{m_J}{
    parameter value \eqn{m=m_J} of the SGC(\eqn{m}) prior,
    which induces a marignal posterior for the heterogeneity standard deviation tau   close to Jeffreys' reference posterior (wrt the Hellinger distance).
    If set to \code{NA} (the default), this parameter is computed 
    using the function \code{m_j_sgc} and some other parameters.
}
\item{M_J}{
    parameter value \eqn{M=M_J} of the SIGC(\eqn{M}) prior,
    which induces a marignal posterior for the heterogeneity standard deviation tau   close to Jeffreys' reference posterior (wrt the Hellinger distance).
    If set to \code{NA} (the default), this parameter is computed 
    using the function \code{M_j_sigc} and some other parameters.
}
\item{upper.J}{
    upper bound for the parameters \code{M_J} and \code{m_J}. Real number in \eqn{(1,\infty)}. Is required only if \code{M_J=NA} or \code{m_J=NA}.
}
\item{digits.J}{
    specifies the desired precision of the parameter values \eqn{M_J} and \eqn{m_J}, i.e. to how many digits these two values
   should be determined. Possible values are 1,2,3. Defaults to 2.
   Is required only if \code{M_J=NA} or \code{m_J=NA}.
}
\item{m_inf}{
    parameter value \eqn{m=m_{inf}} of the SGC(\eqn{m}) prior,
such that the median relative latent model complexity (RLMC) is close to 0.
  If set to \code{NA} (the default), this parameter is computed 
    using the function \code{m_inf_sgc}, such that the median RLMC is
    approximately equal to \code{rlmc0}.
}
\item{M_inf}{
    parameter value \eqn{M=M_{inf}} of the SIGC(\eqn{M}) prior,
such that the median relative latent model complexity (RLMC) is close to 1.
  If set to \code{NA} (the default), this parameter is computed 
    using the function \code{M_inf_sigc}, such that the median RLMC is
    approximately equal to \code{rlmc1}.
}
\item{rlmc0}{
    RLMC target value for the SGC(\eqn{m_{inf}}) benchmark prior (typically close to 0).
    Is required only if \code{m_inf=NA}.
}
\item{rlmc1}{
    RLMC target value for the SIGC(\eqn{M_{inf}}) benchmark prior (typically close to 1).
    Is required only if \code{M_inf=NA}.
}
\item{mu.mean}{
  mean of the normal prior for the effect mu.
}
\item{mu.sd}{
  standard deviation of the normal prior for the effect mu.
}
}
\details{
The methodology for a prior reference analysis and the four proper heterogeneity benchmark priors used are
introduced in the Supplementary Material of Ott et al. (2021, Sections 2.5 and 2.6).
There, these four benchmark priors are denoted by 
SGC(\eqn{m_J}), SIGC(\eqn{M_J}), SGC(\eqn{m_{inf}}) and SIGC(\eqn{M_{inf}}).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with two elements:
The first element named "table" is a matrix containing the Hellinger distance estimates and the
second element called "par" is a named vector giving the parameter values
of the benchmark priors.
The vector "par" has the following five components:
\code{m_inf}, \code{M_J}, \code{m_J}, \code{M_inf} and \code{C}.

The matrix "table" contains the Hellinger distance estimates between actual and benchmark
heterogeneity priors
and has \eqn{4} columns and \eqn{n} rows,
where \eqn{n} is the number of actual heterogeneity priors specified,
i.e. the length of the list \code{tau.prior}.
 
The columns of the matrix give the following Hellinger distance estimates
between two heterogeneity priors
(from left to right):
\item{H(SGC(m_inf), pri_act) }{ benchmark prior SGC(m_inf) and actual prior}
\item{H(SIGC(M_J), pri_act) }{ benchmark prior SIGC(M_J) and actual prior}
\item{H(SGC(m_J), pri_act) }{ benchmark prior SGC(m_J) and actual prior}
\item{H(SIGC(M_inf), pri_act) }{ benchmark prior SIGC(M_inf) and actual prior}
Each row corresponds to one actual heterogeneity prior
specified in the list \code{tau.prior}, in the same order
as in that list. Thus, the row names are:
\item{pri_act_1 }{ first actual prior in \code{tau.prior}}
\item{pri_act_2 }{ second actual prior in \code{tau.prior}}
\item{...}{}
\item{pri_act_n }{ nth (last) actual prior in \code{tau.prior}}
}
\references{
Ott, M., Plummer, M., Roos, M. Supplementary Material:
How vague is vague? How informative is informative? Reference analysis for 
Bayesian meta-analysis. Revised for \emph{Statistics in Medicine}. 2021.
}
\section{Warning }{
This function may take several minutes to run if the parameter \code{m_J} and/or \code{M_J}
is not specified,
especially if the desired precision is \code{digits.J=2} or even \code{digits.J=3}.
}
\seealso{
\code{\link{post_RA_3bm}},
\code{\link{pri_RA_fits}}
}
\examples{
# for aurigular acupuncture (AA) data set with one
# actual half-normal and the "DuMouchel" heterogeneity prior
data(aa)
pri_RA_5bm(df=aa, tau.prior=list(function(t)dhalfnormal(t, scale=1),
                                 "DuMouchel"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.