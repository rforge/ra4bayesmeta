\name{aom}
\alias{aom}
\docType{data}
\title{
Acute otitis media data
%%   ~~ data name/kind ... ~~
}
\description{
Meta-analysis data set including 5 studies on treatment failure of short course (less
than 7 days) versus long course (> 7 days) antibiotics for acute otitis media in children.
The outcome considered is treatment failure at 8-19 days.

%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(aom)}
\format{
  A data frame with 5 observations on the following 7 variables:
  \describe{
    \item{\code{study}}{character string, label of the study}
    %TODO: specify the direction of the OR
    \item{\code{f.t}}{integer, number of patients for whom the short course antibiotics treatment failed}
    \item{\code{n.t}}{integer, number of patients receiving short course antibiotics}
    \item{\code{f.c}}{integer, number of patients for whom the long course antibiotics treatment failed}
    \item{\code{n.c}}{integer, number of patients receiving long course antibiotics}
    \item{\code{y}}{numeric, log odds ratio estimates for the individual studies}
    \item{\code{sigma}}{numeric, standard error of the log odds ratio estimate for the individual studies}
  }
}
\details{
This data set has been analyzed in Kozyrskyj et al. (2000) and 
reanalyzed by Lambert et al. (2005).
Senn (2007, Section 2) discusses problems in this data set.
In particular, one arm of one study has been included twice in the data set.

The counts \code{f.t}, \code{n.t}, \code{f.c} and \code{n.c} and the
names of the studies are taken from Senn (2007),
who reproduces the data set given in Kozyrskyj et al. (2000, Analysis 2.2).
The log odds ratio estimates and standard errors were computed using the \code{escalc} function in the
package \code{metafor} and 
are identical to the estimates used in Lambert et al. (2005).
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%Lambert, P., Sutton, A., Burton, P., Abrams, K., Jones, D. (2005). How vague is
%vague? A simulation study of the impact of the use of vague prior distributions in
%MCMC using WinBUGS. \emph{Statistics in Medicine} \bold{24(15)}, 2401--2428.
%\url{https://doi.org/10.1002/sim.2112}
Senn, S. (2007). Trying to be precise about vagueness. \emph{Statistics in Medicine}, \bold{26}, 1417--1430. \url{https://doi.org/10.1002/sim.2639}
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Kozyrskyj, A., Klassen, T. P., Moffatt, M., Harvey K. (2000).
Short-course antibiotics for acute otitis media.
\emph{Cochrane Database of Systematic Reviews}, Issue 2, Art. No.: CD001095.
\url{https://doi.org/10.1002/14651858.CD001095}

Lambert, P., Sutton, A., Burton, P., Abrams, K., Jones, D. (2005). How vague is
vague? A simulation study of the impact of the use of vague prior distributions in
MCMC using WinBUGS. \emph{Statistics in Medicine} \bold{24(15)}, 2401--2428.
\url{https://doi.org/10.1002/sim.2112}
%Senn, S. (2007). Trying to be precise about vagueness. \emph{Statistics in Medicine}, \bold{26}, 1417--1430. \url{https://doi.org/10.1002/sim.2639}
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(aom)
str(aom)
# forest plot
forest(x=aom$y, sei=aom$sigma, xlab="log odds ratio")
}
\keyword{datasets}
