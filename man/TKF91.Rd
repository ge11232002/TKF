\name{TKF91}
\alias{TKF91}
\alias{TKF91Pair}

\title{
  Evolutionary distance estimation with TKF91 model
}
\description{
  This function implements the TKF91 model to estimate the pairwise distance
  from protein sequences.
}
\usage{
TKF91(fasta, mu=NULL, 
      method=c("gsl", "nlopt", "NM", "Sbplx", "COBYLA","BOBYQA", "PRAXIS"),
      expectedLength=362, 
      substModel, substModelBF, skipFailure=FALSE)
TKF91Pair(seq1, seq2, mu=NULL, distance=NULL,
          method=c("gsl", "nlopt", "NM", "Sbplx", "COBYLA","BOBYQA", "PRAXIS"),
          expectedLength=362, substModel, substModelBF)
}

\arguments{
  \item{fasta}{
    A named list of sequences in vector of characters format.
  \code{read.fasta} from package \code{seqinr} outputs this format 
  when reading from a fasta file.
  }
  \item{mu}{
    A \link{numeric} value between 0 and 1 or NULL.  
    It is the death rate per normal link in TKF91 model.  
    When it is NULL, a joint estimation of \code{mu} and \code{r} will be done.
    When it is given, only the distance will be estimated.
  }
  \item{distance}{
    A \link{numeric} value: the PAM distance between two protein sequences.
    When it is given, TKF91Pair only calculates the negative log-likelihood.
  }
  \item{method}{
    When mu and distance are co-estimated,
    the optimisation method can be one of 
    "gsl", "nlopt", "NM", "Sbplx", "COBYLA","BOBYQA", "PRAXIS".
    When the method is "gsl", 
    the implementation of "nmsimplex2" from gsl is used.
    When it is "nlopt", it will choose the best algorithm from nlopt package 
    based on the minimal negative log-likelihood. 
    
    This argument is ignored when only distance is being estimated.
  }
  \item{expectedLength}{
    A \link{numeric} object: the expected length of input protein sequences.
    By default, the average sequence length, 362, from OMA browser is used.
  }
  \item{substModel}{
    A \link{numeric} \link{matrix}: the mutation probability from one AA
    to another AA at PAM distance 1. 
    The order of AA in the matrix should be identical to \link{AACharacterSet}.
  }
  \item{substModelBF}{
    A vector of \link{numeric}: the backrgound frequency of AAs.
    The order of AA in the vector should also be identical to 
    \link{AACharacterSet}.
  }
  \item{seq1, seq2}{
    A vector of \link{character}: the sequences of two proteins to compare.
  }
  \item{skipFailure}{
    If TRUE, it will skip the failed optimisation of the paired sequences and
    continue to the next pair.
    If FALSE, an error will be raised.
  }
}

\details{
  Currently this implementation only supports the normal 20 AAs.
  Missing or Ambiguous characters are not supported.
  
  The default multidimentional optimisation is "nmsimplex2" from the 
  gsl library \url{http://www.gnu.org/software/gsl/}.
  
  The rest of the multidimentional optimisation implementations 
  are based on the library "nlopt".
  More details regarding the algorithms, please go to 
  \url{http://ab-initio.mit.edu/wiki/index.php/NLopt_Reference#Bound_constraints}.
  
  The one dimentional optimisation implmentation is the "brent" from gsl library.
}

\value{
  A list of matrices are returned: the matrix of estimated distances, 
  the matrix of estimated distance variances, 
  the matrix of negative log-likelihood between the sequences.
}

\references{
Thorne, J.L., Kishino, H., and Felsenstein, J. (1991). An evolutionary model for maximum likelihood alignment of DNA sequences. J. Mol. Evol. 33, 114-124.

Gonnet, G.H., Cohen, M.A., and Benner, S.A. (1992). Exhaustive matching of the entire protein sequence database. Science 256, 1443-1445.
}
\author{
  Ge Tan
}


\seealso{
  \code{\link{AACharacterSet}},
  \code{\link{GONNET}}, \code{\link{GONNETBF}}
}
\examples{
  \donttest{
    ## This example is not tested due to running time > 5s
  data(GONNET)
  data(GONNETBF)
  library(seqinr)
  fasta <- read.fasta(file.path(system.file("extdata", package="TKF"),
                      "pair1.fasta"),
                      seqtype="AA", set.attributes=FALSE)
  ## 1D estimation: only distance
  TKF91(fasta, mu=5.920655e-04, 
        substModel=GONNET, substModelBF=GONNETBF, skipFailure=FALSE)
  ## 2D estimation: joint estimation of distance and mu
  TKF91(fasta, substModel=GONNET, substModelBF=GONNETBF, method="NM")
  ## only apply to a pair of sequences
  seq1 <- fasta[[1]]
  seq2 <- fasta[[2]]
  TKF91Pair(seq1, seq2, mu=5.920655e-04, 
            substModel=GONNET, substModelBF=GONNETBF)
  }
}

