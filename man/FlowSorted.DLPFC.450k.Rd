\name{FlowSorted.DLPFC.450k}
\alias{FlowSorted.DLPFC.450k}
\docType{data}
\title{
  Illumina Human Methylation data from 450k on sorted frontal cortex
  cell populations
}
\description{
  This \code{RGChannelSet} contains Illumina 450k DNA methylation
  measurements  on 58 flow-sorted dorsolateral prefrontal cortex samples
  from non-psychiatric controls from Guintivano et al. 2013. These
  samples were separated into neuronal (NeuN+) and non-neuronal (NeuN-)
  cell types. These data can be used by the \link{minfi} package to
  estimate cellular composition from bulk frontal cortex samples.  This
  data may also be useful to individuals as example Illumina 450k data
  for trying preprocessing methods across a variety of Bioconductor
  packages.
}

\usage{
data(FlowSorted.DLPFC.450k)
}

\details{
  The \code{FlowSorted.DLPFC.450k} objects is based an samples assayed
  as part of Guintivano et al (2013).  Please cite this paper, if the
  data is used.  If you're using this data together with the
  \link{minfi} package, please see the package vignette for details on
  how to cite that package.
}
\format{
  An object of class \code{RGChannelSet}.
}
\references{
  Guintivano J., Aryee M.J., Kaminsky Z.A.
  \emph{A cell epigenotype specific model for the correction of 
    brain cellular heterogeneity bias and its application to age, 
    brain region and major depression.}
  Epigenetics (2013), 8(3):290-302.
  \doi{10.4161/epi.23924}.
}

\examples{
   data(FlowSorted.DLPFC.450k)
}
\keyword{datasets}
