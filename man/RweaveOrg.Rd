% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RweaveAscii.r, R/weaverAscii.r
\name{RweaveOrg}
\alias{RweaveOrg}
\alias{RtangleOrg}
\alias{RweaveOrgOptions}
\alias{RweaveOrgFinish}
\alias{RweaveOrgWritedoc}
\alias{RweaveOrgSetup}
\alias{RweaveOrgRuncode}
\alias{weaverOrg}
\alias{weaverOrgSetup}
\title{A driver to parse org noweb files with Sweave tool
This driver parses org files containing R code and replace pieces of code
with their output.}
\usage{
RweaveOrg()

RweaveOrgSetup(
  file,
  syntax,
  output = NULL,
  quiet = FALSE,
  debug = FALSE,
  extension = "org",
  backend = "org-mode",
  openSchunk = "#+BEGIN_example",
  closeSchunk = "\\n#+END_example\\n",
  openSinput = "",
  closeSinput = "",
  openSoutput = "\\n",
  closeSoutput = "",
  indent = "",
  openInclude = "#+INCLUDE: \\"",
  closeInclude = ".org\\"",
  openFig = "[[file:",
  closeFig = "]]",
  ...
)

weaverOrgSetup(
  file,
  syntax,
  output = NULL,
  quiet = FALSE,
  debug = FALSE,
  extension = "org",
  backend = "org-mode",
  openSchunk = "#+BEGIN_example",
  closeSchunk = "\\n#+END_example\\n",
  openSinput = "",
  closeSinput = "",
  openSoutput = "\\n",
  closeSoutput = "",
  indent = "",
  openInclude = "#+INCLUDE: \\"",
  closeInclude = ".org\\"",
  openFig = "[[file:",
  closeFig = "]]",
  use.cache = TRUE,
  ...
)

weaverOrg()
}
\arguments{
\item{file}{file}

\item{syntax}{syntax}

\item{output}{output}

\item{quiet}{quite}

\item{debug}{debug}

\item{...}{...}

\item{use.cache}{use.cache}
}
\value{
None value is returned. From a .Rnw noweb file, the corresponding
  .org is produced (as eventuals files for graphs).
}
\description{
A driver to parse org noweb files with Sweave tool
This driver parses org files containing R code and replace pieces of code
with their output.

RweaveOrgSetup

weaverOrgSetup

weaverOrg
}
\note{
In order to work properly, noweb codes have to be located at the
  beginning of a line (no indentation).

Compare with RweaveLatex driver, RweaveOrg provides one new option :
  \code{format} to choose the format of figure that will be inserted in the
  final document.

In addition, \code{cache} option from \code{weaver}
  package is also available with \code{weaverOrg} driver.

A wrapper for \code{Sweave} can be used, named \code{Org}.
}
\examples{
  \dontrun{
Org(system.file("examples/Org-test-1.nw",package="ascii"))
Org(system.file("examples/Org-test-1.nw",package="ascii"),weaverOrg())
  }

}
\seealso{
\code{\link[utils]{Sweave}}, \code{\link{Org}}
}
\author{
David Hajage \email{dhajage@gmail.com}
}
\keyword{IO}
\keyword{file}
\keyword{internal}
