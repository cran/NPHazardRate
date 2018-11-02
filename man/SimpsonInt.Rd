\name{SimpsonInt}
\title{Simpson numerical integration}
\alias{SimpsonInt}

\description{
Implements Simpson's extended numerical integration rule
}

\usage{
SimpsonInt(xin, h)
}

\arguments{
\item{xin}{ A vector of data points}
\item{h}{ grid length}
}

\details{
The extended numerical integration rule is given by
\deqn{
\int_0^{x_{2n}} f(x)\,dx = \frac{h}{3}(f(x_0) + 4\{f(x_1) + \dots f(x_{2n-1}) \}
      +2 \{f(x_2) + f(x_4) + \dots f(x_{2n-2})\} + f(x_{2n})) -R_n
}
}

\value{
returns the approximate integral value
}

\references{
  \href{http://mathworld.wolfram.com/SimpsonsRule.html}{Weisstein, Eric W. "Simpson's Rule." From MathWorld--A Wolfram Web Resource}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}




%}
