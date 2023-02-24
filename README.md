# RepresentationStabilityComputations
Calculates limiting multiplicites of families of irreducible representations in the cohomology of complex configuration space. 
Includes algorithm to multiply infinite dimensional symmetric polynomials and represent the Frobenius formula for a Young Tableau as a character polynomial.
In "MainMethod.cs", there are a number of useful computations contained in different code blocks. Each has a variable "k" which can be modified to run
the desired computation. 

The primary function of the code is provided by the function YoungToPoly(part, smallestDegree), which given a partition of n (a list of decreasing positive integers
with sum n) returns an infinite power series in q^{-1} representing the stable multiplicity of the family of irreducible partitions given by "part" in the
cohomology of complex configuration space. The infinite power series is expanded out to the "smallestDegree" term, so it will be cut off at q^{-smallestDegree}.

Abstract:
Representation stability was introduced to study mathematical structures which stabilize when viewed from a representation theoretic framework. 
The instance of representation stability studied in this project is that of ordered complex configuration space, denoted $\pconf_n(\C)$:

$$\pconf_n(\C) := \{ (x_1, x_2, \dots, x_n) \in \C^n \ | \ x_i \neq x_j \}$$

$\pconf_n(\C)$ has a natural $S_n$ action by permuting its coordinates which gives the cohomology groups $H^i(\pconf_n(\C);\Q)$ the structure of an $S_n$ representation. 
The cohomology of $\pconf_n(\C)$ \textit{stabilizes} as $n$ tends toward infinity when viewed as a family of $S_n$ representations. From previous work, there is an 
explicit description for $H^i(\pconf_n(\C);\Q)$ as a direct sum of induced representations for any $i, n$, but this description does not explain the behavior of 
families of irreducible representations as $n\to\infty$. We implement an algorithm which, given a Young Tableau, computes the cohomological degrees where the 
corresponding family of irreducible representations appears stably as $n\to\infty$. Previously, these values were known for only a few Young Tableaus and cohomological 
degrees. Using this algorithm, results have been found for all Young Tableau with up to 8 boxes and certain Tableau with more, which has led us to conjectures based on 
the data collected.


