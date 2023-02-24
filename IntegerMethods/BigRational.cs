using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace IntegerMethods
{
    /// <summary>
    /// Class that represents rational numbers, with BigInteger support, as a/b
    /// </summary>
    public class BigRational
    {
        //fraction is a/b
        public BigInteger a { get; private set; }
        public BigInteger b { get; private set; }

        /// <summary>
        /// Default initializer as a/b
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        public BigRational(BigInteger a, BigInteger b)
        {
            if (b == 0)
                throw new ArgumentException("nooo zero denom");
            if (b < 0)
            {
                b *= -1;
                a *= -1;
            }
            this.a = a;
            this.b = b;
            if (a == 0)
                this.b = 1;

            Reduce();
        }

        /// <summary>
        /// Initalizes as zero
        /// </summary>
        public BigRational() : this(0, 1)
        {
        }

        /// <summary>
        /// Cast from a BigRational to BigInteger
        /// </summary>
        /// <param name="b"></param>
        public static implicit operator BigRational(BigInteger q) => new BigRational(q, 1);

        /// <summary>
        /// Cast from a BigRational to integer
        /// </summary>
        /// <param name="b"></param>
        public static implicit operator BigRational(int q) => new BigRational(q, 1);

        /// <summary>
        /// Absolute value
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static BigRational Abs(BigRational m)
        {
            return new BigRational(BigInteger.Abs(m.a), BigInteger.Abs(m.b));
        }

        /// <summary>
        /// Adds two BigRationals
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static BigRational operator +(BigRational m, BigRational n)
        {
            return new BigRational(m.a * n.b + m.b * n.a, n.b * m.b);
        }

        /// <summary>
        /// Subtracts two BigRationals
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static BigRational operator -(BigRational m, BigRational n)
        {
            return new BigRational(m.a * n.b - m.b * n.a, n.b * m.b);
        }

        /// <summary>
        /// Multiplies two BigRationals
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static BigRational operator *(BigRational m, BigRational n)
        {
            return new BigRational(m.a * n.a, m.b * n.b);
        }

        /// <summary>
        /// Multiplies an int and BigRational
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static BigRational operator *(int m, BigRational n)
        {
            return new BigRational(m * n.a, n.b);
        }

        /// <summary>
        /// Multiplies an int and BigRational
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static BigRational operator *(BigRational m, int n)
        {
            return new BigRational(n * m.a, m.b);
        }

        /// <summary>
        /// Divides two BigRationals
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static BigRational operator /(BigRational m, BigRational n)
        {
            return new BigRational(m.a * n.b, m.b * n.a);
        }

        /// <summary>
        /// Returns the string version of this rational number
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            if (b == 1)
                return "" + a;
            return "" + a + "/" + b;
        }

        /// <summary>
        /// Less than operator
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static bool operator <(BigRational m, BigRational n)
        {
            return m.a * n.b < n.a * m.b;
        }

        /// <summary>
        /// Greater than operator
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static bool operator >(BigRational m, BigRational n)
        {
            return !(m < n);
        }

        /// <summary>
        /// Equals operator
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static bool operator ==(BigRational m, BigRational n)
        {
            return m.a * n.b == n.a * m.b;
        }

        /// <summary>
        /// Not equal operator
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static bool operator !=(BigRational m, BigRational n)
        {
            return m.a * n.b != n.a * m.b;
        }

        /// <summary>
        /// Default HashCode
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }

        /// <summary>
        /// Returns this fraction to reduced form
        /// </summary>
        private void Reduce()
        {
            List<BigInteger> aPrimes = IntegerFunctions.PrimeFactorize(BigInteger.Abs(a));
            List<BigInteger> bPrimes = IntegerFunctions.PrimeFactorize(b);

            int i = 0;
            int j = 0;
            while (i < aPrimes.Count && j < bPrimes.Count)
            {
                if (aPrimes[i] < bPrimes[j])
                    i++;
                else if (aPrimes[i] > bPrimes[j])
                    j++;
                else if (aPrimes[i] == bPrimes[j])
                {
                    a = a / aPrimes[i];
                    b = b / bPrimes[j];
                    i++;
                    j++;
                }
            }
        }

        /// <summary>
        /// Returns a decimal representation of this fraction, to n digits.
        /// The first value in the digits list is the place of the first digit:
        /// i.e., if the decimal representation is .0045623, the first value of the list
        /// is -3, since the first digit is .004 = 4*10^{-3}
        /// If the decimal rep. is 78.432, the first value is 1.
        /// 
        /// Use 
        /// </summary>
        /// <returns></returns>
        public List<int> toDecimal(int n)
        {
            List<int> digitList = new List<int>();
            //the first value of digitlist - see function description
            int exponentOffset = 0;
            BigInteger offset = 1;

            BigInteger r = BigInteger.Abs(a);
            if (r < b)
            {
                while (r * offset <= b)
                {
                    offset *= 10;
                    exponentOffset--;
                }
            }
            else
            {
                while (b * offset * 10 < r)
                {
                    offset *= 10;
                    exponentOffset++;
                }
            }

            digitList.Add(exponentOffset);

            BigInteger offsetB = b;
            if (exponentOffset < 0)
            {
                r = offset * r;
            }
            else
            {
                offsetB = offset * b;
            }
            for (int i = 0; i < n; i++)
            {
                int q = (int)(r / offsetB);
                r = r - offsetB * q;

                digitList.Add(q);
                r *= 10;
                exponentOffset--;
            }

            if (a < 0)
                digitList[1] *= -1;
            return digitList;
        }

        /// <summary>
        /// Returns whether the object is an equal BigRational
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        /// <exception cref="NotImplementedException"></exception>
        public override bool Equals(object? obj)
        {
            if (!(obj is BigRational))
                return false;

            return ((BigRational)obj) == this;
        }
    }
}
