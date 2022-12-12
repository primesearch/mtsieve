Complete documentation for mtsieve and the programs built using the framework can be found
at http://www.mersenneforum.org/rogue/mtsieve.html.

Here is a list of the programs current built on the framework.  The "cl" versions of the programs run on GPUs. 

   afsieve/afsievecl      Find factors of Alternating Factorials.
   cksieve/cksievecl      Find factors of Carol / Kynea numbers. These numbers a form of Near Square numbers with the form (b^n-1)^2-2 and (b^n+1)^2-2.  
   dmdsieve               Find factors of number of the form 2*k*(2^p-1)+1.  Numbers that are not removed from the sieve are potential divisors of Double Mersenne numbers.  
   fbncsieve              Find factors of numbers in the form k*b^n+1 and k*b^n-1.
   fkbnsieve              Find factors of the form k*b^n+c for fixed k, b, and n and variable c.   
   gcwsieve/gcwsievecl    Find factors of Cullen and Woodall numbers.  
   gfndsieve/gfndsievecl  Find factors of numbers in the form k*2^n+1 for a range of k and a range of k.  The output from this sieve should be used with pfgw and the -gxo switch to find GFN divisors. 
   k1b2sieve              Find factors of numbers of the form b^n+c for fixed b and variable c.  
   kbbsieve               Find factors of numbers of the form k*b^b+1 or k*b^b-1 for fixed k and variable b.  
   mfsieve/mfsievecl      Find factors of MultiFactorials.  
   pixsieve/pixsievecl    Find factors of numbers that are a substring of a long decimal string where each successive term adds on decimal digit to the end of the previous decimal term.  
   psieve/psievecl        Find factors of primorials.  
   sgsieve                Find factors of Sophie-Germain numbers of the form k*b^n-1 with variable k, fixed b, and fixed n.
   smsieve/smsievecl      Find factors of Smarandache numbers.
   srsieve2/srsieve2cl    Find factors of Sierpinski/Riesel sequences of the form k*b^n+1 or k*b^n-1 for fixed b, variable n, and multiple k.
   twinsieve              Find factors of twin numbers of the form k*b^n+1 and k*b^n-1.  
   xyyxsieve/xyyxsievecl  Find factors of x^y+y^x and x^y-y^x numbers.  

