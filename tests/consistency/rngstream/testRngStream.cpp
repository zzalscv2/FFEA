//  Program to test the random number streams file:    RngStream.cpp

#include <iostream>
#include <sstream>
#include "RngStream.h"
using namespace std;

int main ()
{
   double sum;
   int  i;
   RngStream g1 ("g1");
   RngStream g2 ("g2");
   RngStream g3 ("g3");

   sum = g2.RandU01 () + g3.RandU01 ();

   g1.AdvanceState (5, 3);   
   sum += g1.RandU01 ();

   g1.ResetStartStream ();
   for (i = 0;  i < 35; i++)
      g1.AdvanceState (0,1);
   sum += g1.RandU01 ();

   g1.ResetStartStream ();
   long sumi = 0;
   for (i = 0;  i < 35; i++)
      sumi += g1.RandInt (1, 10);
   sum += sumi / 100.0;

   double sum3 = 0.0;
   for (i = 0;  i < 100;  i++) {
      sum3 += g3.RandU01 ();
   }
   sum += sum3 / 10.0;

   g3.ResetStartStream ();
   for (i=1; i<=5; i++)
      sum += g3.RandU01 ();

   for (i=0; i<4; i++)
      g3.ResetNextSubstream ();
   for (i=0; i<5; i++)
      sum += g3.RandU01 ();

   g3.ResetStartSubstream ();
   for (i=0; i<5; i++)
      sum += g3.RandU01 ();

   g2.ResetNextSubstream ();
   sum3 = 0.0;
   for (i=1; i<=100000; i++)
      sum3 += g2.RandU01 ();
   sum += sum3 / 10000.0;

   g3.SetAntithetic (true);
   sum3 = 0.0;
   for (i=1; i<=100000; i++)
      sum3 += g3.RandU01 ();
   sum += sum3 / 10000.0;

   unsigned long germe[6] = { 1, 1, 1, 1, 1, 1 };
   RngStream::SetPackageSeed (germe);

   RngStream gar[4] = { "Poisson", "Laplace", "Galois", "Cantor" };
   for  (i = 0; i < 4; i++)
      sum += gar[i].RandU01 ();

   gar[2].AdvanceState (-127, 0);
   sum += gar[2].RandU01 ();

   gar[2].ResetNextSubstream ();
   gar[2].IncreasedPrecis (true);
   sum3 = 0.0;
   for  (i = 0; i < 100000; i++)
      sum3 += gar[2].RandU01 ();
   sum += sum3 / 10000.0;

   gar[2].SetAntithetic (true);
   sum3 = 0.0;
   for  (i = 0; i < 100000; i++)
      sum3 += gar[2].RandU01 ();
   sum += sum3 / 10000.0;
   gar[2].SetAntithetic (false);

   gar[2].IncreasedPrecis (false);
   for  (i = 0; i < 4; i++)
      sum += gar[i].RandU01 ();

   cout.precision(14);
   string result = "39.697547445251";
   stringstream ss; 
   ss.precision(14);
   ss << sum;
   if ( result.compare(ss.str()) ) {
     cout << "wrong result: " << sum << " " << ss.str() << endl;
     return 1;
   } else {
     // cout << "right result: " << sum << " " << ss.str() << endl;
     return 0; 
   }

   //   cin >> i; // for Visual C++
   return 0;
}
