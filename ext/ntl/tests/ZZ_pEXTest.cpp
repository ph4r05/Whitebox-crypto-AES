
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEXFactoring.h>

NTL_CLIENT

int main()
{
   ZZ_p::init(to_ZZ(17));

   ZZ_pX P;
   BuildIrred(P, 10);

   ZZ_pE::init(P);

   ZZ_pEX f, g, h;

   random(f, 20);
   SetCoeff(f, 20);

   random(h, 20);

   g = MinPolyMod(h, f);

   if (deg(g) < 0) Error("bad ZZ_pEXTest (1)");
   if (CompMod(g, h, f) != 0)
      Error("bad ZZ_pEXTest (2)");


   
   vec_pair_ZZ_pEX_long v;

   long i;
   for (i = 0; i < 5; i++) {
      long n = RandomBnd(20)+1;
      cerr << n << " ";

      random(f, n);
      SetCoeff(f, n);

      v = CanZass(f);

      g = mul(v);
      if (f != g) cerr << "oops1\n";

      long i;
      for (i = 0; i < v.length(); i++)
         if (!DetIrredTest(v[i].a))
            Error("bad ZZ_pEXTest (3)");


   }

   cerr << "\n";

   cerr << "ZZ_pEXTest OK\n";
}
