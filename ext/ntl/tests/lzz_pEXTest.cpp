
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pEXFactoring.h>

NTL_CLIENT

int main()
{
   zz_p::init(17);

   zz_pX P;
   BuildIrred(P, 10);

   zz_pE::init(P);

   zz_pEX f, g, h;

   random(f, 20);
   SetCoeff(f, 20);

   random(h, 20);

   g = MinPolyMod(h, f);

   if (deg(g) < 0) Error("bad zz_pEXTest (1)");
   if (CompMod(g, h, f) != 0)
      Error("bad zz_pEXTest (2)");


   
   vec_pair_zz_pEX_long v;

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
            Error("bad zz_pEXTest (3)");


   }

   cerr << "\n";

   cerr << "zz_pEXTest OK\n";
}
