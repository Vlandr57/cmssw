double FitFunVer02(double P0, double P1, double P2, double P3, 
                                         double P4, double P5) 
{

// Loop over histo's bins and calculate chi2_sum:
// ----------------------------------------------

   double chi2_sum = 0.0;
   double depth, sigma;
   double fx1, fx2;
   double zz, delta, xx1, xx2;
   double aa1, bb1, aa2, bb2;

   Int_t maxBins = hev[iev]->GetSize()-1;
   for( Int_t bin=1; bin<maxBins; ++bin ) 
   {
      depth = hev[iev]->GetBinCenter(bin);
      zz    = hev[iev]->GetBinContent(bin);
      if ( zz <=0.0 ) continue;

// 1-st gamma-functon
//-------------------
      xx1 = P1*depth;
      if(xx1 <=0.0) {fx1 = 0.0;}
      else 
      {
        aa1 = pow(xx1,P0-1.)*exp(-xx1);
        bb1 = ROOT::Math::tgamma(P0);
        if (bb1 > 0.0) fx1 = (P2)*(aa1/bb1);
        else fx1 = 0.0;
      }

// 2-nd gamma-function
//--------------------
      xx2 = P4*depth;
      if(xx2 <=0.0) {fx2 = 0.0;}
      else 
      {
        aa2 = pow(xx2,P3-1.)*exp(-xx2);
        bb2 = ROOT::Math::tgamma(P3);
        if (bb2 > 0.0) fx2 = (1.-P2)*(aa2/bb2);
        else fx2 = 0.0;
      }

      sigma = hev[iev]->GetBinError(bin);
      if ( sigma > 0.0) delta = (zz - P5*(fx1+fx2))/sigma;
      else continue;

      chi2_sum += delta*delta;
   }
    return (chi2_sum);
}

//================================================ 
//== Function which is used in MINUIT program == 
//================================================

void myFCN(int& nDim, double* gout, double& result,
                                    double par[], int flg) 
{
    result = FitFunVer02(par[0], par[1], par[2], par[3], par[4], par[5]);
}


