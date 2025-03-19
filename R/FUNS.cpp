#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//



// [[Rcpp::export]]
NumericVector countC(NumericVector cub,NumericVector cua) {
  int n = cub.size();
  int m = cua.size();
  NumericVector cuaP(m+2);
  cuaP[0]=R_NegInf;
  int i,j;
  for (i=0;i<m;i++){
    cuaP[i+1]=cua[i];
  }
  cuaP[m+1]=R_PosInf;

  int m2 = cuaP.size();
  NumericVector counts(m2-1);
  for (i=0;i<(m2-1);i++){
    counts[i]=1;
    for(j=0;j<n;j++){
      if((cub[j]<=cuaP[i+1])&(cub[j]>cuaP[i])) counts[i]=counts[i]+1;
    }
  }

  return counts;
}


// [[Rcpp::export]]
NumericVector stl_sort(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}


// [[Rcpp::export]]
double DW2_discr_C(NumericVector a, NumericVector b,
                   NumericVector wa, NumericVector wb){
  int m=a.size();
  int n=b.size();
  NumericVector ua=wa/sum(wa);
  NumericVector ub=wb/sum(wb);
  NumericVector cua=cumsum(ua[Range(0, m-2)]);
  NumericVector cub=cumsum(ub[Range(0, n-2)]);
  NumericVector arep=countC(cub,cua);
  NumericVector brep=countC(cua,cub);
  NumericVector aa(sum(arep));
  NumericVector bb(sum(brep));
  int i,j,C;
  //     aa <- rep(a, times = arep)
  C=-1;
  for(i=0;i<arep.size();i++){
    for(j=0;j<arep(i);j++){
      C=C+1;
      aa(C)=a(i);
    }
  }

  C=-1;

  //     bb <- rep(b, times = brep)
  for(i=0;i<brep.size();i++){
    for(j=0;j<brep(i);j++){
      C=C+1;
      bb(C)=b(i);
    }
  }
  //     uu <- sort(c(cua, cub),method="quick")
  NumericVector uu(cua.size()+cub.size());
  uu[Range(0, cua.size()-1)]=cua;
  uu[Range(cua.size(),uu.size()-1)]=cub;
  uu=stl_sort(uu);
  NumericVector uu0(uu.size()+1);
  NumericVector uu1(uu.size()+1);
  //     uu0 <- c(0, uu)
  //     uu1 <- c(uu, 1)
  uu0[Range(1,uu0.size()-1)]=uu;
  uu1[Range(0,uu1.size()-2)]=uu;
  uu1[uu1.size()-1]=1;
  double areap;
  areap=sqrt(sum((uu1-uu0)*(bb-aa)*(bb-aa)));

  //     areap <- sqrt(sum((uu1 - uu0) * (bb - aa)^2))
  //double res=2;
  return areap;
}




// [[Rcpp::export]]
NumericMatrix dmatC(DataFrame T) {
  Function WD2 = Environment::global_env()["Wass1d_discrete_2"];
  int cols=T.size();
  // Rcout << "Cols : " << cols << "\n";
  List Temp = T[0];//la colonna
  int rows=Temp.size();
  // Rcout << "Rows : " << rows << "\n";
  NumericMatrix DM(rows);

  int i1,i2,j;
  for (i1=0;i1<(rows-1);i1++){
    for (i2=(i1+1);i2<rows;i2++){
      for (j=0;j<cols;j++){
        List COL = T[j];
        DataFrame ROW1=COL[i1];
        DataFrame ROW2=COL[i2];
        NumericVector a = ROW1[0];
        NumericVector wa = ROW1[2];
        NumericVector b = ROW2[0];
        NumericVector wb = ROW2[2];
        // Rcout << "a : " << a << "\n";
        // Rcout << "wa : " << wa << "\n";
        // Rcout << "b : " << b << "\n";
        // Rcout << "wb : " << wb << "\n";
        SEXP val=WD2(a,b,wa,wb);
        double val2=*REAL(val);
        // Rcout << "DD : " << *REAL(val) << "\n";
        DM(i1,i2)=DM(i1,i2)+(val2*val2);
      }
      DM(i1,i2)=sqrt(DM(i1,i2));
    DM(i2,i1)=DM(i1,i2);
    }
  }
 // DataFrame Temp2 = Temp[0];//la riga
 //  // List Temp3 = Temp2[0];
 //  NumericVector T2 = Temp2[0];

  return DM;
}

// [[Rcpp::export]]
NumericMatrix dmatC2(DataFrame T) {
  int cols=T.size();
  // Rcout << "Cols : " << cols << "\n";
  List Temp = T[0];//la colonna
  int rows=Temp.size();
  // Rcout << "Rows : " << rows << "\n";
  NumericMatrix DM(rows);

  int i1,i2,j;
  for (i1=0;i1<(rows-1);i1++){
    for (i2=(i1+1);i2<rows;i2++){
      for (j=0;j<cols;j++){
        List COL = T[j];
        DataFrame ROW1=COL[i1];
        DataFrame ROW2=COL[i2];
        NumericVector a = ROW1[0];
        NumericVector wa = ROW1["freq"];
        NumericVector b = ROW2[0];
        NumericVector wb = ROW2["freq"];
        // Rcout << "a : " << a << "\n";
        // Rcout << "wa : " << wa << "\n";
        // Rcout << "b : " << b << "\n";
        // Rcout << "wb : " << wb << "\n";

        double val2=DW2_discr_C(a,b,wa,wb);
        // Rcout << "DD : " << *REAL(val) << "\n";
        DM(i1,i2)=DM(i1,i2)+(val2*val2);
      }
      DM(i1,i2)=sqrt(DM(i1,i2));
      DM(i2,i1)=DM(i1,i2);
    }
  }
  // DataFrame Temp2 = Temp[0];//la riga
  //  // List Temp3 = Temp2[0];
  //  NumericVector T2 = Temp2[0];

  return DM;
}

// [[Rcpp::export]]
NumericMatrix DiToCen_C(DataFrame Tib,
                        DataFrame Cen,
                        int n, int p, int k) {
  NumericMatrix DTP(n,k);
  int i,j,clu;
  for(i=0;i<n;i++){
    for(clu=0;clu<k;clu++){
      for(j=0;j<p;j++){
        List COL = Tib[j];
        List COLc = Cen[j];
        DataFrame ROW1=COL[i];
        DataFrame ROW2=COLc[clu];
        NumericVector a = ROW1[0];
        NumericVector wa = ROW1["freq"];
        NumericVector b = ROW2[0];
        NumericVector wb = ROW2["freq"];
        double val2=DW2_discr_C(a,b,wa,wb);
        // Rcout << "DD : " << *REAL(val) << "\n";
        DTP(i,clu)=DTP(i,clu)+(val2*val2);
      }
    }
  }
  return DTP;
}
// Dist_to_prot<-matrix(0,n,k)
//   for(i in 1:n){
//     for (clu in 1:k){
//       for (j in 1:p){
//         Dist_to_prot[i,clu]=Dist_to_prot[i,clu]+
//           DW2_discr_C(a=Tib[i,j][[1]][[1]],
//                       b=centers[clu,j][[1]][[1]],
//                                            wa=Tib[i,j][[1]]$freq,
//                                            wb=centers[clu,j][[1]]$freq)^2
//       }
//     }
//   }



// SEXP Trans_D(SEXP in) {
//   Rcpp::DataFrame dfin(in);
//   Rcpp::DataFrame dfout;
//   Rcpp::CharacterVector namevec;
//   std::string namestem = "Column Heading ";
//   for (int i=0;i<2;i++) {
//     dfout.push_back(dfin(i));
//     namevec.push_back(namestem+std::string(1,(char)(((int)'a') + i)));
//   }
//   dfout.attr("names") = namevec;
//   Rcpp::DataFrame x;
//   Rcpp::Language call("as.data.frame",dfout);
//   x = call.eval();
//   return x;
// }

// [[Rcpp::export]]
List TransposeTib_C(DataFrame Tib){
  List Temp = Tib[0];//la colonna
  int cols = Tib.size();
  int rows=Temp.size();

  List TT;
  Rcpp::CharacterVector namevec;
  std::string namestem = "CH ";

   int i,j;

   for(i=0;i<rows;i++){
      List tmp(cols);
     namevec.push_back(namestem+std::string(1,(char)(((int)'a') + i)));
      for(j=0;j<cols;j++){

        List COL = Tib[j];
        tmp[j]=COL[i];
      }
//      Rcout << "DD : " << rows<<"-"<<cols << "\n";

      TT.push_back(tmp);
    }

  return TT;
}
