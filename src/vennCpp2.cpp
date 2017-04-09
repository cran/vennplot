#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// for new R-devel
void R_init_vennplot(DllInfo* info) {
	R_registerRoutines(info, NULL, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
}

void loss(double& gammar, NumericVector& detaloss, 
          NumericMatrix& xy, int i, 
          NumericVector& radius, NumericMatrix& ED, bool ThreeD,
          NumericVector& tmp) {
  double detaloss1 = 0.0;
  double detaloss2 = 0.0;
  double detaloss3 = 0.0;
  double X2 = 0.0;
  double total = 0.0;
  double totalminus = 0.0;
  int m = xy.nrow();
  for(int j = 0; j<m; j++){
    tmp = xy.row(i) - xy.row(j);
    X2 = sum(tmp*tmp);
    total = radius[i]+radius[j];
    totalminus = fabs(radius[i]-radius[j]);
    if(ED(i,j)>= total  && sqrt(X2)>= total){
      continue;
    }
    else if(ED(i,j)<= totalminus && sqrt(X2)<= totalminus){
      continue;
    }
    else{
      gammar += pow((X2 -pow(ED(i,j),2)),2);
      detaloss1 += 4.0*(X2 - pow(ED(i,j),2))*(xy.row(i)[0] - xy.row(j)[0]);
      detaloss2 += 4.0*(X2 - pow(ED(i,j),2))*(xy.row(i)[1] - xy.row(j)[1]);
      if(ThreeD == true){
        detaloss3 += 4.0*(X2 - pow(ED(i,j),2))*(xy.row(i)[2] - xy.row(j)[2]);}
    }
    
  }
  detaloss[0] = detaloss1;
  detaloss[1] = detaloss2;
  if(ThreeD == true){
    detaloss[2] = detaloss3;
  }
  return;
}

// [[Rcpp::export("lossR")]]
List loss_R(NumericMatrix xy, int iRow, 
            NumericVector radius, NumericMatrix ED, bool ThreeD) {
  double gammar = 0.0;
  int i = iRow-1;
  NumericVector detaloss = xy.row(i);
  NumericVector tmp(xy.ncol());
  loss(gammar, detaloss, xy, i, radius, ED, ThreeD, tmp);
  return List::create(Named("x") = gammar,
                      Named("y") = detaloss);
}

void Loss(double& loss2,
          NumericMatrix& xy,NumericVector& radius, 
          NumericMatrix& ED, bool ThreeD){
  int m = xy.nrow();
  NumericVector loss1(m);
  NumericVector detaloss(xy.ncol());
  NumericVector tmp(xy.ncol());
  for(int j = 0; j<m; j++){
    //loss1[j] = loss(m = m, xy = xy ,j , radius = radius, ED = ED , ThreeD = ThreeD);
    detaloss = xy.row(j);
    loss(loss1[j], detaloss,xy, j, radius, ED, ThreeD, tmp);
    
  }
  loss2 = sum(loss1);
  return;
}


// [[Rcpp::export("LossR")]]
double Loss_R(NumericMatrix xy,
              NumericVector radius, NumericMatrix ED, bool ThreeD) {
  double loss2 = 0.0;
  Loss(loss2, xy, radius,ED, ThreeD);
  return loss2;
}

// [[Rcpp::export("LoopR")]]
List loop_R(NumericMatrix xy, double ALPHA, 
            NumericVector radius, NumericMatrix ED, bool ThreeD) {
  int m = xy.nrow();
  NumericMatrix xy1(clone(xy));
  NumericMatrix xyn(clone(xy));
  NumericVector detax0(xy.ncol());
  NumericVector tmp(xy.ncol());
  NumericMatrix sn(clone(xy));
  int xigma = 0;
  double gammar = 0.0;
  for(int i=0;i<m;i++){
    loss(gammar, detax0, xy1, i, radius, ED, ThreeD, tmp);
    detax0 = -detax0;
    xy1.row(i) = xy.row(i) + ALPHA*detax0; 
  }
  
  double f1 = 0.0;
  double f2 = 0.0;
  Loss(f1, xy1, radius,ED, ThreeD);
  NumericVector detaxn(xy.ncol());
  NumericVector detaxn_1(xy.ncol());
  double betan = 0.0;
  double betan1 = 0.0;
  double tt = 0.0;
  double ts = 0.0;
  while(f1>pow(10,-10)){
    if(xigma == 0){
      for (int i=0;i<m;i++){
        loss(gammar, detaxn,xy1, i, radius, ED, ThreeD, tmp);
        loss(gammar, detaxn_1,xy, i, radius, ED, ThreeD, tmp);
        detaxn = -detaxn;
        detaxn_1 = -detaxn_1;
        ts = sum(detaxn*(detaxn-detaxn_1));
        tt = sum(detaxn_1*detaxn_1);
        betan = ts/tt;
        if(std::isnan(betan) || std::isinf(betan)){betan1 = 1;}
        else if(betan>0){betan1 = betan;}
        else{betan1 = 0;}
        sn.row(i) = detaxn + betan1*detaxn_1;
        xyn.row(i) = xy1.row(i) + ALPHA*sn.row(i);
      }
      xy = xy1;
      xy1 = xyn;
    }
    else{
      NumericMatrix xy11(clone(xyn));
      for (int i=0;i<m;i++){
        loss(gammar, detaxn, xy11, i, radius, ED, ThreeD, tmp);
        loss(gammar, detaxn_1, xy, i, radius, ED, ThreeD, tmp);
        detaxn = -detaxn;
        detaxn_1 = -detaxn_1;
        ts = sum(detaxn*detaxn);
        tt = sum(detaxn_1*detaxn_1);
        betan = ts/tt;
        if(std::isnan(betan) || std::isinf(betan)){betan1 = 1;}
        else if(betan>0){betan1 = betan;}
        else{betan1 = 0;}
        sn.row(i) = detaxn + betan1*sn.row(i);
        xyn.row(i) = xy11.row(i) + ALPHA*sn.row(i);
      }
      Loss(f1, xyn, radius,ED, ThreeD);
      Loss(f2, xy11, radius,ED, ThreeD);
      if(f1>f2){break;}
      NumericMatrix xy1(clone(xy11));
      xy = xy1;
    }
    xigma++;
  }
  return List::create(Named("xy") = xyn,
                      Named("f1") = f1);
}

void transR(NumericVector& xyvec,NumericMatrix& xy, NumericVector& radius, 
            double& radiusvec, NumericVector& radiusall) {
  int out = 1;
  while (out==1) {
    xyvec[0] = R::runif(min(xy.column(0))-3*max(radiusall),max(xy.column(0))+3*max(radiusall));
    xyvec[1] = R::runif(min(xy.column(1))-3*max(radiusall),max(xy.column(1))+3*max(radiusall));
    if(xy.ncol()==3){
      xyvec[2] = R::runif(min(xy.column(2))-3*max(radiusall),max(xy.column(2))+3*max(radiusall));
    }
    NumericVector ed(xy.nrow());
    NumericVector ra(xy.nrow());
    for(int i = 0; i < xy.nrow(); i++){ 
      if(xy.ncol()==2){
        ed[i] = pow((xy.row(i)[0] - xyvec[0]),2)+pow((xy.row(i)[1] - xyvec[1]),2);
      }else{
        ed[i] = pow((xy.row(i)[1] - xyvec[1]),2)+pow((xy.row(i)[2] - xyvec[2]),2)+pow((xy.row(i)[0] - xyvec[0]),2);
      }
      ra[i] = radius[i] + radiusvec;
      
    }
    if(is_true(all(ed>ra))){out = 0;}else{out = 1;} 
  }
  return;
}

// [[Rcpp::export("transR")]]
NumericVector trans_R(NumericMatrix xy, NumericVector radius, double radiusvec,NumericVector radiusall){
  NumericVector xyvec(xy.ncol());
  transR(xyvec,xy,radius,radiusvec,radiusall);
  return xyvec;
}


void alldis(int& out, NumericMatrix& xy1,NumericMatrix& xy2, NumericVector& radius1,  NumericVector& radius2, double& delta){
  NumericVector ed(xy2.nrow());
  NumericVector r1(xy2.nrow());
  NumericVector r2(xy2.nrow());
  NumericMatrix xy(clone(xy1));
  NumericMatrix transxy(clone(xy2));
  out = 1;
  for(int i=0;i<xy.nrow();i++){
    for(int j=0; j<transxy.nrow();j++){
      if(xy.ncol()==2){
        ed[j] = sqrt(pow((xy.row(i)[0] - transxy.row(j)[0]),2)+pow((xy.row(i)[1] - transxy.row(j)[1]),2));
      }else{
        ed[j] = sqrt(pow((xy.row(i)[0] - transxy.row(j)[0]),2)+pow((xy.row(i)[1] - transxy.row(j)[1]),2)+pow((xy.row(i)[2] - transxy.row(j)[2]),2));
      }
      NumericMatrix xy(clone(xy1));
      NumericMatrix transxy(clone(xy2));
      r1[j] = radius1[i] + radius2[j]+delta;
      r2[j] = radius1[i] + radius2[j];
    }
    if(is_true(all(ed > r1))){continue;}
    else if(is_true(any(ed <= r1))){
      if(is_true(all(ed > r2))){out = 2;}
      else{out = 0;break;}
    }
  }
  return;
}

// [[Rcpp::export("alldisR")]]
int alldis_R(NumericMatrix xy1,NumericMatrix xy2, NumericVector radius1, NumericVector radius2, double delta){
  int out;
  alldis(out,xy1,xy2,radius1,radius2,delta);
  return out;
}

// [[Rcpp::export("closeR")]]
List close_R(NumericMatrix xy1,NumericMatrix xy2, NumericVector radius1, 
            NumericVector radius2, double delta, NumericVector direc){
  int out;
  NumericMatrix xy3(clone(xy2));
  alldis(out, xy1,xy3,radius1,radius2,delta);
  while(out!=2){
    for(int i=0; i<xy3.nrow();i++){
      xy3.row(i)[0] = xy3.row(i)[0]+direc[0];
      xy3.row(i)[1] = xy3.row(i)[1]+direc[1];
      if(xy3.ncol()==3){
        xy3.row(i)[2] = xy3.row(i)[2]+direc[2];
      }
    }
    alldis(out, xy1,xy3,radius1,radius2,delta);
  }
  return List::create(Named("out")=out,
                      Named("xy") = xy3);
}

// [[Rcpp::export("listR1")]]
NumericMatrix list_R1(NumericMatrix M, NumericMatrix xy, NumericVector radius,
                      int k ,double yuan, double xuan, int num) {
  k = k-1;
  for(int i=0 ; i<num; i++){
    for(int j=0 ; j<num; j++){
      double detect = ((i+1)*xuan+min(xy.column(0))-max(radius) - xy.row(k)[0]) *
        ((i+1)*xuan+min(xy.column(0))-max(radius) - xy.row(k)[0]) +
        ((j+1)*yuan+min(xy.column(1))-max(radius) - xy.row(k)[1]) *
        ((j+1)*yuan+min(xy.column(1))-max(radius) - xy.row(k)[1]);
      
      if(detect <= radius[k]*radius[k]){
        M(i,j) = 1;
      }
    }
  }
  return M;
}



// [[Rcpp::export("listR2")]]
NumericMatrix list_R2(Rcpp::List myList, int m, int num) {
  NumericMatrix M(num*num,m);
  for(int i=0;i<num;i++){
    for(int j=0;j<num;j++){
      NumericVector L(m);
      for(int k=0; k<m; k++){
        NumericMatrix N = myList[k];
        L[k] = N(i,j);
      }
      M.row(num*i+j) = L;
    }
  }
  return M;
}


// [[Rcpp::export("listR3")]]
NumericVector list_R3(NumericMatrix M, NumericMatrix Me){
  NumericVector Len(Me.nrow());
  for(int i=0;i<Me.nrow();i++){
    int k = 0;
    for(int j = 0;j<M.nrow();j++){
      if(is_true(all(Me.row(i)==M.row(j)))){
        k++;
      }
    }
    Len[i] = k;
  }
  return Len;
}  

// [[Rcpp::export("listR4")]]
NumericVector list_R4(NumericMatrix M) {
  int m = M.nrow();
  NumericVector n(m);
  for(int i=0;i<m;i++){
    if(sum(M.row(i))!=0){
      n[i] = i+1;
    }
  }
  return n;
}
