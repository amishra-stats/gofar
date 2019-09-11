#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>

using namespace arma;
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
double getKappaC0zv(arma::mat X, arma::uvec t1, arma::uvec t2, arma::uvec t3){
  //arma::uvec t1, arma::uvec t2, arma::uvec t3
  arma::vec s = svd( X);
  // double s = conv_to< double >::from(svds( X, 1) );
  arma::vec tem(3);
  // arma::mat x(size(X));
  // arma::uvec t1 = find(familygroup==1),t2 = find(familygroup==2);
  // arma::uvec t3 = find(familygroup==3);
  if(t1.n_elem>0) tem(0) = s(0);
  if(t2.n_elem>0) tem(1) = 0.5*s(0);
  if(t3.n_elem>0) tem(2) = 10*s(0);
  return(tem.max());
}











// [[Rcpp::export]]
arma::mat familyLinkinv3(arma::mat X, arma::uvec t1, arma::uvec t2, arma::uvec t3){
  arma::mat x(size(X));
  // arma::uvec t1 = find(familygroup==1),t2 = find(familygroup==2);
  // arma::uvec t3 = find(familygroup==3);
  if(t1.n_elem>0) x.cols(t1) = X.cols(t1);
  if(t2.n_elem>0) {
    x.cols(t2) = 1/(1+exp(-1*X.cols(t2)));
  }
  if(t3.n_elem>0) {
    x.cols(t3) = exp(X.cols(t3));
  }
  x.elem( find_nonfinite(x) ).zeros();
  return(x);
}


// [[Rcpp::export]]
double absT(double Xm){
  if(Xm<0)
    return(-1.0*Xm);
  else
    return((Xm));
}



// [[Rcpp::export]]
int nzcount(arma::vec x) {
  arma::vec y = nonzeros(x) ;
  return y.n_elem;
}



// [[Rcpp::export]]
double softThres(double x, double lambda) {
  return((x > lambda) ? x - lambda :
           (x < -lambda) ? x + lambda : 0.);
}

// [[Rcpp::export]]
arma::vec softT(arma::vec x, arma::vec lambda) {
  arma::vec y; y.zeros(x.n_elem);
  for(int i=0; i < (int) x.n_elem; i++)
    y(i) = softThres(x(i),lambda(i));
  return(y);
}





// [[Rcpp::export]]
double logisticlik(arma::mat Y, arma::mat MU,arma::mat nind) {
  // vec a = conv_to< vec >::from(sum(nind,0) );
  arma::mat b = (Y+MU);
  uvec c = find((b == 0) || (b == 2));
  b = Y%log(MU) + (1-Y)%log(1-MU); b=-1*(b%nind);
  b.elem(c).zeros();
  b.elem( find_nonfinite(b) ).zeros();
  // if((y+mu) == 0 || (y+mu) == 2) out=0; else out=-1*y*log(mu)-(1-y)*log(1-mu);
  // 0.5*(log(pi2*sigma2)  + ((y-mu)*(y-mu))/sigma2 )
  return(accu(b));
}




// [[Rcpp::export]]
double gausslikehood(arma::mat Y, arma::mat MU, arma::vec Phi,arma::mat nind) {
  double pi2 = datum::pi;
  arma::vec a = conv_to< arma::vec >::from(sum(nind,0) );
  arma::mat b = (Y-MU)%nind;
  b.elem( find_nonfinite(b) ).zeros();
  double temp = 0;
  temp = 0.5*accu(a%log(2*pi2*Phi) + (conv_to< arma::vec >::from(sum(square(b),0) ))/Phi );
  // 0.5*(log(pi2*sigma2)  + ((y-mu)*(y-mu))/sigma2 )
  return(temp);
}






// [[Rcpp::export]]
double poissonlik(arma::mat Y, arma::mat MU,arma::mat nind) {
  // constant part in the negative likelihood is being ignored:
  // vec a = conv_to< vec >::from(sum(nind,0) );
  // arma::mat b = (Y+MU);
  // uvec c = ;
  // b.elem(find(b < 0)).zeros();
  // b = logfactM2(Y); // MU; //
  // b= b+ MU - Y%log(MU);
  arma::mat b= log(MU);
  b.elem( find_nonfinite(b) ).zeros();
  b = MU - Y%b;
  // arma::mat b= MU - Y%log(MU);
  b= b%nind;
  // b.elem(find((Y+MU) < 0)).zeros();
  // b.elem( find_nonfinite(b) ).zeros();
  // if((y+mu) > 0)    out = logfactM((int) y)  + mu - y*log(mu);
  return(accu(b));
}

// [[Rcpp::export]]
double poissondev(arma::mat Y, arma::mat MU,arma::mat nind) {
  // constant part in the negative likelihood is being ignored:
  // vec a = conv_to< vec >::from(sum(nind,0) );
  // arma::mat b = (Y+MU);
  // uvec c = ;
  // b.elem(find(b < 0)).zeros();
  // b = logfactM2(Y); // MU; //
  // b= b+ MU - Y%log(MU);
  arma::mat b= (log(Y));
  b.elem( find_nonfinite(b) ).zeros();
  b= Y%(b -log(MU))+MU -Y;
  b= b%nind;
  // b.elem(find((Y+MU) <= 0)).zeros();
  // b.elem( find_nonfinite(b) ).zeros();
  // if((y+mu) > 0)    out = logfactM((int) y)  + mu - y*log(mu);
  return(2*accu( b ));
}







// [[Rcpp::export]]
double objfun3(arma::mat Y, arma::mat MU, arma::vec Phi,
               arma::uvec t1, arma::uvec t2, arma::uvec t3,
               int msind,arma::mat naind){
  // arma::uvec t1 = find(familygroup==1),t2 = find(familygroup==2);
  // arma::uvec t3 = find(familygroup==3);
  double temp = 0;
  // Gaussian Likelihood
  if(t1.n_elem>0){
    temp = temp + gausslikehood(Y.cols(t1), MU.cols(t1), Phi(t1), naind.cols(t1));
    // cout << 'a';
  }
  if(t2.n_elem>0) {     // Logistic Likelihood
    temp = temp + logisticlik(Y.cols(t2), MU.cols(t2), naind.cols(t2));
    // cout << 'a';
    // x.cols(t2) = 1/(1+exp(-1*X.cols(t2)));
  }
  if(t3.n_elem>0) {     // Poisson Likelihood
    // cout<<t3;
    // cout << size(Y.cols(t3)) << size(MU.cols(t3)) << size(naind.cols(t3)) << endl;
    temp = temp + poissonlik(Y.cols(t3), MU.cols(t3), naind.cols(t3));
  }
  return(temp);
}


// [[Rcpp::export]]
double getdeviance(arma::mat Y, arma::mat MU, arma::vec Phi,
               arma::uvec t1, arma::uvec t2, arma::uvec t3,
               int msind,arma::mat naind){
  // arma::uvec t1 = find(familygroup==1),t2 = find(familygroup==2);
  // arma::uvec t3 = find(familygroup==3);
  double temp = 0;
  arma::mat b;
  // Gaussian Likelihood
  if(t1.n_elem>0){
    temp = temp + accu((conv_to< arma::vec >::from(sum(square((Y.cols(t1)-MU.cols(t1))%naind.cols(t1)),0) ))/Phi(t1));
    // cout << 'a';
  }
  if(t2.n_elem>0) {     // Logistic Likelihood
    temp = temp + 2*logisticlik(Y.cols(t2), MU.cols(t2), naind.cols(t2));
    // cout << 'a';
    // x.cols(t2) = 1/(1+exp(-1*X.cols(t2)));
  }
  if(t3.n_elem>0) {     // Poisson Likelihood
    // cout<<t3;
    // cout << size(Y.cols(t3)) << size(MU.cols(t3)) << size(naind.cols(t3)) << endl;
    temp = temp + poissondev(Y.cols(t3), MU.cols(t3), naind.cols(t3));
  }
  return(temp);
}






// [[Rcpp::export]]
arma::uvec mySdiff(arma::uvec x, arma::uvec y){
  // cout<<zeros<mat>(4,5);
  for (int j = 0; j < (int) y.n_elem; j++)
    x = x.elem(find(x != y(j)));
  return(x);
}





















// [[Rcpp::export]]
Rcpp::List gcure_cpp_init2(arma::mat Y, arma::mat X0, int rnk, arma::vec cindex,
                            arma::mat ofset, arma::vec familygroup,
                            arma::mat  Zini, arma::vec PhiIni,double kappaC0,
                            Rcpp::List control, int msind, arma::mat naind,
                            double ndev){
  bool converged=false; // , equalphi = true ; cObj = false,
  int pt = X0.n_cols, q = Y.n_cols,  maxit = control["initmaxit"]; // n = Y.n_rows,
  int equalphi = control["equalphi"];
  int cObj = control["objI"];
  // int p = pt - cindex.n_elem;
  Rcpp::List out;
  // int msind = miss['msind'];
  // arma::mat naind; naind.ones(n,q);
  // if(msind == 1) naind = miss['naind'];
  // here is columns of X0
  //   arma::mat Au= zeros<mat>(1,p),Av = zeros<mat>(1,q);
  //   arma::vec bu = zeros<vec>(1),bv= zeros<vec>(1);

  double epsilon = control["initepsilon"];  epsilon = epsilon*ndev;

  arma::vec cfamily = unique(familygroup);
  arma::uvec cIndex =  arma::conv_to< uvec >::from(cindex-1);
  arma::uvec cIndexC = mySdiff(linspace<uvec>(0,pt-1,pt), cIndex);
  arma::uvec t4=find(naind==0);
  // arma::uvec dInd= linspace<uvec>(0,n-1,1);
  //double kappaC0 = getKappaCo(X0,cfamily);
  //double kappaC0_ini = kappaC0;

  arma::mat MU(size(Y)), MU0(size(Y)),X2 = X0.cols(cIndexC)/kappaC0;
  arma::mat Res,Ct,X1 = X0.cols(cIndex)/kappaC0;
  // cout << 'a';

  arma::mat C = zeros<mat>(pt,q), X = X0.cols(cIndexC), X3 =X0/kappaC0;

  C.rows(cIndex) = Zini; C = kappaC0*C;
  arma::vec Phi = PhiIni, Phi2 = 1/PhiIni;

  rnk = rnk-1; // rank to extract

  // cout << 'a';
  wall_clock timer;




  // defining required variable used in the loop
  int j,iter;
  double elp,m1=1;
  arma::mat C_temp1=C.rows(cIndexC),C_temp2=C.rows(cIndex);
  arma::mat C_temp,Ut,Vt;
  arma::vec diffobj, obj,dt;
  arma::uvec t1 = find(familygroup==1),t2 = find(familygroup==2);
  arma::uvec t3 = find(familygroup==3);
  arma::vec qc =   conv_to< arma::vec >::from(sum( naind, 0 ) );

  double qg = accu(qc(t1));

  // cout << 'a';
  diffobj.zeros(maxit);obj.zeros(maxit+1);

  // MU = familyLinkinv(ofset + X0*(C/kappaC0),familygroup);
  // MU = familyLinkinv(ofset + X3*C,familygroup);
  MU = familyLinkinv3(ofset + X3*C,t1,t2,t3);
  Res = (Y - MU);
  if(msind == 1) Res.elem(t4).zeros();
  // if(msind == 1) Res = (Y - MU)%naind; else Res = (Y - MU);
  if(cObj!=0) // ********************
    obj(0) = objfun3(Y, MU, Phi,t1,t2,t3,msind,naind);

  // cout << 'a';

  timer.tic();
  for(iter = 1; iter < maxit; iter++){
    C_temp = C;

    // MU = familyLinkinv(ofset + X0*(C/kappaC0),familygroup);
    // if(msind == 1) Res = (Y - MU)%naind; else Res = (Y - MU); // **********
    C_temp1 = X2.t()*Res; C_temp1.each_row() %= Phi2.t();
    Ct = C.rows(cIndexC) +   C_temp1; // (X2.t()*Res)*diagmat(1/Phi);
    svd(Ut,dt,Vt,Ct);
    // cout<<size(Ct)<< std::endl;
    // cout<<Ut<< std::endl;
    // cout<<rnk<< std::endl;
    // cout<<dt<< std::endl;
    // cout<<size(Vt)<< std::endl;
    for(j = 0; j <= rnk; j++)
      Ut.col(j) = Ut.col(j)*dt(j);
    // cout<< dt << std::endl;
    C.rows(cIndexC) = Ut.cols(0,rnk)*Vt.cols(0,rnk).t();
    // cout << 'a';

    // Update Intercept and nonpenalized cofficients
    // MU = familyLinkinv(ofset + X0*(C/kappaC0),familygroup);
    // MU = familyLinkinv(ofset + X3*C,familygroup);
    MU = familyLinkinv3(ofset + X3*C,t1,t2,t3);
    Res = (Y - MU);
    if(msind == 1) Res.elem(t4).zeros();
    // if(msind == 1) Res = (Y - MU)%naind; else Res = (Y - MU);
    C_temp2 = X1.t()*Res; C_temp2.each_row() %= Phi2.t();
    C.rows(cIndex) = C.rows(cIndex) + C_temp2;   //X1.t()*Res*diagmat(1/Phi); // **** //t(crossprod(Res, X1)/Phi)  ##  crossprod(X0[, cIndex], Res)%*%diag(1/Phi) / kappaC0


    // update Phi
    // MU = familyLinkinv(ofset + X3*C,familygroup);
    MU = familyLinkinv3(ofset + X3*C,t1,t2,t3);
    Res = (Y - MU);
    if(msind == 1) Res.elem(t4).zeros();
    // if(msind == 1) Res = (Y - MU)%naind; else Res = (Y - MU);
    Phi.ones();
    if(equalphi == 1){
      if(t1.n_elem>0)  {
        m1 = accu(square(Res.cols(t1)))/qg; // (n*t1.n_elem);
        Phi.elem(t1) = m1*ones<vec>(t1.n_elem);
      }
    } else {
      if(t1.n_elem>0)  {
        Phi.elem(t1) = (conv_to< arma::vec >::from(sum(square(Res.cols(t1)),0))) /( qc(t1) );
      }
    }
    Phi2 = 1/Phi;


    obj(iter) = objfun3(Y, MU, Phi,t1,t2,t3,msind,naind);


    // obj(iter + 1) = objfun2(Y, MU, Phi,familygroup,msind,naind) + as_scalar(lam*wd*(wu.t()*abs(ue))*(wv.t()*abs(ve)));
    diffobj(iter) = absT((obj(iter) -obj(iter- 1)));
    // cout<<iter<< " " << diffobj(iter) << " " << epsilon << std::endl;
    // diffobj(iter) = square(norm( C-C_temp,"fro")/norm(C_temp,"fro")); //sum((C-C_temp)^2)/sum(C_temp^2)
    // diffobj(iter) = accu(square( C-C_temp))/accu(square(C_temp)); //sum((C-C_temp)^2)/sum(C_temp^2)
    // diffobj(iter) = accu(abs( MU-MU_temp))/accu(abs(MU_temp)); //sum((C-C_temp)^2)/sum(C_temp^2)
    // if(t1.n_elem>0) diffobj(iter) = diffobj(iter) + accu(square( C.cols(t1)-C_temp.cols(t1)))/accu(square(C_temp.cols(t1)));
    // if(t2.n_elem>0) diffobj(iter) = diffobj(iter) + accu(square( C.cols(t2)-C_temp.cols(t2)))/accu(square(C_temp.cols(t2)));
    // if(t3.n_elem>0) diffobj(iter) = diffobj(iter) + accu(square( C.cols(t3)-C_temp.cols(t3)))/accu(square(C_temp.cols(t3)));

    if (diffobj(iter) < epsilon ) {
      converged = true;
      break;
    }
  }

  elp = timer.toc();

  out["C"] = C/kappaC0;
  out["PHI"] = Phi;
  out["kappaC0"] = kappaC0;
  out["objval"] = obj;
  out["converged"] = converged;
  out["ExecTimekpath"] = elp;
  out["maxit"] = iter;
  out["converge"] = 0;
  if(converged) out["converge"] = 1;
  return(out);
}











// [[Rcpp::export]]
Rcpp::List gcure_cpp_miss(arma::mat Y, arma::mat Xm,int nlam, arma::vec cindex,
                          arma::mat ofset, arma::vec familygroup,
                          Rcpp::List initw, double Dini,
                          arma::mat  Zini, arma::vec PhiIni,
                          arma::mat  Uini, arma::vec Vini,
                          double kappaC0, double lmax, Rcpp::List control,
                          int msind, arma::mat naind, double ndev,int eea,
                          Rcpp::List zerosol,int maxit, double epsilon){
  arma::mat X0 = Xm, sd0 = stddev(Xm);
  sd0.elem( find(sd0 == 0) ).ones(); sd0.ones();
  sd0 = 1/sd0;
  X0.each_row() %= sd0;
  bool converged=false; // equalphi = true,  // **** cObj = true,
  int pt = X0.n_cols, q = Y.n_cols, n = Y.n_rows; // maxit = control["maxit"];
  int p = pt - cindex.n_elem;//maxit2 = maxit/10; //cout << 'mit ='<<maxit2;
  int cObj = control["objI"];
  Rcpp::List out;
  int equalphi = control["equalphi"];
  arma::mat Z00 = zerosol["Z"];
  arma::vec Phi00 = zerosol["PHI"];


  // double epsilon = control["epsilon"];
  epsilon = epsilon*ndev;
  // cout << epsilon << std::endl;
  // double gammaC0 = control["gammaC0"];
  double lmif = control["lamMinFac"];
  double lmaf = control["lamMaxFac"];
  double spu = control["spU"], spv = control["spV"];
  double wd =  initw["wd"], gamma0 = control["gamma0"];

  arma::vec wu = initw["wu"], wv = initw["wv"],cfamily = unique(familygroup);
  arma::uvec cIndex =  arma::conv_to< uvec >::from(cindex-1);
  arma::uvec cIndexC = mySdiff(linspace<uvec>(0,pt-1,pt), cIndex);
  arma::vec sd1  = conv_to< arma::vec >::from(sd0);
  // arma::uvec dInd= linspace<uvec>(0,n-1,1);

  //double kappaC0 = getKappaCo(X0,cfamily);
  double kappaC0_ini = kappaC0;
  Uini = Uini%(1/sd0(cIndexC)); wu = pow(abs(Uini),-1.0*gamma0);
  Zini.each_col() %= (1/sd0(cIndex));
  // Z00.each_col() %= (1/sd0(cIndex));

  arma::mat MU(size(Y)), MU0(size(Y)),X2 = X0.cols(cIndexC); //kappaC0;
  arma::mat Res,Ct,X1 = X0.cols(cIndex), X3 =X0/kappaC0; //X1 = X0.cols(cIndex)/kappaC0
  arma::mat C = zeros<mat>(pt,q), X = X0.cols(cIndexC);
  arma::mat X2X2 = X.t()*X/n;
  arma::uvec t1 = find(familygroup==1),t2 = find(familygroup==2);
  arma::uvec t3 = find(familygroup==3);
  if(eea==1)
    X2X2 = diagmat(ones<vec>(X.n_cols));


  C.rows(cIndex) = Z00; C.each_col() %= (1/sd0.t());
  C = kappaC0*C;
  arma::vec Phi = Phi00;
  arma::mat Phi2 = trans(1/Phi00);
  // arma::mat C_temp1=C.rows(cIndexC),C_temp2=C.rows(cIndex);

  // generate sequence of lambda
  double lmx = (lmaf*lmax);//elalpha;
  double lmn = lmax*lmif;
  arma::vec lamSeq =  exp(linspace<vec>(log(lmx),  log(lmn), nlam));



  // ------------
  // Initialization of ue and ve resoectively
  // MU0 = familyLinkinv(ofset + X0*(C/kappaC0),familygroup);
  // MU = familyLinkinv3(ofset + X3*C,t1,t2,t3);
  // MU0 = familyLinkinv3(ofset + X3*C,t1,t2,t3);  // X3 =X0/kappaC0;
  MU0 = familyLinkinv3(ofset + X0*(C/kappaC0),t1,t2,t3);
  arma::mat lamMat;
  if(msind == 1) lamMat  = X.t()*((Y-MU0)%naind); else lamMat  = X.t()*(Y-MU0);  // ***
  arma::mat ablamMat = abs(lamMat);
  arma::uvec indxx = find(ablamMat == ablamMat.max());indxx =indxx(0);


  double ob0 = objfun3(Y, MU0, Phi,t1,t2,t3,msind,naind);

  arma::vec ue,ve;
  ue.zeros(p);ve.zeros(q);
  // ue(as_scalar(indxx)%p) = 1; ve(as_scalar(indxx)/p) = 1; ve = ve*sign(lamMat(indxx));
  double de = 0; //kappaC0;
  C.rows(cIndexC) = de*(ue*ve.t());

  int chuk = 0, chvk = 0,ik = 0;
  // MU = familyLinkinv3(ofset + X0*(C/kappaC0),familygroup);
  MU = familyLinkinv3(ofset + X0*(C/kappaC0),t1,t2,t3);
  if(msind == 1) Res = (Y - MU)%naind; else Res = (Y - MU);
  // Res[is.na(Res)] = 0
  Ct = C.rows(cIndexC) +   X2.t()*Res*diagmat(1/(kappaC0*Phi)); //X2.t()*Res*diagmat(1/Phi);
  double xtxv = accu(square(ue)); //#*diag(q)
  arma::vec  xtyv = Ct.t()*ue; //crossprod(Ct,ue)
  double  xtxu = accu(square(ve)); //#*diag(p)
  arma::vec  xtyu = Ct*ve;
  // ------------

  arma::mat uklam,vklam,philam, BIClam,objval;
  arma::vec dklam,lselectSeq,indlam,execTime;
  arma::cube zpath,Ckpath,MUkpath;

  arma::mat objval1;objval1.zeros(5*(maxit+1),nlam+1);

  uklam.zeros(p,nlam+1);vklam.zeros(q,nlam+1);
  BIClam.zeros(4,nlam+1);philam.zeros(q,nlam+1);objval.zeros(maxit+1,nlam+1);
  zpath.zeros(pt-p,q,nlam+1);Ckpath.zeros(pt,q,nlam+1);
  MUkpath.zeros(n,q,nlam+1);
  dklam.zeros(nlam+1);  lselectSeq.zeros(nlam+1);indlam.zeros(nlam+1);
  execTime.zeros(nlam+1);

  arma::mat zb =ofset + X0.cols(cIndex)*(C.rows(cIndex)/kappaC0);
  MU = familyLinkinv3( zb,t1,t2,t3);
  double df = (pt-p)*q, qn = accu(naind);
  // double SSE = getsse2(Y, MU, Phi,familygroup,msind,naind); // ********************
  double SSE = getdeviance(Y, MU, Phi, t1,t2,t3,  msind,naind)/qn;



  BIClam(0,ik) = SSE + (df*log((double)qn))/(qn); //BIC
  BIClam(1,ik) = SSE + 2*df*log((double) pt*q)/(qn); //BICP
  BIClam(2,ik) = SSE + log(log( (double) qn))*df*log((double) pt*q)/(qn); //GIC
  BIClam(3,ik) = SSE + (2/qn)*(df); //AIC

  Ct = C/kappaC0; Ct.each_col() %= sd0.t();
  Ckpath.slice(ik) = Ct; MUkpath.slice(ik) = MU; //C/kappaC0
  zpath.slice(ik) = Ct.rows(cIndex); //C.rows(cIndex)/kappaC0;
  philam.col(ik) = Phi;

  wall_clock timer;


  // defining required variable used in the loop
  int ii=0,iter,mj;
  double lam,elp,fac;//,kappaC0_temp;
  double svk=0,suk=0,m1=1,de_temp;
  arma::mat C_temp,MU_temp;
  arma::vec diffobj, obj, obj2, plfacv, plfacu,Phi_temp,PhiI, ue_temp, ve_temp;
  arma::vec vest, uest, maxitc,convval;maxitc.zeros(nlam);convval.zeros(nlam);
  // arma::uvec t1 = find(familygroup==1),t2 = find(familygroup==2);
  arma::uvec t4=find(naind==0);
  double alp = control["elnetAlpha"];
  arma::vec qc =   conv_to< arma::vec >::from(sum( naind, 0 ) );
  double qg = accu(qc(t1));
  arma::mat facW = alp*wd*(wu*wv.t()),facL;
  arma::mat X2Y = X2.t()*Y,X1Y = X1.t()*Y,xuv(size(Y)),ResN;
  double kappaC0u, kappaC0v,kappaC0z = getKappaC0zv(X0.cols(cIndex), t1, t2, t3);

  // double ab = getKappaC0zv(X2, t1, t2, t3);
  // kappaC0z = getKappaC0zv(X1, t1, t2, t3);
  arma::vec tv2;
  kappaC0z = 1*getKappaC0zv(X1, t1, t2, t3);
  // kappaC0z = kappaC0;
  kappaC0z = 1*kappaC0z*kappaC0z;

  kappaC0u = 1*getKappaC0zv(X2, t1, t2, t3);
  // kappaC0u = kappaC0;
  kappaC0u = 1/(1*kappaC0u*kappaC0u);

  kappaC0v = 1*getKappaC0zv(X0.cols(cIndexC)*(Uini), t1, t2, t3);
  // kappaC0v = kappaC0;
  // kappaC0v = 1/(0.5* (double) n);
  kappaC0v = 1/(kappaC0v * kappaC0v);
  double k0v; // mfc = control["gammaC0"],


  // cout << 1/kappaC0u<< " " << 1/kappaC0v<< " "  << n/PhiIni.min()<< " " <<kappaC0*kappaC0 <<std::endl;

  for(ii=0; ii < nlam; ii++){
    lam = lamSeq(ii); //cout << ii<< std::endl;
    kappaC0 =  kappaC0_ini;
    facL = lam*facW;
    fac=(1-alp)*lam;
    diffobj.zeros(maxit);obj.zeros(maxit+1);obj2.zeros(5*(maxit+1));

    // set initialization for the optimization;
    ue = Uini; ve = Vini; de = Dini; //kappaC0*Dini;
    C.rows(cIndexC) = de*(ue*ve.t());
    C.rows(cIndex) = Zini;  // kappaC0*Zini;
    Phi = PhiIni;
    zb =ofset + X0.cols(cIndex)*(C.rows(cIndex));


    plfacv =  facL.t()*abs(ue);//alp*wv*(wu.t()*abs(ue))*wd ;
    plfacu = facL*abs(ve); //alp*wu*(wv.t()*abs(ve))*wd ;

    // MU = familyLinkinv(ofset + X0*(C/kappaC0),familygroup);
    // MU = familyLinkinv3(ofset + X3*C,t1,t2,t3); //X3 =X0/kappaC0;
    // MU = familyLinkinv3(ofset + X0*(C/kappaC0),t1,t2,t3);
    MU = familyLinkinv3(ofset + X0*C,t1,t2,t3);
    obj(0) = ob0;
    mj = 1;
    // if(cObj!=0)    // cout<<ue<< ve<<de <<std::endl;
    timer.tic();
    for(iter = 1; iter < maxit; iter++){
      C_temp = C;
      ue_temp = ue; ve_temp = ve; de_temp = de;
      MU_temp = MU;
      Phi_temp = Phi;
      Phi2 = trans(1/Phi);
      k0v = 1;
      // if(t1.n_elem > 0)
      k0v = 1/Phi.min();

      // kappaC0_temp = kappaC0;


      // obj2(5*(iter-1)+0) = objfun3(Y, MU, Phi,t1,t2,t3,msind,naind) +
      //   alp*as_scalar(de*lam*wd*(wu.t()*abs(ue))*(wv.t()*abs(ve))) +
      //   (1-alp)*lam*(de*de)*as_scalar((ue.t()*ue)*(ve.t()*ve)); // /(kappaC0*kappaC0);


      // kappaC0v = getKappaC0zv(X0.cols(cIndexC)*(ue/sqrt(Phi.min())), t1, t2, t3);
      // kappaC0v = kappaC0v*kappaC0v;
      // MU = familyLinkinv(ofset + X0*(C/kappaC0),familygroup);
      if(msind == 1) MU.elem(t4).zeros();
      Ct =(X2Y-X2.t()*MU); if(t1.n_elem>0) Ct.each_row() %= Phi2;
      // Ct = C.rows(cIndexC) +   Res;

      // update  ve
      // if(chuk==1){
      tv2 = ue*(kappaC0v/k0v);
      xtxv = accu(square(ue))*(kappaC0v/k0v);
      xtyv = de*ve + Ct.t()*tv2;
      plfacv =  facL.t()*abs(tv2);//alp*wv*(wu.t()*abs(ue))*wd ; alp*wd*(wu*wv.t())
      // }
      vest = softT(xtyv,plfacv)/(1+fac*xtxv);
      svk = norm( vest,2); //accu(abs(vest));
      if(svk==0){
        de = 0;ue.zeros(p);ve.zeros(q);
        // C.rows(cIndexC) = 0*C.rows(cIndexC);
        chvk = 0;break; // lselectSeq(0) = lam;
      } else {
        de = svk; //norm( vest,2);
        ve = vest/svk; //#print(de)
        chvk = 1;
      }
      // C.rows(cIndexC) = de*(ue*ve.t());
      MU = familyLinkinv3(zb + (X2*(de*ue))*(ve.t()),t1,t2,t3); //*****


      // calcObj  after ve update
      // obj2(5*(iter-1)+1) = objfun3(Y, MU, Phi,t1,t2,t3,msind,naind) +
      //   alp*as_scalar(de*lam*wd*(wu.t()*abs(ue))*(wv.t()*abs(ve))) +
      //   (1-alp)*lam*(de*de)*as_scalar((ue.t()*ue)*(ve.t()*ve));  //(kappaC0*kappaC0);

      // update ue
      // kappaC0u = getKappaC0zv(X2/sqrt(Phi.min()), t1, t2, t3);
      // kappaC0u = kappaC0u*kappaC0u;
      if(msind == 1) MU.elem(t4).zeros();
      // Ct =(X2Y-X2.t()*MU)/kappaC0u; if(t1.n_elem>0) Ct.each_row() %= Phi2;
      Ct =(X2Y-X2.t()*MU); if(t1.n_elem>0) Ct.each_row() %= Phi2;
      // Ct = C.rows(cIndexC) +   Res;
      // if(chvk==1){
      tv2 = ve*(kappaC0u/k0v);
      xtxu = kappaC0u/k0v;//accu(square(ve)); //#*diag(p) 1/kappaC0u
      xtyu = de*ue + Ct*(tv2);
      plfacu = facL*abs(tv2); //alp*wu*(wv.t()*abs(ve))*wd ; alp*wd*(wu*wv.t())
      // }
      uest = softT(xtyu,plfacu)/(1+fac*xtxu);
      suk = as_scalar(sqrt((uest.t()*X2X2*uest))); //accu(abs(uest));
      if(suk==0){
        de = 0;ue.zeros(p);ve.zeros(q);
        C.rows(cIndexC) = 0*C.rows(cIndexC);
        chuk = 0;break; // lselectSeq(0) = lam;
      } else {
        de = suk; //as_scalar(sqrt((uest.t()*X2X2*uest)));
        ue = uest/suk; //#print(de)
        chuk = 1;
      }

      C.rows(cIndexC) = de*(ue*ve.t());
      xuv = (X2*(de*ue))*(ve.t());
      MU = familyLinkinv3(zb +xuv,t1,t2,t3);


      // calcObj  after ue update
      // obj2(5*(iter-1)+2) = objfun3(Y, MU, Phi,t1,t2,t3,msind,naind) +
      //   alp*as_scalar(de*lam*wd*(wu.t()*abs(ue))*(wv.t()*abs(ve))) +
      //   (1-alp)*lam*(de*de)*as_scalar((ue.t()*ue)*(ve.t()*ve));

      // kappaC0z = getKappaC0zv(X1/sqrt(Phi.min()), t1, t2, t3);
      // kappaC0z = kappaC0z*kappaC0z;
      if(msind == 1) MU.elem(t4).zeros();
      Res =(X1Y-X1.t()*MU)/(kappaC0z*k0v); if(t1.n_elem>0) Res.each_row() %= Phi2;
      C.rows(cIndex) = C.rows(cIndex) +   Res;
      zb = ofset + X1*C.rows(cIndex);
      MU = familyLinkinv3(zb + xuv,t1,t2,t3); //***********


      // calcObj  after Z update
      // obj2(5*(iter-1)+3) = objfun3(Y, MU, Phi,t1,t2,t3,msind,naind) +
      //   alp*as_scalar(de*lam*wd*(wu.t()*abs(ue))*(wv.t()*abs(ve))) +
      //   (1-alp)*lam*(de*de)*as_scalar((ue.t()*ue)*(ve.t()*ve));//(kappaC0*kappaC0);

      if(msind == 1) MU.elem(t4).zeros();
      Res = (Y - MU);
      Phi.ones();
      if(equalphi == 1){
        if(t1.n_elem>0)  {
          m1 = accu(square(Res.cols(t1)))/qg; // (n*t1.n_elem);
          Phi.elem(t1) = m1*ones<vec>(t1.n_elem);
        }
      } else {
        if(t1.n_elem>0)  {
          Phi.elem(t1) = (conv_to< arma::vec >::from(sum(square(Res.cols(t1)),0))) /( qc(t1) );
        }
      }


      if((cObj!=0) ){ //&&  (iter>2)
        obj(iter) = objfun3(Y, MU, Phi,t1,t2,t3,msind,naind) +
          alp*as_scalar(de*lam*wd*(wu.t()*abs(ue))*(wv.t()*abs(ve))) +
          (1-alp)*lam*(de*de)*as_scalar((ue.t()*ue)*(ve.t()*ve)); //(kappaC0*kappaC0);
        diffobj(iter) = ((obj(iter) -obj(iter- 1)));
        // if(diffobj(iter) < 0){
        //   k0v = k0v/mfc;
        //   mj =1;
        //   // kappaC0 <- kappaC0/gammaC0
        // } else {
        //   ue = ue_temp; ve = ve_temp; de = de_temp;
        //   MU = MU_temp;
        //   Phi = Phi_temp;
        //   k0v = 1;
        //   iter = iter - 1;
        //   mj = mj +1;
        //   // if(mj > 50){
        //   //   cout << "Increase initial kappaC0"<< std::endl;
        //   //   break;
        //   // }
        // }
      } else {
        if(t1.n_elem>0) diffobj(iter) = diffobj(iter) + accu(square( C.cols(t1)-C_temp.cols(t1)))/accu(square(C_temp.cols(t1)));
        if(t2.n_elem>0) diffobj(iter) = diffobj(iter) + accu(square( C.cols(t2)-C_temp.cols(t2)))/accu(square(C_temp.cols(t2)));
        if(t3.n_elem>0) diffobj(iter) = diffobj(iter) + accu(square( C.cols(t3)-C_temp.cols(t3)))/accu(square(C_temp.cols(t3)));
      }

      // obj2(5*(iter-1)+4) = obj(iter);

      // obj(iter + 1) = objfun2(Y, MU, Phi,familygroup,msind,naind) + as_scalar(lam*wd*(wu.t()*abs(ue))*(wv.t()*abs(ve)));
      // cout<<iter<< std::endl;
      // diffobj(iter) = square(norm( C-C_temp,"fro")/norm(C_temp,"fro")); //sum((C-C_temp)^2)/sum(C_temp^2)
      // diffobj(iter) = accu(square( C-C_temp))/accu(square(C_temp)); //sum((C-C_temp)^2)/sum(C_temp^2)
      // diffobj(iter) = accu(abs( MU-MU_temp))/accu(abs(MU_temp)); //sum((C-C_temp)^2)/sum(C_temp^2)
      // if(t1.n_elem>0) diffobj(iter) = diffobj(iter) + accu(square( C.cols(t1)-C_temp.cols(t1)))/accu(square(C_temp.cols(t1)));
      // if(t2.n_elem>0) diffobj(iter) = diffobj(iter) + accu(square( C.cols(t2)-C_temp.cols(t2)))/accu(square(C_temp.cols(t2)));
      // if(t3.n_elem>0) diffobj(iter) = diffobj(iter) + accu(square( C.cols(t3)-C_temp.cols(t3)))/accu(square(C_temp.cols(t3)));

      if (absT(diffobj(iter)) < epsilon ) {
        converged = true;
        convval(ii) = 1;
        break;
      }
    }
    maxitc(ii) = iter;



    elp = timer.toc();
    //  add solution path information
    // if (svk > 0)
    //   if( suk > 0){
    indlam(ii) = ii;
    uklam.col(ik) = ue%sd1(cIndexC);
    vklam.col(ik) = ve;
    philam.col(ik) = Phi;
    dklam(ik) = de; //kappaC0;
    C.rows(cIndexC) = de*(ue*ve.t());
    if(eea==1){
      m1 = norm(uklam.col(ik),2);
      uklam.col(ik) = uklam.col(ik)/m1;
      dklam(ik) = m1*dklam(ik);
    }
    lselectSeq(ik) = lam;
    objval.col(ik) = obj;//
    objval1.col(ik) = obj2;//
    execTime(ik) = elp;


    Ct = C; //kappaC0;
    Ct.each_col() %= sd0.t();
    zpath.slice(ik) = Ct.rows(cIndex); //C.rows(cIndex)/kappaC0;
    Ckpath.slice(ik) = Ct; // C/kappaC0;

    MU = familyLinkinv3(ofset + X0*C,t1,t2,t3);
    // MU = familyLinkinv(ofset + X0*(C/kappaC0),familygroup);
    // SSE = getsse2(Y, MU, Phi,familygroup,msind,naind); // ********************
    SSE = getdeviance(Y, MU, Phi, t1,t2,t3,  msind,naind)/qn;
    df = nzcount(ue) + nzcount(ve)  +(pt-p)*q -1;

    BIClam(0,ik) = SSE + (df*log((double)qn))/(qn); //BIC
    BIClam(1,ik) = SSE + 2*df*log((double) pt*q)/(qn); //BICP
    BIClam(2,ik) = SSE + log(log( (double) qn))*df*log((double) pt*q)/(qn); //GIC
    BIClam(3,ik) = SSE + (2/qn)*(df); //AIC
    MUkpath.slice(ik) = MU;
    ik = ik+1;
    // }

    if( (iter>1) && ((nzcount(ue) > (p*spu)) || (nzcount(ve) > (q*spv)))  ) {break;}
  }
  ik = ik-1;

  out["ukpath"] = uklam.cols(0,ik); // ----
  out["vkpath"] = vklam.cols(0,ik);
  out["dkpath"] = arma::conv_to<arma::vec>::from(dklam.head(ik+1)) ;  // ----
  out["phipath"] = philam.cols(0,ik);
  out["zpath"] = arma::conv_to<arma::cube>::from(zpath.head_slices(ik+1));  //.subcube(0,0,1,  pt-p,q,ik);//    //
  out["mukpath"] = arma::conv_to<arma::cube>::from(MUkpath.head_slices(ik+1));  //.subcube(0,0,1,  pt-p,q,ik);//    //
  out["ICKpath"] = BIClam.cols(0,ik);
  out["nkpath"] = ik+1;
  out["kappaC0"] = kappaC0;
  out["objkval"] = objval.cols(0,ik);
  out["objkval1"] = objval1.cols(0,ik);
  out["Ckpath"] = arma::conv_to<arma::cube>::from(Ckpath.head_slices(ik+1)); //.subcube(0,0,1,  pt,q,ik);  //
  // out["converged"] = converged;
  // cout << lselectSeq << std::endl;
  out["lamKpath"] = arma::conv_to<arma::vec>::from(lselectSeq.head(ik+1));
  out["ExecTimekpath"] = arma::conv_to<arma::vec>::from(execTime.head(ik+1));
  out["lamseq"] = lamSeq;
  out["maxit"] = maxitc;
  out["converge"] = convval;
  return(out);
}





















