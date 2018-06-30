#include <Rcpp.h>
using namespace Rcpp;

#include "code/baseClasses.cpp"
#include "code/utilities.cpp"
#include "code/alphaBetaNode.cpp"
#include "code/gmusNode.cpp"
#include "code/gammaNode.cpp"
#include "code/dataNode.cpp"
#include "code/model.cpp"
#include "code/etaTools.cpp"


// [[Rcpp::export]]
SEXP cMakeModel(List initValues, 
                List initData, 
                List initPriors){
  XPtr<Model> ptr(new Model(initValues, 
                            initData, 
                            initPriors));
  return ptr;
}

// [[Rcpp::export]]
NumericVector c_get_ab_vals(SEXP model){
  XPtr<Model> ptr(model);
  NumericVector ans = ptr->ab_node->getVals();
  return(ans);
}

// [[Rcpp::export]]
NumericVector c_get_gmus_vals(SEXP model, int ind){
  XPtr<Model> ptr(model);
  int n_vec = ptr->gmus_vec.size();
  if(ind < 1 || ind > n_vec){
    stop("ind must be in 1:(number of random effects)");
  }
  NumericVector ans = ptr->gmus_vec[ind - 1]->getVals();
  return(ans);
}

// [[Rcpp::export]]
NumericVector c_get_gamma_vals(SEXP model, int ind1, int ind2){
  XPtr<Model> ptr(model);
  int n_vec1 = ptr->gamma_vec.size();
  if(ind1 < 1 || ind1 > n_vec1){
    stop("ind1 must be in 1:(number of random effects)");
  }
  int n_vec2 = ptr->gamma_vec[ind1-1].size();
  if(ind2 < 1 || ind2 > n_vec2){
    stop("ind must be in 1:(number of unique ids)");
  }
  NumericVector ans = ptr->gamma_vec[ind1-1][ind2-1]->getVals();
  return(ans);
}

// [[Rcpp::export]]
void c_set_ab(SEXP model, NumericVector new_vals){
  XPtr<Model> ptr(model);
  int n = ptr->ab_node->getSize();
  int new_n = new_vals.size();
  if(n != new_n){
    stop("incorrect length for new_vals");
  }
  ptr->set_ab(new_vals);
}

// [[Rcpp::export]]
double c_setAndCalc_ab(SEXP model, NumericVector new_vals){
  XPtr<Model> ptr(model);
  int n = ptr->ab_node->getSize();
  int new_n = new_vals.size();
  if(n != new_n){
    stop("incorrect length for new_vals");
  }
  ptr->set_ab(new_vals);
  
  double ans = ptr->ab_node->calc_lp();
  ans += ptr->computeLLK();
  
  return(ans);
}

// [[Rcpp::export]]
void c_set_gmus(SEXP model, NumericVector new_vals, int ind){
  XPtr<Model> ptr(model);
  int ind_use = ind - 1;
  int n_gmus = ptr->gmus_vec.size();
  if(ind < 1 || ind > n_gmus){ stop("ind out of range"); }
  gmusNode* this_node = ptr->gmus_vec[ind_use];
  int n = ptr->gmus_vec[ind_use]->vals.size();
  int new_n = new_vals.size();
  if(n != new_n){ stop("incorrect length for new_vals");}
  this_node->setVals(new_vals);
}

// [[Rcpp::export]]
double c_setAndCalc_gmus(SEXP model, NumericVector new_vals, int ind){
  XPtr<Model> ptr(model);
  int ind_use = ind - 1;
  int n_gmus = ptr->gmus_vec.size();
  if(ind < 1 || ind > n_gmus){ stop("ind out of range"); }
  gmusNode* this_node = ptr->gmus_vec[ind_use];
  int n = ptr->gmus_vec[ind_use]->vals.size();
  int new_n = new_vals.size();
  if(n != new_n){ stop("incorrect length for new_vals");}
  this_node->setVals(new_vals);
  double ans = this_node->calc_lp();
  int n_gamma = ptr->gamma_vec[ind_use].size();
  for(int i = 0; i < n_gamma; i++){
    ans += ptr->gamma_vec[ind_use][i]->calc_lp();
  }
  return(ans);
}

// [[Rcpp::export]]
void c_set_gamma(SEXP model, NumericVector new_vals, int ind1, int ind2){
  XPtr<Model> ptr(model);
  int ind1_use = ind1 - 1;
  int ind2_use = ind2 - 1;
  int n_gamma1 = ptr->get_n_gmus();
  if(ind1 < 1 || ind1 > n_gamma1){ stop("ind1 out of range"); }
  int n_gamma2 = ptr->get_n_gamma(ind1_use);
  if(ind2 < 1 || ind2 > n_gamma2){ stop("ind2 out of range"); }
  ptr->set_gamma(new_vals, ind1_use, ind2_use);
}


// [[Rcpp::export]]
double c_setAndCalc_gamma(SEXP model, NumericVector new_vals, int ind1, int ind2){
  XPtr<Model> ptr(model);
  int ind1_use = ind1 - 1;
  int ind2_use = ind2 - 1;
  int n_gamma1 = ptr->get_n_gmus();
  if(ind1 < 1 || ind1 > n_gamma1){ stop("ind1 out of range"); }
  int n_gamma2 = ptr->get_n_gamma(ind1_use);
  if(ind2 < 1 || ind2 > n_gamma2){ stop("ind2 out of range"); }
  ptr->set_gamma(new_vals, ind1_use, ind2_use);
  gammaNode* this_node = ptr->gamma_vec[ind1_use][ind2_use];
  double ans = this_node->calc_lp();
  ans += ptr->computeLLK(this_node->rows);
  return(ans);
}

// [[Rcpp::export]]
int n_gmus(SEXP model){
  XPtr<Model> ptr(model);
  int ans = ptr->gmus_vec.size();
  return(ans);
}

// [[Rcpp::export]]
int n_gamma(SEXP model, int gamma_l1){
  int n_l1 = n_gmus(model);
  if(gamma_l1 < 1 || n_l1 < gamma_l1){ stop("gamma_l1 index out of range"); }
  XPtr<Model> ptr(model);
  int ans = ptr->gamma_vec[gamma_l1 - 1].size();
  return(ans);
}

// [[Rcpp::export]]
double c_full_lp(SEXP model){
  XPtr<Model> ptr(model);
  double ans = ptr->ab_node->calc_lp();
  int n_gmus = ptr->gmus_vec.size();
  int n_gamma;
  for(int i = 0; i < n_gmus; i++){
    ans += ptr->gmus_vec[i]->calc_lp();
    n_gamma = ptr->gamma_vec[i].size();
    for(int j = 0; j < n_gamma; j++){
      ans += ptr->gamma_vec[i][j]->calc_lp();
    }
  }
  ans += ptr->computeLLK();
  return(ans);
}