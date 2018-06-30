class alphaBetaNode : public Node{
public:
  int n_alpha;
  int n_beta;
  void setVals(NumericVector new_vals);
  double calc_lp();
  NumericVector aVals;
  NumericVector bVals;
  Rcpp::Function alpha_logDens;
  Rcpp::Function beta_logDens;
  alphaBetaNode(NumericVector r_aVals, NumericVector r_bVals, 
                NumericMatrix Rdata, Function a_prior, 
                Function b_prior);
  
  NumericMatrix data;
  
  ~alphaBetaNode(){}
};



alphaBetaNode::alphaBetaNode(NumericVector r_aVals, NumericVector r_bVals, 
                             NumericMatrix Rdata, Function a_prior, 
                             Function b_prior) :
  alpha_logDens(a_prior), 
  beta_logDens(b_prior){
  aVals = clone(r_aVals);
  bVals = clone(r_bVals);
  n_alpha = aVals.length();
  n_beta = bVals.length();
  data = clone(Rdata);
  vals = NumericVector(n_alpha + n_beta);
  for(int i = 0; i < n_alpha; i++){
    vals[i] = aVals[i];
  }
  for(int i = 0; i < n_beta; i++){
    vals[i + n_alpha] = bVals[i];
  }
}

void alphaBetaNode::setVals(NumericVector new_vals){
  vals = new_vals;
  int input_length = new_vals.length();
  if(input_length != (n_alpha + n_beta)){
    Rcpp::stop("wrong dimensions for setting alphaBeta node\n");
  }
  for(int i = 0; i < n_alpha; i++){
    aVals[i] = new_vals[i];
  }
  for(int i = 0; i < n_beta; i++){
    bVals[i] = new_vals[n_alpha + i];
  }
}

double alphaBetaNode::calc_lp(){
  NumericVector aDens = alpha_logDens(aVals);
  NumericVector bDens = beta_logDens(bVals);
  double ans = aDens[0] + bDens[0];
  return(ans);
}


