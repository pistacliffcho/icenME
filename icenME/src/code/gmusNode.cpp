class gmusNode : public Node{
public:
  int n_gamma;
  void setVals(NumericVector new_vals);
  void setVals(NumericVector gmu, NumericMatrix sqrt_S);
  double calc_lp();
  
  NumericVector mu_vals;
  NumericMatrix sqrt_S_vals;
  
  Function mu_prior;
  Function sqrt_S_prior;
  
  gmusNode(Function r_gmu_prior, Function r_sqrt_S, 
           NumericVector r_gmu_vals, NumericMatrix r_sqrt_S_vals);
}; 

gmusNode::gmusNode(Function r_gmu_prior, Function r_sqrt_S, 
                   NumericVector r_gmu_vals, NumericMatrix r_sqrt_S_vals) :
  mu_prior(r_gmu_prior),
  sqrt_S_prior(r_sqrt_S){
  setVals(r_gmu_vals, r_sqrt_S_vals);
  n_gamma = r_gmu_vals.length();
}

double gmusNode::calc_lp(){
  NumericVector mu_lp = mu_prior(mu_vals);
  NumericVector sqrt_S_lp = sqrt_S_prior(sqrt_S_vals);
  double ans = mu_lp[0] + sqrt_S_lp[0];
  return(ans);
}

void gmusNode::setVals(NumericVector gmu, NumericMatrix sqrt_S){
  mu_vals = clone(gmu);
  sqrt_S_vals = clone(sqrt_S);
  n_gamma = mu_vals.size();
  int vals_size = n_gamma + n_gamma * n_gamma;
  if(vals.length()!= vals_size){ vals = NumericVector(vals_size); }  
  for(int i = 0; i < n_gamma; i++){
    vals[i] = mu_vals[i];
    for(int j = 0; j < n_gamma; j++){
      vals[i + (j+1) * n_gamma] = sqrt_S(i,j);
    }
  }
}

void gmusNode::setVals(NumericVector new_vals){
  for(int i = 0; i < n_gamma; i++){
    mu_vals[i] = new_vals[i];
    for(int j = 0; j < n_gamma; j++){
      sqrt_S_vals(i,j) = new_vals[i + (j+1) * n_gamma];
    }
  }
  int n = new_vals.size();
  for(int i = 0; i < n; i++){
    vals[i] = new_vals[i];
  }
}