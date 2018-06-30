class gammaNode : public Node{
public:
  void setVals(NumericVector new_vals);
  double calc_lp();
  gmusNode* gmus;
  Function logDens;
  
  IntegerVector rows;
  NumericMatrix data;
  
  gammaNode(NumericVector r_vals, Function logDens, List data_info, 
            gmusNode* gmusPtr);
};

gammaNode::gammaNode(NumericVector r_vals, Function logDens, List data_info, 
                     gmusNode* gmusPtr) : 
  logDens(logDens) {
  gmus = gmusPtr;
  setVals(r_vals);
  NumericMatrix r_data = data_info["data"];
  NumericVector r_rows = data_info["rows"];
  data = clone(r_data);
  rows = clone(r_rows);
}


double gammaNode::calc_lp(){
  NumericVector ans = logDens(vals, 
                              gmus->mu_vals,
                              gmus->sqrt_S_vals);
  return(ans[0]);
}

void gammaNode::setVals(NumericVector new_vals){
  int n_newVals = new_vals.length();
  int n_oldVals = vals.length();
  if(n_newVals != n_oldVals){
    if(n_oldVals > 0){
      stop("vector size for gammaNode does not match current size\n");
    }
    vals = NumericVector(n_newVals);
  }
  for(int i = 0; i < n_newVals; i++){ vals[i] = new_vals[i]; }
}