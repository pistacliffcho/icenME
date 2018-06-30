class Node{
public: 
  NumericVector vals;
  int getSize(){return(vals.size());}
  virtual void setVals(NumericVector new_vals) = 0;
  NumericVector getVals(){return(clone(vals));}
  virtual double calc_lp() = 0;
/*
    double calc_lp(NumericVector new_vals){
    NumericVector old_vals = getVals();
    setVals(new_vals);
    NumericVector ans = calc_lp();
    setVals(old_vals);
    return(ans);
  }
 */
};

