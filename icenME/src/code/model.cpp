class Model{
public:
  alphaBetaNode* ab_node;
  void ab_setVals(NumericVector new_vals){ ab_node->setVals(new_vals); }
  NumericVector ab_getVals(){return(ab_node->getVals());}
  std::vector<gmusNode*> gmus_vec;
  std::vector<std::vector<gammaNode*> > gamma_vec;
  
  NumericVector etas;
  NumericVector combined_lps;
  NumericVector fixed_lps;
  NumericMatrix rand_lps;
  
  void initialize_eta();
  void update_eta_fixed();
  void update_eta_rand(int gamma_num, int gamma_ind);
  
  int n_gmus;
  
  dataNode* data_llk;
  LogicalVector isCen;
  
  Model(List initValues, List initData, List initPriors);
  double computeLLK();
  double computeLLK(IntegerVector inds);
  
  void set_ab(NumericVector new_vals){
    ab_node->setVals(new_vals);
    update_eta_fixed();
  }
  void set_gmus(NumericVector new_vals, int ind){
    gmus_vec[ind]->setVals(ind);
  }
  void set_gamma(NumericVector new_vals, int ind1, int ind2){
    gammaNode* this_gamma = gamma_vec[ind1][ind2];
    this_gamma->setVals(new_vals);
    update_eta_rand(ind1, ind2);
  }
  
  int get_n_gmus(){
    int ans = gmus_vec.size();
    return(ans);
  }
  
  int get_n_gamma(int l1){
    int ans = gamma_vec[l1].size();
    return(ans);
  }
  
  ~Model(){}
};

double Model::computeLLK(){
  NumericVector ans = data_llk->computeLLK(data_llk->respMat, 
                                      etas, ab_node->aVals, 
                                      isCen);
  return(ans[0]);
}

double Model::computeLLK(IntegerVector inds){
  NumericVector ans = data_llk->computeLLK(data_llk->respMat, 
                                      etas, ab_node->aVals, 
                                      isCen, inds);
  return(ans[0]);
}


Model::Model(List initValues, List initData, List initPriors){
  NumericVector a_inits = initValues["alpha"];
  NumericVector b_inits = initValues["beta"];
  NumericMatrix fixed_covs = initData["fixed_cov"];
  NumericMatrix respMat = initData["resp"];
  LogicalVector r_isCen = initData["isCen"];
  Function a_prior = initPriors["alpha"];
  Function b_prior = initPriors["beta"];
  Function r_llk = initPriors["data_llk"];

  isCen = clone(r_isCen);
  
  ab_node = new alphaBetaNode(a_inits, b_inits, fixed_covs, a_prior, b_prior);
  List this_list;
  
  List gamma_mu_priors = initPriors["gamma_mu"];
  List gamma_mu_vals = initValues["gamma_mu"];
  List sqrt_S_priors = initPriors["sqrt_S"];
  List sqrt_S_vals = initValues["sqrt_S"];

  n_gmus = gamma_mu_vals.length();
  gmus_vec.resize(n_gmus);
  gamma_vec.resize(n_gmus);
  
  List gamma_priors = initPriors["gamma"];
  List gamma_vals = initValues["gamma"];
  List gamma_data = initData["rand_cov"];  
  
  int n_gammas;
  for(int i = 0; i < n_gmus; i++){
    Function this_gmu_prior = gamma_mu_priors[i];
    Function this_S_prior = sqrt_S_priors[i];
    NumericVector this_gmu_vals = gamma_mu_vals[i];
    NumericMatrix this_S_vals = sqrt_S_vals[i];
    gmus_vec[i] = new gmusNode(this_gmu_prior, this_S_prior, 
                               this_gmu_vals, this_S_vals);

    List these_gamma_vals = gamma_vals[i];
    List these_gamma_data = gamma_data[i];
    n_gammas = these_gamma_vals.length();
    gamma_vec[i].resize(n_gammas);
    List these_gamma_fxns = gamma_priors[i];
    for(int j = 0; j < n_gammas; j++){
      Function this_fxn = these_gamma_fxns[j];
      NumericVector these_values = these_gamma_vals[j];
      List these_data_info = these_gamma_data[j];
      gamma_vec[i][j] = new gammaNode(these_values, 
                                      this_fxn, 
                                      these_data_info, 
                                      gmus_vec[i]);
    }
  } 
  initialize_eta();
  data_llk = new dataNode(respMat, r_llk);
  isCen = isCen;
}
