void Model::initialize_eta(){
  fixed_lps = mm(ab_node->data, ab_node->bVals);
  int nRow = fixed_lps.length();
  int nCol = n_gmus;
  rand_lps = NumericMatrix(nRow, nCol);
  
  std::vector<gammaNode*>* this_gamma_vec;
  gammaNode* this_gamma;
  NumericVector this_lp;
  for(int i = 0; i < n_gmus; i++){
    this_gamma_vec = &gamma_vec[i];
    int n_gammas = this_gamma_vec->size();
    for(int j = 0; j < n_gammas; j++){
      this_gamma = (*this_gamma_vec)[j];
      this_lp = mm(this_gamma->data, this_gamma->vals);
      insertVecByIndex(this_lp, rand_lps, this_gamma->rows, j);
    }
  }
  combined_lps = NumericVector(nRow);
  for(int i = 0; i < nRow; i++){
    combined_lps[i] = fixed_lps[i];
    for(int j = 0; j < nCol; j++){
      combined_lps[i] += rand_lps(i,j);
    }
  }
  etas = NumericVector(nRow);
  for(int i = 0; i < nRow; i++){ etas[i] = exp(combined_lps[i]); }
}

void Model::update_eta_fixed(){
  int nRow = combined_lps.size();
  for(int i = 0; i < nRow; i++){ combined_lps[i] -= fixed_lps[i]; }
  fixed_lps = mm(ab_node->data, ab_node->bVals);
  for(int i = 0; i < nRow; i++){ 
    combined_lps[i] += fixed_lps[i];
    etas[i] = exp(combined_lps[i]);
  }
}

void Model::update_eta_rand(int gamma_num, int gamma_ind){
//  Rcout << "gamma l1 = " << gamma_num << " gamma l2 = " << gamma_ind << "\n";
  gammaNode* this_gamma = gamma_vec[gamma_num][gamma_ind];
  int k = this_gamma->rows.size();
  int this_row;
  
  for(int i = 0; i < k; i++){
    this_row = this_gamma->rows[i] - 1;
    combined_lps[this_row] -= rand_lps(this_row, gamma_num);
  }
  NumericVector this_lp = mm(this_gamma->data, this_gamma->vals);
  for(int i = 0; i < k; i++){
    this_row = this_gamma->rows[i] - 1;
    combined_lps[this_row] += this_lp[i];
    rand_lps(this_row, gamma_num) = this_lp[i];
    etas[this_row] = exp(combined_lps[this_row]);
  }
}