void insertVecByIndex(NumericVector vals, 
                      NumericMatrix mat_into,
                      IntegerVector row_index, 
                      int col){
  int k = row_index.length();
  for(int i = 0; i < k; i++){
    mat_into(row_index[i] - 1, col) = vals[i];
  }
}

/***
void addByIndex(NumericMatrix vals, 
              NumericVector vec_into,
              IntegerVector row_index, 
              int col){
  int k = row_index.length();
  int this_row;
  for(int i = 0; i < k; i++){
    this_row = row_index[i] - 1;
    vec_into[this_row] += vals(this_row, col);
  }
}

void subtractByIndex(NumericMatrix vals, 
                     NumericVector vec_into,
                     IntegerVector row_index, 
                     int col){
  int k = row_index.length();
  int this_row;
  for(int i = 0; i < k; i++){
    this_row = row_index[i] - 1;
    vec_into[this_row] -= vals(row_index[i], col);
  }
}
***/

NumericVector mm(NumericMatrix m, NumericVector v){
  int nRow = m.rows();
  int nCol = m.cols();
  int v_len = v.length();
  if(v_len != nCol){ stop("matrix dimensions do not align"); }
  NumericVector ans(nRow);
  double v_i;
  for(int i = 0; i < nCol; i++){
    v_i = v[i];
    for(int j = 0; j < nRow; j++){
      ans[j] += m(j,i) * v_i;
    }
  }
  return(ans);
}

Function get_a_fxn(List abList){
  Function ans = abList["alpha_fxn"];
  return(ans);
}

Function get_b_fxn(List abList){
  Function ans = abList["beta_fxn"];
  return(ans);
}

Function get_gmu_fxn(List abList){
  Function ans = abList["g_fxn"];
  return(ans);
}

Function get_S_fxn(List abList){
  Function ans = abList["s_fxn"];
  return(ans);
}

Function get_g_fxn(List gammaList){
  Function ans = gammaList["g_fxn"];
  return(ans);
}