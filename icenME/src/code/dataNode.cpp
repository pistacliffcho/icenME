
class dataNode{
public:
  NumericMatrix respMat;
  Function computeLLK; 
  /**
   * computeLLK should take the following arguments:
   *    y, eta, alpha, isCen, inds = (with default all)
   */
  dataNode(NumericMatrix resp, 
           Function computeLLK):
    computeLLK(computeLLK){ 
    respMat = clone(resp);
  }
    
};

