
// This is the main class for a neural network. Generally, it is best to use this class through the Network interface.

#ifndef __NEURAL_NETWORK
#define __NEURAL_NETWORK

#include <iostream>

#include "Error.h"
#include "Network.h"

using namespace std;



class Neural_network : public Network {

  void train();
  double evaluate();
  
  virtual void load();
  virtual void save();

 public:
  
  Neural_network();
  ~Neural_network();
  
  
  

  
 private:
  
  



};


#endif
