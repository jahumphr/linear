#include <iostream>
#include <vector>
#include <stdlib.h>
#include "Matrix.h"
#include <time.h>
#include <unistd.h>
using namespace std;

Matrix<int> randomMatrix(int row, int col){
  srand ((unsigned) time(0));
  Matrix<int> m(row, col);
  for (int i = 0; i < row; ++i){
    for (int j = 0; j < col; ++j){
      
      m(i, j) = rand() % 5;
    }
  }

  return m;
}

int main(){
  // Matrix<int> test = randomMatrix(2,3);
  // Matrix<int> test2 = randomMatrix(3,1);
  // cout << " -------- test" << endl;
  // test.display();
  // cout << " -------- test 2" << endl;
  // test2.display();
  // Matrix<int> result = test * test2;
  // cout << " ----------------- result " << endl;
  // result.display();

  
  // Matrix<int> test(2,3);
  // cout << test.sizeOf() << endl;
  // test.at(99,999) = 5;
  Matrix<int> matrix1 = randomMatrix(2, 3);
  cout << "m1" << endl;
  matrix1.display();
  cout << endl;

  //makes system pause for a second before randomizing next matrix, ensuring time seed changes from matrix 1 to matrix 2
  usleep(1234000); 
  Matrix<int> scalar = 2 - matrix1;
  scalar.display();

  Matrix<int> matrix2 = randomMatrix(2, 3);
  cout << "m2" << endl;
  matrix2.display();

  Matrix<int> result = matrix1 + matrix2;

  result.display();




  return 0;

}




