#include <stdlib.h>
#include <vector>
#include <stdexcept>
#include <iostream>

namespace lina{
template<typename T>
  class Matrix{
  public:
    Matrix(const int rows, const int cols){
      _numRows = rows;
      _numCols = cols;
      for(int i = 0; i < rows; ++i)
        _matrix.push_back(std::vector<T>(cols));
    }

    Matrix(){
        _numRows = 0;
        _numCols = 0;
    }

    bool operator==(const Matrix& other) const{
      if(_numRows != other._numRows || _numCols != other._numCols)
        return false;

      for(int i = 0; i < _numRows; ++i)
        for(int j = 0; j < _numCols; ++j)
          if(_matrix(i,j) != other(i,j))
            return false;
            
      return true;
    }

    T& at(const int row, const int col){
      if (row < 0 || row >= _numRows)
        throw std::invalid_argument( "Invalid row number" );
      
      if(col < 0 || col >= _numCols)
        throw std::invalid_argument( "Invalid column number" );

      return _matrix[row][col];
    }

    T& operator()(const int row, const int col){
      return _matrix[row][col];
    }

    std::vector<T> getRow(const int row) const{
      return _matrix[row];
    }

    std::vector<T> getCol(const int colNum) const{
      std::vector<T> col;
      for (int i = 0; i < _numRows; ++i){
        col.push_back(_matrix[i][colNum]);
      }
      return col;
    }

    //insert row at position. If no position entered, inserts as last row
    void insertRow(std::vector<T> newRow, int pos = -1){
      if (pos < -1 || pos > _numRows)
        throw std::invalid_argument( "Invalid row position" );
      if (_matrix.size() == 0){
          _matrix.push_back(newRow);
          _numRows++;
          _numCols = newRow.size();
          return;
      }
      if (newRow.size() != _matrix[0].size())
        throw std::invalid_argument("New row invalid length");

      if (pos == -1){
        _matrix.push_back(newRow);
        _numRows++;
      }
  //      else{
  //        _matrix.insert(pos, newRow);
  //        _numRows++;
  //      }
    }

    //insert column at position. If no position entered, inserts as last column.
    //Worst case runtime is quadratic here.
    void insertCol(std::vector<T> newCol, int pos = -1){
      if (pos < 0 || pos > _numCols)
        throw std::invalid_argument("Invalid column position");
      if (newCol.size() != _matrix.size())
        throw std::invalid_argument("New column invalid length");

      if (pos == -1){
        for (int i = 0; i < _numRows; ++i)
          _matrix[i].push_back(newCol[i]);
        _numCols++;
      }
      else{
        for (int i = 0; i < _numRows; ++i)
          _matrix[i].insert(pos, newCol[i]);
        _numCols++;
      }
    }

    Matrix<T> Transpose(){
      Matrix<T> product(_numCols, _numRows);
      for(int i=0; i<_numRows; i++){
        for(int j=0; j<_numCols; j++){
          product(j,i) = _matrix[i][j];
        }
      }
      return product;
    }

    //return the determinant of a matrix
    double det(){
      if (_numCols != _numCols)
        throw std::invalid_argument("Cannot find the determinant of a non-square matrix.");
      
      double intermediate1=1;
      double intermediate2=1;
      double result;

      for(int i=0; i<_numCols; i++){
        intermediate1 *= _matrix[i][i];
      }

      for(int j=0; j<_numCols; j++){
        intermediate2 *= _matrix[j][_numCols-j-1];
      }
      result = intermediate1-intermediate2;
      return result;
    }
    //return the inverse of a matrix
    // Matrix<T> Inverse(){
    //    Matrix<T>product(_numRows, _numCols);
    //   double determinant = _matrix.det();
    //   if (determinant == 0 || _numCols != _numRows)
    //     throw std::invalid_argument("Matrix is not invertible.");

    //   if(_numCols == 2){
    //     product(0,0) = _matrix[1][1]/determinant;
    //     product(1,0) = -1*_matrix[1][0]/determinant;
    //     product(0,1) = -1*_matrix[0][1]/determinant;
    //     product(1,1) = _matrix[0][0]/determinant;
    //   }
    //   else{
    //     for(i = 0; i < _numCols; i++){
    //       for(j = 0; j < _numRows; j++){
    //         double intermediate1 = 1;
    //         double intermediate2 = 1;
    //         for(int k=1; k<_numCols; k++){
    //           intermediate1 *= _matrix[(i+k)%_numRows][(j+k)%_numCols];
    //           intermediate2 *= _matrix[(i+((k+1)%(_numCols-1))+1)%_numCols][(j+((k+2)%(_numCols-1))+1)%_numCols];
    //           /*((mat[(j+1)%3][(i+1)%3] * mat[(j+2)%3][(i+2)%3]) - (mat[(j+1)%3][(i+2)%3] * mat[(j+2)%3][(i+1)%3]))/ determinant<<"\t";*/
    //         if((i+j)%2==0)
    //           product(i,j) = (intermediate1-intermediate2)/determinant;
    //         else
    //           product(i,j) = -1*(intermediate1-intermediate2)/determinant;
    //         }
    //       }
    //   }
      
    // }
    // }
    /**
      Returns the sum of two matrices, throws error if matrices are of different dimensions.
    */
    Matrix<T> operator+(const Matrix& other){
      if (_numRows != other._numRows || _numCols != other._numCols)
        throw std::invalid_argument("Cannot add matrices of different dimensions.");

      Matrix<T> sum(_numRows, _numCols);
      for (int i = 0; i < _numRows; ++i){
        for (int j = 0; j < _numCols; ++j){
          sum(i,j) = _matrix[i][j] + other._matrix[i][j];
        }
      }

      return sum;   
    }

    /**
      Returns the difference between two matrices, throws error if matrices are of different dimensions.
    */
    Matrix<T> operator-(const Matrix& other){
      if (_numRows != other._numRows || _numCols != other._numCols)
        throw std::invalid_argument("Cannot subtract matrices of different dimensions.");
      
      Matrix<T> difference(_numRows, _numCols);
      for (int i = 0; i < _numRows; ++i){
        for (int j = 0; j < _numCols; ++j){
          difference(i,j) = _matrix[i][j] - other._matrix[i][j];
        }
      }

      return difference;
    }

    /**
      Returns the product of two matrices (are we doing auto transpose later or nah? Also might need to throw error if dimensions don't work, ie. 2x2 * 2x3 ok but 2x2 * 3x2 not ok.
    */
    Matrix<T> operator*(const Matrix& other){
      Matrix<T> product(_numRows, other._numCols);
      for (int i = 0; i < _numRows; i++){
        for (int j = 0; j < _numCols; j++){
          for(int k = 0; k < other._numRows; k++){
              product(i,j) = product(i,j) + (_matrix[i][k] * other._matrix[k][j]);
          }
        }
      }
      return product;
    }

    /* SCALAR OPERATIONS */

    //num + matrix
    /**
      Don't think we can add scalars to matrices like this... if we had to the scalar would be added to the matrix after being multiplied by the identity matrix (same goes for subtraction)
    */
    friend Matrix<T> operator+(const T num, Matrix<T> m){

      for (int i = 0; i < m._numRows; ++i){
        for (int j = 0; j < m._numCols; ++j){
          m._matrix[i][j] += num;
        }
      }

      return m;
    }

    //matrix + num
    friend Matrix<T> operator+(Matrix<T> m, const T num){

      for (int i = 0; i < m._numRows; ++i){
        for (int j = 0; j < m._numCols; ++j){
          m._matrix[i][j] += num;
        }
      }

      return m;
    }

    //num - matrix
    friend Matrix<T> operator-(const T num, Matrix<T> m){

      for (int i = 0; i < m._numRows; ++i){
        for (int j = 0; j < m._numCols; ++j){
          m._matrix[i][j] = num - m._matrix[i][j];
        }
      }

      return m;
    }

    //matrix - num
    friend Matrix<T> operator-(Matrix<T> m, const T num){

      for (int i = 0; i < m._numRows; ++i){
        for (int j = 0; j < m._numCols; ++j){
          m._matrix[i][j] -= num;
        }
      }

      return m;
    }

    //num * matrix
    /**
      Returns the scalar multiplication of a constant number and a matrix.
      */
    friend Matrix<T> operator*(const T num, Matrix<T> m){

      for (int i = 0; i < m._numRows; ++i){
        for (int j = 0; j < m._numCols; ++j){
          m._matrix[i][j] *= num;
        }
      }

      return m;
    }

    //matrix * num
    /**
      Returns the scalar multiplication of a matrix and a constant number (scalar multiplication is commutative).
      */
    friend Matrix<T> operator*(Matrix<T> m, const T num){

      for (int i = 0; i < m._numRows; ++i){
        for (int j = 0; j < m._numCols; ++j){
          m._matrix[i][j] *= num;
        }
      }

      return m;
    }
    
    //reduced row echelon form
    void rref() {
      int swapWith; //stores index of row x > i, where m[x][j] != 0, to switch with current row where m[i][j] == 0.
      for (int j = 0; j < _numCols; ++j) {
          for (int i = 0; i < _numRows; ++i) {
              if ( j >= _numCols)
                  return;

              while (_matrix[i][j] == 0) {
                  //if column is all zeros, move to next column
                  if (_allZeros(i, j, swapWith)) {
                      j++;
                      if (j >= _numCols)
                          return;
                      continue;
                  }
                  //if element on current row is 0, attempt interchange with row > current row where element in column j is non-zero.
                  if (swapWith != -1)
                      _matrix[i].swap(_matrix[swapWith]);
                  else {
                      j++;
                      if (j >= _numCols)
                          return;
                  }
              }

              scaleRow(i, j);
              subtractRows(i, j);
              j++;

          }
      }
  }

    private:
    //private helper function for rref, checks if current column is already
    //all zeros. Returns true if all zeros, false otherwise.
    bool _allZeros(int row, int col, int& swapWith){
      bool allZero = true;
      for (int i = 0; i < _numRows; ++i){
        if (_matrix[i][col] != 0){
          allZero = false;
          //if i > row and element in i is 0, swapWith assigned with row to interchange
          if (i > row){
            swapWith = i;
            return allZero;
          }

        }
      }
      swapWith = -1;
      return allZero;
    }

    void scaleRow(int row, int col){
      double scalar = (double)_matrix[row][col];
      
      for (int i = col; i < _numCols; ++i)
        _matrix[row][i] /= scalar;

    }

    void subtractRows(int row, int col) {
      double scalar;
      for (int i = 0; i < _numRows; ++i) {
          if (i != row) {
              scalar = _matrix[i][col];
              for (int j = 0; j < _numCols; ++j) {
                  _matrix[i][j] = _matrix[i][j] - _matrix[row][j] * scalar;
              }
          }
      }
  }

    public:

    void display() const{
      //display
      for(int i=0; i<_numRows; i++){
        std::cout << "( ";

        for(int j=0; j<_numCols; j++)
          std::cout << _matrix[i][j] << " ";
          
        std::cout << ")" << std::endl;
      }
      }
    
    int sizeOf() const{
      return _numRows * _numCols;
    }

  private:
    std::vector<std::vector<T>> _matrix;
    int _numRows, _numCols;


  };
}