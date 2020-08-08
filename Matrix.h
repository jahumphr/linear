#include <stdlib.h>
#include <vector>
#include <stdexcept>
#include <iostream>

template<typename T>
class Matrix{
  public:
    Matrix(const int rows, const int cols){
      _numRows = rows;
      _numCols = cols;
      for(int i = 0; i < rows; ++i)
        _matrix.push_back(std::vector<T>(cols));
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

    Matrix<T> operator*(const Matrix& other){
      //return the product of two matrices
      
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
    friend Matrix<T> operator*(const T num, Matrix<T> m){

      for (int i = 0; i < m._numRows; ++i){
        for (int j = 0; j < m._numCols; ++j){
          m._matrix[i][j] *= num;
        }
      }

      return m;
    }

    //matrix * num
    friend Matrix<T> operator*(Matrix<T> m, const T num){

      for (int i = 0; i < m._numRows; ++i){
        for (int j = 0; j < m._numCols; ++j){
          m._matrix[i][j] *= num;
        }
      }

      return m;
    }

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