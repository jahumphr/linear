# linear
C++ Linear Algebra Library

CREATING A MATRIX

The Matrix class currently only supports doubles.

	Matrix<double> matrix;

Matrix can be populated in multiple ways. Simplest method is using a string. Syntax 
requires each matrix value to be separated by a space, and rows are delineated by 
commas. Space after comma optional.

	matrix = "1 2 3, 4 5 6, 7 8 9";
	
Also supports negative numbers:

	matrix = "1 -2 3, -4 -5 6, 7 8 -9";
	
If rows of different lengths are entered, exception thrown:
	
	matrix = "1 2 3, 4 5 6 7, 8 9 10";
	
	terminate called after throwing an instance of 'std::invalid_argument'
		what():  Number of columns don't match.
  

Alternatively, matrix can be populated by creating individual vector rows, and
combining using insertRow().

Syntax:

insertRow(std::vector<T> newRow,
			int position
			);
			
If no position argument, row is appended as last row. Otherwise, row inserted at 
position argument:

	//create vectors
	vector<double> one{3, 2, 0, 1};
    vector<double> two{4, 0, 1, 2};
    vector<double> three{3, 0, 2, 1};
	
	//build matrix
	matrix.insertRow(one);
    matrix.insertRow(two);
    matrix.insertRow(three, 0); //insert as first row


First vector inserted sets row length. Any subsequent insertions that don't match
initial length will throw an exception.

insertCol() can also be used, and works identically.

Identity matrix can be quickly created with identity(size):
	
	matrix = matrix.identity(4);
	matrix.display()
	
Output:

	( 1 0 0 0 )
	( 0 1 0 0 )
	( 0 0 1 0 )
	( 0 0 0 1 )



DISPLAYING MATRIX

matrix.display() will display the matrix.

	Matrix<double> matrix;
	matrix = "1 2 3, 4 5 6, 7 8 9";
	matrix.display():
	
Output:

	( 1 2 3 )
	( 4 5 6 )
	( 7 8 9 )
	


MATRIX OPERATIONS

==

Returns true if matrix1 == matrix2.

		  ( 1 2 3 )
matrix1 = ( 4 5 6 )
		  ( 7 8 9 )
		  
		  ( 1 2 3 )
matrix2 = ( 4 5 6 )
		  ( 7 8 9 )

	return matrix1 == matrix2;

Output:

	true
	

=

Sets matrix1 = matrix2.


*

When used with a matrix and int or double, performs scalar multiplication:

	scalarM = 3 * matrix;

Or

	scalarM = matrix * 3;
	
When used with two matrices, returns product of matrices:

	matrix3 = matrix1 * matrix2;


matrix.det() -- returns determinant of matrix. Throws exception if matrix is non-square.

matrix.inverse() -- returns inverse of a matrix. Throws exception if matrix is not invertible.

matrix.ref() -- returns row echelon form of a matrix.

matrix.rref() -- returns reduced row echelon form of a matrix.

matrix.rank() -- returns rank of a matrix.
	
