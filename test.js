class Matrix {
    //remember Array[row][column]
    body = [];
    swaps = 0;
    constructor(matrix) {
        for(var i = 0; i < matrix.length; i++) {
            let temp = []
            for(var j = 0; j < matrix[i].length; j++) {
                temp.push(matrix[i][j]);
            }
            this.body.push(temp);
        }
    }

    print() {
        //console.log(this.body);
        for(let i = 0; i < this.body.length; i++) {
            let row = "";
            for(let j = 0; j < this.body[i].length; j++) {
                row += this.body[i][j] + " ";
            }
            console.log(row);
        }
    }

    rowSwap(a,b) {
        var temp = this.body[a];
        this.body[a] = this.body[b];
        this.body[b] = temp;
    }

    rowToEnd(row) {
        this.body.push(this.body.splice(row,1)[0]);
    }

    rowMultiply(row, scalar) {
        for(var i = 0; i < this.body[row].length; i++) {
            this.body[row][i] *= scalar;
        }
    }

    colMultiply(col, scalar) {
        let A = Matrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {
            A.body[i][col] *= scalar;
        }
        return A;
    }

    rowAdd(tRow,sRow,scalar) {
        for(var i = 0; i < this.body[tRow].length; i++) {
            var temp = this.body[sRow][i] * scalar;
            this.body[tRow][i] += temp;
        }
    }

    getColumnVector(column) {
        let c = [];
        let A = Matrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {
            let c1 = [];
            c1.push(A.body[i][column]);
            c.push(c1);
        }
        //yes this returns a matrix. a column vector isn't actually a vector you fuckhead
        return new Matrix(c);
    }

    getColumnVectorAsVector(column) {
        let c = [];
        let A = Matrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {
            c.push(A.body[i][column]);
        }
        return new Vector(c);
    }

    getRowVector(row) {
        //yea this ones a vector :)
        return new Vector(Matrix.copy(this).body[row]);
    }

    getRowVectorAsMatrix(row) {
        let arr = Matrix.copy(this).body[row];
        return new Matrix([arr]);
    }

    rref() {
        let A = Matrix.copy(this);
        //A.print();
        let pivotPositions = [0,0];
        let pivotCount = 0;
        //go column by column
        //remember that there can be up to m pivots in a m*n matrix
        for(var c = 0; c < A.body[0].length; c++) {
            //initial variables
            let length = A.body.length; //row count
            //make an array from the column
            let columnArray = [];
            for(var i = 0; i < length; i++) {
                columnArray.push(A.body[i][c]);
            }
            //console.log(columnArray);
            //ensure the column is not all zeros
            //remove duplicate indicies
            let columnSet = new Set(columnArray);
            let shortColumnArray = [...columnSet];
            //check if the column is all zeros
            if (!(shortColumnArray[0] == 0 && shortColumnArray.length == 1)) {
                //now we identify specific points (pivot points) where we want to scale the row so this point becomes 1
                //pivotPositions[0] refers to the row position, pivotPositions[1] refers to the column
                //if pivot count is the size of rows do not apply the transformation
                if (pivotCount < length) {
                    //first identify a point of interest
                    let row = pivotPositions[0];
                    let column = pivotPositions[1];
                    let num = A.body[row][column];
                    if (num != 0) {
                        A.rowMultiply(row, 1/num);
                        pivotPositions[0]++;
                        pivotPositions[1]++;
                        pivotCount++;
                    } else {
                        //this generally assumes we will find a non-negative (existing) index
                        let index = -1;
                        let searchArray = columnArray.slice(row); //get the array starting at the target row
                        for(var i = searchArray.length-1; i > -1; i--) { //search backwards
                            //if a number is nonzero, set it to the index
                            index = searchArray[i] != 0 ? i : index;
                        }
                        if (index != -1) {
                            let targetRow = columnArray.indexOf(searchArray[index], row);
                            A.rowSwap(targetRow, row);
                            //update numbers
                            num = A.body[row][column];
                            A.rowMultiply(row, 1/num);
                            pivotPositions[0]++;
                            pivotPositions[1]++;
                            pivotCount++;
                        } else {
                            pivotPositions[1]++;
                        }
                    }
                    //zero out other elements in the column
                    //console.log(row);
                    for(var i = 0; i < length; i++) {
                        //dont zero out the prime index
                        if (i != row) {
                            A.rowAdd(i, row, -1*A.body[i][c]);
                        }
                    }
                    //A.print();
                }
            } else {
                pivotPositions[1]++;
            }
        }
        return A;
    }

    kernel() {
        let A = Matrix.copy(this).transpose();
        let I = Matrix.identity(A.body.length);
        //A.print();
        let pivotPositions = [0,0];
        let pivotCount = 0;
        //go column by column
        //remember that there can be up to m pivots in a m*n matrix
        for(var c = 0; c < A.body[0].length; c++) {
            //initial variables
            let length = A.body.length; //row count
            //make an array from the column
            let columnArray = [];
            for(var i = 0; i < length; i++) {
                columnArray.push(A.body[i][c]);
            }
            //console.log(columnArray);
            //ensure the column is not all zeros
            //remove duplicate indicies
            let columnSet = new Set(columnArray);
            let shortColumnArray = [...columnSet];
            //check if the column is all zeros
            if (!(shortColumnArray[0] == 0 && shortColumnArray.length == 1)) {
                //now we identify specific points (pivot points) where we want to scale the row so this point becomes 1
                //pivotPositions[0] refers to the row position, pivotPositions[1] refers to the column
                //if pivot count is the size of rows do not apply the transformation
                if (pivotCount < length) {
                    //first identify a point of interest
                    let row = pivotPositions[0];
                    let column = pivotPositions[1];
                    let num = A.body[row][column];
                    if (num != 0) {
                        I.rowMultiply(row, 1/num);
                        A.rowMultiply(row, 1/num);
                        pivotPositions[0]++;
                        pivotPositions[1]++;
                        pivotCount++;
                    } else {
                        //this generally assumes we will find a non-negative (existing) index
                        let index = -1;
                        let searchArray = columnArray.slice(row); //get the array starting at the target row
                        for(var i = searchArray.length-1; i > -1; i--) { //search backwards
                            //if a number is nonzero, set it to the index
                            index = searchArray[i] != 0 ? i : index;
                        }
                        if (index != -1) {
                            let targetRow = columnArray.indexOf(searchArray[index], row);
                            I.rowSwap(targetRow, row);
                            A.rowSwap(targetRow, row);
                            //update numbers
                            num = A.body[row][column];
                            I.rowMultiply(row, 1/num);
                            A.rowMultiply(row, 1/num);
                            pivotPositions[0]++;
                            pivotPositions[1]++;
                            pivotCount++;
                        } else {
                            pivotPositions[1]++;
                        }
                    }
                    //zero out other elements in the column
                    //console.log(row);
                    for(var i = 0; i < length; i++) {
                        //dont zero out the prime index
                        if (i != row) {
                            I.rowAdd(i, row, -1*A.body[i][c]);
                            A.rowAdd(i, row, -1*A.body[i][c]);
                        }
                    }
                    //A.print();
                }
            } else {
                pivotPositions[1]++;
            }
        }
        //remove rows that correspond to nonzero
        let vectors = [];
        for(var i = 0; i < A.body.length; i++) {
            let set = new Set();
            for(var j = 0; j < A.body[i].length; j++) {
                set.add(A.body[i][j]);
            }
            if (set.has(0) && set.size == 1) {
                vectors.push(new Vector(I.body[i]));
            }
        } 
        return vectors;
    }

    add(matrix) {
        let A = Matrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {
            for(var j = 0; j < A.body[i].length; j++) {
                A.body[i][j] += matrix.body[i][j];
            }
        }
        return A;
    }

    subtract(matrix) {
        let A = Matrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {
            for(var j = 0; j < A.body[i].length; j++) {
                A.body[i][j] -= matrix.body[i][j];
            }
        }
        return A;
    }

    multiplyConstant(c) {
        let A = Matrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {
            for(var j = 0; j < A.body[i].length; j++) {
                A.body[i][j] *= c;
            }
        }
        return A;
    }

    multiplyMatrix(matrix) {
        //dot the row of this matrix with the column of the specified matrix
        let A = Matrix.copy(this);
        let B = Matrix.copy(matrix);
        //we execute as many multiplications as the amount of rows in the first matrix times the amount of columns in the second matrix
        //when finding a column vector times its transpose we should get a 3x3 matrix
        let m = [];
        for(var i = 0; i < A.body.length; i++) {
            //row
            let newRow = [];
            for(var j = 0; j < B.body[0].length; j++) {
                //column
                //now get the row and column vectors
                let r = A.getRowVector(i);
                let c = B.getColumnVectorAsVector(j);
                //now multiply them together and push to the new row array
                let d = r.dot(c);
                newRow.push(d);
            }
            //newRow has the amount of elements it needs
            //push to the matrix
            m.push(newRow);
        }
        return new Matrix(m);
    }

    static copy(matrix) {
        let arr = [];
        for(var i = 0; i < matrix.body.length; i++) {
            let row = [];
            for(var j = 0; j < matrix.body[i].length; j++) {
                row.push(matrix.body[i][j]);
            }
            arr.push(row);
        }
        return new Matrix(arr);
    }

    transpose() {
        let cols = [];
        for(var i = 0; i < this.body[0].length; i++) { //column
            let col = [];
            for(var j = 0; j < this.body.length; j++) { //row
                col.push(this.body[j][i]);
            }
            cols.push(col);
        }
        return new Matrix(cols);
    }

    removeRow(row) {
        let A = Matrix.copy(this);
        let matrix = [];
        for(var i = 0; i < A.body.length; i++) {
            if (i != row) {
                matrix.push(A.body[i]);
            }
        }
        return new Matrix(matrix);
    }

    removeCol(col) {
        let A = Matrix.copy(this);
        let matrix = [];
        for(var i = 0; i < A.body.length; i++) {
            let row = [];
            for(var j = 0; j < A.body[i].length; j++) {
                if (j != col) {
                    row.push(A.body[i][j]);
                }
            }
            matrix.push(row);
        }
        return new Matrix(matrix);
    }

    setRow(row, newRow) {
        let A = Matrix.copy(this);
        A.body[row] = newRow;
        return A;
    }

    setCol(col, newCol) {
        let A = Matrix.copy(this);
        A = A.transpose();
        A = A.setRow(col, newCol);
        A = A.transpose();
        return A;
    }

    charpoly() {
        let A = Matrix.copy(this);
        let length = A.body.length;
        let M = [];
        let c = [];
        M[0] = Matrix.copy(A).multiplyConstant(0); //zero matrix
        c[length] = 1;
        for(var i = 1; i <= length; i++) {
            M[i] = A.multiplyMatrix(M[i-1]).add(Matrix.identity(length).multiplyConstant(c[length-i+1]));
            c[length-i] = (-1/i)*(A.multiplyMatrix(M[i]).trace());
        }
        return c.reverse();
    }

    trace() {
        let tr = 0;
        for(var i = 0; i < this.body.length; i++) {
            tr += this.body[i][i];
        }
        return tr;
    }

    nullspace() {
        let A = Matrix.copy(this);
        let vectors = [];
        let R = A.rref();
        if (A.rref().equals(Matrix.identity(A.body.length))) {
            let arr = [];
            for(var i = 0; i < A.body.length; i++) {
                arr.push(0);
            }
            return new Vector(arr);
        } else {
            A = A.transpose();
            let length = A.body.length; //row count
            //we now have two matricies, now we have to append the identity matrix to the right of the original matrix
            let mArray = [];
            for(var i = 0; i < A.body.length; i++) { //row iterator
                let row = [];
                for(var j = 0; j < (A.body[0].length+length); j++) { //columns
                    if (j < A.body[0].length) {
                        row.push(A.body[i][j]);
                    } else {
                        let k = Math.abs(A.body[0].length - j);
                        //console.log(k);
                        row.push(k == i ? 1 : 0);
                    }
                }
                mArray.push(row);
            }
            let M = new Matrix(mArray);
            let R = M.rref();
            //identify zero rows, these are the nullspace vectors
            for(var i = 0; i < R.body.length; i++) {
                let zero = R.body[i].slice(0,length);
                let set = new Set(zero);
                let arr = [...set];
                if (arr.length == 1 && arr[0] == 0) {
                    //remember, this is in zyx instead of xyz, so reverse the vector
                    let vector = R.body[i].splice(length);
                    vector.reverse();
                    vectors.push(new Vector(vector));
                }
            }
        }
        //if we have no vectors, we have a problem
        //remember, if a vector has more pivots than columns there will be residual columns that correspond to the
        //null space, nth empty columns within the reduced echelon form represent the nth variable
        if (vectors.length == 0) {
            //find the last pivot point
            let lastPivotColumn = -1;
            for(var i = R.body[0].length-1; i > -1; i--) {
                let row = R.body.length-1; //start from zero, not one
                if (R.body[row][i] == 1) {
                    lastPivotColumn = i;
                } 
            }

            //now find empty columns and add them to the the vectors
            for(var i = 0; i < R.body[0].length; i++) {//column iteration
                //create column and convert it to a set
                let column = [];
                for(var j = 0; j < R.body.length; j++) { //rows
                    column.push(R.body[j][i]);
                }
                let set = new Set(column);
                let array = [...set];
                if (array.length == 1 && array[0] == 0) {
                    //all zeros
                    let index = i; //nth column = nth variable in the vector
                    let vector = [];
                    for(var k = 0; k < R.body[0].length; k++) { //vector size will always be of the amount of columns
                        vector.push(k == i ? 1 : 0);
                    }
                    vectors.push(new Vector(vector));
                } else {
                    //see if we're past the last pivot
                    if (i > lastPivotColumn) {
                        //if we are, get the index we're on and negate the whole column
                        let vector = [];
                        for (var j = 0; j < R.body.length; j++) {//rows
                            R.body[j][i] *= -1;
                            //now add it to our vector
                            //for the index, the first occurence in the row is the vector index
                            let index = R.body[j].indexOf(1);
                            vector[index] = R.body[j][i];
                        }
                        //i represents the index we need to make one, so we push column size - vector size zeroes to the vector
                        for(var j = vector.length; j < R.body[0].length; j++) {
                            vector.push(0);
                        }
                        //now we set the ith index of the vector to 1
                        vector[i] = 1;
                        vectors.push(new Vector(vector));
                    }
                }
            }
        }
        return vectors;
    }

    inverse() {
        if (this.toComplexMatrix().isSingular()) {
            return null;
        }
        if (this.body.length == this.body[0].length) {
            let length = this.body.length;
            let matrix = [];
            for(var i = 0; i < length; i++) {
                let row = [];
                for(var j = 0; j < length*2; j++) {
                    if (j < length) {
                        row.push(this.body[i][j]);
                    } else {
                        let k = j-length;
                        row.push(k == i ? 1 : 0);
                    }
                }
                matrix.push(row);
            }
            let n = new Matrix(matrix);
            //console.log(n);
            n = n.rref();
            //console.log(n);
            let m = [];
            for(var i = 0; i < length; i++) {
                let r = [];
                for(var j = length; j < length*2; j++) {
                    r.push(n.body[i][j]);
                }
                m.push(r);
            }
            return new Matrix(m);
        }
        return null;
    }

    ref() {
        let m = this;
        if (m.body.length == m.body[0].length) {
            for(var c = 0; c < m.body[0].length; c++) {
                let column = [];
                for(var i = 0; i < m.body.length; i++) {
                    //create column list
                    //console.log(m.body[i][c]);
                    column.push(m.body[i][c]);
                }
                //console.log(column);
                for(var i = 0; i < m.body.length; i++) { //iterate through each row, and subtract the number to get rid of with the reciprocal of the number above
                    if (i > c && m.body[i][c] != 0) { //lower triangular matrix portion, check if not equal to zero to shorten computation time
                        m.rowAdd(i, c, -1*m.body[i][c]*(1/m.body[c][c]));
                    } else {
                        //if our row we want is zero, swap with nonzero row
                        if ((i >= c && m.body[i][c] == 0) && (i != this.body.length-1)) {
                            m.rowSwap(i,i+1);
                            m.swaps++;
                        }
                    }
                }
            }
        }

        return m;
    }

    det() {
        let A = Matrix.copy(this);
        let eigen = A.eigen()[0];
        let c = eigen[0];
        for(var i = 1; i < eigen.length; i++) {
            c = c.multiply(eigen[i]);
        }
        return c;
    }

    power(n) {
        let A = Matrix.copy(this);
        let B = A;
        for(var i = 0; i < n-1; i++) {
            B = B.multiplyMatrix(A);
        }
        return B;
    }

    static randomize(low,high, z) {
        let m = Matrix.null(z);
        for(var i = 0; i < m.body.length; i++) {
            for(var j = 0; j < m.body[i].length; j++) {
                let n = Math.random() * (high-low)+low;
                m.body[i][j] = n;
            }
        }
        return m;
    }

    static null(n) {
        let ar = [];
        for(var i = 0; i < n; i++) {
            let r = [];
            for(var j = 0; j < n; j++) {
                r.push(0);
            }
            ar.push(r);
        }
        return new Matrix(ar);
    }
 
    static identity(n) {
        let ar = [];
        for(var i = 0; i < n; i++) {
            let r = [];
            for(var j = 0; j < n; j++) {
                r.push(j != i ? 0 : 1);
            }
            ar.push(r);
        }
        return new Matrix(ar);
    }

    equals(matrix) {
        for(var i = 0; i < this.body.length; i++) {
            for(var j = 0; j < this.body[i].length; j++) {
                if (this.body[i][j] != matrix.body[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    columnVectorMagnitude() {
        //only works if its a column vector
        let A = Matrix.copy(this).transpose();
        let v = new Vector(A.body[0]);
        return v.magnitude();
    }

    identityReplace(n) {
        //n dictates how much of the identity matrix to preserve
        let A = Matrix.copy(this);
        for(var i = 0; i < n; i++) {
            console.log(Vector.identity(A.body.length,i).vector);
            A = A.setCol(i,Vector.identity(A.body.length,i).vector);
            A = A.setRow(i,Vector.identity(A.body.length,i).vector);
        }
        return A;
    }

    identityMerge(n) {
        //n dictates how many rows and columns to add to the matrix
        let A = Matrix.copy(this);
        let I = Matrix.identity(n + A.body.length);
        for(var i = n; i < I.body.length; i++) {
            for(var j = n; j < I.body.length; j++) {
                I.body[i][j] = A.body[i-n][j-n];
            }
        }
        return I;
    }

    minor(n) {
        let A = Matrix.copy(this);
        for (var i = n - 1; i >= 0; i--) {
            A = A.removeRow(i);
            A = A.removeCol(i);
        }
        return A;
    }

    householder() {
        /*
        let m = new Matrix([[12,-51,4],[6,167,-68],[-4,24,-41]]);
        let x = m.getColumnVector(0); //the arbitrary column vector
        let a = Math.abs(x.columnVectorMagnitude()); //alpha
        let e1 = Vector.identity(3,0); //e1 must be represented as a column vector as well
        e1 = e1.toMatrix().transpose(); //column vector transformation
        let u = x.subtract(e1.multiplyConstant(a)); //u
        let v = u.multiplyConstant(1/u.columnVectorMagnitude()); //unit vector version of u
        let vt = v.transpose(); //transpose v into a row vector vt
        let Q = [];
        Q.push(Matrix.identity(3).subtract(v.multiplyMatrix(vt).multiplyConstant(2)));
        let qa = Q[0].multiplyMatrix(m);
        let m11 = qa.minor(1);
        x = m11.getColumnVector(0);
        a = Math.abs(x.columnVectorMagnitude());
        e1 = Vector.identity(2,0);
        e1 = e1.toMatrix().transpose();
        u = x.subtract(e1.multiplyConstant(a));
        v = u.multiplyConstant(1/u.columnVectorMagnitude());
        vt = v.transpose();
        Q.push(Matrix.identity(2).subtract(v.multiplyMatrix(vt).multiplyConstant(2)).identityMerge(1));
        let qt = Q[0].transpose().multiplyMatrix(Q[1].transpose()).transpose();
        qt.transpose().print();
        let R = qt.multiplyMatrix(m);
        R.print();
        qt.transpose().multiplyMatrix(R).print();
        */
        let A = Matrix.copy(this);
        if (A.isTriangular()) {
            return [Matrix.identity(A.body.length), A];
        }
        let length = A.body.length;
        let x = A.getColumnVector(0);
        let a = Math.abs(x.columnVectorMagnitude());
        let e1 = Vector.identity(length, 0);
        e1 = e1.toMatrix().transpose();
        let u = x.subtract(e1.multiplyConstant(a));
        let v = u.multiplyConstant(1/u.columnVectorMagnitude());
        let vt = v.transpose();
        let Q = [];
        Q.push(Matrix.identity(length).subtract(v.multiplyMatrix(vt).multiplyConstant(2)));
        for (var i = 1; i < length-1; i++) {
            let qa = Q[i-1].multiplyMatrix(A);
            let m = qa.minor(i);
            x = m.getColumnVector(0);
            a = Math.abs(x.columnVectorMagnitude());
            e1 = Vector.identity(length-i,0);
            e1 = e1.toMatrix().transpose();
            u = x.subtract(e1.multiplyConstant(a));
            v = u.multiplyConstant(1/u.columnVectorMagnitude());
            vt = v.transpose();
            Q.push(Matrix.identity(length-i).subtract(v.multiplyMatrix(vt).multiplyConstant(2)).identityMerge(i));
        }
        let q = Matrix.identity(length);
        for(var i = 0; i < Q.length; i++) {
            q = q.multiplyMatrix(Q[i].transpose());
        }
        let qt = q.transpose();
        let R = qt.multiplyMatrix(A);
        return [q,R];
    }

    isTriangular() {
        let A = Matrix.copy(this);
        for(var i = 0; i < A.body.length-1; i++) {
            if (A.body[i+1][i] != 0) {
                return false;
            }
        }
        return true;
    }

    diagonal() {
        let A = Matrix.copy(this);
        let data = [];
        for(var i = 0; i < A.body.length; i++) {
            data.push(A.body[i][i]);
        }
        return new Vector(data);
    }

    deflate() {
        let A = Matrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {
            for(var j = 0; j < A.body.length; j++) {
                if (Math.abs(A.body[i][j]) < 0.001) {
                    A.body[i][j] = 0;
                }
            }
        }
        return A;
    }

    realEigenvalues() {
        let A = Matrix.copy(this);
        for(var i = 0; i < A.body.length-1; i++) {
            if (A.body[i+1][i] != 0 && A.body[i+1][i]) {
                return false;
            }
        }
        return true;
    }

    round() {
        let A = Matrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {    
            for(var j = 0; j < A.body.length; j++) {
                A.body[i][j] = Math.round(A.body[i][j]);
            }
        }
        return A;
    }

    schurEigenMatrices() {
        let A = Matrix.copy(this);
        //we need to decompose this matrix into a list of matrices
        //this assumes there are complex eigenvalues in this matrix
        let list = [];
        let boolList = [];
        let prevVector = Vector.null(A.body.length);
        for(var i = 0; i < A.body.length; i++) {
            //get the eigenvalues, so check if it's a real eigenvalue and push the column to the list of Matrices
            let col = A.getColumnVectorAsVector(i);
            //figure out if the eigenvalue is real
            if (col.vector[i+1] == 0 || col.vector[i+1] == undefined) {
                if (boolList[i] != false && boolList[i] != true) {
                    boolList.push(true);
                }
            } else {
                boolList.push(false);
                boolList.push(false);
            }
        }
        //now we know according to the list of booleans which are real and which are complex
        for(var i = 0; i < boolList.length; i++) {
            if (boolList[i]) {
                list.push(new Matrix([[A.getColumnVectorAsVector(i).vector[i]]]));
            } else {
                let a = A.getColumnVectorAsVector(i);
                let b = A.getColumnVectorAsVector(i+1);
                let m = new Matrix([[a.vector[i],b.vector[i]],[a.vector[i+1],b.vector[i+1]]]);
                list.push(m);
                i++;
            }
        }
        return list;
    }

    schurEigenvalues() {
        let A = Matrix.copy(this);
        let list = [];
        if(A.realEigenvalues()) {
            for(var i = 0; i < A.body.length; i++) {
                list.push(new Complex(A.body[i][i], 0));
            }
        } else {
            let matrixList = A.schurEigenMatrices();
            for(var i = 0; i < matrixList.length; i++) {
                if (matrixList[i].body.length == 1) {
                    list.push(new Complex(matrixList[i].body[0][0], 0));
                } else {
                    //finding eigenvalues and eigenvectors is pretty easy but we do need to account for special cases: this is where we use our complex variable class
                    //go for the eigenvalues first
                    //figure out how tf i do eigenvectors later
                    //we can just use the quadratic formula lol
                    let coefficients = matrixList[i].charpoly().reverse();
                    let a = coefficients[0];
                    let b = coefficients[1];
                    let c = coefficients[2];
                    let real = -1*b/(2*a);
                    //we know b^2-4ac is negative
                    let complex = Math.sqrt(-1*(b*b-4*a*c))/(2*a);
                    let e1 = new Complex(real, complex);
                    let e2 = e1.conjugate();
                    list.push(e1);
                    list.push(e2);
                }
            }
        }
        return list;
    }

    qr(n) {
        let A = Matrix.copy(this);
        let Q = Matrix.identity(A.body.length);
        for(var i = 0; i < n; i++) {
            let S = A.householder();
            A = S[1].multiplyMatrix(S[0]);
            Q = Q.multiplyMatrix(S[0]);
        }
        return [A, Q]; //A is the schur matrix, Q is the orthogonal matrix
    }

    eigen() {
        //apply qr algorithm to get the product of all orthogonal matrices Q and the schur matrix S
        let A = Matrix.copy(this);
        let QR = A.qr(100);
        let S = QR[0];
        let realEigenvectorMatrix = QR[1];
        let eigenvalues = [];
        let eigenvectors = [];
        eigenvalues = S.deflate().schurEigenvalues();
        // eigenvectors.forEach(e => {
        //     e.print();
        // })
        // if (S.deflate().realEigenvalues()) {
        //     // for(var i = 0; i < realEigenvectorMatrix.body.length; i++) {
        //     //     let col = realEigenvectorMatrix.getColumnVectorAsVector(i);
        //     //     col = col.multiply(1/col.vector[col.vector.length-1]);
        //     //     eigenvectors.push(col);
        //     // }
        //     for(var i = 0; i < eigenvalues.length; i++) {
        //         let vectors = A.subtract(Matrix.identity(A.body.length).multiplyConstant(eigenvalues[i])).kernel();
        //         for(var j = 0; j < vectors.length; j++) {
        //             eigenvectors.push(vectors[j].toMatrix().toComplexMatrix().body[0]);
        //         }
        //     }
        // } else {
            for(var i = 0; i < eigenvalues.length; i++) {
                let L = A.toComplexMatrix().subtract(ComplexMatrix.identity(A.body.length).multiplyConstant(eigenvalues[i]));
                let list = L.kernel();
                eigenvectors.push(list);
            }
        // }
        return [eigenvalues, eigenvectors];
    }

    toComplexMatrix() {
        let body = [];
        let A = Matrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {
            let list = [];
            for(var j = 0; j < A.body[i].length; j++) {
                list.push(new Complex(A.body[i][j], 0));
            }
            body.push(list);
        }
        return new ComplexMatrix(body);
    }

    diagonalize() {
        let A = Matrix.copy(this);
        let eigen = A.eigen();
        let eigenvectors = eigen[1]; //list
        let m = ComplexMatrix.null(eigenvectors.length,eigenvectors.length);
        //eigenvectors is a list of complex numbers
        for(var i = 0; i < eigenvectors.length; i++) {
            m.body[i] = eigenvectors[i];
        }
        if (m.isSingular()) {
            return null;
        }
        m = m.transpose();
        let mi = m.inverse();
        A = A.toComplexMatrix();
        return mi.multiplyMatrix(A).multiplyMatrix(m).deflate();
    }

    split() {
        //assuming we have even rows
        let A = Matrix.copy(this);
        let b1 = [];
        let b2 = [];
        for(var i = 0; i < A.body.length; i++) {
            if (i < A.body.length/2) {
                b1.push(A.body[i]);
            } else {
                b2.push(A.body[i]);
            }
        }
        return [new Matrix(b1), new Matrix(b2)];
    }
}

class ComplexMatrix {
    body = [];
    constructor(matrix) {
        for(var i = 0; i < matrix.length; i++) {
            let temp = []
            for(var j = 0; j < matrix[i].length; j++) {
                temp.push(matrix[i][j]);
            }
            this.body.push(temp);
        }
    }
    //know that this is a 2d array of Complex objects
    static copy(matrix) {
        let arr = [];
        for(var i = 0; i < matrix.body.length; i++) {
            let row = [];
            for(var j = 0; j < matrix.body[i].length; j++) {
                row.push(matrix.body[i][j]);
            }
            arr.push(row);
        }
        return new ComplexMatrix(arr);
    }

    static identity(n) {
        let vectors = [];
        for(var i = 0; i < n; i++) {
            let vector = [];
            for(var j = 0; j < n; j++) {
                if (i == j) {
                    vector.push(new Complex(1,0))
                } else {
                    vector.push(new Complex(0,0))
                }
            }
            vectors.push(vector);
        }
        return new ComplexMatrix(vectors);
    }

    static null(r,c) {
        let vectors = [];
        for(var i = 0; i < r; i++) {
            let vector = [];
            for(var j = 0; j < c; j++) {
                vector.push(new Complex(0,0))
            }
            vectors.push(vector);
        }
        return new ComplexMatrix(vectors);
    }

    print() {
        for(var i = 0; i < this.body.length; i++) {
            let str = "";
            for(var j = 0; j < this.body[i].length; j++) {
                str += (this.body[i][j].toString() + " ");
            }
            console.log(str);
        }
    }

    add(matrix) {
        let A = ComplexMatrix.copy(this)
        for(var i = 0; i < A.body.length; i++) {
            for(var j = 0; j < A.body[i].length; j++) {
                A.body[i][j] = A.body[i][j].add(matrix.body[i][j]);
            }
        }
        return A;
    }

    subtract(matrix) {
        let A = ComplexMatrix.copy(this)
        for(var i = 0; i < A.body.length; i++) {
            for(var j = 0; j < A.body[i].length; j++) {
                A.body[i][j] = A.body[i][j].subtract(matrix.body[i][j]);
            }
        }
        return A;
    }

    multiplyConstant(c) {
        //c is a complex number!
        let A = ComplexMatrix.copy(this)
        for(var i = 0; i < A.body.length; i++) {
            for(var j = 0; j < A.body[i].length; j++) {
                A.body[i][j] = A.body[i][j].multiply(c);
            }
        }
        return A;
    }

    rowSwap(r1, r2) {
        let A = ComplexMatrix.copy(this);
        let temp = A.body[r1];
        A.body[r1] = A.body[r2];
        A.body[r2] = temp;
        return A;
    }

    rowMultiply(row, c) {
        let A = ComplexMatrix.copy(this);
        for(var i = 0; i < A.body[0].length; i++) {
            A.body[row][i] = A.body[row][i].multiply(c);
        }
        return A;
    }

    rowAdd(orig, dest) {
        //origin is the row you are adding
        //destination is the row to add to
        let A = ComplexMatrix.copy(this);
        for(var i = 0; i < A.body[dest].length; i++) {
            A.body[dest][i] = A.body[dest][i].add(A.body[orig][i]);
        }
        return A;
    }

    rowMultiplyAdd(orig, dest, c) {
        //origin is the row you are adding
        //destination is the row to add to
        //c is a complex number you are multiplying origin with
        let A = ComplexMatrix.copy(this);
        for(var i = 0; i < A.body[dest].length; i++) {
            A.body[dest][i] = A.body[dest][i].add(A.body[orig][i].multiply(c));
        }
        return A;
    }

    deflate() {
        let A = ComplexMatrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {
            for(var j = 0; j < A.body.length; j++) {
                if (Math.abs(A.body[i][j].re) < 0.001) {
                    A.body[i][j] = new Complex(0,A.body[i][j].im);
                }
                if (Math.abs(A.body[i][j].im) < 0.001) {
                    A.body[i][j] = new Complex(A.body[i][j].re,0);
                }
            }
        }
        return A;
    }

    rref() {
        let A = ComplexMatrix.copy(this);
        //ensure [0][0] is nonzero
        let zero = new Complex(0,0);
        if (A.body[0][0].equals(zero)) {
            for(var i = 0; i < A.body.length; i++) {
                if (!A.body[i][0].equals(zero)) {
                    A = A.rowSwap(0,i);
                    break;
                }
            }
        }
        let nullMatrix = ComplexMatrix.null(A.body.length, A.body[0].length);
        //now we go through each column and ensure [i][i'] is one
        //we have to check each column before applying transforms because if the whole column is zero we move to the next
        //target row and target column
        //add 1 to the target column if the column is zero
        //add 1 to the target row and column if the column has been reduced correctly (0 entries aside from [i][i'])
        let targetCol = 0;
        let targetRow = 0;
        for(var i = 0; i < A.body[0].length; i++) {
            if (targetCol > A.body[0].length-1 || targetRow > A.body.length-1) {
                break;
            }
            if (A.body[targetRow][targetCol].equals(zero)) {
                targetRow++;
                targetCol++;
            }
            //unify column
            if (targetCol > A.body.length-1 || targetRow > A.body.length-1) {
                break;
            }
            if (!A.body[targetRow][targetCol].equals(new Complex(1,0))) {
                //console.log(targetRow);
                if (A.body[targetRow][targetCol].equals(new Complex(0,0))) {
                    //console.log("sqitch");
                    //switch rows, otherwise do nothing
                    for(var j = targetRow; j < A.body.length; j++) {
                        if (!A.body[j][targetCol].equals(new Complex(0,0))) {
                            A = A.rowSwap(j, targetRow);
                            break;
                        }
                    }
                } else {
                    let constant = A.body[targetRow][targetCol];
                    A = A.rowMultiply(targetRow, constant.reciprocal());
                }
            }
            //check if column is zero
            let nullColumn = true;
            for(var i = 0; i < A.body.length; i++) {
                if (!A.body[i][targetCol].equals(zero)) {
                    nullColumn = false;
                    break;
                }
            }
            if (!nullColumn) {
                //zero rows
                for(var j = 0; j < A.body.length; j++) {
                    if (j != targetRow) {
                        A = A.rowMultiplyAdd(targetRow, j, A.body[targetRow][targetCol].multiply(A.body[j][targetCol]).negate());
                    }
                }
                targetCol++;
                targetRow++;
            } else {
                targetCol++;
            }
            if (A.body.length == A.body[0].length) {
                A = A.deflate().subtract(ComplexMatrix.identity(A.body.length)).deflate().add(ComplexMatrix.identity(A.body.length));
            }
        }
        return A;
    }
    
    kernel() {
        let A = ComplexMatrix.copy(this).rref();
        let vector = [];
        for(var i = 0; i < A.body.length; i++) {
            vector.push(A.body[i][A.body.length-1].negate())
        }
        vector[vector.length-1] = new Complex(1,0);
        return vector;
    }
    
    nullspace() {
        let A = ComplexMatrix.copy(this).rref();
        //last column has our null space times negative 1
        let vector = [];
        for(var i = 0; i < A.body.length; i++) {
            vector.push([A.body[i][A.body.length-1].negate()])
        }
        let n = new ComplexMatrix(vector);
        n.body[A.body.length-1][0] = new Complex(1,0);
        return n;
    }

    appendIdentity() {
        let A = ComplexMatrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {
            let row = A.body[i];
            for(var j = 0; j < A.body.length; j++) {
                if (j == i) {
                    row.push(new Complex(1,0))
                } else {
                    row.push(new Complex(0,0))
                }
            }
            A.body[i] = row;
        }
        return A;
    }

    isSingular() {
        return this.det().equals(new Complex(0,0));
    }

    det() {
        let A = ComplexMatrix.copy(this);
        let eigen = A.eigenvalues();
        let c = eigen[0];
        for(var i = 1; i < eigen.length; i++) {
            c = c.multiply(eigen[i]);
        }
        if (Math.abs(c.re) < 0.001) {
            c = new Complex(0,c.im);
        }
        if (Math.abs(c.im) < 0.001) {
            c = new Complex(c.re,0);
        }
        return c;
    }

    eigenvalues() {
        let A = ComplexMatrix.copy(this);
        let q = A.qr(100);
        //q[0] is the schur matrix
        let eigenvalues = []
        for(var i = 0; i < q[0].body.length; i++) {
            eigenvalues.push(q[0].body[i][i])
        }
        return eigenvalues;
    }

    qr(n) {
        let A = ComplexMatrix.copy(this);
        let Q = ComplexMatrix.identity(A.body.length);
        for(var i = 0; i < n; i++) {
            let S = A.householder();
            A = S[1].multiplyMatrix(S[0]);
            Q = Q.multiplyMatrix(S[0]);
        }
        return [A, Q]; //A is the schur matrix, Q is the orthogonal matrix
    }

    isTriangular() {
        let A = ComplexMatrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {
            for(var j = 0; j < i; j++) {
                if (!A.body[i][j].equals(new Complex(0,0))) {
                    return false;
                }
            }
        }
        return true;
    }

    getColumnVector(n) {
        let list = [];
        let A = ComplexMatrix.copy(this);
        for(var i = 0; i < A.body.length; i++) {
            list.push([A.body[i][n]])
        }
        return new ComplexMatrix(list);
    }

    isComplex() {
        let A = ComplexMatrix.copy(this)
        for(var i = 0; i < A.body.length; i++) {
            for(var j = 0; j < A.body[i].length; j++) {
                if (A.body[i][j].im != 0) {
                    return true;
                }
            }
        }
        return false;
    }

    columnVectorMagnitude() {
        //assuming 1 column
        let A = ComplexMatrix.copy(this);
        let sum = new Complex(0,0);
        for(var i = 0; i < A.body.length; i++) {
            sum = sum.add(A.body[i][0].multiply(A.body[i][0]));
        }
        return new Complex(Math.sqrt(sum.re), 0);
    }

    static e(index, n) {
        let list = [];
        for(var i = 0; i < n; i++) {
            if (i == index) {
                list.push([new Complex(1,0)])
            } else {
                list.push([new Complex(0,0)])
            }
        }
        return new ComplexMatrix(list);
    }

    multiplyMatrix(matrix) {
        let A = ComplexMatrix.copy(this);
        let B = ComplexMatrix.copy(matrix);
        let m = [];
        for(var i = 0; i < A.body.length; i++) {
            //row
            let newRow = [];
            for(var j = 0; j < B.body[0].length; j++) {
                //column
                //now get the row and column vectors
                let r = A.body[i];
                let c = B.getColumnVector(j);
                let C = [];
                for(var k = 0; k < B.body.length; k++) {
                    C.push(c.body[k][0]);
                }
                //now multiply them together and push to the new row array
                let d = new Complex(0,0);
                for(var k = 0; k < Math.min(r.length, C.length); k++) {
                    d = d.add(r[k].multiply(C[k]));
                }
                newRow.push(d);
            }
            //newRow has the amount of elements it needs
            //push to the matrix
            m.push(newRow);
        }
        return new ComplexMatrix(m);
        /*
        //dot the row of this matrix with the column of the specified matrix
        let A = Matrix.copy(this);
        let B = Matrix.copy(matrix);
        //we execute as many multiplications as the amount of rows in the first matrix times the amount of columns in the second matrix
        //when finding a column vector times its transpose we should get a 3x3 matrix
        let m = [];
        for(var i = 0; i < A.body.length; i++) {
            //row
            let newRow = [];
            for(var j = 0; j < B.body[0].length; j++) {
                //column
                //now get the row and column vectors
                let r = A.getRowVector(i);
                let c = B.getColumnVectorAsVector(j);
                //now multiply them together and push to the new row array
                let d = r.dot(c);
                newRow.push(d);
            }
            //newRow has the amount of elements it needs
            //push to the matrix
            m.push(newRow);
        }
        return new Matrix(m);
        */
    }

    identityMerge(n) {
        //n dictates how many rows and columns to add to the matrix
        let A = ComplexMatrix.copy(this);
        let I = ComplexMatrix.identity(n + A.body.length);
        for(var i = n; i < I.body.length; i++) {
            for(var j = n; j < I.body.length; j++) {
                I.body[i][j] = A.body[i-n][j-n];
            }
        }
        return I;
    }

    minor(n) {
        let A = ComplexMatrix.copy(this);
        for(var i = -1; i < n; i++) {
            A.body.reverse().pop();
            for(var j = 0; j < A.body.length; j++) {
                A.body[j].reverse().pop();
                A.body[j] = A.body[j].reverse();
            }
            A.body = A.body.reverse()
        }
        return A;
    }

    householder() {
        let A = ComplexMatrix.copy(this);
        if (A.isTriangular()) {
            return [ComplexMatrix.identity(A.body.length), A];
        }
        let length = A.body.length;
        let x = A.getColumnVector(0);
        let m = x.columnVectorMagnitude();
        let e = ComplexMatrix.e(0,length);
        let a = new Complex(0,0);
        if (!x.isComplex()) {
            a = new Complex(m, 0);
        } else {
            a = new Complex(Math.cos(x.body[0][0].arg()),Math.sin(x.body[0][0].arg())).multiply(m).negate();
        }
        let u = x.subtract(e.multiplyConstant(m));
        let v = u.multiplyConstant(u.columnVectorMagnitude().reciprocal());
        let I = ComplexMatrix.identity(length);
        let vvt = v.multiplyMatrix(v.conjugateTranspose());
        let Q = [];
        Q.push(I.subtract(vvt.multiplyConstant(new Complex(2,0))));
        for(var i = 1; i < length-1; i++) {
            let qa = Q[i-1].multiplyMatrix(A);
            let m = qa.minor(i-1);
            x = m.getColumnVector(0);
            m = x.columnVectorMagnitude();
            e = ComplexMatrix.e(0,length-i);
            a = new Complex(0,0);
            if (!x.isComplex()) {
                a = new Complex(m, 0);
            } else {
                a = new Complex(Math.cos(x.body[0][0].arg()),Math.sin(x.body[0][0].arg())).multiply(m).negate();
            }
            u = x.subtract(e.multiplyConstant(m));
            v = u.multiplyConstant(u.columnVectorMagnitude().reciprocal());
            I = ComplexMatrix.identity(length-i);
            vvt = v.multiplyMatrix(v.conjugateTranspose());

            Q.push(I.subtract(vvt.multiplyConstant(new Complex(2,0))).identityMerge(i));
        }
        let q = ComplexMatrix.identity(length);
        for(var i = 0; i < Q.length; i++) {
            q = q.multiplyMatrix(Q[i].transpose());
        }
        let qt = q.transpose();
        let R = qt.multiplyMatrix(A);
        return [q,R];


        // let A = Matrix.copy(this);
        // if (A.isTriangular()) {
        //     return [Matrix.identity(A.body.length), A];
        // }
        // let length = A.body.length;
        // let x = A.getColumnVector(0);
        // let a = Math.abs(x.columnVectorMagnitude());
        // let e1 = Vector.identity(length, 0);
        // e1 = e1.toMatrix().transpose();
        // let u = x.subtract(e1.multiplyConstant(a));
        // let v = u.multiplyConstant(1/u.columnVectorMagnitude());
        // let vt = v.transpose();
        // let Q = [];
        // Q.push(Matrix.identity(length).subtract(v.multiplyMatrix(vt).multiplyConstant(2)));
        // for (var i = 1; i < length-1; i++) {
        //     let qa = Q[i-1].multiplyMatrix(A);
        //     let m = qa.minor(i);
        //     x = m.getColumnVector(0);
        //     a = Math.abs(x.columnVectorMagnitude());
        //     e1 = Vector.identity(length-i,0);
        //     e1 = e1.toMatrix().transpose();
        //     u = x.subtract(e1.multiplyConstant(a));
        //     v = u.multiplyConstant(1/u.columnVectorMagnitude());
        //     vt = v.transpose();
        //     Q.push(Matrix.identity(length-i).subtract(v.multiplyMatrix(vt).multiplyConstant(2)).identityMerge(i));
        // }
        // let q = Matrix.identity(length);
        // for(var i = 0; i < Q.length; i++) {
        //     q = q.multiplyMatrix(Q[i].transpose());
        // }
        // let qt = q.transpose();
        // let R = qt.multiplyMatrix(A);
        // return [q,R];
    }

    inverse() {
        
        let A = ComplexMatrix.copy(this);
        if (A.isSingular()) {
            return null;
        }
        let inverse = A.appendIdentity().rref();
        for(var i = 0; i < inverse.body.length; i++) {
            let I = inverse.body[i].reverse();
            for(var j = 0; j < inverse.body.length; j++) {
                I.pop();
            }
            A.body[i] = I.reverse();
        }
        return A;
    }

    transpose() {
        let A = ComplexMatrix.copy(this);
        let B = ComplexMatrix.null(A.body[0].length, A.body.length);
        for(var i = 0; i < A.body.length; i++) {
            for(var j = 0; j < A.body[i].length; j++) {
                B.body[j][i] = A.body[i][j];
            }
        }
        return B;
    }

    conjugateTranspose() {
        let A = ComplexMatrix.copy(this);
        let B = ComplexMatrix.null(A.body[0].length, A.body.length);
        for(var i = 0; i < A.body.length; i++) {
            for(var j = 0; j < A.body[i].length; j++) {
                B.body[j][i] = A.body[i][j].conjugate();
            }
        }
        return B;
    }

    static dotVectors(v1,v2) {
        //v1 and v2 are both lists of complex numbers
        let dot = new Complex(0,0);
        for(var i = 0; i < v1.length; i++) {
            dot = dot.add(v1[i].multiply(v2[i]))
        }
        return dot;
    }
}

class Complex {
    re = 0;
    im = 0;
    constructor(a, b) {
        this.re = a;
        this.im = b;
    }

    equals(z) {
        let c = this;
        return z.re == c.re && z.im == c.im;
    }

    add(z2) {
        let z1 = this;
        return new Complex(z1.re + z2.re, z1.im + z2.im);
    }

    subtract(z2) {
        let z1 = this;
        return new Complex(z1.re - z2.re, z1.im - z2.im);
    }
    
    multiply(z2) {
        let z1 = this;
        let a = z1.re;
        let b = z1.im;
        let c = z2.re;
        let d = z2.im;
        return new Complex(a*c-b*d,a*d+b*c);
    }

    divide(z2) {
        let z1 = this;
        let a = z1.re;
        let b = z1.im;
        let c = z2.re;
        let d = z2.im;
        return new Complex((a*c+b*d)/(Math.pow(c,2)+Math.pow(d,2)),(b*c-a*d)/(Math.pow(c,2)+Math.pow(d,2)));
    }

    arg() {
        return Math.atan2(this.im,this.re);
    }

    reciprocal() {
        let z = this;
        let a = z.re;
        let b = z.im;
        return new Complex(a/(a*a+b*b),(b*-1)/(a*a+b*b));
    }

    negate() {
        let z = this;
        let a = z.re;
        let b = z.im;
        return new Complex(a*-1,b*-1);
    }

    conjugate() {
        let z1 = this;
        return new Complex(z1.re, -1*z1.im);
    }

    sqrt() {
        //square rooting a real number. we should not have a complex value
        let z = this;
        if (z < 0) {
            return new Complex(0,Math.sqrt(Math.abs(z.re)));
        } else {
            return new Complex(Math.sqrt(z.re), 0);
        }
    }

    toString() {
        let r = "";
        let i = "";
        if (this.re != 0) {
            r += this.re;
        }
        if (this.im != 0) {
            if (this.im == 1) {
                i = "+i";
            } else if (this.im == -1) {
                i = "-i"
            } else {
                if (this.im >= 0) {
                    i += "+"
                }
                i += this.im;
                i += "i";
            }
        }
        if (this.re == 0 && this.im == 0) {
            r = "0";
            i = "";
        }
        if (this.re == 0 && this.im > 0) {
            r = "";
            i = this.im + "i"
        }
        if (this.re == 0 && this.im == 1) {
            i = "i";
        }
        return r + i;
    }

    print() {
        console.log(this.toString());
    }
}

class Vector {
    vector = [];
    constructor(vec) {
        for(var i = 0; i < vec.length; i++) {
            this.vector.push(vec[i]);
        }
    }

    static identity(num, place) {
        let list = [];
        for(var i = 0; i < num; i++) {
            list.push(0);
        }
        list[place] = 1;
        return new Vector(list);
    }

    static null(num) {
        let list = [];
        for(var i = 0; i < num; i++) {
            list.push(0);
        }
        return new Vector(list);
    }

    dot(vector) {
        //make dimensions equal
        while(Math.abs(this.vector.length - vector.vector.length) != 0) {
            if (this.vector.length < vector.vector.length) {
                this.vector.push(0);
            } else if (this.vector.length > vector.vector.length) {
                vector.vector.push(0);
            }
        }
        //apply dot
        let dot = 0;
        for(var i = 0; i < this.vector.length; i++) {
            dot += (vector.vector[i] * this.vector[i]) 
        }
        return dot;
    }

    multiply(c) {
        //use disposable variable by making a new vector and operating on it instead of the core vector
        let arr = [];
        for(var i = 0; i < this.vector.length; i++) {
            arr.push(this.vector[i]);
        }
        let vec = new Vector(arr);
        for(var i = 0; i < vec.vector.length; i++) {
            vec.vector[i] *= c;
        }
        return vec;
    }

    add(v) {
        //use disposable variable by making a new vector and operating on it instead of the core vector
        let arr = [];
        for(var i = 0; i < this.vector.length; i++) {
            arr.push(this.vector[i]);
        }
        let vec = new Vector(arr);
        //make dimensions equal
        while(Math.abs(vec.vector.length - v.vector.length) != 0) {
            if (vec.vector.length < v.vector.length) {
                vec.vector.push(0);
            } else if (vec.vector.length > v.vector.length) {
                v.vector.push(0);
            }
        }
        //add
        for(var i = 0; i < vec.vector.length; i++) {
            vec.vector[i] += v.vector[i];
        }
        return vec;
    }

    subtract(v) {
        //use disposable variable by making a new vector and operating on it instead of the core vector
        let arr = [];
        for(var i = 0; i < this.vector.length; i++) {
            arr.push(this.vector[i]);
        }
        let vec = new Vector(arr);
        //make dimensions equal
        while(Math.abs(vec.vector.length - v.vector.length) != 0) {
            if (vec.vector.length < v.vector.length) {
                vec.vector.push(0);
            } else if (vec.vector.length > v.vector.length) {
                v.vector.push(0);
            }
        }
        //subtract
        for(var i = 0; i < vec.vector.length; i++) {
            vec.vector[i] -= v.vector[i];
        }
        return vec;
    }

    magnitude() {
        let sum = 0;
        for(var i = 0; i < this.vector.length; i++) {
            sum += Math.pow(this.vector[i], 2);
        }
        return Math.sqrt(sum);
    }

    print() {
        console.log(this.vector);
    }

    toMatrix() {
        return new Matrix([this.vector]);
    }
}


let selectedItem = "";
let requireSquare = false;
let itemsSquare = ["Inverse", "Trace", "Determinant", "Characteristic Polynomial", "Eigenvalues", "Eigenvectors", "Diagonalize"];
let dualMatrix = ["Add", "Subtract", "Multiply by Matrix"];
let constant = false;
let twoMatrices = false;
function storeSelection() {
    var selectedOption = document.getElementById('dropdown').value;
    document.getElementById('selectedValue').textContent = selectedOption;
    console.log(selectedOption);

    // For demonstration, storing the selected value in a variable
    var storedSelection = selectedOption;
    if (storedSelection == "None") {
        alert("Select something else you fucking idiot")
    }
    console.log("Stored Selection:", storedSelection);
    selectedItem = selectedOption;
    requireSquare = itemsSquare.includes(selectedItem);
    constant = (selectedItem == "Multiply by Constant" || selectedItem == "Power");

    // You can also store the selection in the browser's local storage
    // localStorage.setItem('selectedDropdownValue', selectedOption);
}
//setup for grid generation
const container = document.getElementById('setup-boxes');
container.innerHTML = ''; // Clear existing content

for (let i = 0; i < 1; i++) {
    //create square or nonsquare matrix depending on operation
    const row = document.createElement('div');
    row.classList.add('row');

    for (let j = 0; j < 2; j++) {
        const input = document.createElement('input');
        input.type = 'text';
        input.classList.add('box-input');
        input.name = `${i}-${j}`; // Assign unique name to each input

        // You can add additional attributes or event listeners to the input element here

        row.appendChild(input);
    }

    container.appendChild(row);
}
let rows = 0;
let cols = 0;
document.getElementById('matrix-setup').addEventListener('submit', function (event) {
    twoMatrices = dualMatrix.includes(selectedItem);
    console.log(twoMatrices);
    event.preventDefault(); // Prevent default form submission

    // Retrieve form data
    const formData = new FormData(this);

    // Process form data (example: log values to console)
    let rc = [];
    for (const entry of formData.entries()) {
        rc.push(Number.parseInt(entry[1]));
    }
    rows = rc[0];
    cols = rc[1];
    if (requireSquare) {
        if (cols != rows) {
            alert(`A square matrix is required for this operation, so a matrix with ${rows} rows and columns will be created.`);
        }
        cols = rows;
    }
    generateGrid(rows, cols, twoMatrices);
});
// Function to generate the grid
function generateGrid(rows, cols, two) {
    const container = document.getElementById('matrix-container');
    container.innerHTML = ''; // Clear existing content

    for (let i = 0; i < rows; i++) {
        const row = document.createElement('div');
        row.classList.add('row');

        for (let j = 0; j < cols; j++) {
            const input = document.createElement('input');
            input.type = 'text';
            input.classList.add('matrix-input');
            input.name = `${i}-${j}`; // Assign unique name to each input

            // You can add additional attributes or event listeners to the input element here

            row.appendChild(input);
        }

        container.appendChild(row);
    }
    if (two) {
        let output = document.createElement('body');
        output.innerHTML = `<body>Second Matrix:</body><div></div>`;
        container.appendChild(output);
        for (let i = 0; i < rows; i++) {
            const row = document.createElement('div');
            row.classList.add('row');
    
            for (let j = 0; j < cols; j++) {
                const input = document.createElement('input');
                input.type = 'text';
                input.classList.add('matrix-input-two');
                input.name = `${i}-${j}`; // Assign unique name to each input
    
                // You can add additional attributes or event listeners to the input element here
    
                row.appendChild(input);
            }
    
            container.appendChild(row);
        }
    } else if (constant) {
        let output = document.createElement('body');
        output.innerHTML = `<body>Constant:</body><div></div>`;
        container.appendChild(output);
        for (let i = 0; i < 1; i++) {
            const row = document.createElement('div');
            row.classList.add('row');
    
            for (let j = 0; j < 1; j++) {
                const input = document.createElement('input');
                input.type = 'text';
                input.classList.add('matrix-input-two');
                input.name = `${i}-${j}`; // Assign unique name to each input
    
                // You can add additional attributes or event listeners to the input element here
    
                row.appendChild(input);
            }
    
            container.appendChild(row);
        }
    }
}

function convertTo2DArray(array1D, numRows, numCols) {
    const array2D = [];
    let index = 0;
  
    for (let i = 0; i < numRows; i++) {
      const row = [];
      for (let j = 0; j < numCols; j++) {
        if (index < array1D.length) {
          row.push(array1D[index]);
          index++;
        } else {
          // If there are not enough elements in the 1D array, fill the remaining cells with null or any other default value
          row.push(null); // You can change this to any default value you prefer
        }
      }
      array2D.push(row);
    }
  
    return array2D;
  }

let matrix = Matrix.null(rows, cols);
let matrix2 = Matrix.copy(matrix);
let entries = [];
let eigen = [];

// Event listener for form submission
document.getElementById('matrix-form').addEventListener('submit', function (event) {
    
    console.log(requireSquare);
    event.preventDefault(); // Prevent default form submission

    // Retrieve form data
    const formData = new FormData(this);

    // Process form data (example: log values to console)
    for (const entry of formData.entries()) {
        const regex = /[^0-9+\-i]/g;
        if (regex.test(entry[1])) {
            alert("invalid character detected, please try again");
        }
    }
    //go through again since its guaranteed to have no invalid boxes
    for (const entry of formData.entries()) {
        entries.push(Number.parseFloat(entry[1]))
    }
    //now parse the matrix/matrices
    let cons = null;
    if (!constant) {
        arr = twoMatrices ? convertTo2DArray(entries, rows*2, cols) : convertTo2DArray(entries, rows, cols);
    } else {
        cons = entries.pop()
        console.log(cons);
        arr = convertTo2DArray(entries, rows, cols);
    }
    console.log(arr);
    if (twoMatrices) {
        matrix.body = arr;
        matrix2 = matrix.split()[1];
        matrix = matrix.split()[0];
    } else {
        matrix.body = arr;
    }
    matrix.print();
    // matrix2.print();
    if (constant) {
        console.log(cons);
    }
    //setup
    let container = document.getElementById('matrix-form');
    container.lastChild.remove();
    //computation/logic
    let out = null;
    //output
    let output = document.createElement('body');
    output.innerHTML = ``;//`<body>f</body><div></div><body>f</body><div></div><body>f</body>`;
    switch(selectedItem) {
        case "Add":
            out = matrix.add(matrix2);
            for(var i = 0; i < out.body.length; i++) {
                let string = ``
                for(var j = 0; j < out.body[0].length; j++) {
                    string += `${out.body[i][j]} `;
                }
                output.innerHTML += `<body>${string}</body><div></div>`;
            }
            break;
        case "Subtract":
            out = matrix.subtract(matrix2);
            for(var i = 0; i < out.body.length; i++) {
                let string = ``
                for(var j = 0; j < out.body[0].length; j++) {
                    string += `${out.body[i][j]} `;
                }
                output.innerHTML += `<body>${string}</body><div></div>`;
            }
            break;
        case "Multiply by Constant":
            out = matrix.multiplyConstant(cons);
            for(var i = 0; i < out.body.length; i++) {
                let string = ``
                for(var j = 0; j < out.body[0].length; j++) {
                    string += `${out.body[i][j]} `;
                }
                output.innerHTML += `<body>${string}</body><div></div>`;
            }
            break;
        case "Power":
            out = matrix.power(cons);
            for(var i = 0; i < out.body.length; i++) {
                let string = ``
                for(var j = 0; j < out.body[0].length; j++) {
                    string += `${out.body[i][j]} `;
                }
                output.innerHTML += `<body>${string}</body><div></div>`;
            }
            break;
        case "Multiply by Matrix":
            out = matrix.multiplyMatrix(matrix2);
            for(var i = 0; i < out.body.length; i++) {
                let string = ``
                for(var j = 0; j < out.body[0].length; j++) {
                    string += `${out.body[i][j]} `;
                }
                output.innerHTML += `<body>${string}</body><div></div>`;
            }
            break;
        case "Transpose":
            out = matrix.transpose();
            for(var i = 0; i < out.body.length; i++) {
                let string = ``
                for(var j = 0; j < out.body[0].length; j++) {
                    string += `${out.body[i][j]} `;
                }
                output.innerHTML += `<body>${string}</body><div></div>`;
            }
            break;
        case "Reduced Row-Echelon Form":
            out = matrix.rref();
            for(var i = 0; i < out.body.length; i++) {
                let string = ``
                for(var j = 0; j < out.body[0].length; j++) {
                    string += `${out.body[i][j]} `;
                }
                output.innerHTML += `<body>${string}</body><div></div>`;
            }
            break;
        case "Kernel":
            out = matrix.kernel();
            for(var i = 0; i < out.length; i++) {
                let string = ``
                string += `${out[i].vector}`;
                output.innerHTML += `<body>${string}</body><div></div>`;
            }
            break;    
        case "Trace":
            out = matrix.trace();
            output.innerHTML += `<body>${out}</body><div></div>`;
            break;
        case "Inverse":
            out = matrix.inverse();
            if (out == null) {
                output.innerHTML += `<body>Impossible Operation</body><div></div>`;
            } else {
                for(var i = 0; i < out.body.length; i++) {
                    let string = ``
                    for(var j = 0; j < out.body[0].length; j++) {
                        string += `${out.body[i][j]} `;
                    }
                    output.innerHTML += `<body>${string}</body><div></div>`;
                }
            }
            break;
        case "Determinant":
            out = matrix.det();
            output.innerHTML += `<body>${out}</body><div></div>`;
            break;
        case "Characteristic Polynomial":
            out = matrix.charpoly();
            string = ``;
            for(var i = 0; i < out.length; i++) {
                if (i != out.length-1) {
                    if (out[i] < 0) {
                        if (out.length-i-1 != 1) {
                            if (out[i] != -1) {
                                string += `${out[i]}x^${out.length-i-1}`
                            } else {
                                string += `-x^${out.length-i-1}`
                            }
                        } else {
                            if (out[i] != -1) {
                                string += `${out[i]}x`
                            } else {
                                string += `-x`
                            }
                        }
                    } else if (out[i] == 1) {
                        if (out.length-i-1 != 1) {
                            if (i == 0) {
                                string += `x^${out.length-i-1}`
                            } else {
                                string += `+x^${out.length-i-1}`
                            }
                        } else {
                            if (i == 0) {
                                string += `x`
                            } else {
                                string += `+x`
                            }
                        }
                    } else if (out[i] != 1) {
                        if (out.length-i-1 != 1) {
                            string += `+${out[i]}x^${out.length-i-1}`
                        } else {
                            string += `+${out[i]}x`
                        }
                    }
                } else {
                    if (out[i] < 0) {
                        string += `${out[i]}`
                    } else if (out[i] != 0) {
                        string += `+${out[i]}`
                    }
                }
            }
            output.innerHTML += `<body>${string}</body><div></div>`;
            break;
        case "Eigenvalues":
            out = matrix.eigen()[0];
            for(var i = 0; i < out.length; i++) {
                let string = ``
                //for(var j = 0; j < out.body[0].length; j++) {
                    string += `${out[i]} `;
                //}
                output.innerHTML += `<body>${string}</body><div></div>`;
            }
            break;
        case "Eigenvectors":
            out = matrix.eigen()[1];
            for(var i = 0; i < out.length; i++) {
                let string = ``
                for(var j = 0; j < out[i].length; j++) {
                    string += `${out[j][i]} `;
                }
                output.innerHTML += `<body>${string}</body><div></div>`;
            }
            break;
        case "Diagonalize":
            out = matrix.diagonalize();
            if (out == null) {
                output.innerHTML += `<body>Impossible Operation</body><div></div>`;
            } else {
                for(var i = 0; i < out.body.length; i++) {
                    let string = ``
                    for(var j = 0; j < out.body[0].length; j++) {
                        string += `${out.body[i][j]} `;
                    }
                    output.innerHTML += `<body>${string}</body><div></div>`;
                }
            }
            break;
    }
    entries = [];
    container.appendChild(output);
});