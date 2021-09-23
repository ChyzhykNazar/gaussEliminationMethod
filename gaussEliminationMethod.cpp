/*Знайти розв'язок СЛАР методом Гауса*/
/*To find solution for system of linear equation using Gaussian elimination method*/

#include <iostream>
#include <fstream>

using namespace std;

void printFullMatrix(double** mas, double* b, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << mas[i][j] << " ";
        }
        cout << "|" << b[i] << endl;
    }
}

void printMatrix(double** mas, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << mas[i][j] << " ";
        }
        cout << endl;
    }
}

void methodGaus(double** a, double* b, double* x, int n) {
    cout << "___________________________Gauss Method________________________________" << endl;
    double maxValue, s, d;
    int maxIndex;
    int p = 0;
    while (p < n)
    {
        // Search for the string with the maximum value
        maxValue = abs(a[p][p]);
        maxIndex = p;
        for (int i = p + 1; i < n; i++)
        {
            if (abs(a[i][p]) > maxValue)
            {
                maxValue = abs(a[i][p]);        //remember the greatest value
                maxIndex = i;                   //remember the index
            }
        }
        //permutation of lines
        if (maxValue < 0.00001)
        {
            //there are no non-zero diagonal elements
            cout << "The solution cannot be obtained because of the zero column ";
            cout << maxIndex << " matrix A" << endl;
            return;
        }
        for (int j = 0; j < n; j++)
        {
            swap(a[maxIndex][j], a[p][j]);      //rearrange the rows in the matrix
        }
        swap(b[maxIndex], b[p]);                //rearrange the coefficient
        p++;
    }
    //printFullMatrix(a, b, n);
    // Applying Gauss Elimination
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            s = a[j][i] / a[i][i];              // виконую для того, щоб отримати у скільки разів число яке має вийти 0 
            for (int q = i; q < n; q++) {
                a[j][q] = a[j][q] - s * a[i][q];        //subtraction in the array
            }
            b[j] = b[j] - s * b[i];                     //subtraction in the result
            //printFullMatrix(a, b, n);
            //cout << endl;
        }
    }
    // Obtaining Solution by Back Substitution Method
    for (int i = n - 1; i >= 0; i--) {
        d = 0;
        s = 0;
        for (int j = i + 1; j < n; j++) {
            s = a[i][j] * x[j];     //find each multiplier  
            d = d + s;              //sum all the members of the matrix
        }
        x[i] = (b[i] - d) / a[i][i];    //finding roots
    }

    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "]=" << x[i] << " " << endl << endl;       //displaying solution
    }
}

int main()
{
    int n = 0;
    double eps = 0.01;
    int p = 0;
    ifstream file;
    file.open("gaussEliminationMethod.txt");
    file >> n;
    double** a;
    double** aa;
    double* b = new double[n];
    double* bb = new double[n];
    double* x = new double[n];
    a = new double* [n];
    aa = new double* [n];
    for (int i = 0; i < n; i++) {
        a[i] = new double[n];
        aa[i] = new double[n];
    }
    cout << "___________________________Input matrix________________________________" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file >> a[i][j];
            cout << a[i][j] << " ";
        }
        file >> b[i];
        cout << "|" << b[i] << endl;
    }
    file.close();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aa[i][j] = a[i][j];
        }
        bb[i] = b[i];
    }
    methodGaus(a, b, x, n);
}