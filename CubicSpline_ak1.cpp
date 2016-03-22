// CubicSpline_ak1.cpp - Program for calculating a cubic spline to input data.
// Written in Microsoft Visual Studio Express 2013 for Windows Desktop
// 21 March 2016
//
// The sub-routines included in this program are translations of the FORTRAN routines PCHEV and PCHEZ 
// written by David K.Kahaner, National Bureau of Standards.
// From the book "Numerical Methods and Software"
// D. Kahaner, C. Moler, and S. Nash
// Prentice Hall, 1988
//
// These are top-level programs that control several sub-routines from the SLATEC collection:
//
// http://www.netlib.org/slatec/
//
// To distinguish the routines posted below from others, an _ak1 suffix has been appended to them.
//
// A small main program is included also, to provide an example of how to use CubicSpline. In this 
// example, data is input from a file to eliminate the need for a user to type data in via
// the console.

#include <iostream>
#include <fstream>
#include <cctype>
#include <cmath>
#include <vector>
#include <cfloat>

using namespace std;

typedef vector<double> C1DArray;
typedef vector<vector<double> > C2DArray;

int main()
{
	char rflag = 0;	//Readiness flag 

	cout << "                     CubicSpline_ak1   (21 March 2016)\n";
	cout << "=========================================================================== \n";
	cout << "This program calculates the cubic spline for input data.\n";
	cout << "The (x,y) data pairs should have been saved beforehand in a file named\n";
	cout << "splinedata.txt, which should be in the same folder as the CubicSpline \n";
	cout << "executable.\n";
	cout << "\nThe first entry in this file should be N, the number of (x, y)\n";
	cout << "data pairs.\n";
	cout << "\nThe second entry in this file should be NVAL, the number of x-values at which\n";
	cout << "interpolating values are to be calculated.\n";
	cout << "\nNext, the (x,y) data pairs themselves should be listed.\n";
	cout << "THEY MUST BE LISTED IN ASCENDING VALUES OF X.\n";
	cout << "\nFinally, the x-values at which interpolating values are to be calculated should\n";
	cout << "be given in order of ascending value because the program execution is most\n";
	cout << "efficient when they are given in this order.\n";
	cout << "\nThe data is assumed to be of type double. Variables used within this program\n";
	cout << "are type double.\n";
	cout << "\nOutput is written to the file splineout.txt\n";
	cout << "\nx-Val, Y-val, and derivative-Val.\n";

	cout << "\nIMPORTANT-- Note the Error Codes: \n";
	cout << "\nierr >  0 : Extrapolation was performed at ierr points\n";
	cout << "\nierr =  0 : Normal Completion (No errors)\n";
	cout << "ierr = -1 : N < 2 \n";
	cout << "ierr = -3 : the x-values are NOT given in order of increasing value\n";

	cout << "\nIs everything ready (are you ready to continue?)? If yes, Enter y. \n";
	cout << "Otherwise Enter any other key. \n";
	cin >> rflag;

	if (toupper(rflag) == 'Y') {
		
		int mDim, NVAL;
		
		// NVAL is the number of points at which the interpolating function will be used to provide values

		cout << "Appear to be ready. \n";

		ifstream in("splinedata.txt", ios::in);

		if (!in) {
			cout << "Cannot open the input file.\n";
			return 0;
		}

		ofstream out("splineout.txt", ios::out);

		if (!out) {
			in.close();  //Close the input file
			cout << "Cannot open the output file. \n";
			return 0;
		}

		in >> mDim;  //Input the number of known data pairs from the file
		in >> NVAL;  //Input the number of x-values at which interpolating y-values are to be computed

		if (mDim < 2)  {
			in.close();  //Close the opened files
			out.close();
			cout << "ierr = -1. Fewer than two data pairs input. Program terminated. \n";
			return 0;
		}

		if (NVAL < 1){
			in.close();  //Close the opened files
			out.close();
			cout << "Fewer than one evaluation point specified.\n";
			return 0;
		}

		int i, ierr = 0, ir = 1, j, jfirst = 0, nj;
		int next[2];

		double dummy, temp;		// Dummy variables
		double h, xmi, xma;		// variables for cubic evaluation
		double del1, del2, delta;	// variables for cubic evaluation
		double c2, c2t2, c3, c3t3;// variables for cubic evaluation

		C1DArray dVec, xVec, yVec;	// Arrays for d-, x-, and y-values.
		C1DArray dval, fval, xval;
		C2DArray wk;				// wk array, a scratch (work) array

		// xval is the array of x-values, the NVAL points at which interpolating function values will be calculated
		// fval is the array of y-values, calculated by the interpolating function, calculated at each xval point
		// dval is the array of derivative values, calculated by the interpolating function, calculated at each xval point

		try { // Resize the x, y, and d arrays to their required sizes
			xVec.resize(mDim);
			yVec.resize(mDim);
			dVec.resize(mDim);
		} // End of try block

		catch (bad_alloc& xa) { // Catch block, for exceptions
			in.close();
			out.close();
			cerr << "In catch block for resizing d, x, and y arrays: " << xa.what() << "\n";
			cout << "\nEnter any key to continue. \n";
			cin >> rflag;
			return 0;
		} // End of catch block

		try { // Resize the xval array to its required size
			xval.resize(NVAL);
			dval.resize(NVAL);
			fval.resize(NVAL);
		} // End of try block

		catch (bad_alloc& xa) { // Catch block, for exceptions
			in.close();
			out.close();
			cerr << "In catch block for resizing xval array: " << xa.what() << "\n";
			cout << "\nEnter any key to continue. \n";
			cin >> rflag;
			return 0;
		} // End of catch block

		try { // Resize the wk array to its required size
			wk.resize(2);
			wk[0].resize(mDim);
			wk[1].resize(mDim);
		} // End of try block

		catch (bad_alloc& xa) { // Catch block, for exceptions
			in.close();
			out.close();
			cerr << "In catch block for resizing wk array: " << xa.what() << "\n";
			cout << "\nEnter any key to continue. \n";
			cin >> rflag;
			return 0;
		} // End of catch block

		for (i = 0; i < mDim; ++i) //Input the data pairs from the file
			in >> xVec[i] >> yVec[i];

		for (i = 0; i < NVAL; ++i) //Input the x-values at which interpolants are to be computed
			in >> xval[i];

		in.close();  //Close the input file

		// Echo the data pairs to the console, to ensure it was input correctly
		cout << "\n";
		cout << mDim << " data pairs were input.\n";
		cout << "The data pairs follow:\n";
		for (i = 0; i < mDim; ++i)
			cout << xVec[i] << "     " << yVec[i] << "\n";

		// Echo the xInt values to the console, to ensure they were input correctly
		cout << "\n";
		cout << NVAL << " xInt values were input.\n";
		cout << "The xInt values follow:\n";
		for (i = 0; i < NVAL; ++i)
			cout << xval[i] << "\n";

//		cout << "\nEnter any key to continue. \n";
//		cin >> rflag;
//		return 0;

		// The main program follows.
		// Confirm that the data pairs are in order of increasing x
		for (i = 1; i < mDim; ++i){
			if (xVec[i] <= xVec[i - 1]){
				cout << "ierr = -3. x array not strictly increasing. Program terminated.\n";
				return 0;
			}  // End if (xVec[i] <= xVec[i-1])
		} //End for i

		// Differences between the x-values are now stored in wk[0][. .], starting with wk[0][1]
		// Divided differences between y-values are stored in wk[1][. .], starting with wk[1][1]
		// Note that wk[0][0] and wk[1][0] are not filled in this first loop; they are presently left unassigned

		j = 1;
		for (i = 0; i < (mDim - 1); ++i){
			wk[0][j] = xVec[j] - xVec[i];
			wk[1][j] = (yVec[j] - yVec[i]) / wk[0][j];
			++j;
		}//End for i

		if (mDim == 2){
			wk[1][0] = wk[0][0] = 1.0;
			dVec[0] = 2.0*wk[1][1];
		}// End if (mDim == 2)
		else { // else mDim > 2
			temp = dummy = wk[0][1];
			wk[1][0] = wk[0][2];
			wk[0][0] = temp + wk[0][2];
			dummy *= dummy*wk[1][2];
			dVec[0] = ((temp + 2.0*wk[0][0])*wk[1][1] * wk[0][2] + dummy) / wk[0][0];
		} // End else mDim > 2

		nj = mDim - 1;
		for (i = 1; i < nj; ++i){
			if (wk[1][i - 1] == 0){
				cout << "Error Code 5008.\n";
				return 0;
			}
			temp = -(wk[0][i + 1] / wk[1][i - 1]);
			dVec[i] = temp*dVec[i - 1] + 3.0*(wk[0][i] * wk[1][i + 1] + wk[0][i + 1] * wk[1][i]);
			wk[1][i] = temp*wk[0][i - 1] + 2.0*(wk[0][i] + wk[0][i + 1]);
		}//End for i

		if (mDim == 2){
			dVec[1] = wk[1][1];
		} // End if mDim == 2
		else { //else mDim != 2
			if (mDim == 3){
				dVec[mDim - 1] = 2.0*wk[1][mDim - 1];
				wk[1][mDim - 1] = 1.0;
				if (wk[1][mDim - 2] == 0){
					cout << "Error Code 5008.\n";
					return 0;
				}
				temp = -(1.0 / wk[1][mDim - 2]);
			}// End if (mDim == 3)
			else {
				temp = wk[0][mDim - 2] + wk[0][mDim - 1];
				dummy = wk[0][mDim - 1] * wk[0][mDim - 1] * (yVec[mDim - 2] - yVec[mDim - 3]);
				dummy /= wk[0][mDim - 2];
				dVec[mDim - 1] = ((wk[0][mDim - 1] + 2.0*temp)*wk[1][mDim - 1] * wk[0][mDim - 2] + dummy) / temp;
				if (wk[1][mDim - 2] == 0){
					cout << "Error Code 5008.\n";
					return 0;
				}
				temp = -(temp / wk[1][mDim - 2]);
				wk[1][mDim - 1] = wk[0][mDim - 2];
			}//End else

			// Complete forward pass of Gauss Elimination

			wk[1][mDim - 1] = temp*wk[0][mDim - 2] + wk[1][mDim - 1];
			if (wk[1][mDim - 1] == 0){
				cout << "Error Code 5008.\n";
				return 0;
			}
			dVec[mDim - 1] = (temp*dVec[mDim - 2] + dVec[mDim - 1]) / wk[1][mDim - 1];

		} // End else mDim != 2

		//Carry out back substitution

		for (i = mDim - 2; i >= 0; --i){
			if (wk[1][i] == 0){
				cout << "Error Code 5008.\n";
				return 0;
			}
			dVec[i] = (dVec[i] - wk[0][i] * dVec[i + 1]) / wk[1][i];
		}//End for i

		// End of PCHEZ

		// Start of PCHEV

		// ===============================================================
		// Main loop; go through and calculate interpolant at each xval value

		while (jfirst < NVAL){

			// Locate all points in interval. 
			for (i = jfirst; i < NVAL; i++){
				if (xval[i] >= xVec[ir]) break;
			} // End for loop

			if (i < NVAL){
				if (ir == (mDim - 1))
					i = NVAL;
			} // End if (i < NVAL)

			nj = i - jfirst;

			// Skip evaluation if no points in interval

			if (nj > 0){

				// Evaluate Cubic at xval[i], j = jfirst (1) to i-1
				// ===============================================================
				// Begin CHFDV

				xma = h = xVec[ir] - xVec[ir - 1];

				next[1] = next[0] = 0;
				xmi = 0.0;

				// Compute Cubic Coefficients (expanded about x1)

				delta = (yVec[ir] - yVec[ir - 1]) / h;
				del1 = (dVec[ir - 1] - delta) / h;
				del2 = (dVec[ir] - delta) / h;

				//delta is no longer needed

				c2 = -(del1 + del1 + del2);
				c2t2 = c2 + c2;
				c3 = (del1 + del2) / h;

				// h, del1, and del2 are no longer needed

				c3t3 = c3 + c3 + c3;

				// Evaluation loop

				for (j = 0; j < nj; ++j){
					temp = xval[jfirst + j] - xVec[ir - 1];
					fval[jfirst + j] = yVec[ir - 1] + temp*(dVec[ir - 1] + temp*(c2 + temp*c3));
					dval[jfirst + j] = dVec[ir - 1] + temp*(c2t2 + temp*c3t3);
					if (temp < xmi) next[0] = next[0] + 1;
					if (temp > xma) next[1] = next[1] + 1;
					// Note the redundancy: if either condition is true, other is false
				} // End for j loop

				// End CHFDV

				// ===============================================================
				if ((next[1] > 0) && (ir != (mDim - 1))) {
					cout << "Error Code 5005, Number of evaluation points less than 1.\n";  // This option should never happen
					return 0;
				}  // End if ((next[1] > 0) && (ir != (mDim - 1)))

				if ((next[0] > 0) && (ir != 1)) {
					// xval is not ordered relative to xVec, so must adjust evaluation interval
					// First, locate first point to left of xVec[ir - 1]
					for (j = jfirst; j < i; ++j){
						if (xval[j] < xVec[ir - 1])
							break;
					}
					if (j == i){
						cout << "Error Code 5005, Should not have made it all the way through the loop unless there was an error in CHFDV.\n";  // This option should never happen(?)
						return 0;
					}
					i = j; //Reset i. This will be the first new jfirst

					// Now find out how far to back up in the xVec array

					for (j = 0; j < ir; ++j){
						if (xval[i] < xVec[j])
							break;
					} // End for j

					// The above loop should NEVER run to completion because xval[i] < xVec[ir - 1]

					// At this point, either xval[i] < xVec[0] or
					//                    xVec[j-1] <= xval[i] < xVec[j]
					// Reset ir, recognizing that it will be incremented before cycling

					ir = (((j - 1) > 0) ? (j - 1) : 0);

				} // End if ((next[0] > 0) && (ir != 1))

				jfirst = i;

			} // End if (nj > 0)

			++ir;

			if (ir >= mDim) break;

		}// End while (jfirst < NVAL)


		// End of PCHEV
		// =================================================================

		// Count the number of points which were extrapolated
		for (i = 0; i < NVAL; ++i){
			if ((xval[i] > xVec[mDim - 1]) || (xval[i] < xVec[0]))    ++ierr;
		} // End for loop

		out.precision(DBL_DIG);

		out << "ierr = " << "  " << ierr << ".\n";
		out << "The interpolating data follows (xval, fval, dval):\n";
		out << "\n";
		for (i = 0; i < NVAL; ++i)
			out << xval[i] << "  " << fval[i] << "  " << dval[i] << " \n";

		out << "\n";

		out.close();

		cout << "\nDone! The solution is in the text file splineout.txt \n";

	} //End if rflag = 'Y'
	else cout << "\nNot ready. Try again when ready with information. \n";
	cout << "\nEnter any key to continue. \n";
	cin >> rflag;
	return 0;
}                      // End main program.