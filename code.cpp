#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

// --- MATRIX MATH UTILITIES ---
vector<double> identityMatrix() {
    return {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
}

vector<double> addMat(const vector<double>& A, const vector<double>& B) {
    vector<double> C(9);
    for(int i = 0; i < 9; i++) C[i] = A[i] + B[i];
    return C;
}

vector<double> subMat(const vector<double>& A, const vector<double>& B) {
    vector<double> C(9);
    for(int i = 0; i < 9; i++) C[i] = A[i] - B[i];
    return C;
}

vector<double> multMat(const vector<double>& A, const vector<double>& B) {
    vector<double> C(9, 0.0);
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            for(int k = 0; k < 3; k++) {
                C[i * 3 + j] += A[i * 3 + k] * B[k * 3 + j];
            }
        }
    }
    return C;
}

vector<double> scaleMat(double scalar, const vector<double>& A) {
    vector<double> C(9);
    for(int i = 0; i < 9; i++) C[i] = scalar * A[i];
    return C;
}

vector<double> transposeMat(const vector<double>& A) {
    vector<double> C(9);
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            C[j * 3 + i] = A[i * 3 + j];
        }
    }
    return C;
}

double traceMat(const vector<double>& A) {
    return A[0] + A[4] + A[8];
}

double detMat(const vector<double>& A) {
    return A[0] * (A[4]*A[8] - A[5]*A[7]) -
           A[1] * (A[3]*A[8] - A[5]*A[6]) +
           A[2] * (A[3]*A[7] - A[4]*A[6]);
}

// Computes the Frobenius norm of a matrix (useful for error calculation)
double frobeniusNorm(const vector<double>& A) {
    double sum = 0.0;
    for (double val : A) {
        sum += val * val;
    }
    return sqrt(sum);
}

// --- CONTINUUM MECHANICS CORE ---
vector<double> applyDeformationMap(double X, double Y, double Z, double alpha) {
    double x = X * cos(alpha * Z) - Y * sin(alpha * Z);
    double y = X * sin(alpha * Z) + Y * cos(alpha * Z);
    return {x, y, Z};
}

vector<double> getAnalyticalF(double alpha, double Z, double x, double y) {
    double c = cos(alpha * Z);
    double s = sin(alpha * Z);
    return {
        c, -s, -alpha * y,
        s,  c,  alpha * x,
        0.0, 0.0, 1.0
    };
}

// 7. COMPUTE NUMERICAL DEFORMATION GRADIENT (Central Finite Difference)
vector<double> getNumericalF(double X, double Y, double Z, double alpha, double delta = 1e-6) {
    // Forward perturbations
    vector<double> pX_fwd = applyDeformationMap(X + delta, Y, Z, alpha);
    vector<double> pY_fwd = applyDeformationMap(X, Y + delta, Z, alpha);
    vector<double> pZ_fwd = applyDeformationMap(X, Y, Z + delta, alpha);

    // Backward perturbations
    vector<double> pX_bwd = applyDeformationMap(X - delta, Y, Z, alpha);
    vector<double> pY_bwd = applyDeformationMap(X, Y - delta, Z, alpha);
    vector<double> pZ_bwd = applyDeformationMap(X, Y, Z - delta, alpha);

    vector<double> F_num(9);
    double denom = 2.0 * delta;
   
    // Column 1 (Derivative with respect to X)
    F_num[0] = (pX_fwd[0] - pX_bwd[0]) / denom;
    F_num[3] = (pX_fwd[1] - pX_bwd[1]) / denom;
    F_num[6] = (pX_fwd[2] - pX_bwd[2]) / denom;
   
    // Column 2 (Derivative with respect to Y)
    F_num[1] = (pY_fwd[0] - pY_bwd[0]) / denom;
    F_num[4] = (pY_fwd[1] - pY_bwd[1]) / denom;
    F_num[7] = (pY_fwd[2] - pY_bwd[2]) / denom;
   
    // Column 3 (Derivative with respect to Z)
    F_num[2] = (pZ_fwd[0] - pZ_bwd[0]) / denom;
    F_num[5] = (pZ_fwd[1] - pZ_bwd[1]) / denom;
    F_num[8] = (pZ_fwd[2] - pZ_bwd[2]) / denom;

    return F_num;
}

vector<double> getInfinitesimalStrain(const vector<double>& grad_u) {
    vector<double> grad_u_T = transposeMat(grad_u);
    return scaleMat(0.5, addMat(grad_u, grad_u_T));
}

vector<double> computeSmallStrainStress(const vector<double>& eps, double lambda, double mu) {
    double trEps = traceMat(eps);
    return addMat(scaleMat(lambda * trEps, identityMatrix()), scaleMat(2.0 * mu, eps));
}

vector<double> computeFiniteStrainStress(const vector<double>& F, const vector<double>& E_strain, double lambda, double mu) {
    double J = detMat(F);
    double trE = traceMat(E_strain);
    vector<double> S = addMat(scaleMat(lambda * trE, identityMatrix()), scaleMat(2.0 * mu, E_strain)); // Second Piola-Kirchhoff
    vector<double> FSFt = multMat(multMat(F, S), transposeMat(F)); // Push forward to Cauchy
    return scaleMat(1.0 / J, FSFt);
}

double extractStressComponent(const vector<double>& sigma, const vector<double>& vec_left, const vector<double>& vec_right) {
    double result = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result += vec_left[i] * sigma[i * 3 + j] * vec_right[j];
        }
    }
    return result;
}

// --- CLASSICAL THEORY ---
double getClassicalShearStress(double G, double alpha, double r) {
    return G * alpha * r;
}

double getClassicalTorque(double G, double alpha, double a) {
    double Jp = (M_PI * pow(a, 4)) / 2.0;
    return G * Jp * alpha;
}

// --- MAIN EXECUTION: CSV GENERATOR ---
int main() {
    // 1. Define Default Parameters
    double E_modulus = 200e9;
    double nu = 0.3;
    double radius = 0.05;
    double length = 0.5;
    double alpha = 5.0;
    int nR = 20;
    int nTheta = 36;

    // 2. Read Parameters directly from Text File
    ifstream paramFile("input.txt");
    if (paramFile.is_open()) {
        string key;
        while (paramFile >> key) {
            if (key == "E_modulus") paramFile >> E_modulus;
            else if (key == "nu") paramFile >> nu;
            else if (key == "radius") paramFile >> radius;
            else if (key == "length") paramFile >> length;
            else if (key == "alpha") paramFile >> alpha;
            else if (key == "nR") paramFile >> nR;
            else if (key == "nTheta") paramFile >> nTheta;
        }
        paramFile.close();
        cout << "Loaded parameters from input.txt successfully.\n";
    } else {
        cout << "Warning: 'input.txt' not found. Using default values.\n";
    }

    // Calculate dependent material properties
    double mu = E_modulus / (2.0 * (1.0 + nu));
    double lambda = (E_modulus * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double G = mu;
   
    vector<double> I = identityMatrix();

    // 3. Setup CSV File
    ofstream outFile("shaft_torsion_data.csv");
    if (!outFile.is_open()) {
        cerr << "Error: Could not open file for writing." << endl;
        return 1;
    }

    // Write CSV Header
    outFile << "R,Theta,Z,"
            << "X_Ref,Y_Ref,"
            << "x_Def,y_Def,"
            << "u_x,u_y,u_z,";
           
    auto writeMatrixHeaders = [&](const string& prefix) {
        outFile << prefix<<"11,"<<prefix<<"12,"<<prefix<<"13,"
                << prefix<<"21,"<<prefix<<"22,"<<prefix<<"23,"
                << prefix<<"31,"<<prefix<<"32,"<<prefix<<"33,";
    };
   
    writeMatrixHeaders("F_");
    writeMatrixHeaders("F_Num_"); // Headers for numerical F
    writeMatrixHeaders("E_");
    writeMatrixHeaders("Eps_");
    writeMatrixHeaders("SigmaSmall_");
    writeMatrixHeaders("SigmaFinite_");

    outFile << "Tau_Small,Tau_Finite,Tau_Classical,"
            << "Sigma_ZZ_Small,Sigma_ZZ_Finite,"
            << "dTorque_Small,dTorque_Finite,"
            << "F_Error_Frobenius\n"; // Header for error tracking

    // Integration & Discretization Setup
    double dR = radius / nR;
    double dTheta = (2.0 * M_PI) / nTheta;
    double Z = length;
   
    double total_torque_small = 0.0;
    double total_torque_finite = 0.0;
    double max_F_error = 0.0;

    cout << "Generating CSV data..." << endl;

    auto writeMatrixData = [&](const vector<double>& mat) {
        for(int k = 0; k < 9; k++) outFile << mat[k] << ",";
    };

    // Iterate over the cross-section
    for (int i = 0; i < nR; i++) {
        double R = (i + 0.5) * dR;
        for (int j = 0; j < nTheta; j++) {
            double Theta = (j + 0.5) * dTheta;
           
            // Generate material points in the shaft (Reference)
            double X = R * cos(Theta);
            double Y = R * sin(Theta);
           
            // Apply the deformation map
            vector<double> pos = applyDeformationMap(X, Y, Z, alpha);
            double x = pos[0];
            double y = pos[1];
           
            // Compute displacement
            double u_x = x - X;
            double u_y = y - Y;
            double u_z = 0.0;
           
            // Compute Analytical deformation gradient
            vector<double> F = getAnalyticalF(alpha, Z, x, y);
           
            // Compute Numerical deformation gradient and Error (Central Difference)
            vector<double> F_num = getNumericalF(X, Y, Z, alpha);
            double f_err = frobeniusNorm(subMat(F, F_num));
            if (f_err > max_F_error) max_F_error = f_err;
           
            // Strain measures
            vector<double> C = multMat(transposeMat(F), F);
            vector<double> E_strain = scaleMat(0.5, subMat(C, I));
            vector<double> eps_strain = getInfinitesimalStrain(subMat(F, I));
           
            // Compute stress tensors
            vector<double> sigma_small = computeSmallStrainStress(eps_strain, lambda, mu);
            vector<double> sigma_finite = computeFiniteStrainStress(F, E_strain, lambda, mu);
           
            // Extract torsional shear stress and axial normal stress
            double theta_deformed = Theta + alpha * Z;
            vector<double> e_theta = {-sin(theta_deformed), cos(theta_deformed), 0.0};
            vector<double> e_z = {0.0, 0.0, 1.0};
           
            double tau_small = extractStressComponent(sigma_small, e_z, e_theta);
            double tau_finite = extractStressComponent(sigma_finite, e_z, e_theta);
            double tau_class = getClassicalShearStress(G, alpha, R);
           
            double szz_small = extractStressComponent(sigma_small, e_z, e_z);
            double szz_finite = extractStressComponent(sigma_finite, e_z, e_z);
           
            // Compute torque by numerical integration
            double dA = R * dR * dTheta;
            double dT_small = R * tau_small * dA;
            double dT_finite = R * tau_finite * dA;
           
            total_torque_small += dT_small;
            total_torque_finite += dT_finite;

            // --- Write row to CSV ---
            outFile << R << "," << Theta << "," << Z << ","
                    << X << "," << Y << ","
                    << x << "," << y << ","
                    << u_x << "," << u_y << "," << u_z << ",";
           
            writeMatrixData(F);
            writeMatrixData(F_num);
            writeMatrixData(E_strain);
            writeMatrixData(eps_strain);
            writeMatrixData(sigma_small);
            writeMatrixData(sigma_finite);
           
            outFile << tau_small << "," << tau_finite << "," << tau_class << ","
                    << szz_small << "," << szz_finite << ","
                    << dT_small << "," << dT_finite << ","
                    << f_err << "\n";
        }
    }

    outFile.close();
   
    double total_torque_classical = getClassicalTorque(G, alpha, radius);
   
    cout << "Success! Data saved to 'shaft_torsion_data.csv'" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "Maximum F_error (Numerical vs Analytical): " << max_F_error << endl;
    cout << "-------------------------------------------" << endl;
    cout << "INTEGRATION RESULTS (Z = " << length << " m, Alpha = " << alpha << " rad/m):" << endl;
    cout << "Total Torque (Small Strain):  " << total_torque_small << " N*m" << endl;
    cout << "Total Torque (Finite Strain): " << total_torque_finite << " N*m" << endl;
    cout << "Total Torque (Classical):     " << total_torque_classical << " N*m" << endl;

    return 0;
}