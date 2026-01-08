#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TAxis.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <set>
#include <vector>
#include <TMath.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TPaveStatsEditor.h>
#include "TLorentzVector.h"
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>
#include "LHAPDF/LHAPDF.h"
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/rgh/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgh/output_dir

// PER LE AUTOGENERAZIONI
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/rgh_simulation/analysis/computed_first_gen_pipkp/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgh/agen_1_output_dir
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/rgh_simulation/analysis/computed_second_gen_pipkp/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgh/agen_2_output_dir
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/rgh_simulation/analysis/computed_third_gen_pipkp/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgh/agen_3_output_dir
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/rgc/summer22/NH3/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgc/sum22_NH3_data
using namespace std;
using namespace LHAPDF;
namespace fs = std::filesystem;
gROOT->SetBatch(kTRUE);
// To generate bin in log scale, necessary for the canvas
std::vector<float> CreateLogBinning(int nbins, float xmin, float xmax) {
    std::vector<float> bin_edges(nbins + 1);
    float logxmin = std::log10(xmin);
    float logxmax = std::log10(xmax);
    float bin_width = (logxmax - logxmin) / nbins;
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = std::pow(10, logxmin + i * bin_width);
    }
    return bin_edges;
}

TRandom3 randGen(0);  
int random_helicity(float phi){
    double r = randGen.Uniform(0, 1); 
    //double epsilon = 0.01 + 0.01 * cos(phi);  
    //double epsilon = randGen.Uniform(0.01, 0.02);
    //double epsilon = std::min(0.04, std::max(0.02, randGen.Gaus(0.03, 0.004)));
    //float p_up = (1 + epsilon * sin(phi)) / 2;  
    float p_up = 0.51;
    float prob = (r < p_up) ? 1 : -1;  
    return prob;
}

void SetStatsBox2(TH1F* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.75);  // Posizione pannello (sinistra)
        stats0->SetX2NDC(0.88);  // Posizione pannello (destra)
        stats0->SetY1NDC(0.68);  // Posizione pannello (basso)
        stats0->SetY2NDC(0.88);  // Posizione pannello (alto)
    }
}

double AUT_loglike(const double* Aut, const vector<double>& phi_h, const vector<double>& phi_s, const vector<double>& depol, const vector<double>& helicity) {
    double logLike = 0.0;
    double P_t = 1.0;
    double A_sivers = Aut[0];
    double A_collins = Aut[1];
    double A_pretz = Aut[2];
    for(size_t i = 0; i < phi_h.size(); i++) {   
        double modulation = (A_collins * sin(phi_h[i] + phi_s[i])) + (A_sivers * sin(phi_h[i] - phi_s[i])) + A_pretz * sin(3*phi_h[i] - phi_s[i]);
        //double modulation = ((A_collins) * sin(phi_h[i] + phi_s[i])) + ((A_sivers )* sin(phi_h[i] - phi_s[i]));
        double prob = 0.5*(1 + P_t * depol[i] * modulation);
        if (prob <= 0) return 1e10; // Avoid log(negative)
        //logLike += std::log(prob);
        
        if (helicity[i] == +1)
            logLike += log(prob);
        else
            logLike += log(1 - prob);
            
    }
    return -logLike; // Minuit minimizes, so we minimize -logL
}

int getBinIndex_xQ2(double xB, double Q2) {
    double Clas_Q2Bins[][2] = {{1,3}, {3, 5}, {5, 7}, {7,11}};
    double Clas_xBins_1[][2] = {{0, 0.142}, {0.142, 0.172}, {0.172, 0.2}, {0.2, 0.232}, {0.232, 0.269}, {0.269, 0.315}, {0.315, 1}};
    double Clas_xBins_2[][2] = {{0, 0.269}, {0.269, 0.315}, {0.315, 0.384}, {0.384, 0.45}, {0.45, 1}};
    double Clas_xBins_3[][2] = {{0, 0.45}, {0.45, 0.55}, {0.55, 1}};
    double Clas_xBins_4[][2] = {{0, 0.6}, {0.6, 1}};
    int bin_xB = -1;
    int bin_Q2 = -1;
    for (int i = 0; i < 4; i++) {
        if(Q2 >= Clas_Q2Bins[i][0] && Q2 < Clas_Q2Bins[i][1]) {
            if(i == 0){
                for (int j = 0; j < 7; j++){
                    if(xB >= Clas_xBins_1[j][0] && xB < Clas_xBins_1[j][1]){
                        bin_Q2 = i;
                        bin_xB = j;
                        break;
                    }
                }
            } else if(i == 1){
                for (int j = 0; j < 5; j++){
                    if(xB >= Clas_xBins_2[j][0] && xB < Clas_xBins_2[j][1]){
                        bin_Q2 = i;
                        bin_xB = j + 7 ;
                        break;
                    }
                }
            } else if(i == 2){
                for (int j = 0; j < 3; j++){
                    if(xB >= Clas_xBins_3[j][0] && xB < Clas_xBins_3[j][1]){
                        bin_Q2 = i;
                        bin_xB = j + 12;
                        break;
                    }
                }
            } else if(i == 3){
                for (int j = 0; j < 2; j++){
                    if(xB >= Clas_xBins_4[j][0] && xB < Clas_xBins_4[j][1]){
                        bin_Q2 = i;
                        bin_xB = j + 15;
                        break;
                    }
                }
            }
        }
    }
    // non valid bin
    if (bin_xB == -1 || bin_Q2 == -1) {
        return -1; 
    }

    /*
    COME USARLA
    auto [bin_Q2, bin_xB] = getBinIndex_xQ2(xB, Q2); // C++17
    // oppure
    std::pair<int, int> result = getBinIndex_xQ2(xB, Q2);
    int bin_xB = result.first;
    int bin_Q2 = result.second;
    */
    return bin_xB;  
}

int getBinIndex_zPt(double z, double Pt) {

    double Clas_PtBins[][2] = {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 1.25}};
    double Clas_zBins_12[][2] = {{0, 0.15}, {0.15, 0.25}, {0.25, 0.35}, {0.35, 0.45}, {0.45, 0.55}, {0.55, 0.7}, {0.7, 0.9}};
    double Clas_zBins_34[][2] = {{0, 0.25}, {0.25, 0.35}, {0.35, 0.45}, {0.45, 0.55}, {0.55, 0.7}, {0.7, 0.9}};
    double Clas_zBins_5[][2] = {{0, 0.35}, {0.35, 0.45}, {0.45, 0.55}, {0.55, 0.9}};

    int bin_z = -1;
    int bin_Pt = -1;
    for (int i = 0; i < 5; i++) {
        if(Pt >= Clas_PtBins[i][0] && Pt < Clas_PtBins[i][1]) {
            if(i == 0){
                for (int j = 0; j < 7; j++){
                    if(z >= Clas_zBins_12[j][0] && z < Clas_zBins_12[j][1]){
                        bin_Pt = i;
                        bin_z = j;
                        break;
                    }
                }
            } else if(i == 1){
                for (int j = 0; j < 7; j++){
                    if(z >= Clas_zBins_12[j][0] && z < Clas_zBins_12[j][1]){
                        bin_Pt = i;
                        bin_z = j + 7;
                        break;
                    }
                }
            } else if(i == 2){
                for (int j = 0; j < 6; j++){
                    if(z >= Clas_zBins_34[j][0] && z < Clas_zBins_34[j][1]){
                        bin_Pt = i;
                        bin_z = j + 14 ;
                        break;
                    }
                }
            } else if(i == 3){
                for (int j = 0; j < 6; j++){
                    if(z >= Clas_zBins_34[j][0] && z < Clas_zBins_34[j][1]){
                        bin_Pt = i;
                        bin_z = j + 20 ;
                        break;
                    }
                }
            } else if(i == 4){
                for (int j = 0; j < 4; j++){
                    if(z >= Clas_zBins_5[j][0] && z < Clas_zBins_5[j][1]){
                        bin_Pt = i;
                        bin_z = j + 26;
                        break;
                    }
                }
            } 
        }
    }

    if (bin_z == -1 || bin_Pt == -1) {
        return -1; 
    }

    return bin_z;  
}



struct HermesPoint {
    float Q2;
    float x;
    float z;
    float Pt;
    float Aut;
    float stat_err;
};

std::vector<HermesPoint> ReadHermesData(const std::string& filename) {
    std::ifstream infile(filename);
    std::vector<HermesPoint> data;
    std::string line;

    while (std::getline(infile, line)) {
        // Skip comments or blank lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string dummy1, dummy2, dummy3;
        HermesPoint point;
        float dummy_vals[5]; // for columns we want to skip

        iss >> dummy1 >> dummy2 >> dummy3    // skip bin ranges
            >> point.Q2 >> point.x >> dummy_vals[0]
            >> point.z >> point.Pt >> dummy_vals[1]
            >> point.Aut                     // asymmetry value 
            >> point.stat_err;               // statistical error

        data.push_back(point);
    }

    return data;
}

struct CompassPoint {
    float x;
    float z;
    float Pt;
    float Q2;
    float Aut;
    float stat_err;
};

std::vector<CompassPoint> ReadCompassData(const std::string& filename) {
    std::ifstream infile(filename);
    std::vector<CompassPoint> data;
    std::string line;

    while (std::getline(infile, line)) {
        // Skip comments or blank lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string dummy1, dummy2, dummy3;
        CompassPoint point;
        float dummy_vals[5]; // for columns we want to skip

        iss >> dummy1 >> dummy2                  // skip bin ranges
            >> point.x >> dummy_vals[0] >> point.z 
            >> point.Pt >> dummy_vals[1]
            >> point.Q2 >> point.Aut             // asymmetry value 
            >> point.stat_err;                   // statistical error

        data.push_back(point);
    }

    return data;
}

void AddFilesToChain(TChain& chain, const string& dir) {
    for (const auto &entry : fs::directory_iterator(dir)) {
        if (entry.path().extension() == ".root") {
            string fileName = entry.path().filename().string();
            if (fileName.rfind("output", 0) == 0) {
                string filePath = entry.path().string();
                chain.Add(Form("%s/PionTree", filePath.c_str()));
            }
        }
    }
}

double Parametrization_Aut_Sivers_pip(double x, double Q2, double z, double Pt){
  double par_Aut_siv_pip = (0.045228 + 0.347852*x - 0.010543*Q2 + 0.070264*z + 0.0100019*Pt + 0.089687*x*x + 0.007352*z*z - 0.126305*Pt*Pt + 0.241113*x*z + 0.424116*x*Pt + 0.296046*z*Pt);
  if(par_Aut_siv_pip > 0.8) par_Aut_siv_pip = 0.8;
  if(par_Aut_siv_pip < -0.8) par_Aut_siv_pip = -0.8;
  return par_Aut_siv_pip;
}
double Parametrization_Aut_Collins_pip(double x, double Q2, double z, double Pt){
  double par_Aut_col_pip = (0.022848 + 0.388968*x - 0.003056*Q2 - 0.018206*z + 0.046010*Pt - 1.374359*x*x + 0.164153*z*z - 0.037786*Pt*Pt + 0.881125*x*z - 0.314619*x*Pt - 0.030364*z*Pt);
  if(par_Aut_col_pip > 0.8) par_Aut_col_pip = 0.8;
  if(par_Aut_col_pip < -0.8) par_Aut_col_pip = -0.8;
  return par_Aut_col_pip;
}
double PolarizFunction_pip(double Phi_coll, double Phi_siv, double A_coll, double A_siv, double Dnn, double A_pretz, double Phi_pretz, double eps){
  double funct = 0.5*(1 + Dnn*(A_coll*TMath::Sin(Phi_coll) + A_siv*TMath::Sin(Phi_siv) + A_pretz*TMath::Sin(Phi_pretz)));
  return funct;
}

/*
double fq_GRV98(double x, std::string flavour) {
    if (x <= 0.0001 || x >= 0.999) return 0.0;
    if (flavour == "u") {
        // u valence (approssimazione GRV98LO a Q² ~ 2.5 GeV²)
        return 2 * pow(x, -0.5) * pow(1 - x, 3); 
    } else if (flavour == "d") {
        // d valence
        return 1 * pow(x, -0.5) * pow(1 - x, 4);
    } else if (flavour == "sea") {
        // sea
        return 0.5 * pow(x, -0.2) * pow(1 - x, 7);
    }
    return 0.0;
}
*/
// SE USIAMO LHAPDF POSSIAMO USARE LA PDF DEL PROTONE DA:
//cteq6l1, NNPDF31_lo_as_0130, HERAPDF15LO_EIG
// per scaricarle:
// cd /Users/lorenzopolizzi/Software/lhapdf/share/LHAPDF
// wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/NNPDF31_lo_as_0130.tar.gz
// tar -xzf NNPDF31_lo_as_0130.tar.gz

// PDF collineare non_TMD
const PDF* pdf = LHAPDF::mkPDF("NNPDF31_lo_as_0130", 0); // provo con la replica centrale
double fq_GRV98(double x, const std::string& flavour, const PDF* pdf) {
    if (x <= 0.00001 || x >= 0.999) return 0.0;
    int pid = 0;
    double Q = 2.5;
    if (flavour == "u")      pid = 2;   // u quark
    else if (flavour == "d") pid = 1;   // d quark
    else if (flavour == "sea") {
        // gluon has 21 as pdg
        // Somma delle antiquark: ū (−2), d̄ (−1), s̄ (−3)
        double qbar_u = pdf->xfxQ(-2, x, Q)/x; // x*ū
        double qbar_d = pdf->xfxQ(-1, x, Q)/x; // x*d̄
        double qbar_s = pdf->xfxQ(-3, x, Q)/x; // x*s̄
        double q_s = pdf->xfxQ(3, x, Q)/x;
        return qbar_u + qbar_d + qbar_s + q_s;
    } else {
        return 0.0;
    }
    return pdf->xfxQ(pid, x, Q)/x; // x*f(x,Q)
}
// Questa è la funzione di Sivers Delta^N f_{q/P}
double SiversModel(double x, double Pt, const double* params) { 
    // params = {Nu, Nd, Nsea, alpha, beta, M1}
    double Nu   = params[0];
    double Nd   = params[1];
    double Nsea = params[2];
    double alpha_u   = params[3];
    double alpha_d   = params[4];
    double alpha_sea = params[5];
    double beta  = params[6];
    double M1    = params[7];
    double kT2_avg = 0.25;
    if (x <= 0.0 || x >= 1.0 || M1 <= 0.0) return 1e10;
    double norm_u   = pow(alpha_u + beta, alpha_u + beta) / (pow(alpha_u, alpha_u) * pow(beta, beta));
    double norm_d   = pow(alpha_d + beta, alpha_d + beta) / (pow(alpha_d, alpha_d) * pow(beta, beta));
    double norm_sea = pow(alpha_sea + beta, alpha_sea + beta) / (pow(alpha_sea, alpha_sea) * pow(beta, beta));
    double Nx_u   = Nu   * pow(x, alpha_u)   * pow(1 - x, beta) * norm_u;
    double Nx_d   = Nd   * pow(x, alpha_d)   * pow(1 - x, beta) * norm_d;
    double Nx_sea = Nsea * pow(x, alpha_sea) * pow(1 - x, beta) * norm_sea;

    double fq_u   = fq_GRV98(x, "u", pdf);
    double fq_d   = fq_GRV98(x, "d", pdf);
    double fq_sea = fq_GRV98(x, "sea", pdf);

    double gauss_kt = exp(-Pt * Pt / kT2_avg) / (M_PI * kT2_avg);
    double hkt = sqrt(2 * M_E) * Pt / M1 * exp(-Pt * Pt / (M1 * M1));
    double f_sivers = Nx_u * fq_u + Nx_d * fq_d + 2.0 * Nx_sea * fq_sea;

    return 2.0 * hkt * f_sivers * gauss_kt;
}


// === CHI2 ===
double Chi2Sivers(const double* params, const std::vector<double>& x_vals, const std::vector<double>& Pt_vals, const std::vector<double>& A_vals, const std::vector<double>& A_err) {
    double chi2 = 0.0;
    for (size_t i = 0; i < x_vals.size(); ++i) {
        double model = SiversModel(x_vals[i], Pt_vals[i], params);
        double diff = A_vals[i] - model;
        double err = A_err[i];
        if (err < 0.001) err = 0.001; // protezione
        chi2 += diff * diff / (err * err);
        //if (std::abs(params[0]) > 1) chi2 += 1e4;
        //if (std::abs(params[0]) < 0.01) chi2 += 1e4; // Nq troppo vicino a zero? inutile
        //if (params[3] < 0.03 ) chi2 += 1e4;   // penalizza M1 troppo largo
        //if (params[2] < 1 || params[2] > 7) chi2 += 1e2;   // penalizza beta troppo alta
        //if (params[1] < 0.2 || params[1] > 4) chi2 += 1e3;
    }
    return chi2;
}
randGen.SetSeed(98732); 
// 98763 bello per down
// 98762 bello per up
// 98758 forse buono per tutto, errore da controllare, non da girare
// 98753 ricorda che devi girare mi sa
// 98742 bello se non mettiamo il -
// 98740 down troppo piccolo
// 98739 limite 0.2, interessante
// 98732 a me piace, forse picchi a x alti 
// === FIT CON MINUIT2 ===
void FitSiversFunction(const std::vector<double>& x_vals, const std::vector<double>& Pt_vals,const std::vector<double>& A_vals, 
                       const std::vector<double>& A_err, double* best_params_out,double cov_matrix[8][8])  {
    ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
    minimizer.SetMaxFunctionCalls(20000);
    minimizer.SetMaxIterations(10000);
    minimizer.SetTolerance(0.0005);

    auto chi2_lambda = [&](const double* p) {
        return Chi2Sivers(p, x_vals, Pt_vals, A_vals, A_err);
    };

    ROOT::Math::Functor functor(chi2_lambda, 8);
    minimizer.SetFunction(functor);

    minimizer.SetLimitedVariable(0, "Nu",       -0.3,  0.01,  -1.0, 1.0); // -0.3
    minimizer.SetLimitedVariable(1, "Nd",        0.6,  0.01,  -1.0, 1.0); //  0.6
    minimizer.SetLimitedVariable(2, "Nsea",      0.3,  0.01,  -1.0, 1.0); //  0.3 nice
    minimizer.SetLimitedVariable(3, "alpha_u",   1.0,  0.02,   0.2, 2);
    minimizer.SetLimitedVariable(4, "alpha_d",   1.0,  0.02,   0.2, 2);
    minimizer.SetLimitedVariable(5, "alpha_sea", 1.0,  0.02,   0.2, 2); // 0.5-1.5 o 0.1-1.5
    minimizer.SetLimitedVariable(6, "beta",      3.0,  0.02,   1, 5.0);
    minimizer.SetLimitedVariable(7, "M1",        0.6,  0.01,   0.2, 1.0);

    bool success = minimizer.Minimize();
    int N_params = 8;
    int N_points = x_vals.size();
    double chi2_dof = minimizer.MinValue() / (N_points - N_params);
    std::cout << "Fit status: " << minimizer.Status() << std::endl;
    std::cout << "Final chi2/dof: " << chi2_dof << std::endl;

    if (!success)
        std::cerr << "Warning: Fit did not converge!" << std::endl;

    for (int i = 0; i < 8; ++i)
        best_params_out[i] = minimizer.X()[i];

    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
            cov_matrix[i][j] = minimizer.CovMatrix(i, j);
}

double SiversModel_u(double x, double Pt, const double* p) {
    double Nq = p[0];
    double alpha = p[3];
    double beta = p[6];
    double M1 = p[7];
    double kT2_avg = 0.25;

    if (x <= 0 || x >= 1 || M1 <= 0) return 1e10;

    double norm = pow(alpha + beta, alpha + beta) / (pow(alpha, alpha) * pow(beta, beta));
    double Nx = Nq * pow(x, alpha) * pow(1 - x, beta) * norm;

    double hkt = sqrt(2 * M_E) * Pt / M1 * exp(-Pt * Pt / (M1 * M1));
    double fq = fq_GRV98(x, "u", pdf);
    double fqkT = fq * exp(-Pt * Pt / kT2_avg) / (M_PI * kT2_avg);

    return 2.0 * Nx * hkt * fqkT;
}

double SiversModel_d(double x, double Pt, const double* p) {
    double Nq = p[1];
    double alpha = p[4];
    double beta = p[6];
    double M1 = p[7];
    double kT2_avg = 0.25;

    if (x <= 0 || x >= 1 || M1 <= 0) return 1e10;

    double norm = pow(alpha + beta, alpha + beta) / (pow(alpha, alpha) * pow(beta, beta));
    double Nx = Nq * pow(x, alpha) * pow(1 - x, beta) * norm;

    double hkt = sqrt(2 * M_E) * Pt / M1 * exp(-Pt * Pt / (M1 * M1));
    double fq = fq_GRV98(x, "d", pdf);
    double fqkT = fq * exp(-Pt * Pt / kT2_avg) / (M_PI * kT2_avg);

    return 2.0 * Nx * hkt * fqkT;
}

double SiversModel_sea(double x, double Pt, const double* p) {
    double Nq = p[2];     // normalmente si assume 0 o fit separato
    double alpha = p[5];  // dummy params
    double beta = p[6];
    double M1 = p[7];
    double kT2_avg = 0.25;

    if (x <= 0 || x >= 1 || M1 <= 0) return 1e10;

    double norm = pow(alpha + beta, alpha + beta) / (pow(alpha, alpha) * pow(beta, beta));
    double Nx = Nq * pow(x, alpha) * pow(1 - x, beta) * norm;

    double hkt = sqrt(2 * M_E) * Pt / M1 * exp(-Pt * Pt / (M1 * M1));
    double fq = 2.0 * fq_GRV98(x, "sea", pdf); // somma anti-up + anti-down
    double fqkT = fq * exp(-Pt * Pt / kT2_avg) / (M_PI * kT2_avg);

    return 2.0 * Nx * hkt * fqkT;
}


double DerivativeSiversModel(double x, double Pt, const double* params, int i_param, double h = 1e-5) {
    double p_forward[8], p_backward[8];
    for (int i = 0; i < 8; ++i) {
        p_forward[i] = params[i];
        p_backward[i] = params[i];
    }
    p_forward[i_param]  += h;
    p_backward[i_param] -= h;

    double f_plus  = SiversModel(x, Pt, p_forward);
    double f_minus = SiversModel(x, Pt, p_backward);
    return (f_plus - f_minus) / (2.0 * h);
}
TGraphErrors* BuildSiversErrorBand(double Pt_fixed,const double* best_params,double cov_matrix[8][8], int up) {
    const int n_points = 1000;
    TGraphErrors* error_band = new TGraphErrors(n_points);
    for (int i = 0; i < n_points; ++i) {
        if (up == 1){
            double x = 0.001 + i * (0.7 - 0.001) / (n_points - 1);
            double fx = SiversModel(x, Pt_fixed, best_params);

            double grad[8];
            for (int j = 0; j < 8; ++j)
                grad[j] = DerivativeSiversModel(x, Pt_fixed, best_params, j);

            double variance = 0.0;
            for (int j = 0; j < 8; ++j)
                for (int k = 0; k < 8; ++k)
                    variance += grad[j] * grad[k] * cov_matrix[j][k];

            error_band->SetPoint(i, x, fx);
            error_band->SetPointError(i, 0.0, 1.96*sqrt(std::abs(variance))); //1.96 to be at 95% CL
        } else if (up == -1) {
            double x = 0.001 + i * (0.999 - 0.001) / (n_points - 1);
            double fx = SiversModel(x, Pt_fixed, best_params);

            double grad[8];
            for (int j = 0; j < 8; ++j)
                grad[j] = DerivativeSiversModel(x, Pt_fixed, best_params, j);

            double variance = 0.0;
            for (int j = 0; j < 8; ++j)
                for (int k = 0; k < 8; ++k)
                    variance += grad[j] * grad[k] * cov_matrix[j][k];
            error_band->SetPoint(i, x, -fx);
            error_band->SetPointError(i, 0.0, 1.96*sqrt(std::abs(variance))); //1.96 to be at 95% CL
        }
    }

    error_band->SetFillColorAlpha(kRed, 0.25);
    error_band->SetLineColor(kRed);
    error_band->SetLineWidth(1);
    return error_band;
}

TGraphErrors* BuildSiversErrorBand_flavour(double Pt_fixed, const double* params, double cov[8][8], double flavour, double confidence) {
    const int n_points = 5000;
    TGraphErrors* graph = new TGraphErrors(n_points);

    for (int i = 0; i < n_points; ++i) {
        double x = 0.002 + i * (0.85 - 0.002) / (n_points - 1);
        double f;
        if (flavour == 0) f = SiversModel_u(x, Pt_fixed, params);
        else if (flavour == 1) f = SiversModel_d(x, Pt_fixed, params);
        else if (flavour == 2) f = SiversModel_sea(x, Pt_fixed, params);

        // Derivata rispetto a Nu, α, β, M1 = p[0], p[3], p[4], p[5]
        double grad[4];
        double h = 1e-5;
        for (int j = 0; j < 4; ++j) {
            double p_forward[8], p_backward[8];
            for (int k = 0; k < 8; ++k) {
                p_forward[k] = params[k];
                p_backward[k] = params[k];
            }
            int index;
            if (flavour == 0) index = (j == 0) ? 0 : (j == 1) ? 3 : j + 4;
            else if (flavour == 1) index = (j == 0) ? 1 : (j == 1) ? 4 : j + 4;
            else if (flavour == 2) index = (j == 0) ? 2 : (j == 1) ? 5 : j + 4;
            p_forward[index] += h;
            p_backward[index] -= h;
            if (flavour == 0) grad[j] = (SiversModel_u(x, Pt_fixed, p_forward) - SiversModel_u(x, Pt_fixed, p_backward)) / (2 * h);
            else if (flavour == 1) grad[j] = (SiversModel_d(x, Pt_fixed, p_forward) - SiversModel_d(x, Pt_fixed, p_backward)) / (2 * h);
            else if (flavour == 2) grad[j] = (SiversModel_sea(x, Pt_fixed, p_forward) - SiversModel_sea(x, Pt_fixed, p_backward)) / (2 * h);
        }

        double variance = 0.0;
        for (int j = 0; j < 4; ++j) {
            int idx_j;
            if (j == 0) idx_j = flavour;
            else if (j == 1) idx_j = 3 + flavour;
            else if (j == 2) idx_j = 6;
            else if (j == 3) idx_j = 7;
            for (int k = 0; k < 4; ++k) {
                int idx_k;
                if (k == 0) idx_k = flavour;
                else if (k == 1) idx_k = 3 + flavour;
                else if (k == 2) idx_k = 6;
                else if (k == 3) idx_k = 7;
                variance += grad[j] * grad[k] * cov[idx_j][idx_k];
            }
        }

        graph->SetPoint(i, x, -f*x);
        graph->SetPointError(i, 0.0, confidence * sqrt(std::abs(variance))*x);
    }

    return graph;
}

TGraphErrors* BuildSiversErrorBand_flavour_Pt(double xB_fixed, const double* params, double cov[8][8], double flavour, double confidence) {
    const int n_points = 5000;
    TGraphErrors* graph = new TGraphErrors(n_points);

    for (int i = 0; i < n_points; ++i) {
        double Pt = 0.002 + i * (0.95 - 0.002) / (n_points - 1);
        double f;
        if (flavour == 0) f = SiversModel_u(xB_fixed, Pt, params);
        else if (flavour == 1) f = SiversModel_d(xB_fixed, Pt, params);
        else if (flavour == 2) f = SiversModel_sea(xB_fixed, Pt, params);

        // Derivata rispetto a Nu, α, β, M1 = p[0], p[3], p[4], p[5]
        double grad[4];
        double h = 1e-5;
        for (int j = 0; j < 4; ++j) {
            double p_forward[8], p_backward[8];
            for (int k = 0; k < 8; ++k) {
                p_forward[k] = params[k];
                p_backward[k] = params[k];
            }
            int index;
            if (flavour == 0) index = (j == 0) ? 0 : (j == 1) ? 3 : j + 4;
            else if (flavour == 1) index = (j == 0) ? 1 : (j == 1) ? 4 : j + 4;
            else if (flavour == 2) index = (j == 0) ? 2 : (j == 1) ? 5 : j + 4;
            p_forward[index] += h;
            p_backward[index] -= h;
            if (flavour == 0) grad[j] = (SiversModel_u(xB_fixed, Pt, p_forward) - SiversModel_u(xB_fixed, Pt, p_backward)) / (2 * h);
            else if (flavour == 1) grad[j] = (SiversModel_d(xB_fixed, Pt, p_forward) - SiversModel_d(xB_fixed, Pt, p_backward)) / (2 * h);
            else if (flavour == 2) grad[j] = (SiversModel_sea(xB_fixed, Pt, p_forward) - SiversModel_sea(xB_fixed, Pt, p_backward)) / (2 * h);
        }

        double variance = 0.0;
        for (int j = 0; j < 4; ++j) {
            int idx_j;
            if (j == 0) idx_j = flavour;
            else if (j == 1) idx_j = 3 + flavour;
            else if (j == 2) idx_j = 6;
            else if (j == 3) idx_j = 7;
            for (int k = 0; k < 4; ++k) {
                int idx_k;
                if (k == 0) idx_k = flavour;
                else if (k == 1) idx_k = 3 + flavour;
                else if (k == 2) idx_k = 6;
                else if (k == 3) idx_k = 7;
                variance += grad[j] * grad[k] * cov[idx_j][idx_k];
            }
        }

        graph->SetPoint(i, Pt, -f*xB_fixed);
        graph->SetPointError(i, 0.0, confidence * sqrt(std::abs(variance))*xB_fixed);
    }

    return graph;
}



// _________________________________________________ COLLINS ___________________________________________________________________________


double Transversity(double x, const double* p, std::string flavour) {
    if (x <= 0.0 || x >= 1.0) return 0.0;
    double N_Tu    = p[0];   // N_T^u
    double N_Td    = p[1];   // N_T^d
    double N_Tsea  = p[2];   // N_T^sea
    double alpha   = p[3];   // alpha common
    double beta    = p[4];   // beta common
    double N_T = 0.0;
    if (flavour == 'u') N_T = N_Tu;
    else if (flavour == 'd') N_T = N_Td;
    else N_T = N_Tsea * 0.1; // to reduce the impact
    double norm = pow(alpha + beta, alpha + beta) / (pow(alpha, alpha) * pow(beta, beta));

    return N_T * pow(x, alpha) * pow(1 - x, beta) * norm;
}
const PDF* delta_q_pdf = mkPDF("JAM20-SIDIS_PDF_proton_nlo", 0);
double Delta_q(double x, const std::string &flavour, const PDF* pdf) {
    if (x <= 0.0 || x >= 1.0) return 0.0;
    double Q = 2.5;
    /*
    double N_q, alpha_q, beta_q, gamma_q;
    if (flavour == "u") {
        N_q = 0.7; alpha_q = 0.5; beta_q = 3.0; gamma_q = 1.0;
    } else if (flavour == "d") {
        N_q = -0.5; alpha_q = 1.0; beta_q = 4.0; gamma_q = 0.5;
    } else if (flavour == "sea") {
        N_q = -0.1; alpha_q = 1.5; beta_q = 5.0; gamma_q = 0.0;
    } else {
        return 0.0;
    }
    return N_q * pow(x, alpha_q) * pow(1.0 - x, beta_q) * (1.0 + gamma_q * x);
    */
   if (flavour == "u") return pdf->xfxQ(2, x, Q)/x; // le pdf restituiscono sempre x*PDF
   else if (flavour == "d") return pdf->xfxQ(1, x, Q)/x;
   else if (flavour == "sea"){
        double ubar = pdf->xfxQ(-2, x, Q)/x;
        double dbar = pdf->xfxQ(-1, x, Q)/x;
        double sbar = pdf->xfxQ(-3, x, Q)/x;
        double s = pdf->xfxQ(3, x, Q)/x;
        return ubar + dbar + sbar + s;
   }
   else return 0.0;
}
const PDF* D1_frag_pdf = mkPDF("JAM20-SIDIS_FF_pion_nlo", 0);
double D1_fragmentation(double z, const std::string& flavour, const PDF* pdf) {
    if (z <= 0.0 || z >= 1.0) return 0.0;
    double Q = 2.5;
    if (flavour == "u") {
        return pdf->xfxQ(2, z, Q)/z;
        //return pow(z, -0.8) * pow(1 - z, 1.5);
    } else if (flavour == "d") {
        return pdf->xfxQ(1, z, Q)/z;
        //return 0.5 * pow(z, -0.9) * pow(1 - z, 2.5);
    } else {
        double ubar = pdf->xfxQ(-2, z, Q)/z;
        double dbar = pdf->xfxQ(-1, z, Q)/z;
        double sbar = pdf->xfxQ(-3, z, Q)/z;
        double s = pdf->xfxQ(3, z, Q)/z;
        return ubar + dbar + sbar + s;
        //return 0.3 * pow(z, -1.0) * pow(1 - z, 3.0);
    }
}

double H1perp(double z, double Pt, double N_C, double gamma_C, double delta_C, double MH, std::string flavour) {
    if (z <= 0.0 || z >= 1.0) return 0.0;
    double norm = pow(gamma_C + delta_C, gamma_C + delta_C) / (pow(gamma_C, gamma_C) * pow(delta_C, delta_C));
    double Nz; // only N_fav has gamma and delta, N_dis is just the free parameter -> following the Anselmino paper
    if (flavour == 'u') Nz = N_C * pow(z, gamma_C) * pow(1 - z, delta_C) * norm;
    else Nz = N_C * norm;
    double hpt = sqrt(2 * M_E) * Pt / MH * exp(-Pt * Pt / (MH * MH));
    return Nz * D1_fragmentation(z, flavour, D1_frag_pdf) * hpt;
}

double DeltaTildeD(double z, double N_C, double gamma_C, double delta_C, double MH, std::string flavour) {
    if (z <= 0.0 || z >= 1.0) return 0.0;
    double norm = pow(gamma_C + delta_C, gamma_C + delta_C) / (pow(gamma_C, gamma_C) * pow(delta_C, delta_C));
    double D1 = D1_fragmentation(z, flavour, D1_frag_pdf);
    double Nz = (flavour == "u") ? N_C * pow(z, gamma_C) * pow(1 - z, delta_C) : N_C;
    // Integrale trasversale
    double prefactor = sqrt(2 * M_E) * MH * MH / z;
    return 2 * Nz * norm * D1; // da capire o meno se ci va il prefactor
}


/*
double CollinsModel(double x, double z, double Pt, const double* p) {
    // parametri
    double Nu   = p[0];
    double Nd   = p[1];
    double Nsea = p[2];
    double N_C_fav = p[5];
    double N_C_dis = p[6];
    double gamma_C = p[7];
    double delta_C = p[8];
    double MH = p[9];
    // transversity
    double hu   = Transversity(x, p, "u");
    double hd   = Transversity(x, p, "d");
    double hsea = Transversity(x, p, "sea");
    // frammentazioni
    double D_factor = (2*Pt)/(z*MH);  // pass from Delta D to H^perp
    double Hu   = D_factor*H1perp(z, Pt, N_C_fav, gamma_C, delta_C, MH, "u"); // favoured
    double Hd   = D_factor*H1perp(z, Pt, N_C_dis, gamma_C, delta_C, MH, "d"); // disfavoured
    double Hsea = D_factor*H1perp(z, Pt, N_C_dis, gamma_C, delta_C, MH, "sea"); // di solito anche sea ~ disfavoured

    double gaussian_Pt = exp(-Pt * Pt / 0.181) / (M_PI * 0.181); // remove the Pt dependence
    double H_tot = (hu * Hu + hd * Hd + hsea * Hsea);
    return H_tot*gaussian_Pt;
}
*/

// maybe more solid with everything inside!
double CollinsAsymmetry(double x, double z, double Pt, double y, const double* p) {
    // Parametri trasversità e FF
    double Nu   = p[0];
    double Nd   = p[1];
    double Nsea = p[2];
    double N_C_fav = p[5];
    double N_C_dis = p[6];
    double gamma_C = p[7];
    double delta_C = p[8];
    double MH = p[9];
    // Larghezze trasverse (puoi cambiarle!)
    double pT2_avg = 0.181;
    double kT2_avg = 0.25;
    // Larghezza convoluzione di Collins
    double pT2_C = (MH * MH * pT2_avg) / (MH * MH + pT2_avg);
    double PT2_C_total = pT2_C + z * z * kT2_avg;
    double PT2_total = pT2_avg + z * z * kT2_avg;
    // Fattore cinetico
    double prefactor_kin_num = sqrt(2 * M_E) * (Pt / MH) * (pT2_C * pT2_C / pT2_avg) * exp(-Pt * Pt / PT2_C_total) / (PT2_C_total * PT2_C_total);
    double prefactor_kin_den = (exp(-Pt * Pt / PT2_total) / PT2_total);
    double prefactor_kin = prefactor_kin_num / prefactor_kin_den;
    // Fattore y
    double prefactor_y = (1.0 - y) / (1.0 + pow(1.0 - y, 2));
    // Transversità
    double Dq_u = Transversity(x, p, "u") * 0.5 * (fq_GRV98(x, "u", pdf) + Delta_q(x, "u", delta_q_pdf));
    double Dq_d = Transversity(x, p, "d") * 0.5 * (fq_GRV98(x, "d", pdf) + Delta_q(x, "d", delta_q_pdf));
    double Dq_sea = Transversity(x, p, "sea") * 0.5 * (fq_GRV98(x, "sea", pdf) + Delta_q(x, "sea", delta_q_pdf));
    // FF di Collins
    double D_fav = DeltaTildeD(z, N_C_fav, gamma_C, delta_C, MH, "u"); // favoured
    double D_dis = DeltaTildeD(z, N_C_dis, gamma_C, delta_C, MH, "d"); // disfavoured both for down and sea
    // Charge q
    double e_u = 2.0/3.0;
    double e_d = -1.0/3.0; 
    double e_sea_2 = 7.0/36.0; // sea has s, s_bar, u_bar, d_bar -> e^2 = (4/9 + 1/9 + 1/9 + 1/9)/4 = 7/36
    // Denominatore
    double F_uu = fq_GRV98(x, "u", pdf) * D1_fragmentation(z, "u", D1_frag_pdf);
    double F_dd = fq_GRV98(x, "d", pdf) * D1_fragmentation(z, "d", D1_frag_pdf);
    double F_ss = fq_GRV98(x, "sea", pdf) * D1_fragmentation(z, "sea", D1_frag_pdf);
    double denom = pow(e_u, 2) * F_uu + pow(e_d, 2) * F_dd + e_sea_2 * F_ss;
    // Numeratore
    double num = pow(e_u, 2) * (Dq_u * D_fav) + pow(e_d, 2) * (Dq_d * D_dis) + e_sea_2 * (Dq_sea * D_dis);
    
    if(denom == 0){
        cout << "denom = " << denom << '\n';
        cout << "F_uu,dd,ss :" << F_uu << "      " << F_dd << "    " << F_ss << endl;
        cout << "A = " << e_u * e_u * F_uu << "        B = " << e_d * e_d * F_dd << "      C = " << e_sea_2 * F_ss << endl;
        cout << "x = " << x << "    z = " << z << "    Pt = " << Pt << endl;
    }
    
    
    return prefactor_kin * prefactor_y * (num / denom);
}


double Chi2Collins(const double* p, const std::vector<double>& x_vals, const std::vector<double>& z_vals,
                   const std::vector<double>& Pt_vals, const vector<double>& y_vals, const std::vector<double>& A_vals, const std::vector<double>& A_errs) {
    double chi2 = 0.0;
    
    for (size_t i = 0; i < x_vals.size(); ++i) {
        double model = CollinsAsymmetry(x_vals[i], z_vals[i], Pt_vals[i], y_vals[i], p);
        double diff = A_vals[i] - model;
        double err = A_errs[i] < 0.001 ? 0.001 : A_errs[i];
        chi2 += diff * diff / (err * err);
        if (isnan(model) || isinf(model)) {
            double F_uu = fq_GRV98(x_vals[i], "u", pdf) * D1_fragmentation(z_vals[i], "u", D1_frag_pdf);
            double F_dd = fq_GRV98(x_vals[i], "d", pdf) * D1_fragmentation(z_vals[i], "d", D1_frag_pdf);
            double F_ss = fq_GRV98(x_vals[i], "sea", pdf) * D1_fragmentation(z_vals[i], "sea", D1_frag_pdf);
            if (F_uu == 0) cout << " F_uu = " << F_uu << "   x = " << x_vals[i] << "   z = " << z_vals[i] << endl;
        }
    }
    return chi2;
}

void FitCollins(const std::vector<double>& x_vals, const std::vector<double>& z_vals,
                const std::vector<double>& Pt_vals, const std::vector<double>& A_vals,
                const std::vector<double>& A_errs, double* best_params, double cov[10][10], const std::vector<double>& y_vals) {
    ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
    minimizer.SetMaxFunctionCalls(10000);
    minimizer.SetMaxIterations(10000);
    minimizer.SetTolerance(0.0005);

    auto chi2 = [&](const double* p) {
        return Chi2Collins(p, x_vals, z_vals, Pt_vals, y_vals,A_vals, A_errs);
    };
    ROOT::Math::Functor functor(chi2, 10);
    minimizer.SetFunction(functor);

    minimizer.SetLimitedVariable(0, "Nu_T",     0.5,  0.01, -1.0,  1.0);
    minimizer.SetLimitedVariable(1, "Nd_T",    -0.3,  0.01, -1.0,  1.0);
    minimizer.SetLimitedVariable(2, "Nsea_T",   0.1,  0.01, -1.0,  1.0);
    minimizer.SetLimitedVariable(3, "alpha_T",  1.0,  0.01,  0.2,  4.0);
    minimizer.SetLimitedVariable(4, "beta_T",   2.0,  0.01,  0.5,  4.0);
    minimizer.SetLimitedVariable(5, "N_C_fav",  0.4,  0.01, -0.5,  1.0);
    minimizer.SetLimitedVariable(6, "N_C_dis",  0.2,  0.01, -1.0,  0.5);
    minimizer.SetLimitedVariable(7, "gamma_C",  2.0,  0.01,  0.0,  5.0);
    minimizer.SetLimitedVariable(8, "delta_C",  1.0,  0.01,  0.0,  5.0);
    minimizer.SetLimitedVariable(9, "MH",       0.6,  0.01,  0.2,  1.0);

    bool success = minimizer.Minimize();
    int N_params = 10;
    int N_points = x_vals.size();
    double chi2_dof = minimizer.MinValue() / (N_points - N_params);
    std::cout << "Fit status: " << minimizer.Status() << std::endl;
    std::cout << "Final chi2/dof: " << chi2_dof << std::endl;
    if (!success) std::cerr << "Warning: Fit did not converge!" << std::endl;

    for (int i = 0; i < 10; ++i) best_params[i] = minimizer.X()[i];
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            cov[i][j] = minimizer.CovMatrix(i, j);
}

TGraphErrors* BuildCollinsErrorBand(double z_fixed, double Pt_fixed, double y_fixed, const double* best_params, double cov_matrix[10][10], int up) {
    const int n_points = 1000;
    TGraphErrors* error_band = new TGraphErrors(n_points);

    for (int i = 0; i < n_points; i++) {
        double x = 0.001 + i * ((up == 1 ? 0.7 : 0.99) - 0.001) / (n_points - 1);
        double fx = CollinsAsymmetry(x, z_fixed, Pt_fixed, y_fixed, best_params);

        // Derivate numeriche su 10 parametri
        double grad[10]; double h = 1e-5;
        for (int j = 0; j < 10; j++) {
            double p_forward[10], p_backward[10];
            for (int k = 0; k < 10; k++) { p_forward[k]=best_params[k]; p_backward[k]=best_params[k]; }
            p_forward[j]  += h;
            p_backward[j] -= h;
            double f_plus  = CollinsAsymmetry(x, z_fixed, Pt_fixed, y_fixed, p_forward);
            double f_minus = CollinsAsymmetry(x, z_fixed, Pt_fixed, y_fixed, p_backward);
            grad[j] = (f_plus - f_minus) / (2.0 * h);
        }

        double variance = 0.0;
        for (int j = 0; j < 10; j++)
            for (int k = 0; k < 10; k++)
                variance += grad[j] * grad[k] * cov_matrix[j][k];

        double y = (up == 1 ? fx : -fx);
        double yerr = 1.96 * std::sqrt(std::abs(variance));
        error_band->SetPoint(i, x, y);
        error_band->SetPointError(i, 0.0, yerr);
    }
    return error_band;
}


double CollinsModel_u(double x, double z, double Pt, double y, const double* p) {
    double Nu   = p[0];
    double Nd   = p[1];
    double Nsea = p[2];
    double N_C_fav = p[5];
    double N_C_dis = p[6];
    double gamma_C = p[7];
    double delta_C = p[8];
    double MH = p[9];
    // Larghezze trasverse (puoi cambiarle!)
    double pT2_avg = 0.181;
    double kT2_avg = 0.25;
    // Larghezza convoluzione di Collins
    double pT2_C = (MH * MH * pT2_avg) / (MH * MH + pT2_avg);
    double PT2_C_total = pT2_C + z * z * kT2_avg;
    double PT2_total = pT2_avg + z * z * kT2_avg;
    // Fattore cinetico
    double prefactor_kin_num = sqrt(2 * M_E) * (Pt / MH) * (pT2_C * pT2_C / pT2_avg) * exp(-Pt * Pt / PT2_C_total) / (PT2_C_total * PT2_C_total);
    double prefactor_kin_den = (exp(-Pt * Pt / PT2_total) / PT2_total);
    double prefactor_kin = prefactor_kin_num / prefactor_kin_den;
    // Fattore y
    double prefactor_y = (1.0 - y) / (1.0 + pow(1.0 - y, 2));
    // Transversità
    double Dq_u = Transversity(x, p, "u") * 0.5 * (fq_GRV98(x, "u", pdf) + Delta_q(x, "u", delta_q_pdf));
    // FF di Collins
    double D_fav = DeltaTildeD(z, N_C_fav, gamma_C, delta_C, MH, "u"); // favoured
    // Charge q
    double e_u = 2.0/3.0;
    // Denominatore
    double F_uu = fq_GRV98(x, "u", pdf) * D1_fragmentation(z, "u", D1_frag_pdf);
    double denom = pow(e_u, 2) * F_uu ;
    // Numeratore
    double num = pow(e_u, 2) * (Dq_u * D_fav);
    
    return prefactor_kin * prefactor_y * (num / denom);
}

double DeltaT_q(double x, double z, double Pt, double y, const double* p, std::string flavour){
    double Nu   = p[0];
    double Nd   = p[1];
    double Nsea = p[2];
    double N_C_fav = p[5];
    double N_C_dis = p[6];
    double gamma_C = p[7];
    double delta_C = p[8];
    double MH = p[9];
    double deltaT_q;
    if (flavour == 'u') deltaT_q = Transversity(x, p, "u") * 0.5 * (fq_GRV98(x, "u", pdf) + Delta_q(x, "u", delta_q_pdf));
    else if (flavour == 'd') deltaT_q = Transversity(x, p, "d") * 0.5 * (fq_GRV98(x, "d", pdf) + Delta_q(x, "d", delta_q_pdf));
    else deltaT_q = Transversity(x, p, "sea") * 0.5 * (fq_GRV98(x, "sea", pdf) + Delta_q(x, "sea", delta_q_pdf));
    return deltaT_q;
}


double CollinsModel_d(double x, double z, double Pt, double y, const double* p) {
    double Nu   = p[0];
    double Nd   = p[1];
    double Nsea = p[2];
    double N_C_fav = p[5];
    double N_C_dis = p[6];
    double gamma_C = p[7];
    double delta_C = p[8];
    double MH = p[9];
    // Larghezze trasverse (puoi cambiarle!)
    double pT2_avg = 0.181;
    double kT2_avg = 0.25;
    // Larghezza convoluzione di Collins
    double pT2_C = (MH * MH * pT2_avg) / (MH * MH + pT2_avg);
    double PT2_C_total = pT2_C + z * z * kT2_avg;
    double PT2_total = pT2_avg + z * z * kT2_avg;
    // Fattore cinetico
    double prefactor_kin_num = sqrt(2 * M_E) * (Pt / MH) * (pT2_C * pT2_C / pT2_avg) * exp(-Pt * Pt / PT2_C_total) / (PT2_C_total * PT2_C_total);
    double prefactor_kin_den = (exp(-Pt * Pt / PT2_total) / PT2_total);
    double prefactor_kin = prefactor_kin_num / prefactor_kin_den;
    // Fattore y
    double prefactor_y = (1.0 - y) / (1.0 + pow(1.0 - y, 2));
    // Transversità
    double Dq_d = Transversity(x, p, "d") * 0.5 * (fq_GRV98(x, "d", pdf) + Delta_q(x, "d", delta_q_pdf));
    // FF di Collins
    double D_dis = DeltaTildeD(z, N_C_dis, gamma_C, delta_C, MH, "d"); // disfavoured both for down and sea
    // Charge q
    double e_d = -1.0/3.0; 
    // Denominatore
    double F_dd = fq_GRV98(x, "d", pdf) * D1_fragmentation(z, "d", D1_frag_pdf);
    double denom = pow(e_d, 2) * F_dd ;
    // Numeratore
    double num = pow(e_d, 2) * (Dq_d * D_dis);

    return prefactor_kin * prefactor_y * (num / denom);
}


double CollinsModel_sea(double x, double z, double Pt, double y, const double* p) {
    double Nu   = p[0];
    double Nd   = p[1];
    double Nsea = p[2];
    double N_C_fav = p[5];
    double N_C_dis = p[6];
    double gamma_C = p[7];
    double delta_C = p[8];
    double MH = p[9];
    // Larghezze trasverse (puoi cambiarle!)
    double pT2_avg = 0.181;
    double kT2_avg = 0.25;
    // Larghezza convoluzione di Collins
    double pT2_C = (MH * MH * pT2_avg) / (MH * MH + pT2_avg);
    double PT2_C_total = pT2_C + z * z * kT2_avg;
    double PT2_total = pT2_avg + z * z * kT2_avg;
    // Fattore cinetico
    double prefactor_kin_num = sqrt(2 * M_E) * (Pt / MH) * (pT2_C * pT2_C / pT2_avg) * exp(-Pt * Pt / PT2_C_total) / (PT2_C_total * PT2_C_total);
    double prefactor_kin_den = (exp(-Pt * Pt / PT2_total) / PT2_total);
    double prefactor_kin = prefactor_kin_num / prefactor_kin_den;
    // Fattore y
    double prefactor_y = (1.0 - y) / (1.0 + pow(1.0 - y, 2));
    // Transversità
    double Dq_sea = Transversity(x, p, "sea") * 0.5 * (fq_GRV98(x, "sea", pdf) + Delta_q(x, "sea", delta_q_pdf));
    // FF di Collins
    double D_dis = DeltaTildeD(z, N_C_dis, gamma_C, delta_C, MH, "sea"); // disfavoured both for down and sea
    // Charge q
    double e_sea_2 = 7.0/36.0; // sea has s, s_bar, u_bar, d_bar -> e^2 = (4/9 + 1/9 + 1/9 + 1/9)/4 = 7/36
    // Denominatore
    double F_ss = fq_GRV98(x, "sea", pdf) * D1_fragmentation(z, "sea", D1_frag_pdf);
    double denom = e_sea_2 * F_ss;
    // Numeratore
    double num = e_sea_2 * (Dq_sea * D_dis);
       
    return prefactor_kin * prefactor_y * (num / denom);
}


TGraphErrors* BuildCollinsErrorBand_flavour(double z_fixed, double Pt_fixed, double y_fixed, const double* params, double cov[10][10], int flavour, double confidence) {
    const int n_points = 5000;
    TGraphErrors* graph = new TGraphErrors(n_points);

    // Indici dei parametri da variare in base al flavour
    int paramIndices[7]; 
    if (flavour == 0) { // up
        paramIndices[0] = 0; // Nu_T
        paramIndices[1] = 3; // alpha_T
        paramIndices[2] = 4; // beta_T
        paramIndices[3] = 5; // N_C_fav
        paramIndices[4] = 7; // gamma_C
        paramIndices[5] = 8; // delta_C
        paramIndices[6] = 9; // M_C
    } else if (flavour == 1) { // down
        paramIndices[0] = 1; // Nd_T
        paramIndices[1] = 3; // alpha_T
        paramIndices[2] = 4; // beta_T
        paramIndices[3] = 6; // N_C_dis
        paramIndices[4] = 7; // gamma_C
        paramIndices[5] = 8; // delta_C
        paramIndices[6] = 9; // M_C
    } else if (flavour == 2) { // sea
        paramIndices[0] = 2; // Nsea_T
        paramIndices[1] = 3; // alpha_T
        paramIndices[2] = 4; // beta_T
        paramIndices[3] = 6; // N_C_dis
        paramIndices[4] = 7; // gamma_C
        paramIndices[5] = 8; // delta_C
        paramIndices[6] = 9; // M_C
    }

    double h = 1e-5;
    for (int i = 0; i < n_points; ++i) {
        double x = 0.002 + i * (0.9 - 0.002) / (n_points - 1);
        double f;
        if (flavour == 0) f = DeltaT_q(x, z_fixed, Pt_fixed, y_fixed, params, "u");
        else if (flavour == 1) f = DeltaT_q(x, z_fixed, Pt_fixed, y_fixed, params, "d");
        else f = DeltaT_q(x, z_fixed, Pt_fixed, y_fixed, params, "sea");

        // Calcolo del gradiente numerico
        double grad[7]; 
        for (int j = 0; j < 7; ++j) {
            int index = paramIndices[j];
            double p_forward[10], p_backward[10];
            for (int k = 0; k < 10; ++k) {
                p_forward[k] = params[k];
                p_backward[k] = params[k];
            }
            p_forward[index] += h;
            p_backward[index] -= h;

            double f_plus, f_minus;
            if (flavour == 0) {
                f_plus = DeltaT_q(x, z_fixed, Pt_fixed, y_fixed, p_forward, "u");
                f_minus = DeltaT_q(x, z_fixed, Pt_fixed, y_fixed, p_backward, "u");
            } else if (flavour == 1) {
                f_plus = f = DeltaT_q(x, z_fixed, Pt_fixed, y_fixed, p_forward, "d");
                f_minus = f = DeltaT_q(x, z_fixed, Pt_fixed, y_fixed, p_backward, "d");
            } else {
                f_plus = DeltaT_q(x, z_fixed, Pt_fixed, y_fixed, p_forward, "sea");
                f_minus = DeltaT_q(x, z_fixed, Pt_fixed, y_fixed, p_backward, "sea");
            }

            grad[j] = (f_plus - f_minus) / (2 * h);
        }

        // Propagazione dell'errore
        double variance = 0.0;
        for (int j = 0; j < 7; ++j) {
            int idx_j = paramIndices[j];
            for (int k = 0; k < 7; ++k) {
                int idx_k = paramIndices[k];
                variance += grad[j] * grad[k] * cov[idx_j][idx_k];
            }
        }

        graph->SetPoint(i, x, f*x); 
        graph->SetPointError(i, 0.0, confidence * std::sqrt(std::abs(variance)) * x); 
    }

    return graph;
}

double compute_weighted_sivers_mean(const std::vector<double>& param_siv,
                                    const std::vector<double>& phi_h, const std::vector<double>& phi_s, const std::vector<double>& depol) {
    double num = 0.0;
    double denom = 0.0;
    for (size_t i = 0; i < param_siv.size(); ++i) {
        double delta_phi = phi_h[i] - phi_s[i];
        double weight = depol[i] * std::sin(delta_phi);
        num += weight * param_siv[i];
        denom += weight;
    }

    if (std::abs(denom) < 1e-8) return 0.0; // Evita divisione per zero
    return num / denom;
}

double compute_weighted_collins_mean(const std::vector<double>& param_siv,
                                    const std::vector<double>& phi_h, const std::vector<double>& phi_s, const std::vector<double>& depol) {
    double num = 0.0;
    double denom = 0.0;
    for (size_t i = 0; i < param_siv.size(); ++i) {
        double delta_phi = phi_h[i] + phi_s[i];
        double weight = depol[i] * std::sin(delta_phi);
        num += weight * param_siv[i];
        denom += weight;
    }

    if (std::abs(denom) < 1e-8) return 0.0; // Evita divisione per zero
    return num / denom;
}

// _____________________________________________________________________________________________________________________________________




void misure_rgh_pion(){

    // quantity
    float event_number, epsilon, gamma, helicity;
    float sector_e, sector_pi;
    float el_px, el_py, el_pz, el_theta, el_En;
    float pion_px, pion_py, pion_pz, pion_eta, pion_Q2, pion_xB, pion_y, pion_z, pion_PhT, pion_Phi_h, pion_Mom, pion_theta;
    float pion_px_mc, pion_py_mc, pion_pz_mc, pion_eta_mc, pion_Q2_mc, pion_xB_mc, pion_y_mc, pion_z_mc, pion_PhT_mc, pion_Phi_h_mc, pion_Mom_mc, pion_theta_mc;
    float pion_xF, pion_Mx, pion_phi_lab, el_phi_lab;
    float pion_eta_res, pion_Q2_res, pion_xB_res, pion_y_res, pion_z_res, pion_PhT_res, pion_Phi_h_res, pion_Mom_res, pion_theta_res;
    float pion_eta_res_ratio, pion_Q2_res_ratio, pion_xB_res_ratio, pion_y_res_ratio, pion_z_res_ratio, pion_PhT_res_ratio, pion_Phi_h_res_ratio, pion_Mom_res_ratio;
    float pion_Q2_new, pion_Q2_new_res, PhT_over_zQ, pion_Phi_s, pion_Phi_s_noA, pion_dSivers;
    float polarization = 0.9;
    float target_pol = 1.0;
    float depolariz;
    TLorentzVector beam(0, 0, 10.6, 10.6);
    TVector3 s_axis;

    //std::vector<double> phiVector, polVector, epsilonVector;
    //std::vector<int> helicityVector;
    const int nBins = 9;
    const int nBins_3d = 5;
    const int nBins_z = 4;
    const int nBins_Pt = 4;
    std::vector<std::vector<double>> phiVector(nBins);
    std::vector<std::vector<double>> phi_s_Vector(nBins);
    std::vector<std::vector<double>> polVector(nBins);
    std::vector<std::vector<double>> TpolVector(nBins);
    std::vector<std::vector<double>> epsilonVector(nBins);
    std::vector<std::vector<float>> helicityVector(nBins);
    std::vector<std::vector<std::vector<double>>> phiVector_z(nBins, std::vector<std::vector<double>>(nBins_z));
    std::vector<std::vector<std::vector<double>>> polVector_z(nBins, std::vector<std::vector<double>>(nBins_z));
    std::vector<std::vector<std::vector<double>>> TpolVector_z(nBins, std::vector<std::vector<double>>(nBins_z));
    std::vector<std::vector<std::vector<double>>> epsilonVector_z(nBins, std::vector<std::vector<double>>(nBins_z));
    std::vector<std::vector<std::vector<float>>> helicityVector_z(nBins, std::vector<std::vector<float>>(nBins_z));   
    vector<float> A_LU_values_unbinned(nBins, 0);
    vector<float> A_LU_errors_unbinned(nBins, 0);
    vector<float> A_UT_values_unbinned(nBins, 0);
    vector<float> A_UT_errors_unbinned(nBins, 0);
    TVector3 electron_beam(0, 0, 10.6);
    vector<vector<float>> A_LU_values_3D_unbinned(nBins, vector<float> (nBins_z, 0));
    vector<vector<float>> A_LU_errors_3D_unbinned(nBins, vector<float> (nBins_z, 0));
    
    const int Q2_nBins = 4;
    const int xB_nBins = 17;
    const int Pt_nBins = 5;
    const int z_nBins = 30;

    double Clas_Q2Bins[][2] = {{1,3}, {3, 5}, {5, 7}, {7,10}};
    double Clas_xBins_1[][2] = {{0, 0.142}, {0.142, 0.172}, {0.172, 0.2}, {0.2, 0.232}, {0.232, 0.269}, {0.269, 0.315}, {0.315, 1}};
    double Clas_xBins_2[][2] = {{0, 0.269}, {0.269, 0.315}, {0.315, 0.384}, {0.384, 0.45}, {0.45, 1}};
    double Clas_xBins_3[][2] = {{0, 0.45}, {0.45, 0.55}, {0.55, 1}};
    double Clas_xBins_4[][2] = {{0, 0.6}, {0.6, 1}};
    double Clas_PtBins[][2] = {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 1.25}};
    double Clas_zBins_12[][2] = {{0, 0.15}, {0.15, 0.25}, {0.25, 0.35}, {0.35, 0.45}, {0.45, 0.55}, {0.55, 0.7}, {0.7, 0.9}};
    double Clas_zBins_34[][2] = {{0, 0.25}, {0.25, 0.35}, {0.35, 0.45}, {0.45, 0.55}, {0.55, 0.7}, {0.7, 0.9}};
    double Clas_zBins_5[][2] = {{0, 0.35}, {0.35, 0.45}, {0.45, 0.55}, {0.55, 0.9}};
    double Clas_zBins_12_edge[] = {0, 0.15, 0.25, 0.35, 0.45, 0.55, 0.7, 0.9};
    double Clas_PtBins_edge[] = {0, 0.2, 0.4, 0.6, 0.8, 1.25};
    // double Q2Bins[][2] = {{1,3}, {3, 5}, {5, 7}, {7,11}};
    //double xBins[][2] = {{0, 0.142}, {0.142, 0.172}, {0.172, 0.2}, {0.2, 0.232}, {0.232, 0.269}, {0.269, 0.315}, {0.315, 1}, {0, 0.269}, {0.269, 0.315}, {0.315, 0.384}, {0.384, 0.45}, {0.45, 1}, {0, 0.45}, {0.45, 0.55}, {0.55, 1}, {0, 0.6}, {0.6, 1}};
    double xBins[][2] = {{0.05, 0.12}, {0.12, 0.2}, {0.12, 0.2}, {0.2, 0.3}, {0.2, 0.3}, {0.3, 0.45}, {0.3, 0.45}, {0.45, 0.85}, {0.45, 0.85}};
    double xBin_3d[][2] = {{0.05, 0.12}, {0.12, 0.2}, {0.2, 0.3}, {0.3, 0.45}, {0.45, 0.85}};
    double Q2Bins[][2] = {{1.0, 2.5}, {1, 2}, {2, 4}, {1.2, 2.75}, {2.75, 5.75}, {1.4, 3.75}, {3.75, 7.5}, {1.75, 5.5}, {5.5, 9.0}};
    double zBins[][2] = {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.55}, {0.55, 0.8}};
    double PhTBins[][2] = {{0.0, 0.2}, {0.2, 0.4}, {0.4, 0.7}, {0.7, 1.2}};
    double HERMES_xBins[][2] = {{0.023, 0.072}, {0.072, 0.092}, {0.092, 0.138}, {0.138, 0.600}};
    double HERMES_zBins[][2] = {{0.2, 0.28}, {0.28, 0.37}, {0.37, 0.49}, {0.49, 0.7}};
    double HERMES_PtBins[][2] = {{0.0, 0.23}, {0.23, 0.36}, {0.36, 0.54}, {0.54, 2}};
    double COMPASS_Q2Bins[][2] = {{1,4}, {4, 6.25}, {6.25, 16}, {16, 81}};
    vector<HermesPoint> hermes_data = ReadHermesData("HERMES_AUT_SSA.txt");
    vector<CompassPoint> compass_data = ReadCompassData("COMPASS_AUT_SSA.txt");
    /*
    for (const auto& point : hermes_data) {
        std::cout << "Q2: " << point.Q2
                  << ", x: " << point.x
                  << ", z: " << point.z
                  << ", Pt: " << point.Pt
                  << ", Aut: " << point.Aut
                  << ", StatErr: " << point.stat_err << std::endl;
    }
    */
    

    // chain formation
    TChain pionTree("PionTree");
    AddFilesToChain(pionTree, "agen_2_output_dir");
    AddFilesToChain(pionTree, "agen_1_output_dir");
    AddFilesToChain(pionTree, "agen_3_output_dir");
    AddFilesToChain(pionTree, "agen_4_output_dir");
    AddFilesToChain(pionTree, "agen_5_output_dir");
    AddFilesToChain(pionTree, "agen_6_output_dir");
    AddFilesToChain(pionTree, "agen_7_output_dir");
    AddFilesToChain(pionTree, "agen_8_output_dir");
    AddFilesToChain(pionTree, "agen_9_output_dir");

    // creo un output root 
    const char* outputFile = "plot_rgh_asymmetries.root";
    //const char* outputFile = "plot_sivers_pip.root"; 
    TFile outFile(outputFile, "RECREATE");  // File di output ROOT
    TDirectory* dir_xQ2 = outFile.mkdir("Binning xB-Q2");
    TDirectory* dir_zPt = outFile.mkdir("Binning z-Pt");
    TDirectory* dir_xQ2_zPt = outFile.mkdir("Binning 4D xB-Q2 & z-Pt");
    TDirectory* dir_comparison = outFile.mkdir("Comparison Hermes-Compass-Clas");
    TDirectory* dir_Sivers = outFile.mkdir("Sivers uncertainty");
    TDirectory* dir_Collins = outFile.mkdir("Collins uncertainty");
    TDirectory* dir_Sivers_ext_xB = outFile.mkdir("Sivers asymmetry extraction vs xB");
    TDirectory* dir_Sivers_ext_z = outFile.mkdir("Sivers asymmetry extraction vs z");
    TDirectory* dir_Collins_ext_xB = outFile.mkdir("Collins asymmetry extraction vs xB");
    TDirectory* dir_Collins_ext_z = outFile.mkdir("Collins asymmetry extraction vs z");
    TDirectory* dir_Pretz_ext_xB = outFile.mkdir("Pretzelosity asymmetry extraction vs xB");
    TDirectory* dir_Pretz_ext_z = outFile.mkdir("Pretzelosity asymmetry extraction vs z");
    TDirectory* dir_Sivers_funct = outFile.mkdir("Sivers function");
    TDirectory* dir_Collins_funct = outFile.mkdir("Collins function");

    pionTree.SetBranchAddress("Mom", &pion_Mom);
    pionTree.SetBranchAddress("Q2", &pion_Q2);
    pionTree.SetBranchAddress("xB", &pion_xB);
    pionTree.SetBranchAddress("xF", &pion_xF);
    pionTree.SetBranchAddress("z", &pion_z);
    pionTree.SetBranchAddress("PhT", &pion_PhT);
    pionTree.SetBranchAddress("Phi_h", &pion_Phi_h);
    pionTree.SetBranchAddress("Phi_s", &pion_Phi_s_noA);
    pionTree.SetBranchAddress("theta", &pion_theta);
    pionTree.SetBranchAddress("eta", &pion_eta);
    pionTree.SetBranchAddress("y", &pion_y);
    pionTree.SetBranchAddress("Mx", &pion_Mx);
    //pionTree.SetBranchAddress("helicity", &helicity);
    pionTree.SetBranchAddress("epsilon", &epsilon);
    pionTree.SetBranchAddress("pion_px", &pion_px);
    pionTree.SetBranchAddress("pion_py", &pion_py);
    pionTree.SetBranchAddress("pion_pz", &pion_pz);
    pionTree.SetBranchAddress("el_px", &el_px);
    pionTree.SetBranchAddress("el_py", &el_py);
    pionTree.SetBranchAddress("el_pz", &el_pz);
    //pionTree.SetBranchAddress("sector_e", &sector_e);
    //pionTree.SetBranchAddress("sector_pi", &sector_pi);
    
    // maybe plot
    int nbin = 150;
    int nbin2 = 150;
    double xmin_xbj = 2e-2;
    double xmax_xbj = 1;
    double xmin_Q2 = 1;
    double xmax_Q2 = 12;
    std::vector<float> log_bins_Q2 = CreateLogBinning(nbin2, xmin_Q2, xmax_Q2);
    std::vector<float> log_bins_xbj = CreateLogBinning(nbin2, xmin_xbj, xmax_xbj);
    TH1F spaciong_4 ("------------------------------", "spacing", nbin, 0, 0);
    TH1F* Hist_pion_Mom = new TH1F("pion_Mom", "Momentum of the pion+ | z>0.2, M_{x}>1.6 GeV; Mom [GeV]; counts", nbin, 0, 7);
    TH1F* Hist_pion_Q2 = new TH1F("pion_Q2", "Q^{2} for the pion+ | z>0.2, M_{x}>1.6 GeV; Q^{2} [GeV^{2}]; counts", nbin, 1, 10);
    TH1F* Hist_pion_Q2_xCut = new TH1F("pion_Q2_xCut", "Q^{2} with x_{B}>0.3, z>0.2, M_{x}>1.6 GeV | #pi+; Q^{2} [GeV^{2}]; counts", nbin, 1, 10);
    TH1F* Hist_pion_xB = new TH1F("pion_xB", "x_{B} for the pion+ | z>0.2, M_{x}>1.6 GeV; x_{B}; counts", nbin, 0, 1);
    TH1F* Hist_pion_xB_Q2cut = new TH1F("pion_xB_Q2cut", "x_{B} for the pion+ | Q^{2}>2 GeV^{2}, z>0.2, M_{x}>1.6 GeV; x_{B}; counts", nbin, 0, 1);
    TH1F* Hist_pion_xF = new TH1F("pion_xF", "x_{F} of the pion+ | z>0.2, M_{x}>1.6 GeV; x_{F}; counts", nbin, 0, 1);
    TH1F* Hist_pion_z = new TH1F("pion_z", "z of the pion+ | z>0.2, M_{x}>1.6 GeV; z; counts", nbin, 0.2, 1);
    TH1F* Hist_pion_PhT = new TH1F("pion_PhT", "P_{hT} of the pion+ | z>0.2, M_{x}>1.6 GeV; P_{hT} [GeV]; counts", nbin, 0, 2);
    TH1F* Hist_pion_PhT2 = new TH1F("pion_PhT_2", "P_{hT}^{2} of the pion+ | z>0.2, M_{x}>1.6 GeV; P_{hT}^{2} [GeV^{2}]; counts", nbin, 0, 3);
    TH1F* Hist_pion_Phi_h = new TH1F("pion_Phi_h", "#Phi_{h} of the pion+ | z>0.2, M_{x}>1.6 GeV; #Phi_{h} [Rad]; counts", nbin, -TMath::Pi(), TMath::Pi());
    TH1F* Hist_pion_Phi_s = new TH1F("pion_Phi_s", "#Phi_{s} of the pion+ | z>0.2, M_{x}>1.6 GeV; #Phi_{s} [Rad]; counts", nbin, -TMath::Pi(), TMath::Pi());
    TH1F* Hist_pion_Phi_s_noA = new TH1F("pion_Phi_s_noA", "#Phi_{s} of the pion+ | no Asymmetry | z>0.2, M_{x}>1.6 GeV; #Phi_{s} [Rad]; counts", nbin, -TMath::Pi(), TMath::Pi());
    TH1F* Hist_pion_dSivers = new TH1F("pion_Phi_dSivers", "#Phi_{h} - #Phi_{s} of the pion+ | z>0.2, M_{x}>1.6 GeV; #Delta#Phi [Rad]; counts", nbin, -TMath::Pi(), TMath::Pi());
    TH1F* Hist_pion_Phi_h_Hp = new TH1F("pion_Phi_h_Hp", "#Phi_{h} of the pion+ with positive helicity | z>0.2, M_{x}>1.6 GeV; #Phi_{h} [Rad]; counts", nbin, -TMath::Pi(), TMath::Pi());
    TH1F* Hist_pion_Phi_h_Hm = new TH1F("pion_Phi_h_Hm", "#Phi_{h} of the pion+ with negative helicity | z>0.2, M_{x}>1.6 GeV; #Phi_{h} [Rad]; counts", nbin, -TMath::Pi(), TMath::Pi());
    TH1F* Hist_pion_theta = new TH1F("pion_theta", "#theta of the pion+ | z>0.2, M_{x}>1.6 GeV; #theta [Rad]; counts", nbin, 0, 0.5);
    TH1F* Hist_pion_eta = new TH1F("pion_eta", "Pseudorapidity of the pion+ | z>0.2, M_{x}>1.6 GeV; #eta; counts", nbin, 0.5, 3.5);
    TH1F* Hist_helicity = new TH1F("target_polarization", "Target polarization | pion+| z>0.2, M_{x}>1.6 GeV; s; counts", 10, -2, 2);
    TH1F* Hist_helicity_ratio = new TH1F("polarization_ratio", "#spin up / #spin tot | pion+| z>0.2, M_{x}>1.6 GeV; s_{up}%; counts", 100, 0, 1);
    TH1F* Hist_pion_y = new TH1F("pion_y", "y for pion+ | z>0.2, M_{x}>1.6 GeV; y; counts", nbin, 0, 1);
    TH1F* Hist_pion_Mx = new TH1F("pion_Mx", "M_{x} of the pion+ | z>0.2, M_{x}>1.6 GeV; M_{x} [GeV]; counts", nbin, 0, 4.5);
    TH1F* Hist_epsilon = new TH1F("epsilon", "epsilon pion+ | z>0.2, M_{x}>1.6 GeV; #epsilon; counts", nbin, 0, 1);
    TH1F* Hist_depolariz = new TH1F("depolarization", "depolarization pion+ | z>0.2, M_{x}>1.6 GeV; D_{NN}; counts", nbin, 0, 1);
    TH1F* Hist_pion_Pt_over_zQ = new TH1F("pion_Pt_over_zQ", "P_{hT} / zQ^{2} for pion+ | z>0.2, M_{x}>1.6 GeV; P_{hT}/zQ; counts", nbin, 0, 3);
    TH2F* Hist_Q2VsXb = new TH2F("_Q2vsXb", "Correlation Q^{2} vs x_{B} | z > 0.2 , M_{x} > 1.6 GeV | Pion; x_{B}; Q^{2} [GeV^{2}]", nbin2, 0, 1, nbin2, 1, 10);
    TH2F* Hist_zVsPt = new TH2F("_zvsPt", "Correlation z vs P_{hT} | z > 0.2 , M_{x} > 1.6 GeV | Pion; z; P_{hT} [GeV]", nbin2, 0, 1, nbin2, 0, 1.6);
    TH2F* Hist_Q2VsXb_log = new TH2F("_Q2vsXb_log", "Correlation Q^{2} vs x_{B} | z > 0.2 , M_{x} > 1.6 GeV | Pion; x_{B}; Q^{2} [GeV^{2}]", nbin2, log_bins_xbj.data(), nbin2, log_bins_Q2.data());
    /*
    TH1F* Hist_pion_Q2_Elsect3 = new TH1F("pion_Q2_el_sect3", "Q^{2} for the pion+ when e is in sect 3 | z>0.2, M_{x}>1.6 GeV; Q^{2} [GeV^{2}]; counts", nbin, 1, 10);
    TH1F* Hist_pion_Q2_Elsect4 = new TH1F("pion_Q2_el_sect4", "Q^{2} for the pion+ when e is in sect 4 | z>0.2, M_{x}>1.6 GeV; Q^{2} [GeV^{2}]; counts", nbin, 1, 10);
    TH1F* Hist_pion_Q2_Pisect4 = new TH1F("pion_Q2_Pi_sect4", "Q^{2} for the pion+ when Pi is in sect 4 | z>0.2, M_{x}>1.6 GeV; Q^{2} [GeV^{2}]; counts", nbin, 1, 10);
    TH1F* Hist_pion_Q2_Elsect5 = new TH1F("pion_Q2_el_sect5", "Q^{2} for the pion+ when e is in sect 5 | z>0.2, M_{x}>1.6 GeV; Q^{2} [GeV^{2}]; counts", nbin, 1, 10);
    TH1F* Hist_pion_xB_Elsect3 = new TH1F("pion_xB_el_sect3", "x_{B} for the pion+ when e is in sect 3 | z>0.2, M_{x}>1.6 GeV; x_{B}; counts", nbin, 0, 1);
    TH1F* Hist_pion_xB_Elsect4 = new TH1F("pion_xB_el_sect4", "x_{B} for the pion+ when e is in sect 4 | z>0.2, M_{x}>1.6 GeV; x_{B}; counts", nbin, 0, 1);
    TH1F* Hist_pion_xB_Pisect4 = new TH1F("pion_xB_Pi_sect4", "x_{B} for the pion+ when Pi is in sect 4 | z>0.2, M_{x}>1.6 GeV; x_{B}; counts", nbin, 0, 1);
    TH1F* Hist_pion_xB_Elsect5 = new TH1F("pion_xB_el_sect5", "x_{B} for the pion+ when e is in sect 5 | z>0.2, M_{x}>1.6 GeV; x_{B}; counts", nbin, 0, 1);
    */
    TH1F* Hist_pion_Phi_lab = new TH1F("pion_Phi_lab", "#Phi_{Lab} of the pion+ in the lab system| z>0.2, M_{x}>1.6 GeV; #Phi_{Lab} [Rad]; counts", nbin, -TMath::Pi(), TMath::Pi());
    //TH1F* Hist_pion_Phi_lab_No4 = new TH1F("pion_Phi_lab_No4", "#Phi of the pion+ in the lab system w/o sect 4| z>0.2, M_{x}>1.6 GeV; #Phi_{h} [Rad]; counts", nbin, -TMath::Pi(), TMath::Pi());
    TH1F* Hist_el_Phi_lab = new TH1F("el_Phi_lab", "#Phi_{Lab} of the electron in the lab system| z>0.2, M_{x}>1.6 GeV; #Phi_{Lab} [Rad]; counts", nbin, -TMath::Pi(), TMath::Pi());
    //TH1F* Hist_el_Phi_lab_No4 = new TH1F("el_Phi_lab_No4", "#Phi of the electron in the lab system w/o sect 4| z>0.2, M_{x}>1.6 GeV; #Phi_{h} [Rad]; counts", nbin, -TMath::Pi(), TMath::Pi());
    TH2F* Hist_pion_PhiTheta_lab = new TH2F("pion_Phi_vs_Theta_lab", "#Phi_{Lab} vs #theta of the pion+ in the lab system| z>0.2, M_{x}>1.6 GeV; #Phi_{Lab} [Rad]; #Theta [Rad]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0, 0.5);
    TH2F* Hist_pion_PhiTheta = new TH2F("pion_Phi_vs_Theta", "#Phi_{h} vs #theta of the pion+ | z>0.2, M_{x}>1.6 GeV; #Phi_{h} [Rad]; #Theta [Rad]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0, 0.5);
    //TH2F* Hist_pion_PhiTheta_lab_No4 = new TH2F("pion_Phi_vs_Theta_lab_No4", "#Phi vs #theta of the pion+ in the lab system| z>0.2, M_{x}>1.6 GeV; #Phi_{h} [Rad]; #Theta [Rad]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0, 0.75);
    TH2F* Hist_el_PhiTheta_lab = new TH2F("el_Phi_vs_Theta_lab", "#Phi_{Lab} vs #theta of the electron in the lab system| z>0.2, M_{x}>1.6 GeV; #Phi_{Lab} [Rad]; #Theta [Rad]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0, 0.7);
    //TH2F* Hist_el_PhiTheta_lab_No4 = new TH2F("el_Phi_vs_Theta_lab_No4", "#Phi vs #theta of the electron in the lab system| z>0.2, M_{x}>1.6 GeV; #Phi_{h} [Rad]; #Theta [Rad]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0, 0.75);

    TH1F spaciong_1 ("----------------------------", "spacing", nbin, 0, 0);
    //
    std::vector<TH2F*> hist_xBvsQ2_binned(xB_nBins);
    std::vector<TH2F*> hist_zvsPt_binned(z_nBins);
    std::vector<TH1F*> hist_Phi_Hp_binned(xB_nBins);
    std::vector<TH1F*> hist_Phi_Hm_binned(xB_nBins);
    std::vector<TH1F*> hist_z_Phi_Hp_binned(z_nBins);
    std::vector<TH1F*> hist_z_Phi_Hm_binned(z_nBins);
    std::vector<TH1F*> hist_Pt_Phi_Hp_binned(z_nBins);
    std::vector<TH1F*> hist_Pt_Phi_Hm_binned(z_nBins);
    // one dimensional plot
    std::vector<TH1F*> hist_x_dist_2D_bin(xB_nBins);
    std::vector<TH1F*> hist_Q2_dist_2D_bin(xB_nBins);
    std::vector<TH1F*> hist_z_dist_2D_bin(z_nBins);
    std::vector<TH1F*> hist_PhT_dist_2D_bin(z_nBins);
    float z_mean_1D[z_nBins], z_meanError_1D[z_nBins];
    float x_mean_2D[xB_nBins], x_meanError_2D[xB_nBins];
    float Q2_mean_2D[xB_nBins], Q2_meanError_2D[xB_nBins];
    float PhT_mean_1D[z_nBins], PhT_meanError_1D[z_nBins];
    // multi-dimensional plot
    //vector<vector<TH1F*>> hist_AUT_Hp_multidim_xQPz(nBins, vector<TH1F*> (z_nBins));
    //vector<vector<TH1F*>> hist_AUT_Hm_multidim_xQPz(nBins, vector<TH1F*> (z_nBins));
    vector<vector<TH1F*>> hist_dist_PhT_4d(xB_nBins, vector<TH1F*> (z_nBins));
    vector<vector<TH1F*>> hist_dist_z_4d(xB_nBins, vector<TH1F*> (z_nBins));
    vector<vector<TH1F*>> hist_dist_Q2_4d(xB_nBins, vector<TH1F*> (z_nBins));
    vector<vector<TH1F*>> hist_dist_xB_4d(xB_nBins, vector<TH1F*> (z_nBins));
    //
    vector<TH2F*> z_vs_Pt_xQ_bin(xB_nBins);
    //
    // VECTOR FOR THE MLE
    vector<double> Param_A_sivers;
    vector<double> Param_A_collins;
    // xB-Q2 binning
    vector<vector<double>> vec_pion_phi_h(xB_nBins);
    vector<vector<double>> vec_pion_phi_s(xB_nBins);
    vector<vector<double>> vec_pion_depol(xB_nBins);
    vector<vector<double>> vec_pion_epsilon(xB_nBins);
    vector<vector<double>> vec_pion_Pt(xB_nBins);
    vector<vector<double>> vec_pion_xB(xB_nBins);
    vector<vector<double>> vec_pion_helicity(xB_nBins);
    vector<double> vec_pion_xB_mean;
    vector<double> vec_pion_Pt_mean;
    vector<double> A_UT_sivers(xB_nBins);
    vector<double> A_UT_sivers_err(xB_nBins);
    vector<double> A_UT_collins(xB_nBins);
    vector<double> A_UT_collins_err(xB_nBins);
    vector<double> A_UT_pretz(xB_nBins);
    vector<double> A_UT_pretz_err(xB_nBins);
    // z-Pt binning
    vector<vector<double>> vec_pion_helicity_z(z_nBins);
    vector<vector<double>> vec_pion_phi_h_z(z_nBins);
    vector<vector<double>> vec_pion_phi_s_z(z_nBins);
    vector<vector<double>> vec_pion_depol_z(z_nBins);
    vector<vector<double>> vec_pion_epsilon_z(z_nBins);
    vector<vector<double>> vec_pion_z(z_nBins);
    vector<double> A_UT_sivers_z(z_nBins);
    vector<double> A_UT_sivers_err_z(z_nBins);
    vector<double> A_UT_collins_z(z_nBins);
    vector<double> A_UT_collins_err_z(z_nBins);
    vector<double> A_UT_pretz_z(z_nBins);
    vector<double> A_UT_pretz_err_z(z_nBins);
    // xB-Q2 + z-Pt bin
    vector<vector<vector<double>>> vec_pion_helicity_4d(xB_nBins, vector<vector<double>> (z_nBins));
    vector<vector<vector<double>>> vec_pion_phi_h_4d(xB_nBins, vector<vector<double>> (z_nBins));
    vector<vector<vector<double>>> vec_pion_phi_s_4d(xB_nBins, vector<vector<double>> (z_nBins));
    vector<vector<vector<double>>> vec_pion_depol_4d(xB_nBins, vector<vector<double>> (z_nBins));
    vector<vector<vector<double>>> vec_pion_epsilon_4d(xB_nBins, vector<vector<double>> (z_nBins));
    vector<vector<vector<double>>> vec_pion_xB_4d(xB_nBins, vector<vector<double>> (z_nBins));
    vector<vector<vector<double>>> vec_pion_y_4d(xB_nBins, vector<vector<double>> (z_nBins));
    vector<vector<vector<double>>> vec_pion_Pt_4d(xB_nBins, vector<vector<double>> (z_nBins));
    vector<vector<vector<double>>> vec_pion_z_4d(xB_nBins, vector<vector<double>> (z_nBins));
    vector<vector<double>> A_UT_sivers_4d(xB_nBins, vector<double> (z_nBins));
    vector<vector<double>> A_UT_sivers_err_4d(xB_nBins, vector<double> (z_nBins));
    vector<vector<double>> A_UT_collins_4d(xB_nBins, vector<double> (z_nBins));
    vector<vector<double>> A_UT_collins_err_4d(xB_nBins, vector<double> (z_nBins));
    vector<vector<double>> A_UT_pretz_4d(xB_nBins, vector<double> (z_nBins));
    vector<vector<double>> A_UT_pretz_err_4d(xB_nBins, vector<double> (z_nBins));
    vector<double> vec_pion_xB_4d_mean;
    vector<double> vec_pion_z_4d_mean;
    vector<double> vec_pion_Pt_4d_mean;
    vector<double> vec_pion_y_4d_mean;
    vector<double> vec_pion_xB_4d_mean_coll;
    vector<double> vec_pion_z_4d_mean_coll;
    vector<double> vec_pion_Pt_4d_mean_coll;
    vector<double> vec_pion_y_4d_mean_coll;
    vector<double> A_UT_sivers_single_4d;
    vector<double> A_UT_sivers_single_err_4d;
    vector<double> A_UT_collins_single_4d;
    vector<double> A_UT_collins_single_err_4d;

    vector<vector<vector<double>>> A_UT_mean_Sivers_4D(xB_nBins, vector<vector<double>> (z_nBins));
    vector<vector<double>> A_UT_mean_Sivers_2D(xB_nBins);
    vector<double> Siv_mean_2D;
    vector<double> Siv_mean_4D;
    vector<vector<vector<double>>> A_UT_mean_Collins_4D(xB_nBins, vector<vector<double>> (z_nBins));
    vector<vector<double>> A_UT_mean_Collins_2D(xB_nBins);
    vector<double> Col_mean_2D;
    vector<double> Col_mean_4D;

    // plots
    for (int i = 0; i < xB_nBins; i++) {
        // modo dim
        hist_x_dist_2D_bin[i] = new TH1F(Form("hist_x_dist_2D_bin%d", i+1),
            Form("x_{B} distribution for bin (%d, x_{B}-Q^{2}) | pion+; x_{B}", i+1), 120, 0, 1); 
        hist_Q2_dist_2D_bin[i] = new TH1F(Form("hist_Q2_dist_2D_bin%d", i+1),
            Form("Q^{2} distribution for bin (%d, x_{B}-Q^{2}) | pion+; Q^{2} [GeV^{2}]", i+1), 120, 1, 10); 
        // FOR AUT UNCERTAINTY
        z_vs_Pt_xQ_bin[i] = new TH2F(Form("hist_PhTvsZ_xQ2_bin%d", i+1),
            Form("Asymmetry statistical uncertainty as P_{hT} vs z for bin (%d, x_{B}-Q^{2}) | pion+ | <0.04; z; P_{hT} [GeV]", i+1), 7, Clas_zBins_12_edge, 5, Clas_PtBins_edge);
        // more dim
        hist_xBvsQ2_binned[i] = new TH2F(Form("hist_xBvsQ2_bin%d", i+1),
                                Form("x_{B} vs Q^{2} for bin (%d, x_{B}-Q^{2}) | pion+; x_{B}; Q^{2} [GeV^{2}]", i+1), 120, 0, 1, 120, 1, 10); 
        hist_Phi_Hp_binned[i] = new TH1F(Form("hist_Phi_Hp_bin%d", i+1),
                                Form("#Phi_{h} with positive helicity for bin (%d, x_{B}-Q^{2}) | pion+; #Phi_{h}; counts", i+1), 10, -TMath::Pi(), TMath::Pi()); 
        hist_Phi_Hm_binned[i] = new TH1F(Form("hist_Phi_Hm_bin%d", i+1),
                                Form("#Phi_{h} with negative helicity for bin (%d, x_{B}-Q^{2}) | pion+; #Phi_{h}; counts", i+1), 10, -TMath::Pi(), TMath::Pi()); 
        // siamo dentro xQ2 apriamo prima Pt e poi z
        for (int z = 0; z < z_nBins; z++){
            hist_dist_xB_4d[i][z] = new TH1F(Form("hist_xB_bin_xQ2%d_Pt_z%d", i+1, z+1),
            Form("dist_xB_4d bin (%d, x_{B}-Q^{2}) and (%d, z-P_{hT}) | pion+; x_{B}; counts", i+1, z+1), 60, 0, 1); 
            hist_dist_Q2_4d[i][z] = new TH1F(Form("hist_Q2_bin_xQ2%d_Pt_z%d", i+1, z+1),
            Form("dist_Q2_4d bin (%d, x_{B}-Q^{2}) and (%d, z-P_{hT}) | pion+; Q^{2}; counts", i+1, z+1), 60, 0, 9);
            hist_dist_PhT_4d[i][z] = new TH1F(Form("hist_PhT_bin_xQ2%d_Pt_z%d", i+1, z+1),
            Form("dist_PhT_4d bin (%d, x_{B}-Q^{2}) and (%d, z-P_{hT}) | pion+; P_{hT}; counts", i+1, z+1), 60, 0, 1.6);
            hist_dist_z_4d[i][z] = new TH1F(Form("hist_z_bin_xQ2%d_Pt_z%d", i+1, z+1),
            Form("dist_z_4d bin (%d, x_{B}-Q^{2}) and (%d, z-P_{hT}) | pion+; z; counts", i+1, z+1), 60, 0, 1);
        }
    }
    for(int j=0; j<z_nBins; j++){
        hist_z_dist_2D_bin[j] = new TH1F(Form("hist_z_dist_2D_bin%d", j+1),
            Form("z distribution for bin (%d, z-P_{hT}) | pion+; z", j+1), 120, 0, 1); 
        hist_PhT_dist_2D_bin[j] = new TH1F(Form("hist_PhT_dist_2D_bin%d", j+1),
            Form("P_{hT} distribution for bin (%d, z-P_{hT}) | pion+; P_{hT} [GeV]", j+1), 120, 0, 1.6); 
        hist_zvsPt_binned[j] = new TH2F(Form("hist_zvsPt_bin%d", j+1),
            Form("z vs P_{hT} for bin (%d, z-P_{hT}) | pion+; z; P_{hT} [GeV]", j+1), 120, 0, 1, 120, 0, 1.6); 
        hist_z_Phi_Hp_binned[j] = new TH1F(Form("hist_z_Phi_Hp_bin%d", j+1),
            Form("#Phi_{h} with positive helicity for bin (%d, z-P_{hT}) | pion+; #Phi_{h}; counts", j+1), 10, -TMath::Pi(), TMath::Pi()); 
        hist_z_Phi_Hm_binned[j] = new TH1F(Form("hist_z_Phi_Hm_bin%d", j+1),
            Form("#Phi_{h} with negative helicity for bin %d z and P_{hT} | pion+; #Phi_{h}; counts", j+1), 10, -TMath::Pi(), TMath::Pi()); 
    }

    double h_up = 0, h_down = 0;
    // analysis
    Long64_t nEntries = pionTree.GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) { 
        pionTree.GetEntry(i);
        //if (sector_e == 4 || sector_pi == 4) continue;
        //if (i == 851392) continue; // Phi_h is Nan here, idk why
        //if (!std::isnan(pion_Phi_h) && !std::isnan(pion_Phi_s)) pion_dSivers = TVector2::Phi_mpi_pi(pion_Phi_h - pion_Phi_s);
        //else std::cerr << "Warning! NaN detected → pion_Phi_h: " << pion_Phi_h << ", pion_Phi_s: " << pion_Phi_s << "event number: " << i << std::endl;
        TVector3 mom_pion(pion_px, pion_py, pion_pz);
        TVector3 mom_el(el_px, el_py, el_pz);
        //
        //helicity = random_helicity(pion_Phi_h);
        //depolariz = (1-epsilon) / (1 + epsilon);
        depolariz = (1-pion_y)/(1-pion_y+(pion_y*pion_y*0.5));
        // parameters centered in their mean value to have more stability
        double Param_Siv = Parametrization_Aut_Sivers_pip(pion_xB - 0.2406, pion_Q2 - 2.709, pion_z - 0.3805, pion_PhT - 0.3736);
        double Param_Col = Parametrization_Aut_Collins_pip(pion_xB - 0.2406, pion_Q2 - 2.709, pion_z - 0.3805, pion_PhT - 0.3736);
        double Param_Pretz = 0.01;
        //double Param_Siv = Parametrization_Aut_Sivers_pip(pion_xB, pion_Q2, pion_z, pion_PhT);
        //double Param_Col = Parametrization_Aut_Collins_pip(pion_xB, pion_Q2, pion_z, pion_PhT);
        double pion_phi_sivers = TVector2::Phi_mpi_pi(pion_Phi_h - pion_Phi_s_noA);
        double pion_phi_collins = TVector2::Phi_mpi_pi(pion_Phi_h + pion_Phi_s_noA);
        double pion_phi_pretz = TVector2::Phi_mpi_pi(3*pion_Phi_h - pion_Phi_s_noA);
        double pol_funct = PolarizFunction_pip(pion_phi_collins, pion_phi_sivers, Param_Col, Param_Siv, depolariz, Param_Pretz, pion_phi_pretz, epsilon);
        Param_A_collins.push_back(Param_Col);
        Param_A_sivers.push_back(Param_Siv);
        double extracted_siv = depolariz*(0.05 + Param_Siv*sin(pion_phi_sivers));
        double extracted_col = depolariz*(0.05 + Param_Col*sin(pion_phi_collins));
        double rho = randGen.Uniform(0, 1); 
        pion_Phi_s = pion_Phi_s_noA;
        if(rho <= pol_funct){
            s_axis.SetXYZ(0,1,0);
            helicity = 1;
            h_up ++;
        } else if(rho > pol_funct){ 
            s_axis.SetXYZ(0,-1,0);
            helicity = -1;
            h_down ++;
        }
        //
        // new Phi_s measurement
        TVector3 el(el_px, el_py, el_pz);
        double philab_el = el.Phi();
        double el_mom = el.Mag();
        el_En = sqrt(el_mom*el_mom + 5e-3*5e-3);
        TLorentzVector electron(el, el_En);
        TLorentzVector q = beam - electron;
        TVector3 q_vect = q.Vect();
        TVector3 z_axis = q_vect.Unit();
        TVector3 y_axis = beam.Vect().Cross(electron.Vect()).Unit();
        TVector3 x_axis = y_axis.Cross(q_vect).Unit();
        TVector3 t1 = x_axis.Cross(s_axis);
        double termine1 = t1 * z_axis;
        double termine2 = x_axis.Dot(s_axis);
        //pion_Phi_s = std::atan2(termine1, termine2);
        pion_dSivers = TVector2::Phi_mpi_pi(pion_Phi_h - pion_Phi_s);
        //
        //
        pion_phi_lab = mom_pion.Phi(), el_phi_lab = mom_el.Phi(), el_theta = mom_el.Theta();
        PhT_over_zQ = pion_PhT/(pion_z*sqrt(pion_Q2));
        if(abs(pion_phi_lab) > 2.6) continue;
        if(abs(philab_el) > 2.6) continue;
        Hist_pion_eta->Fill(pion_eta);
        Hist_pion_Q2->Fill(pion_Q2);
        Hist_pion_PhT->Fill(pion_PhT);
        Hist_pion_PhT2->Fill(pion_PhT*pion_PhT);
        Hist_pion_xB->Fill(pion_xB);
        Hist_pion_Phi_h->Fill(pion_Phi_h), Hist_pion_Phi_h->SetMinimum(0);
        Hist_pion_Phi_s->Fill(pion_Phi_s);
        Hist_pion_Phi_s_noA->Fill(pion_Phi_s_noA);
        Hist_pion_dSivers->Fill(pion_dSivers), Hist_pion_dSivers->SetMinimum(0);
        Hist_pion_z->Fill(pion_z);
        Hist_pion_Mom->Fill(pion_Mom);
        Hist_helicity->Fill(helicity);
        Hist_depolariz->Fill(depolariz);
        Hist_pion_theta->Fill(pion_theta);
        Hist_pion_y->Fill(pion_y);
        Hist_epsilon->Fill(epsilon);
        Hist_pion_Mx->Fill(pion_Mx);
        Hist_pion_xF->Fill(pion_xF);
        Hist_Q2VsXb->Fill(pion_xB, pion_Q2);
        Hist_Q2VsXb_log->Fill(pion_xB, pion_Q2);
        Hist_zVsPt->Fill(pion_z, pion_PhT);
        Hist_pion_Pt_over_zQ->Fill(PhT_over_zQ);
        Hist_pion_Phi_lab->Fill(pion_phi_lab);
        Hist_el_Phi_lab->Fill(el_phi_lab);
        Hist_pion_PhiTheta_lab->Fill(pion_phi_lab, pion_theta);
        Hist_pion_PhiTheta->Fill(pion_Phi_h, pion_theta);
        Hist_el_PhiTheta_lab->Fill(el_phi_lab, el_theta);
        if(helicity == 1) Hist_pion_Phi_h_Hp->Fill(pion_Phi_h);
        else if(helicity == -1) Hist_pion_Phi_h_Hm->Fill(pion_Phi_h);
        /*
        if (sector_pi != 4){
            Hist_pion_Phi_lab_No4->Fill(pion_phi_lab);
            Hist_pion_PhiTheta_lab_No4->Fill(pion_phi_lab, pion_theta);
        } else if (sector_pi == 4){
            Hist_pion_Q2_Pisect4->Fill(pion_Q2);
            Hist_pion_xB_Pisect4->Fill(pion_xB);
        }
        if (sector_e != 4){
            Hist_el_Phi_lab_No4->Fill(el_phi_lab);
            Hist_el_PhiTheta_lab_No4->Fill(el_phi_lab, el_theta);
        }
        if (sector_e == 3){
            Hist_pion_Q2_Elsect3->Fill(pion_Q2);
            Hist_pion_xB_Elsect3->Fill(pion_xB);
        } else if (sector_e == 4){
            Hist_pion_Q2_Elsect4->Fill(pion_Q2);
            Hist_pion_xB_Elsect4->Fill(pion_xB);
        }else if (sector_e == 5){
            Hist_pion_Q2_Elsect5->Fill(pion_Q2);
            Hist_pion_xB_Elsect5->Fill(pion_xB);
        }
        */
        if(pion_Q2 > 2) Hist_pion_xB_Q2cut->Fill(pion_xB);
        if(pion_xB > 0.3) Hist_pion_Q2_xCut->Fill(pion_Q2);
        //int binIndex = getBinIndex(pion_xB, pion_Q2); // Funzione che determina il bin;
        int bin_xQ2 = getBinIndex_xQ2(pion_xB, pion_Q2);
        int bin_zPt = getBinIndex_zPt(pion_z, pion_PhT);
        if(bin_xQ2 >= 0){ // if -1 error
            hist_x_dist_2D_bin[bin_xQ2]->Fill(pion_xB);
            hist_Q2_dist_2D_bin[bin_xQ2]->Fill(pion_Q2);
            hist_xBvsQ2_binned[bin_xQ2]->Fill(pion_xB, pion_Q2);
            // Fill the MLE vector
            vec_pion_phi_h[bin_xQ2].push_back(pion_Phi_h);
            vec_pion_phi_s[bin_xQ2].push_back(pion_Phi_s);
            vec_pion_depol[bin_xQ2].push_back(depolariz);
            vec_pion_epsilon[bin_xQ2].push_back(epsilon);
            vec_pion_xB[bin_xQ2].push_back(pion_xB);
            vec_pion_Pt[bin_xQ2].push_back(pion_PhT);
            A_UT_mean_Sivers_2D[bin_xQ2].push_back(extracted_siv);
            A_UT_mean_Collins_2D[bin_xQ2].push_back(extracted_col);
            vec_pion_helicity[bin_xQ2].push_back(helicity);

            if(bin_zPt >= 0){
                hist_dist_xB_4d[bin_xQ2][bin_zPt]->Fill(pion_xB);
                hist_dist_Q2_4d[bin_xQ2][bin_zPt]->Fill(pion_Q2);
                hist_dist_z_4d[bin_xQ2][bin_zPt]->Fill(pion_z);
                hist_dist_PhT_4d[bin_xQ2][bin_zPt]->Fill(pion_PhT);
                vec_pion_helicity_4d[bin_xQ2][bin_zPt].push_back(helicity);
                // Fill the MLE vector 4D
                vec_pion_phi_h_4d[bin_xQ2][bin_zPt].push_back(pion_Phi_h);
                vec_pion_phi_s_4d[bin_xQ2][bin_zPt].push_back(pion_Phi_s);
                vec_pion_depol_4d[bin_xQ2][bin_zPt].push_back(depolariz);
                vec_pion_epsilon_4d[bin_xQ2][bin_zPt].push_back(epsilon);
                vec_pion_xB_4d[bin_xQ2][bin_zPt].push_back(pion_xB);
                vec_pion_y_4d[bin_xQ2][bin_zPt].push_back(pion_y);
                vec_pion_Pt_4d[bin_xQ2][bin_zPt].push_back(pion_PhT);
                vec_pion_z_4d[bin_xQ2][bin_zPt].push_back(pion_z);
                A_UT_mean_Sivers_4D[bin_xQ2][bin_zPt].push_back(extracted_siv);
                A_UT_mean_Collins_4D[bin_xQ2][bin_zPt].push_back(extracted_col);
            }
        }
        if(bin_zPt >= 0){ // empty bin due to the z>0.2 cut -> 0 and 7
            hist_z_dist_2D_bin[bin_zPt]->Fill(pion_z);
            hist_PhT_dist_2D_bin[bin_zPt]->Fill(pion_PhT);
            hist_zvsPt_binned[bin_zPt]->Fill(pion_z, pion_PhT);
            vec_pion_helicity_z[bin_zPt].push_back(helicity);
            // Fill the MLE vector
            vec_pion_phi_h_z[bin_zPt].push_back(pion_Phi_h);
            vec_pion_phi_s_z[bin_zPt].push_back(pion_Phi_s);
            vec_pion_depol_z[bin_zPt].push_back(depolariz);
            vec_pion_epsilon_z[bin_zPt].push_back(epsilon);
            vec_pion_z[bin_zPt].push_back(pion_z);
        }
        
    }

    double h_ratio = h_up / (h_up+h_down);
    Hist_helicity_ratio->Fill(h_ratio);

    double sum_col = 0, sum_siv = 0;
    for (double val : Param_A_collins) sum_col += val;
    double mean_A_col = sum_col / Param_A_collins.size();
    for (double val : Param_A_sivers) sum_siv += val;
    double mean_A_siv = sum_siv / Param_A_sivers.size();
    cout << "\n======================== A_UT ========================" << endl;
    cout << "" << endl;
    cout << "<A_UT_Collins> = " << mean_A_col << "   <A_UT_Sivers> = " << mean_A_siv << endl;
    cout << "\n======================================================" << endl;
    cout << "" << endl;

    for(int i = 0; i < xB_nBins; i++){
        const auto& vettori = A_UT_mean_Sivers_2D[i];
        double sum = 0;
        for (double val : vettori) sum += val;
        double mean_siv = sum / A_UT_mean_Sivers_2D[i].size();

        const auto& vettors = A_UT_mean_Collins_2D[i];
        double sumc = 0;
        for (double valc : vettors) sumc += valc;
        double mean_col = sumc / A_UT_mean_Collins_2D[i].size();

        Siv_mean_2D.push_back(mean_siv);
        Col_mean_2D.push_back(mean_col);
    }

    TGraphErrors* sivv = new TGraphErrors();
    TGraphErrors* coll = new TGraphErrors();
    for(int i = 0; i < Siv_mean_2D.size(); i++){
        double x = hist_x_dist_2D_bin[i]->GetMean();
        sivv->SetPoint(i, x, Siv_mean_2D[i]);
        sivv->SetPointError(i, 0, 0);
    }
    for(int i = 0; i < Col_mean_2D.size(); i++){
        double x = hist_x_dist_2D_bin[i]->GetMean();
        coll->SetPoint(i, x, Col_mean_2D[i]);
        coll->SetPointError(i, 0, 0);
    }

    TCanvas* c_siv = new TCanvas("c_siv", "");
    c_siv->cd();
    sivv->SetTitle("Siv implemented; x_{B}; A_{UT}");
    sivv->SetMarkerStyle(20);
    sivv->SetMarkerColor(kRed);
    sivv->SetLineColor(kRed);
    sivv->Draw("AP");
    c_siv->Update();
    c_siv->Write();

    TCanvas* c_col = new TCanvas("c_col", "");
    c_col->cd();
    coll->SetTitle("Col implemented; x_{B}; A_{UT}");
    coll->SetMarkerStyle(20);
    coll->SetMarkerColor(kRed);
    coll->SetLineColor(kRed);
    coll->Draw("AP");
    c_col->Update();
    c_col->Write();

    vector<double> depo;
    for(int x = 0; x < xB_nBins; x++){
        for (int z = 0; z < z_nBins; z++){
            const auto& vect_s = A_UT_mean_Sivers_4D[x][z];
            double sums = 0;
            for (double vals : vect_s) sums += vals;
            double mean_s = sums / A_UT_mean_Sivers_4D[x][z].size();
            Siv_mean_4D.push_back(mean_s);

            const auto& vect_c = A_UT_mean_Collins_4D[x][z];
            double sumc = 0;
            for (double valc : vect_c) sumc += valc;
            double mean_c = sumc / A_UT_mean_Collins_4D[x][z].size();
            Col_mean_4D.push_back(mean_c);

            const auto& vect_depo = vec_pion_depol_4d[x][z];
            double sumd = 0;
            for (double vald : vect_depo) sumd += vald;
            double meand = sumd / vec_pion_depol_4d[x][z].size();

            depo.push_back(meand);

        }
    }

        
    spaciong_4.Write();
    Hist_Q2VsXb->Write(), Hist_Q2VsXb_log->Write(), Hist_zVsPt->Write();
    Hist_pion_eta->Write(), Hist_pion_Q2->Write(), Hist_pion_Q2_xCut->Write(), Hist_pion_PhT->Write(), Hist_pion_PhT2->Write(), Hist_pion_xB->Write();
    Hist_pion_xB_Q2cut->Write(), Hist_pion_Phi_h->Write(), Hist_pion_Phi_s->Write(), Hist_pion_Phi_s_noA->Write(), Hist_pion_dSivers->Write(), Hist_pion_z->Write();
    Hist_pion_Mom->Write(), Hist_helicity->Write(), Hist_helicity_ratio->Write(), Hist_pion_theta->Write(), Hist_pion_y->Write();
    Hist_epsilon->Write(), Hist_depolariz->Write();
    Hist_pion_Mx->Write(), Hist_pion_xF->Write(), Hist_pion_Pt_over_zQ->Write(); 
    Hist_pion_Phi_lab->Write(), Hist_el_Phi_lab->Write(), Hist_pion_PhiTheta_lab->Write(), Hist_pion_PhiTheta->Write(), Hist_el_PhiTheta_lab->Write();

    
    // XQ2 DIRECTORY
    dir_xQ2->cd();
    TCanvas *c11 = new TCanvas("Q2_vs_xB_Bin", "Q^{2} vs x_{B} bin", 800, 700);
    c11->SetLogz();
    Hist_Q2VsXb->Draw("COLZ");
    std::vector<TPolyLine*> rectangles2;
    std::vector<TPolyLine*> rectangles3;
    std::vector<TPolyLine*> rectangles4;
    std::vector<TPolyLine*> rectangles5;
    for (int i = 0; i < 7; ++i) {
        double xB[5] = {Clas_xBins_1[i][0], Clas_xBins_1[i][1], Clas_xBins_1[i][1], Clas_xBins_1[i][0], Clas_xBins_1[i][0]};
        double Q2[5] = {Clas_Q2Bins[0][0], Clas_Q2Bins[0][0], Clas_Q2Bins[0][1], Clas_Q2Bins[0][1], Clas_Q2Bins[0][0]};
        TPolyLine *rect = new TPolyLine(5, xB, Q2);
        //rect->SetLineWidth(2);
        rect->Draw("same");
        rectangles2.push_back(rect);
    }
    std::vector<TText*> labels2;
    for (int i = 0; i < 7; ++i) {
        double x_center = 0.5 * (Clas_xBins_1[i][0] + Clas_xBins_1[i][1]);
        double Q2_center = 0.5 * (Clas_Q2Bins[0][0] + Clas_Q2Bins[0][1]);
        TText *label = new TText(x_center, Q2_center, Form("%d", i+1)); // etichetta come numero del bin
        label->SetTextAlign(22); // centrato
        label->SetTextSize(0.025); // regola la dimensione in base al tuo plot
        label->Draw("same");
        labels2.push_back(label);
    }
    for (int i = 0; i < 5; ++i) {
        double xB[5] = {Clas_xBins_2[i][0], Clas_xBins_2[i][1], Clas_xBins_2[i][1], Clas_xBins_2[i][0], Clas_xBins_2[i][0]};
        double Q2[5] = {Clas_Q2Bins[1][0], Clas_Q2Bins[1][0], Clas_Q2Bins[1][1], Clas_Q2Bins[1][1], Clas_Q2Bins[1][0]};
        TPolyLine *rect = new TPolyLine(5, xB, Q2);
        //rect->SetLineWidth(2);
        rect->Draw("same");
        rectangles3.push_back(rect);
    }
    std::vector<TText*> labels3;
    for (int i = 0; i < 5; ++i) {
        double x_center = 0.5 * (Clas_xBins_2[i][0] + Clas_xBins_2[i][1]);
        double Q2_center = 0.5 * (Clas_Q2Bins[1][0] + Clas_Q2Bins[1][1]);
        TText *label = new TText(x_center, Q2_center, Form("%d", i+8)); // etichetta come numero del bin
        label->SetTextAlign(22); // centrato
        label->SetTextSize(0.025); // regola la dimensione in base al tuo plot
        label->Draw("same");
        labels3.push_back(label);
    }
    for (int i = 0; i < 3; ++i) {
        double xB[5] = {Clas_xBins_3[i][0], Clas_xBins_3[i][1], Clas_xBins_3[i][1], Clas_xBins_3[i][0], Clas_xBins_3[i][0]};
        double Q2[5] = {Clas_Q2Bins[2][0], Clas_Q2Bins[2][0], Clas_Q2Bins[2][1], Clas_Q2Bins[2][1], Clas_Q2Bins[2][0]};
        TPolyLine *rect = new TPolyLine(5, xB, Q2);
        //rect->SetLineWidth(2);
        rect->Draw("same");
        rectangles4.push_back(rect);
    }
    std::vector<TText*> labels4;
    for (int i = 0; i < 3; ++i) {
        double x_center = 0.5 * (Clas_xBins_3[i][0] + Clas_xBins_3[i][1]);
        double Q2_center = 0.5 * (Clas_Q2Bins[2][0] + Clas_Q2Bins[2][1]);
        TText *label = new TText(x_center, Q2_center, Form("%d", i+13)); // etichetta come numero del bin
        label->SetTextAlign(22); // centrato
        label->SetTextSize(0.025); // regola la dimensione in base al tuo plot
        label->Draw("same");
        labels4.push_back(label);
    }
    for (int i = 0; i < 2; ++i) {
        double xB[5] = {Clas_xBins_4[i][0], Clas_xBins_4[i][1], Clas_xBins_4[i][1], Clas_xBins_4[i][0], Clas_xBins_4[i][0]};
        double Q2[5] = {Clas_Q2Bins[3][0], Clas_Q2Bins[3][0], Clas_Q2Bins[3][1], Clas_Q2Bins[3][1], Clas_Q2Bins[3][0]};
        TPolyLine *rect = new TPolyLine(5, xB, Q2);
        //rect->SetLineWidth(2);
        rect->Draw("same");
        rectangles5.push_back(rect);
    }
    std::vector<TText*> labels5;
    for (int i = 0; i < 2; ++i) {
        double x_center = 0.5 * (Clas_xBins_4[i][0] + Clas_xBins_4[i][1]);
        double Q2_center = 0.5 * (Clas_Q2Bins[3][0] + Clas_Q2Bins[3][1]);
        TText *label = new TText(x_center, Q2_center, Form("%d", i+16)); // etichetta come numero del bin
        label->SetTextAlign(22); // centrato
        label->SetTextSize(0.025); // regola la dimensione in base al tuo plot
        label->Draw("same");
        labels5.push_back(label);
    }
    c11->Update();
    c11->Write();
    for(int i = 0; i < xB_nBins; i++){
        hist_xBvsQ2_binned[i]->Write();
        hist_x_dist_2D_bin[i]->Write();
        hist_Q2_dist_2D_bin[i]->Write();
    }


    // ZPT DIRECTORY
    dir_zPt->cd();
    TCanvas *c22 = new TCanvas("z_vs_PhT_Bin", "z vs P_{hT} bin", 800, 700);
    c22->SetLogz();
    Hist_zVsPt->Draw("COLZ");
    std::vector<TPolyLine*> gridLines2;
    std::vector<TPolyLine*> gridLines3;
    std::vector<TPolyLine*> gridLines4;
    int binCounter = 1; // Se vuoi numerare i rettangoli da 1 a N
    std::vector<TText*> labels;
    // Disegna i rettangoli per ogni bin
    for (int j = 0; j < 5; ++j) {
        if (j < 2) {
            for (int i = 0; i < 7; ++i) {
                double p_temp[] = {Clas_PtBins[j][0], Clas_PtBins[j][0], Clas_PtBins[j][1], Clas_PtBins[j][1], Clas_PtBins[j][0]};
                double z_temp[] = {Clas_zBins_12[i][0], Clas_zBins_12[i][1], Clas_zBins_12[i][1], Clas_zBins_12[i][0], Clas_zBins_12[i][0]};
                TPolyLine *rect = new TPolyLine(5, z_temp, p_temp);
                rect->SetLineColor(kBlack);
                rect->Draw("same");
                gridLines2.push_back(rect);
                // Calcola il centro
                double p_center = 0.5 * (Clas_PtBins[j][0] + Clas_PtBins[j][1]);
                double z_center = 0.5 * (Clas_zBins_12[i][0] + Clas_zBins_12[i][1]);
                // Disegna il numero
                TText *label = new TText(z_center, p_center, Form("%d", binCounter));
                label->SetTextAlign(22);
                label->SetTextSize(0.025);
                label->Draw("same");
                labels.push_back(label);
                binCounter++;
            }
        } else if (j < 4) {
            for (int i = 0; i < 6; ++i) {
                double p_temp[] = {Clas_PtBins[j][0], Clas_PtBins[j][0], Clas_PtBins[j][1], Clas_PtBins[j][1], Clas_PtBins[j][0]};
                double z_temp[] = {Clas_zBins_34[i][0], Clas_zBins_34[i][1], Clas_zBins_34[i][1], Clas_zBins_34[i][0], Clas_zBins_34[i][0]};
                TPolyLine *rect = new TPolyLine(5, z_temp, p_temp);
                rect->SetLineColor(kBlack);
                rect->Draw("same");
                gridLines3.push_back(rect);
                double p_center = 0.5 * (Clas_PtBins[j][0] + Clas_PtBins[j][1]);
                double z_center = 0.5 * (Clas_zBins_34[i][0] + Clas_zBins_34[i][1]);
                TText *label = new TText(z_center, p_center, Form("%d", binCounter));
                label->SetTextAlign(22);
                label->SetTextSize(0.025);
                label->Draw("same");
                labels.push_back(label);
                binCounter++;
            }
        } else {
            for (int i = 0; i < 4; ++i) {
                double p_temp[] = {Clas_PtBins[j][0], Clas_PtBins[j][0], Clas_PtBins[j][1], Clas_PtBins[j][1], Clas_PtBins[j][0]};
                double z_temp[] = {Clas_zBins_5[i][0], Clas_zBins_5[i][1], Clas_zBins_5[i][1], Clas_zBins_5[i][0], Clas_zBins_5[i][0]};
                TPolyLine *rect = new TPolyLine(5, z_temp, p_temp);
                rect->SetLineColor(kBlack);
                rect->Draw("same");
                gridLines4.push_back(rect);
                double p_center = 0.5 * (Clas_PtBins[j][0] + Clas_PtBins[j][1]);
                double z_center = 0.5 * (Clas_zBins_5[i][0] + Clas_zBins_5[i][1]);
                TText *label = new TText(z_center, p_center, Form("%d", binCounter));
                label->SetTextAlign(22);
                label->SetTextSize(0.025);
                label->Draw("same");
                labels.push_back(label);
                binCounter++;
            }
        }
    }
    c22->Update();
    c22->Write();
    for(int j = 0; j < z_nBins; j++){
        if(j == 0 || j == 7) continue;
        hist_zvsPt_binned[j]->Write();
        hist_z_dist_2D_bin[j]->Write();
        hist_PhT_dist_2D_bin[j]->Write();
    }
    dir_xQ2_zPt->cd();
    for(int x = 0; x < xB_nBins; x++){
        for(int z = 0; z < z_nBins; z++){
            //if(z == 0 || z == 7) continue; AAA
            hist_dist_xB_4d[x][z]->Write();
            hist_dist_Q2_4d[x][z]->Write();
            hist_dist_z_4d[x][z]->Write();
            hist_dist_PhT_4d[x][z]->Write();
        }
    }


    // HERMESS AND COMPASS DIRECTORY
    dir_comparison->cd();
    TGraphErrors* mean_xvsQ2 = new TGraphErrors(xB_nBins+1); // metto più 1 così da avere un punto vuoto sullo zero e modificare l'asse x
    TGraphErrors* mean_HERMES = new TGraphErrors(4);
    TGraphErrors* mean_COMPASS = new TGraphErrors();
    int id = 0;
    int id_comp = 0;
    for(int j = 0; j < xB_nBins; j++){
        double xB_mean = hist_x_dist_2D_bin[j]->GetMean();
        double Q2_mean = hist_Q2_dist_2D_bin[j]->GetMean();
        mean_xvsQ2->SetPoint(j, xB_mean , Q2_mean);
        //mean_xvsQ2->SetPointError(j, x_meanError_2D[j], Q2_meanError_2D[j]);
        double sum_x = 0.0;
        double sum_Q2 = 0.0;
        double sum_x_comp = 0;
        double sum_Q2_comp = 0;
        int count = 0;
        int count_comp = 0;
        if(j <= 3){
            for (const auto& pt : hermes_data) {
                if (pt.x >= HERMES_xBins[j][0] && pt.x < HERMES_xBins[j][1]) {
                    sum_x += pt.x;
                    sum_Q2 += pt.Q2;
                    count++;
                }
            }
            if (count > 0) {
                double mean_x = sum_x / count;
                double mean_Q2 = sum_Q2 / count;
                mean_HERMES->SetPoint(id, mean_x, mean_Q2);
                mean_HERMES->SetPointError(id, 0, 0);  // se non hai errore associato
                id++;
            }
            for(int p = 0; p < nBins_Pt; p++){
                for(int z = 0; z < nBins_z; z++){
                    for (const auto& pt : compass_data) {
                        if (pt.Q2 >= COMPASS_Q2Bins[j][0] && pt.Q2 < COMPASS_Q2Bins[j][1]) {
                            mean_COMPASS->SetPoint(id_comp, pt.x, pt.Q2);
                            mean_COMPASS->SetPointError(id_comp, 0, 0);  // se non hai errore associato
                            id_comp++;
                        }
                    }
                }
            }
            if (count_comp > 0) {
                double mean_x = sum_x_comp / count_comp;
                double mean_Q2 = sum_Q2_comp / count_comp;
                mean_COMPASS->SetPoint(id_comp, mean_x, mean_Q2);
                mean_COMPASS->SetPointError(id_comp, 0, 0);  // se non hai errore associato
                id_comp++;
            }
        }
    }
    TCanvas* can_mean_xvsQ2 = new TCanvas("mean_xB_vs_Q2", "<x_{B}> vs <Q^{2}> distribution | pion+; #LTx_{B}#GT; #LTQ^{2}#GT");
    can_mean_xvsQ2->cd();
    can_mean_xvsQ2->SetLogy();
    mean_COMPASS->SetTitle("<x_{B}> vs <Q^{2}> distribution for multi-dim SSA analysis | pion+; #LTx_{B}#GT; #LTQ^{2}#GT [GeV^{2}]");
    mean_xvsQ2->SetMarkerStyle(20);
    mean_xvsQ2->SetMarkerSize(1.5);
    mean_xvsQ2->SetLineColor(kBlue-2);  
    mean_xvsQ2->SetMarkerColor(kBlue-2);
    mean_xvsQ2->GetYaxis()->SetMaxDigits(3);
    mean_xvsQ2->GetXaxis()->SetRangeUser(0.01, 0.7);
    mean_COMPASS->GetYaxis()->SetRangeUser(1, 37.5);
    mean_HERMES->SetMarkerStyle(24);
    mean_HERMES->SetMarkerSize(1.5);
    mean_HERMES->SetLineColor(kRed-3);  // Different colors per xQ region
    mean_HERMES->SetMarkerColor(kRed-3);
    mean_COMPASS->SetMarkerStyle(25);
    mean_COMPASS->SetMarkerSize(1.5);
    mean_COMPASS->SetLineColor(kGreen-2);  // Different colors per xQ region
    mean_COMPASS->SetMarkerColor(kGreen-2);
    mean_COMPASS->Draw("AP");
    mean_HERMES->Draw("P SAME");
    mean_xvsQ2->Draw("P SAME");
    /*
    latex.SetTextSize(0.025);
    latex.SetTextAlign(12); // left-bottom alignment
    for (int j = 0; j < nBins; j++) {
        double x, y;
        mean_xvsQ2->GetPoint(j, x, y);
        TString label = Form("Bin %d", j+1);  // customize as needed
        latex.DrawLatex(x + 0.005, y + 0.05, label);  // small offset to avoid overlapping the marker
    }
        */
    TLegend *leg_A = new TLegend(0.7, 0.65, 0.88, 0.8);
    leg_A->AddEntry(mean_xvsQ2, "CLAS12 RGH (4D)", "ep");
    leg_A->AddEntry(mean_HERMES, "HERMES JHEP12(2020) (3D)", "epf");
    leg_A->AddEntry(mean_COMPASS, "COMPASS 1512.06590 (2D)", "epf");
    leg_A->Draw();
    can_mean_xvsQ2->Update();
    can_mean_xvsQ2->Write();


    TCanvas* c_single_plot = new TCanvas("Pt_over_zQ_comparison", "", 800, 800);
    TGraphErrors* g_CLAS = new TGraphErrors();
    TGraphErrors* g_COMPASS = new TGraphErrors();
    TGraphErrors* g_HERMES = new TGraphErrors();
    int idxo = 0;
    for (const auto& pt : compass_data) {
        float y_axis = pt.Pt / (sqrt(pt.Q2) * pt.z);
        g_COMPASS->SetPoint(idxo, pt.x, y_axis);
        g_COMPASS->SetPointError(idxo, 0, 0); 
        idxo++;
    }
    int pointIndex = 0;
    for(int z = 0; z < z_nBins; z++){
        int idx = 0;
        for (const auto& pt : hermes_data) {
            float y_axis = pt.Pt / (sqrt(pt.Q2) * pt.z);
            g_HERMES->SetPoint(idx, pt.x, y_axis);
            g_HERMES->SetPointError(idx, 0, 0); 
            idx++;
        }
        for(int j = 0; j < xB_nBins; j++){
            float A = hist_dist_PhT_4d[j][z]->GetMean();
            float B = hist_dist_Q2_4d[j][z]->GetMean();
            float C = hist_dist_z_4d[j][z]->GetMean();
            float sigma_A = hist_dist_PhT_4d[j][z]->GetMeanError();
            float sigma_B = hist_dist_Q2_4d[j][z]->GetMeanError();
            float sigma_C = hist_dist_z_4d[j][z]->GetMeanError();
            float xValue = hist_dist_xB_4d[j][z]->GetMean(); 
            float yValue = A / (sqrt(B) * C);
            float yError = yValue * sqrt( pow(sigma_A / A, 2) + pow(sigma_B / B, 2) + pow(sigma_C / C, 2) );
            g_CLAS->SetPoint(pointIndex, xValue, yValue);
            g_CLAS->SetPointError(pointIndex, 0, yError);
            //g_CLAS->SetPoint(30, 0.6, -0.1);
            pointIndex++;
            //if (i == 1) cout << "i: " << i << " j: " << j << " z: " << z << " p: " << p << " Q2 mean: " << hist_mean_Q2_4d[j][z]->GetMean() << endl;
        }
    }
    g_CLAS->SetPoint(pointIndex+1, 0, -1);
    g_CLAS->SetPoint(pointIndex+2, 0.6, -1);
    g_CLAS->SetTitle("Phase space population for all the avaiable bin");
    g_CLAS->GetXaxis()->SetTitleSize(0.04);
    g_CLAS->GetYaxis()->SetTitle("P_{hT} / (z#upointQ)");
    g_CLAS->GetXaxis()->SetTitle("<x_{B}>");
    g_CLAS->SetMarkerStyle(20);
    g_CLAS->SetMarkerSize(0.8);
    g_CLAS->SetMarkerColor(kBlue - 2);
    g_CLAS->SetLineColor(kBlue - 2);
    g_CLAS->GetXaxis()->SetRangeUser(0.0, 0.6);
    g_CLAS->GetYaxis()->SetRangeUser(0, 2.7);
    g_CLAS->Draw("AP");
    g_HERMES->SetMarkerStyle(24); // Open square
    g_HERMES->SetMarkerSize(1.2);
    g_HERMES->SetMarkerColor(kRed-3);
    g_HERMES->SetLineColor(kRed-3);
    g_HERMES->SetTitle("Phase space population for all the avaiable bin");
    g_HERMES->Draw("P SAME");
    g_COMPASS->SetMarkerStyle(25); 
    g_COMPASS->SetMarkerSize(1.2);
    g_COMPASS->SetTitle("Phase space population for all the avaiable bin");
    g_COMPASS->SetMarkerColor(kGreen-2);
    g_COMPASS->SetLineColor(kGreen-2);
    g_COMPASS->GetXaxis()->SetRangeUser(0., 0.7);
    g_COMPASS->Draw("P SAME");
    TLegend *leg_All = new TLegend(0.7, 0.75, 0.88, 0.88);
    leg_All->AddEntry(g_CLAS, "CLAS12 RGH (4D)", "ep");
    leg_All->AddEntry(g_HERMES, "HERMES (3D)", "epf");
    leg_All->AddEntry(g_COMPASS, "COMPASS (2D)", "epf");
    leg_All->Draw();

    c_single_plot->Update();
    c_single_plot->Write();


// ________________________________________________________________________________________________________________________________________________________

    // MLE calculation
    for (int x = 0; x < xB_nBins; x++){
        ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE([&](const double* p) {return AUT_loglike(p, vec_pion_phi_h[x], vec_pion_phi_s[x], vec_pion_depol[x], vec_pion_helicity[x]);}, 3);
        minimizer.SetFunction(MLE);
        minimizer.SetMaxFunctionCalls(10000);
        minimizer.SetMaxIterations(10000);
        minimizer.SetTolerance(0.001);
        minimizer.SetVariable(0, "A_Sivers", mean_A_siv, 0.01); // initial guess 0, step 0.05
        minimizer.SetVariable(1, "A_Collins", mean_A_col, 0.01); 
        minimizer.SetVariable(2, "A_Pretz", 0.01, 0.005);
        minimizer.Minimize();
        A_UT_sivers[x] = minimizer.X()[0];
        A_UT_sivers_err[x] = minimizer.Errors()[0];
        A_UT_collins[x] = minimizer.X()[1];
        A_UT_collins_err[x] = minimizer.Errors()[1];
        A_UT_pretz[x] = minimizer.X()[2];
        A_UT_pretz_err[x] = minimizer.Errors()[2];
        //cout << " A_UT = " << A_UT_Sivers[x] << " ± " << A_UT_Sivers_err[x] << endl;
    }
    TH1F* Hist_system2d = new TH1F("AUT_Siv_syst_2d", "syst | A_{UT}^{err} < 0.04; #delta; count", 30, -2, 2);
    for (int i = 0; i < A_UT_sivers.size(); i++){
        Hist_system2d->Fill((Siv_mean_2D[i] - A_UT_sivers[i])/Siv_mean_2D[i]) ;
    } // abs(abs(Siv_mean_4D[i]) - abs(aut_4d[i]/(0.85*depo[i]))) / (abs(Siv_mean_4D[i]) + abs(aut_4d[i]/(0.85*depo[i])))
    Hist_system2d->Write();


    dir_Sivers_ext_xB->cd();
    TGraphErrors* graph_AUT_vs_xB = new TGraphErrors();
    int p_idx = 0;
    for (int x = 0; x < xB_nBins; x++) {
        if (vec_pion_xB[x].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_xB = 0.0;
        for (double val : vec_pion_xB[x]) sum_xB += val;
        double mean_xB = sum_xB / vec_pion_xB[x].size();
        vec_pion_xB_mean.push_back(mean_xB);
        // Get A_UT and error from vectors
        double A_UT = A_UT_sivers[x];
        double A_UT_err = A_UT_sivers_err[x];
        if(A_UT > 1 || A_UT < -1) continue;
        // Fill the graph
        graph_AUT_vs_xB->SetPoint(p_idx, mean_xB, A_UT);
        graph_AUT_vs_xB->SetPointError(p_idx, 0.0, A_UT_err); // No x error
        p_idx++;
    }
    graph_AUT_vs_xB->SetTitle("Sivers Asymmetry vs x_{B} (x_{B}-Q^{2} binning); x_{B}; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}");
    graph_AUT_vs_xB->SetMarkerStyle(20);
    graph_AUT_vs_xB->SetLineColor(kAzure-7);
    graph_AUT_vs_xB->SetMarkerColor(kAzure-7);

    TCanvas* c_Aut_xB = new TCanvas("Aut_Sivers_vs_xB", "Sivers Asymmetry vs x_{B}", 800, 600);
    TLine* guideLine = new TLine(graph_AUT_vs_xB->GetXaxis()->GetXmin(), 0, graph_AUT_vs_xB->GetXaxis()->GetXmax(), 0);
    guideLine->SetLineStyle(2);  
    guideLine->SetLineColor(kGray+1);
    graph_AUT_vs_xB->Draw("AP");
    guideLine->Draw();
    c_Aut_xB->Update();
    c_Aut_xB->Write();


    // Q2 color scale
    TGraphErrors* graph_AUT_vs_xB_1 = new TGraphErrors();
    TGraphErrors* graph_AUT_vs_xB_2 = new TGraphErrors();
    TGraphErrors* graph_AUT_vs_xB_3 = new TGraphErrors();
    TGraphErrors* graph_AUT_vs_xB_4 = new TGraphErrors();
    int p_idx_1 = 0, p_idx_2 = 0, p_idx_3 = 0, p_idx_4 = 0;
    for (int x = 0; x < xB_nBins; x++) {
        if (vec_pion_xB[x].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_xB = 0.0;
        for (double val : vec_pion_xB[x]) sum_xB += val;
        double mean_xB = sum_xB / vec_pion_xB[x].size();
        // Get A_UT and error from vectors
        double A_UT = A_UT_sivers[x];
        double A_UT_err = A_UT_sivers_err[x];
        if(A_UT > 1 || A_UT < -1) continue;
        // Fill the graph
        if(x < 7){
            graph_AUT_vs_xB_1->SetPoint(p_idx_1, mean_xB, A_UT);
            graph_AUT_vs_xB_1->SetPointError(p_idx_1, 0.0, A_UT_err); // No x error
            p_idx_1++;
        } else if (x < 12){
            graph_AUT_vs_xB_2->SetPoint(p_idx_2, mean_xB, A_UT);
            graph_AUT_vs_xB_2->SetPointError(p_idx_2, 0.0, A_UT_err);
            p_idx_2++;
        } else if (x < 15){
            graph_AUT_vs_xB_3->SetPoint(p_idx_3, mean_xB, A_UT);
            graph_AUT_vs_xB_3->SetPointError(p_idx_3, 0.0, A_UT_err);
            p_idx_3++;
        } else if (x < 18){
            graph_AUT_vs_xB_4->SetPoint(p_idx_4, mean_xB, A_UT);
            graph_AUT_vs_xB_4->SetPointError(p_idx_4, 0.0, A_UT_err);
            p_idx_4++;
        }
    }
    graph_AUT_vs_xB_1->SetTitle("Sivers Asymmetry vs x_{B} (x_{B}-Q^{2} binning); x_{B}; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}");
    graph_AUT_vs_xB_1->SetMarkerStyle(20), graph_AUT_vs_xB_2->SetMarkerStyle(20), graph_AUT_vs_xB_3->SetMarkerStyle(20), graph_AUT_vs_xB_4->SetMarkerStyle(20);;
    graph_AUT_vs_xB_1->SetLineColor(kAzure+5);
    graph_AUT_vs_xB_1->SetMarkerColor(kAzure+5);
    graph_AUT_vs_xB_2->SetLineColor(kViolet+5);
    graph_AUT_vs_xB_2->SetMarkerColor(kViolet+5);
    graph_AUT_vs_xB_3->SetLineColor(kPink+5);
    graph_AUT_vs_xB_3->SetMarkerColor(kPink+5);
    graph_AUT_vs_xB_4->SetLineColor(kOrange+5);
    graph_AUT_vs_xB_4->SetMarkerColor(kOrange+5);

    TCanvas* c_Aut_xB_2 = new TCanvas("Aut_Sivers_vs_xB_sepQ2", "Sivers Asymmetry vs x_{B}", 800, 600);
    graph_AUT_vs_xB->Draw("A");
    graph_AUT_vs_xB_1->Draw("P SAME");
    graph_AUT_vs_xB_2->Draw("P SAME");
    graph_AUT_vs_xB_3->Draw("P SAME");
    graph_AUT_vs_xB_4->Draw("P SAME");
    guideLine->Draw();

    TLegend* legend = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
    legend->AddEntry(graph_AUT_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
    legend->AddEntry(graph_AUT_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
    legend->AddEntry(graph_AUT_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
    legend->AddEntry(graph_AUT_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 10", "ep");
    legend->SetFillStyle(0);  // Transparent background
    legend->Draw();
    c_Aut_xB_2->Update();
    c_Aut_xB_2->Write();


    // MLE calculation 4D
    for (int x = 0; x < xB_nBins; x++){
        for(int z = 0; z < z_nBins; z++){
            ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
            ROOT::Math::Functor MLE([&](const double* p) {return AUT_loglike(p, vec_pion_phi_h_4d[x][z], vec_pion_phi_s_4d[x][z], vec_pion_depol_4d[x][z], vec_pion_helicity_4d[x][z]);}, 3);
            minimizer.SetFunction(MLE);
            minimizer.SetMaxFunctionCalls(10000);
            minimizer.SetMaxIterations(10000);
            minimizer.SetTolerance(0.001);
            minimizer.SetVariable(0, "A_Sivers", mean_A_siv, 0.01); 
            minimizer.SetVariable(1, "A_Collins", mean_A_col, 0.01);
            minimizer.SetVariable(2, "A_Pretz", 0.01, 0.005);
            minimizer.Minimize();
            A_UT_sivers_4d[x][z] = minimizer.X()[0];
            A_UT_sivers_err_4d[x][z] = minimizer.Errors()[0];
            A_UT_collins_4d[x][z] = minimizer.X()[1];
            A_UT_collins_err_4d[x][z] = minimizer.Errors()[1];
            A_UT_pretz_4d[x][z] = minimizer.X()[2];
            A_UT_pretz_err_4d[x][z] = minimizer.Errors()[2];
        }
    }

    TGraphErrors* graph_AUT_vs_xB_4D = new TGraphErrors();
    int p_mean_idx = 0;

    for (int x = 0; x < xB_nBins; x++) {
        double sum_A_UT = 0.0;
        double sum_A_UT_err2 = 0.0;
        double sum_weights = 0.0;
        double sum_xB = 0.0;
        int count = 0;

        for (int z = 0; z < z_nBins; z++) {
            if (vec_pion_xB_4d[x][z].empty()) continue;

            double A_UT = A_UT_sivers_4d[x][z];
            double A_UT_err = A_UT_sivers_err_4d[x][z];
            if (A_UT_err < 0.1) {
                // Valore medio xB per questo bin z
                double mean_xB = 0.0;
                for (double val : vec_pion_xB_4d[x][z]) mean_xB += val;
                mean_xB /= vec_pion_xB_4d[x][z].size();

                // Media pesata di A_UT
                double weight = 1.0 / (A_UT_err * A_UT_err);
                sum_A_UT += A_UT * weight;
                sum_weights += weight;
                sum_xB += mean_xB;
                count++;
            }
        }

        if (count > 0 && sum_weights > 0) {
            double A_UT_mean = sum_A_UT / sum_weights;
            double A_UT_mean_err = std::sqrt(1.0 / sum_weights);
            double xB_mean = sum_xB / count;

            graph_AUT_vs_xB_4D->SetPoint(p_mean_idx, xB_mean, A_UT_mean);
            graph_AUT_vs_xB_4D->SetPointError(p_mean_idx, 0.0, A_UT_mean_err);
            p_mean_idx++;
        }
    }

    graph_AUT_vs_xB_4D->SetTitle("Sivers Asymmetry vs x_{B} 4D; x_{B}; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}");
    graph_AUT_vs_xB_4D->SetMarkerStyle(20);
    graph_AUT_vs_xB_4D->SetLineColor(kAzure-7);
    graph_AUT_vs_xB_4D->SetMarkerColor(kAzure-7);

    TCanvas* c_Aut_xB_4D = new TCanvas("Aut_Sivers_vs_xB_4D", "Sivers Asymmetry vs x_{B} in 4D", 800, 600);
    TLine* guideLine_4D = new TLine(graph_AUT_vs_xB_4D->GetXaxis()->GetXmin(), 0, graph_AUT_vs_xB_4D->GetXaxis()->GetXmax(), 0);
    guideLine_4D->SetLineStyle(2);  
    guideLine_4D->SetLineColor(kGray+1);
    graph_AUT_vs_xB_4D->Draw("AP");
    //guideLine_4D->Draw();
    c_Aut_xB_4D->Update();
    c_Aut_xB_4D->Write();

    vector<double> aut_4d;
    vector<double> aut_err_4d;
    vector<double> aut_col_4d;
    vector<double> aut_col_err_4d;
    for(int x = 0; x<xB_nBins; x++){
        for(int z = 0; z<z_nBins; z++){
            aut_4d.push_back(A_UT_sivers_4d[x][z]);
            aut_err_4d.push_back(A_UT_sivers_err_4d[x][z]);
            aut_col_4d.push_back(A_UT_collins_4d[x][z]);
            aut_col_err_4d.push_back(A_UT_collins_err_4d[x][z]);
        }
    }
    TH1F* Hist_system = new TH1F("AUT_Siv_syst", "syst | A_{UT}^{err} < 0.04; #delta; count", 60, -2, 2);
    TH1F* Hist_system_col = new TH1F("AUT_Col_syst", "syst | A_{UT}^{err} < 0.04; #delta; count", 60, -2, 2);
    TH2F* Hist_system_2 = new TH2F("AUT_Siv_syst_2", "syst | A_{UT}^{err} < 0.04; MC; RECO", 30, 0, 0.12, 30, -0.2, 0.2);
    for (int i = 0; i < aut_4d.size(); i++){
        if(aut_err_4d[i] < 0.04){
            Hist_system->Fill(abs(Siv_mean_4D[i] - aut_4d[i])/abs(Siv_mean_4D[i])) ;
            Hist_system_2->Fill(abs(Siv_mean_4D[i]), abs(aut_4d[i]));
        }
        if(aut_col_err_4d[i] < 0.04){
            Hist_system_col->Fill(abs(Col_mean_4D[i] - aut_col_4d[i])/abs(Col_mean_4D[i])) ;
        }
    } 
    Hist_system->Write();
    Hist_system_col->Write();
    Hist_system_2->Write();

    // Draw the plot of AUT

    // Q2 color scale vs xB 4D
    for (int z = 0; z < z_nBins; z++) {
        //if(z == 0 || z == 7) continue;
        TGraphErrors* graph_AUT_vs_xB = new TGraphErrors();
        TGraphErrors* graph_AUT_vs_xB_1 = new TGraphErrors();
        TGraphErrors* graph_AUT_vs_xB_2 = new TGraphErrors();
        TGraphErrors* graph_AUT_vs_xB_3 = new TGraphErrors();
        TGraphErrors* graph_AUT_vs_xB_4 = new TGraphErrors();
        int p_idx = 0, p_idx_1 = 0, p_idx_2 = 0, p_idx_3 = 0, p_idx_4 = 0;
        for(int x = 0; x < xB_nBins; x++){
            if (vec_pion_xB_4d[x][z].empty()) continue; // Skip empty bins
            if (vec_pion_z_4d[x][z].empty()) continue;
            // Compute mean xB for this bin
            double sum_xB = 0.0;
            for (double val : vec_pion_xB_4d[x][z]) sum_xB += val;
            double mean_xB = sum_xB / vec_pion_xB_4d[x][z].size();
            // Get A_UT and error from vectors
            double A_UT = A_UT_sivers_4d[x][z];
            double A_UT_err = A_UT_sivers_err_4d[x][z];
            if(A_UT > 1 || A_UT < -1) continue;
            // Fill the graph
            graph_AUT_vs_xB->SetPoint(p_idx, mean_xB, A_UT);
            graph_AUT_vs_xB->SetPointError(p_idx, 0.0, A_UT_err); // No x error
            p_idx++;
            if(x < 7){
                graph_AUT_vs_xB_1->SetPoint(p_idx_1, mean_xB, A_UT);
                graph_AUT_vs_xB_1->SetPointError(p_idx_1, 0.0, A_UT_err); // No x error
                p_idx_1++;
            } else if (x < 12){
                graph_AUT_vs_xB_2->SetPoint(p_idx_2, mean_xB, A_UT);
                graph_AUT_vs_xB_2->SetPointError(p_idx_2, 0.0, A_UT_err);
                p_idx_2++;
            } else if (x < 15){
                graph_AUT_vs_xB_3->SetPoint(p_idx_3, mean_xB, A_UT);
                graph_AUT_vs_xB_3->SetPointError(p_idx_3, 0.0, A_UT_err);
                p_idx_3++;
            } else if (x < 18){
                graph_AUT_vs_xB_4->SetPoint(p_idx_4, mean_xB, A_UT);
                graph_AUT_vs_xB_4->SetPointError(p_idx_4, 0.0, A_UT_err);
                p_idx_4++;
            }
        }

        graph_AUT_vs_xB->SetTitle(Form("Sivers Asymmetry vs x_{B} (x_{B}-Q^{2} bin) for bin (%d, z-P_{hT}) ; x_{B}; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}", z+1));
        graph_AUT_vs_xB_1->SetMarkerStyle(20), graph_AUT_vs_xB_2->SetMarkerStyle(20), graph_AUT_vs_xB_3->SetMarkerStyle(20), graph_AUT_vs_xB_4->SetMarkerStyle(20);;
        graph_AUT_vs_xB_1->SetLineColor(kAzure+5);
        graph_AUT_vs_xB_1->SetMarkerColor(kAzure+5);
        graph_AUT_vs_xB_2->SetLineColor(kViolet+5);
        graph_AUT_vs_xB_2->SetMarkerColor(kViolet+5);
        graph_AUT_vs_xB_3->SetLineColor(kPink+5);
        graph_AUT_vs_xB_3->SetMarkerColor(kPink+5);
        graph_AUT_vs_xB_4->SetLineColor(kOrange+5);
        graph_AUT_vs_xB_4->SetMarkerColor(kOrange+5);

        TCanvas* c_Aut_xB_2 = new TCanvas(Form("Aut_Sivers_vs_xB_4d_zPt_bin%d", z+1), "Sivers Asymmetry vs x_{B} for bin (n, z-P_{hT})", 800, 600);
        graph_AUT_vs_xB->Draw("A");
        graph_AUT_vs_xB_1->Draw("P SAME");
        graph_AUT_vs_xB_2->Draw("P SAME");
        graph_AUT_vs_xB_3->Draw("P SAME");
        graph_AUT_vs_xB_4->Draw("P SAME");
        TLine* guideLine2 = new TLine(graph_AUT_vs_xB->GetXaxis()->GetXmin(), 0, graph_AUT_vs_xB->GetXaxis()->GetXmax(), 0);
        guideLine2->SetLineStyle(2);  
        guideLine2->SetLineColor(kGray+1);
        guideLine2->Draw();

        TLegend* legend = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
        legend->AddEntry(graph_AUT_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
        legend->AddEntry(graph_AUT_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
        legend->AddEntry(graph_AUT_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
        legend->AddEntry(graph_AUT_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 10", "ep");
        legend->SetFillStyle(0);  // Transparent background
        legend->Draw();
        c_Aut_xB_2->Update();
        c_Aut_xB_2->Write();
    }

    // ___________________________________________________________________________________________________________________________________________

    // MLE calculation
    for (int z = 0; z < z_nBins; z++){
        ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE([&](const double* p) {return AUT_loglike(p, vec_pion_phi_h_z[z], vec_pion_phi_s_z[z], vec_pion_depol_z[z], vec_pion_helicity_z[z]);}, 3);
        minimizer.SetFunction(MLE);
        minimizer.SetMaxFunctionCalls(10000);
        minimizer.SetMaxIterations(10000);
        minimizer.SetTolerance(0.001);
        minimizer.SetVariable(0, "A_Sivers", mean_A_siv, 0.01); 
        minimizer.SetVariable(1, "A_Collins", mean_A_col, 0.01);
        minimizer.SetVariable(2, "A_Pretz", 0.01, 0.005);
        minimizer.Minimize();
        A_UT_sivers_z[z] = minimizer.X()[0];
        A_UT_sivers_err_z[z] = minimizer.Errors()[0];
        A_UT_collins_z[z] = minimizer.X()[1];
        A_UT_collins_err_z[z] = minimizer.Errors()[1];
        A_UT_pretz_z[z] = minimizer.X()[2];
        A_UT_pretz_err_z[z] = minimizer.Errors()[2];
    }

    dir_Sivers_ext_z->cd();
    TGraphErrors* graph_Sivers_vs_z = new TGraphErrors();
    int p_idxz = 0;
    for (int z = 0; z < z_nBins; z++) {
        if (vec_pion_z[z].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_z = 0.0;
        for (double val : vec_pion_z[z]) sum_z += val;
        double mean_z = sum_z / vec_pion_z[z].size();
        // Get A_UT and error from vectors
        double A_UT = A_UT_sivers_z[z];
        double A_UT_err = A_UT_sivers_err_z[z];
        if(A_UT > 1 || A_UT < -1) continue;
        // Fill the graph
        graph_Sivers_vs_z->SetPoint(p_idxz, mean_z, A_UT);
        graph_Sivers_vs_z->SetPointError(p_idxz, 0.0, A_UT_err); 
        p_idxz++;
    }
    graph_Sivers_vs_z->SetTitle("Sivers Asymmetry vs z (z-P_{hT} binning); z; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}");
    graph_Sivers_vs_z->SetMarkerStyle(20);
    graph_Sivers_vs_z->SetLineColor(kAzure-7);
    graph_Sivers_vs_z->SetMarkerColor(kAzure-7);

    TCanvas* c_Sivers_z = new TCanvas("Aut_Sivers_vs_z", "Sivers Asymmetry vs z", 800, 600);
    TLine* guideLineZ = new TLine(graph_Sivers_vs_z->GetXaxis()->GetXmin(), 0, graph_Sivers_vs_z->GetXaxis()->GetXmax(), 0);
    guideLineZ->SetLineStyle(2);  
    guideLineZ->SetLineColor(kGray+1);
    graph_Sivers_vs_z->Draw("AP");
    guideLineZ->Draw();
    c_Sivers_z->Update();
    c_Sivers_z->Write();

    TGraphErrors* graph_Sivers_vs_z_1 = new TGraphErrors();
    TGraphErrors* graph_Sivers_vs_z_2 = new TGraphErrors();
    TGraphErrors* graph_Sivers_vs_z_3 = new TGraphErrors();
    TGraphErrors* graph_Sivers_vs_z_4 = new TGraphErrors();
    TGraphErrors* graph_Sivers_vs_z_5 = new TGraphErrors();
    int p_idxz_1 = 0, p_idxz_2 = 0, p_idxz_3 = 0, p_idxz_4 = 0, p_idxz_5 = 0;
    for (int z = 0; z < z_nBins; z++) {
        if (vec_pion_z[z].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_z = 0.0;
        for (double val : vec_pion_z[z]) sum_z += val;
        double mean_z = sum_z / vec_pion_z[z].size();
        // Get A_UT and error from vectors
        double A_UT = A_UT_sivers_z[z];
        double A_UT_err = A_UT_sivers_err_z[z];
        if(A_UT > 1 || A_UT < -1) continue;
        // Fill the graph
        if(z < 7){
            graph_Sivers_vs_z_1->SetPoint(p_idxz_1, mean_z, A_UT);
            graph_Sivers_vs_z_1->SetPointError(p_idxz_1, 0.0, A_UT_err); // No x error
            p_idxz_1++;
        } else if (z < 14){
            graph_Sivers_vs_z_2->SetPoint(p_idxz_2, mean_z, A_UT);
            graph_Sivers_vs_z_2->SetPointError(p_idxz_2, 0.0, A_UT_err);
            p_idxz_2++;
        } else if (z < 20){
            graph_Sivers_vs_z_3->SetPoint(p_idxz_3, mean_z, A_UT);
            graph_Sivers_vs_z_3->SetPointError(p_idxz_3, 0.0, A_UT_err);
            p_idxz_3++;
        } else if (z < 26){
            graph_Sivers_vs_z_4->SetPoint(p_idxz_4, mean_z, A_UT);
            graph_Sivers_vs_z_4->SetPointError(p_idxz_4, 0.0, A_UT_err);
            p_idxz_4++;
        } else if (z < 30){
            graph_Sivers_vs_z_5->SetPoint(p_idxz_5, mean_z, A_UT);
            graph_Sivers_vs_z_5->SetPointError(p_idxz_5, 0.0, A_UT_err);
            p_idxz_5++;
        }
    }
    graph_Sivers_vs_z->SetTitle("Sivers Asymmetry vs z (z-P_{hT} bin) ; z; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}");
    graph_Sivers_vs_z_1->SetMarkerStyle(20), graph_Sivers_vs_z_2->SetMarkerStyle(20), graph_Sivers_vs_z_3->SetMarkerStyle(20), graph_Sivers_vs_z_4->SetMarkerStyle(20), graph_Sivers_vs_z_5->SetMarkerStyle(20);
    graph_Sivers_vs_z_1->SetLineColor(kTeal+5);
    graph_Sivers_vs_z_1->SetMarkerColor(kTeal+5);
    graph_Sivers_vs_z_2->SetLineColor(kAzure+5);
    graph_Sivers_vs_z_2->SetMarkerColor(kAzure+5);
    graph_Sivers_vs_z_3->SetLineColor(kViolet+5);
    graph_Sivers_vs_z_3->SetMarkerColor(kViolet+5);
    graph_Sivers_vs_z_4->SetLineColor(kPink+5);
    graph_Sivers_vs_z_4->SetMarkerColor(kPink+5);
    graph_Sivers_vs_z_5->SetLineColor(kOrange+5);
    graph_Sivers_vs_z_5->SetMarkerColor(kOrange+5);

    TCanvas* c_Sivers_z_2 = new TCanvas("Aut_Sivers_vs_z_sepPt", "Sivers Asymmetry vs z for bin (n, x_{B}-Q^{2})", 800, 600);
    graph_Sivers_vs_z->Draw("A");
    graph_Sivers_vs_z_1->Draw("P SAME");
    graph_Sivers_vs_z_2->Draw("P SAME");
    graph_Sivers_vs_z_3->Draw("P SAME");
    graph_Sivers_vs_z_4->Draw("P SAME");
    graph_Sivers_vs_z_5->Draw("P SAME");
    guideLineZ->Draw();

    TLegend* legend_z = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
    legend_z->AddEntry(graph_Sivers_vs_z_1, "0.0 < P_{hT} [GeV] < 0.2", "ep");
    legend_z->AddEntry(graph_Sivers_vs_z_2, "0.2 < P_{hT} [GeV] < 0.4", "ep");
    legend_z->AddEntry(graph_Sivers_vs_z_3, "0.4 < P_{hT} [GeV] < 0.6", "ep");
    legend_z->AddEntry(graph_Sivers_vs_z_4, "0.6 < P_{hT} [GeV] < 0.8", "ep");
    legend_z->AddEntry(graph_Sivers_vs_z_5, "0.8 < P_{hT} [GeV] < 1.2", "ep");
    legend_z->SetFillStyle(0);  // Transparent background
    legend_z->Draw();
    gPad->Update();
    c_Sivers_z_2->Update();
    c_Sivers_z_2->Write();

    gStyle->SetOptStat(0);
    TH1F* Sivers_z_frame = new TH1F("", "", 1, graph_Sivers_vs_z->GetXaxis()->GetXmin(), graph_Sivers_vs_z->GetXaxis()->GetXmax());
    Sivers_z_frame->SetMinimum(graph_Sivers_vs_z->GetYaxis()->GetXmin());
    Sivers_z_frame->SetMaximum(graph_Sivers_vs_z->GetYaxis()->GetXmax());
    TCanvas* c_Sivers_zBin_1 = new TCanvas("Aut_Sivers_vs_z_Pt_1", "Sivers Asymmetry vs z", 800, 600);
    Sivers_z_frame->SetTitle("Sivers Asymmetry vs z | pion+ | 0.0 < P_{hT} [GeV] < 0.2  ; z; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}");
    Sivers_z_frame->Draw();
    graph_Sivers_vs_z_1->Draw("P SAME");
    guideLineZ->Draw();
    gPad->Update();
    c_Sivers_zBin_1->Update();
    c_Sivers_zBin_1->Write();
    TCanvas* c_Sivers_zBin_2 = new TCanvas("Aut_Sivers_vs_z_Pt_2", "Sivers Asymmetry vs z", 800, 600);
    Sivers_z_frame->SetTitle("Sivers Asymmetry vs z | pion+ | 0.2 < P_{hT} [GeV] < 0.4  ; z; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}");
    Sivers_z_frame->Draw();
    graph_Sivers_vs_z_2->Draw("P SAME");
    guideLineZ->Draw();
    gPad->Update();
    c_Sivers_zBin_2->Update();
    c_Sivers_zBin_2->Write();
    TCanvas* c_Sivers_zBin_3 = new TCanvas("Aut_Sivers_vs_z_Pt_3", "Sivers Asymmetry vs z", 800, 600);
    Sivers_z_frame->SetTitle("Sivers Asymmetry vs z | pion+ | 0.4 < P_{hT} [GeV] < 0.6  ; z; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}");
    Sivers_z_frame->Draw();
    graph_Sivers_vs_z_3->Draw("P SAME");
    guideLineZ->Draw();
    gPad->Update();
    c_Sivers_zBin_3->Update();
    c_Sivers_zBin_3->Write();
    TCanvas* c_Sivers_zBin_4 = new TCanvas("Aut_Sivers_vs_z_Pt_4", "Sivers Asymmetry vs z", 800, 600);
    Sivers_z_frame->SetTitle("Sivers Asymmetry vs z | pion+ | 0.6 < P_{hT} [GeV] < 0.8  ; z; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}");
    Sivers_z_frame->Draw();
    graph_Sivers_vs_z_4->Draw("P SAME");
    guideLineZ->Draw();
    gPad->Update();
    c_Sivers_zBin_4->Update();
    c_Sivers_zBin_4->Write();
    TCanvas* c_Sivers_zBin_5 = new TCanvas("Aut_Sivers_vs_z_Pt_5", "Sivers Asymmetry vs z", 800, 600);
    Sivers_z_frame->SetTitle("Sivers Asymmetry vs z | pion+ | 0.8 < P_{hT} [GeV] < 1.2  ; z; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}");
    Sivers_z_frame->Draw();
    graph_Sivers_vs_z_5->Draw("P SAME");
    guideLineZ->Draw();
    gPad->Update();
    c_Sivers_zBin_5->Update();
    c_Sivers_zBin_5->Write();

    // Pt color scale vs z 4D
    for (int x = 0; x < xB_nBins; x++) {
        TGraphErrors* graph_Sivers_vs_z = new TGraphErrors();
        TGraphErrors* graph_Sivers_vs_z_1 = new TGraphErrors();
        TGraphErrors* graph_Sivers_vs_z_2 = new TGraphErrors();
        TGraphErrors* graph_Sivers_vs_z_3 = new TGraphErrors();
        TGraphErrors* graph_Sivers_vs_z_4 = new TGraphErrors();
        TGraphErrors* graph_Sivers_vs_z_5 = new TGraphErrors();
        int p_idx = 0, p_idx_1 = 0, p_idx_2 = 0, p_idx_3 = 0, p_idx_4 = 0, p_idx_5 = 0;
        for(int z = 0; z < z_nBins; z++){
            //if(z == 0 || z == 7) continue; AAA
            if (vec_pion_xB_4d[x][z].empty()) continue; // Skip empty bins
            if (vec_pion_z_4d[x][z].empty()) continue;
            // Compute mean xB for this bin
            double sum_z = 0.0;
            for (double val : vec_pion_z_4d[x][z]) sum_z += val;
            double mean_z = sum_z / vec_pion_z_4d[x][z].size();
            // Get A_UT and error from vectors
            double A_UT = A_UT_sivers_4d[x][z];
            double A_UT_err = A_UT_sivers_err_4d[x][z];
            if(A_UT > 1 || A_UT < -1) continue;
            // Fill the graph
            graph_Sivers_vs_z->SetPoint(p_idx, mean_z, A_UT);
            graph_Sivers_vs_z->SetPointError(p_idx, 0.0, A_UT_err); // No x error
            p_idx++;
            if(z < 7){
                graph_Sivers_vs_z_1->SetPoint(p_idx_1, mean_z, A_UT);
                graph_Sivers_vs_z_1->SetPointError(p_idx_1, 0.0, A_UT_err); // No x error
                p_idx_1++;
            } else if (z < 14){
                graph_Sivers_vs_z_2->SetPoint(p_idx_2, mean_z, A_UT);
                graph_Sivers_vs_z_2->SetPointError(p_idx_2, 0.0, A_UT_err);
                p_idx_2++;
            } else if (z < 20){
                graph_Sivers_vs_z_3->SetPoint(p_idx_3, mean_z, A_UT);
                graph_Sivers_vs_z_3->SetPointError(p_idx_3, 0.0, A_UT_err);
                p_idx_3++;
            } else if (z < 26){
                graph_Sivers_vs_z_4->SetPoint(p_idx_4, mean_z, A_UT);
                graph_Sivers_vs_z_4->SetPointError(p_idx_4, 0.0, A_UT_err);
                p_idx_4++;
            } else if (z < 30){
                graph_Sivers_vs_z_5->SetPoint(p_idx_5, mean_z, A_UT);
                graph_Sivers_vs_z_5->SetPointError(p_idx_5, 0.0, A_UT_err);
                p_idx_5++;
            }
        }

        graph_Sivers_vs_z->SetTitle(Form("Sivers Asymmetry vs z (z-P_{hT} bin) for bin (%d, x_{B}-Q^{2}}) ; z; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}", x+1));
        graph_Sivers_vs_z_1->SetMarkerStyle(20), graph_Sivers_vs_z_2->SetMarkerStyle(20), graph_Sivers_vs_z_3->SetMarkerStyle(20), graph_Sivers_vs_z_4->SetMarkerStyle(20), graph_Sivers_vs_z_5->SetMarkerStyle(20);
        graph_Sivers_vs_z_1->SetLineColor(kTeal+5);
        graph_Sivers_vs_z_1->SetMarkerColor(kTeal+5);
        graph_Sivers_vs_z_2->SetLineColor(kAzure+5);
        graph_Sivers_vs_z_2->SetMarkerColor(kAzure+5);
        graph_Sivers_vs_z_3->SetLineColor(kViolet+5);
        graph_Sivers_vs_z_3->SetMarkerColor(kViolet+5);
        graph_Sivers_vs_z_4->SetLineColor(kPink+5);
        graph_Sivers_vs_z_4->SetMarkerColor(kPink+5);
        graph_Sivers_vs_z_5->SetLineColor(kOrange+5);
        graph_Sivers_vs_z_5->SetMarkerColor(kOrange+5);

        TCanvas* c_Sivers_z_2 = new TCanvas(Form("Aut_Sivers_vs_z_4d_xQ2_bin%d", x+1), "Sivers Asymmetry vs z for bin (n, x_{B}-Q^{2})", 800, 600);
        graph_Sivers_vs_z->Draw("A");
        graph_Sivers_vs_z_1->Draw("P SAME");
        graph_Sivers_vs_z_2->Draw("P SAME");
        graph_Sivers_vs_z_3->Draw("P SAME");
        graph_Sivers_vs_z_4->Draw("P SAME");
        graph_Sivers_vs_z_5->Draw("P SAME");
        TLine* guideLine2 = new TLine(graph_Sivers_vs_z->GetXaxis()->GetXmin(), 0, graph_Sivers_vs_z->GetXaxis()->GetXmax(), 0);
        guideLine2->SetLineStyle(2);  
        guideLine2->SetLineColor(kGray+1);
        guideLine2->Draw();

        TLegend* legend = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
        legend->AddEntry(graph_Sivers_vs_z_1, "0 < P_{hT} [GeV] < 0.2", "ep");
        legend->AddEntry(graph_Sivers_vs_z_2, "0.2 < P_{hT} [GeV] < 0.4", "ep");
        legend->AddEntry(graph_Sivers_vs_z_3, "0.4 < P_{hT} [GeV] < 0.6", "ep");
        legend->AddEntry(graph_Sivers_vs_z_4, "0.6 < P_{hT} [GeV] < 0.8", "ep");
        legend->AddEntry(graph_Sivers_vs_z_5, "0.8 < P_{hT} [GeV] < 1.2", "ep");
        legend->SetFillStyle(0);  // Transparent background
        legend->Draw();
        c_Sivers_z_2->Update();
        c_Sivers_z_2->Write();
    }

    // ________________________________________________________________________________________________________________________________________________________

    
    
    
    
    
    
    dir_Collins_ext_xB->cd();
    TGraphErrors* graph_Collins_vs_xB = new TGraphErrors();
    int P_idx_col = 0;
    for (int x = 0; x < xB_nBins; x++) {
        if (vec_pion_xB[x].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_xB = 0.0, sum_eps = 0;
        for (double val : vec_pion_xB[x]) sum_xB += val;
        double mean_xB = sum_xB / vec_pion_xB[x].size();
        for (double valeps : vec_pion_epsilon[x]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_pion_epsilon[x].size();
        // Get A_UT and error from vectors
        double A_UT = A_UT_collins[x];
        double A_UT_err = A_UT_collins_err[x];
        if(A_UT > 1 || A_UT < -1) continue;
        // Fill the graph
        graph_Collins_vs_xB->SetPoint(P_idx_col, mean_xB, A_UT);
        graph_Collins_vs_xB->SetPointError(P_idx_col, 0.0, A_UT_err); // No x error
        P_idx_col++;
    }
    graph_Collins_vs_xB->SetTitle("Collins Asymmetry vs x_{B} (x_{B}-Q^{2} binning); x_{B}; A_{UT}^{sin(#Phi_{h}+#Phi_{s})}");
    graph_Collins_vs_xB->SetMarkerStyle(20);
    graph_Collins_vs_xB->SetLineColor(kAzure-7);
    graph_Collins_vs_xB->SetMarkerColor(kAzure-7);

    TCanvas* c_Collins_xB = new TCanvas("Aut_Collins_vs_xB", "Collins Asymmetry vs x_{B}", 800, 600);
    TLine* guideLine_col = new TLine(graph_Collins_vs_xB->GetXaxis()->GetXmin(), 0, graph_Collins_vs_xB->GetXaxis()->GetXmax(), 0);
    guideLine_col->SetLineStyle(2);  
    guideLine_col->SetLineColor(kGray+1);
    graph_Collins_vs_xB->Draw("AP");
    guideLine_col->Draw();
    c_Collins_xB->Update();
    c_Collins_xB->Write();


    // Q2 color scale
    TGraphErrors* graph_Collins_vs_xB_1 = new TGraphErrors();
    TGraphErrors* graph_Collins_vs_xB_2 = new TGraphErrors();
    TGraphErrors* graph_Collins_vs_xB_3 = new TGraphErrors();
    TGraphErrors* graph_Collins_vs_xB_4 = new TGraphErrors();
    int P_idx_col_1 = 0, P_idx_col_2 = 0, P_idx_col_3 = 0, P_idx_col_4 = 0;
    for (int x = 0; x < xB_nBins; x++) {
        if (vec_pion_xB[x].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_xB = 0.0, sum_eps = 0;
        for (double val : vec_pion_xB[x]) sum_xB += val;
        double mean_xB = sum_xB / vec_pion_xB[x].size();
        for (double valeps : vec_pion_epsilon[x]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_pion_epsilon[x].size();
        // Get A_UT and error from vectors
        double A_UT = A_UT_collins[x];
        double A_UT_err = A_UT_collins_err[x];
        if(A_UT > 1 || A_UT < -1) continue;
        // Fill the graph
        if(x < 7){
            graph_Collins_vs_xB_1->SetPoint(P_idx_col_1, mean_xB, A_UT);
            graph_Collins_vs_xB_1->SetPointError(P_idx_col_1, 0.0, A_UT_err); // No x error
            P_idx_col_1++;
        } else if (x < 12){
            graph_Collins_vs_xB_2->SetPoint(P_idx_col_2, mean_xB, A_UT);
            graph_Collins_vs_xB_2->SetPointError(P_idx_col_2, 0.0, A_UT_err);
            P_idx_col_2++;
        } else if (x < 15){
            graph_Collins_vs_xB_3->SetPoint(P_idx_col_3, mean_xB, A_UT);
            graph_Collins_vs_xB_3->SetPointError(P_idx_col_3, 0.0, A_UT_err);
            P_idx_col_3++;
        } else if (x < 18){
            graph_Collins_vs_xB_4->SetPoint(P_idx_col_4, mean_xB, A_UT);
            graph_Collins_vs_xB_4->SetPointError(P_idx_col_4, 0.0, A_UT_err);
            P_idx_col_4++;
        }
    }
    graph_Collins_vs_xB_1->SetTitle("Collins Asymmetry vs x_{B} (x_{B}-Q^{2} binning); x_{B}; A_{UT}^{sin(#Phi_{h}+#Phi_{s})}");
    graph_Collins_vs_xB_1->SetMarkerStyle(20), graph_Collins_vs_xB_2->SetMarkerStyle(20), graph_Collins_vs_xB_3->SetMarkerStyle(20), graph_Collins_vs_xB_4->SetMarkerStyle(20);;
    graph_Collins_vs_xB_1->SetLineColor(kAzure+5);
    graph_Collins_vs_xB_1->SetMarkerColor(kAzure+5);
    graph_Collins_vs_xB_2->SetLineColor(kViolet+5);
    graph_Collins_vs_xB_2->SetMarkerColor(kViolet+5);
    graph_Collins_vs_xB_3->SetLineColor(kPink+5);
    graph_Collins_vs_xB_3->SetMarkerColor(kPink+5);
    graph_Collins_vs_xB_4->SetLineColor(kOrange+5);
    graph_Collins_vs_xB_4->SetMarkerColor(kOrange+5);

    TCanvas* c_Collins_xB_2 = new TCanvas("Aut_Collins_vs_xB_sepQ2", "Collins Asymmetry vs x_{B}", 800, 600);
    graph_Collins_vs_xB->Draw("A");
    graph_Collins_vs_xB_1->Draw("P SAME");
    graph_Collins_vs_xB_2->Draw("P SAME");
    graph_Collins_vs_xB_3->Draw("P SAME");
    graph_Collins_vs_xB_4->Draw("P SAME");
    guideLine_col->Draw();

    TLegend* legend_col = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
    legend_col->AddEntry(graph_Collins_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
    legend_col->AddEntry(graph_Collins_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
    legend_col->AddEntry(graph_Collins_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
    legend_col->AddEntry(graph_Collins_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 10", "ep");
    legend_col->SetFillStyle(0);  // Transparent background
    legend_col->Draw();
    c_Collins_xB_2->Update();
    c_Collins_xB_2->Write();



    // Draw the plot of Collins
    // Q2 color scale vs xB 4D
    for (int z = 0; z < z_nBins; z++) {
        //if(z == 0 || z == 7) continue;
        TGraphErrors* graph_Collins_vs_xB = new TGraphErrors();
        TGraphErrors* graph_Collins_vs_xB_1 = new TGraphErrors();
        TGraphErrors* graph_Collins_vs_xB_2 = new TGraphErrors();
        TGraphErrors* graph_Collins_vs_xB_3 = new TGraphErrors();
        TGraphErrors* graph_Collins_vs_xB_4 = new TGraphErrors();
        int P_idx_col = 0, P_idx_col_1 = 0, P_idx_col_2 = 0, P_idx_col_3 = 0, P_idx_col_4 = 0;
        for(int x = 0; x < xB_nBins; x++){
            if (vec_pion_xB_4d[x][z].empty()) continue; // Skip empty bins
            if (vec_pion_z_4d[x][z].empty()) continue;
            // Compute mean xB for this bin
            double sum_xB = 0.0, sum_eps = 0;
            for (double val : vec_pion_xB_4d[x][z]) sum_xB += val;
            double mean_xB = sum_xB / vec_pion_xB_4d[x][z].size();
            for (double valeps : vec_pion_epsilon_4d[x][z]) sum_eps += valeps;
            double mean_eps = sum_eps / vec_pion_epsilon_4d[x][z].size();
            // Get A_UT and error from vectors
            double A_UT = A_UT_collins_4d[x][z];
            double A_UT_err = A_UT_collins_err_4d[x][z];
            if(A_UT > 1 || A_UT < -1) continue;
            // Fill the graph
            graph_Collins_vs_xB->SetPoint(P_idx_col, mean_xB, A_UT);
            graph_Collins_vs_xB->SetPointError(P_idx_col, 0.0, A_UT_err); // No x error
            P_idx_col++;
            if(x < 7){
                graph_Collins_vs_xB_1->SetPoint(P_idx_col_1, mean_xB, A_UT);
                graph_Collins_vs_xB_1->SetPointError(P_idx_col_1, 0.0, A_UT_err); // No x error
                P_idx_col_1++;
            } else if (x < 12){
                graph_Collins_vs_xB_2->SetPoint(P_idx_col_2, mean_xB, A_UT);
                graph_Collins_vs_xB_2->SetPointError(P_idx_col_2, 0.0, A_UT_err);
                P_idx_col_2++;
            } else if (x < 15){
                graph_Collins_vs_xB_3->SetPoint(P_idx_col_3, mean_xB, A_UT);
                graph_Collins_vs_xB_3->SetPointError(P_idx_col_3, 0.0, A_UT_err);
                P_idx_col_3++;
            } else if (x < 18){
                graph_Collins_vs_xB_4->SetPoint(P_idx_col_4, mean_xB, A_UT);
                graph_Collins_vs_xB_4->SetPointError(P_idx_col_4, 0.0, A_UT_err);
                P_idx_col_4++;
            }
        }

        graph_Collins_vs_xB->SetTitle(Form("Collins Asymmetry vs x_{B} (x_{B}-Q^{2} bin) for bin (%d, z-P_{hT}) ; x_{B}; A_{UT}^{sin(#Phi_{h}+#Phi_{s})}", z+1));
        graph_Collins_vs_xB_1->SetMarkerStyle(20), graph_Collins_vs_xB_2->SetMarkerStyle(20), graph_Collins_vs_xB_3->SetMarkerStyle(20), graph_Collins_vs_xB_4->SetMarkerStyle(20);;
        graph_Collins_vs_xB_1->SetLineColor(kAzure+5);
        graph_Collins_vs_xB_1->SetMarkerColor(kAzure+5);
        graph_Collins_vs_xB_2->SetLineColor(kViolet+5);
        graph_Collins_vs_xB_2->SetMarkerColor(kViolet+5);
        graph_Collins_vs_xB_3->SetLineColor(kPink+5);
        graph_Collins_vs_xB_3->SetMarkerColor(kPink+5);
        graph_Collins_vs_xB_4->SetLineColor(kOrange+5);
        graph_Collins_vs_xB_4->SetMarkerColor(kOrange+5);

        TCanvas* c_Collins_xB_2 = new TCanvas(Form("Aut_Collins_vs_xB_4d_zPt_bin%d", z+1), "Collins Asymmetry vs x_{B} for bin (n, z-P_{hT})", 800, 600);
        graph_Collins_vs_xB->Draw("A");
        graph_Collins_vs_xB_1->Draw("P SAME");
        graph_Collins_vs_xB_2->Draw("P SAME");
        graph_Collins_vs_xB_3->Draw("P SAME");
        graph_Collins_vs_xB_4->Draw("P SAME");
        TLine* guideLine_col2 = new TLine(graph_Collins_vs_xB->GetXaxis()->GetXmin(), 0, graph_Collins_vs_xB->GetXaxis()->GetXmax(), 0);
        guideLine_col2->SetLineStyle(2);  
        guideLine_col2->SetLineColor(kGray+1);
        guideLine_col2->Draw();

        TLegend* legend = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
        legend->AddEntry(graph_Collins_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
        legend->AddEntry(graph_Collins_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
        legend->AddEntry(graph_Collins_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
        legend->AddEntry(graph_Collins_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 10", "ep");
        legend->SetFillStyle(0);  // Transparent background
        legend->Draw();
        c_Collins_xB_2->Update();
        c_Collins_xB_2->Write();
    }






    // ________________________________________________________________________________________________________________________________________________________

    dir_Collins_ext_z->cd();
    TGraphErrors* graph_Collins_vs_z = new TGraphErrors();
    int p_idxz_col = 0;
    for (int z = 0; z < z_nBins; z++) {
        if (vec_pion_z[z].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_z = 0.0, sum_eps = 0;
        for (double val : vec_pion_z[z]) sum_z += val;
        double mean_z = sum_z / vec_pion_z[z].size();
        for (double valeps : vec_pion_epsilon_z[z]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_pion_epsilon_z[z].size();
        // Get A_UT and error from vectors
        double A_UT = A_UT_collins_z[z];
        double A_UT_err = A_UT_collins_err_z[z];
        if(A_UT > 1 || A_UT < -1) continue;
        // Fill the graph
        graph_Collins_vs_z->SetPoint(p_idxz_col, mean_z, A_UT);
        graph_Collins_vs_z->SetPointError(p_idxz_col, 0.0, A_UT_err); 
        p_idxz_col++;
    }
    graph_Collins_vs_z->SetTitle("Collins Asymmetry vs z (z-P_{hT} binning); z; A_{UT}^{sin(#Phi_{h}+#Phi_{s})}");
    graph_Collins_vs_z->SetMarkerStyle(20);
    graph_Collins_vs_z->SetLineColor(kAzure-7);
    graph_Collins_vs_z->SetMarkerColor(kAzure-7);

    TCanvas* c_Collins_z = new TCanvas("Aut_Collins_vs_z", "Collins Asymmetry vs z", 800, 600);
    TLine* guideLineZ_col = new TLine(graph_Collins_vs_z->GetXaxis()->GetXmin(), 0, graph_Collins_vs_z->GetXaxis()->GetXmax(), 0);
    guideLineZ_col->SetLineStyle(2);  
    guideLineZ_col->SetLineColor(kGray+1);
    graph_Collins_vs_z->Draw("AP");
    guideLineZ_col->Draw();
    c_Collins_z->Update();
    c_Collins_z->Write();

    TGraphErrors* graph_Collins_vs_z_1 = new TGraphErrors();
    TGraphErrors* graph_Collins_vs_z_2 = new TGraphErrors();
    TGraphErrors* graph_Collins_vs_z_3 = new TGraphErrors();
    TGraphErrors* graph_Collins_vs_z_4 = new TGraphErrors();
    TGraphErrors* graph_Collins_vs_z_5 = new TGraphErrors();
    int p_idxz_col_1 = 0, p_idxz_col_2 = 0, p_idxz_col_3 = 0, p_idxz_col_4 = 0, p_idxz_col_5 = 0;
    for (int z = 0; z < z_nBins; z++) {
        if (vec_pion_z[z].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_z = 0.0, sum_eps = 0;
        for (double val : vec_pion_z[z]) sum_z += val;
        double mean_z = sum_z / vec_pion_z[z].size();
        for (double valeps : vec_pion_epsilon_z[z]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_pion_epsilon_z[z].size();
        // Get A_UT and error from vectors
        double A_UT = A_UT_collins_z[z];
        double A_UT_err = A_UT_collins_err_z[z];
        if(A_UT > 1 || A_UT < -1) continue;
        // Fill the graph
        if(z < 7){
            graph_Collins_vs_z_1->SetPoint(p_idxz_col_1, mean_z, A_UT);
            graph_Collins_vs_z_1->SetPointError(p_idxz_col_1, 0.0, A_UT_err); // No x error
            p_idxz_col_1++;
        } else if (z < 14){
            graph_Collins_vs_z_2->SetPoint(p_idxz_col_2, mean_z, A_UT);
            graph_Collins_vs_z_2->SetPointError(p_idxz_col_2, 0.0, A_UT_err);
            p_idxz_col_2++;
        } else if (z < 20){
            graph_Collins_vs_z_3->SetPoint(p_idxz_col_3, mean_z, A_UT);
            graph_Collins_vs_z_3->SetPointError(p_idxz_col_3, 0.0, A_UT_err);
            p_idxz_col_3++;
        } else if (z < 26){
            graph_Collins_vs_z_4->SetPoint(p_idxz_col_4, mean_z, A_UT);
            graph_Collins_vs_z_4->SetPointError(p_idxz_col_4, 0.0, A_UT_err);
            p_idxz_col_4++;
        } else if (z < 30){
            graph_Collins_vs_z_5->SetPoint(p_idxz_col_5, mean_z, A_UT);
            graph_Collins_vs_z_5->SetPointError(p_idxz_col_5, 0.0, A_UT_err);
            p_idxz_col_5++;
        }
    }
    graph_Collins_vs_z->SetTitle("Collins Asymmetry vs z (z-P_{hT} bin) ; z; A_{UT}^{sin(#Phi_{h}+#Phi_{s})}");
    graph_Collins_vs_z_1->SetMarkerStyle(20), graph_Collins_vs_z_2->SetMarkerStyle(20), graph_Collins_vs_z_3->SetMarkerStyle(20), graph_Collins_vs_z_4->SetMarkerStyle(20), graph_Collins_vs_z_5->SetMarkerStyle(20);
    graph_Collins_vs_z_1->SetLineColor(kTeal+5);
    graph_Collins_vs_z_1->SetMarkerColor(kTeal+5);
    graph_Collins_vs_z_2->SetLineColor(kAzure+5);
    graph_Collins_vs_z_2->SetMarkerColor(kAzure+5);
    graph_Collins_vs_z_3->SetLineColor(kViolet+5);
    graph_Collins_vs_z_3->SetMarkerColor(kViolet+5);
    graph_Collins_vs_z_4->SetLineColor(kPink+5);
    graph_Collins_vs_z_4->SetMarkerColor(kPink+5);
    graph_Collins_vs_z_5->SetLineColor(kOrange+5);
    graph_Collins_vs_z_5->SetMarkerColor(kOrange+5);

    TCanvas* c_Collins_z_2 = new TCanvas("Aut_Collins_vs_z_sepPt", "Collins Asymmetry vs z for bin (n, x_{B}-Q^{2})", 800, 600);
    graph_Collins_vs_z->Draw("A");
    graph_Collins_vs_z_1->Draw("P SAME");
    graph_Collins_vs_z_2->Draw("P SAME");
    graph_Collins_vs_z_3->Draw("P SAME");
    graph_Collins_vs_z_4->Draw("P SAME");
    graph_Collins_vs_z_5->Draw("P SAME");
    guideLineZ_col->Draw();

    TLegend* legend_z_col = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
    legend_z_col->AddEntry(graph_Collins_vs_z_1, "0.0 < P_{hT} [GeV] < 0.2", "ep");
    legend_z_col->AddEntry(graph_Collins_vs_z_2, "0.2 < P_{hT} [GeV] < 0.4", "ep");
    legend_z_col->AddEntry(graph_Collins_vs_z_3, "0.4 < P_{hT} [GeV] < 0.6", "ep");
    legend_z_col->AddEntry(graph_Collins_vs_z_4, "0.6 < P_{hT} [GeV] < 0.8", "ep");
    legend_z_col->AddEntry(graph_Collins_vs_z_5, "0.8 < P_{hT} [GeV] < 1.2", "ep");
    legend_z_col->SetFillStyle(0);  // Transparent background
    legend_z_col->Draw();
    gPad->Update();
    c_Collins_z_2->Update();
    c_Collins_z_2->Write();

    gStyle->SetOptStat(0);
    TH1F* Collins_z_frame = new TH1F("", "", 1, graph_Collins_vs_z->GetXaxis()->GetXmin(), graph_Collins_vs_z->GetXaxis()->GetXmax());
    Collins_z_frame->SetMinimum(graph_Collins_vs_z->GetYaxis()->GetXmin());
    Collins_z_frame->SetMaximum(graph_Collins_vs_z->GetYaxis()->GetXmax());
    TCanvas* c_Collins_zBin_1 = new TCanvas("Aut_Collins_vs_z_Pt_1", "Collins Asymmetry vs z", 800, 600);
    Collins_z_frame->SetTitle("Collins Asymmetry vs z | pion+ | 0.0 < P_{hT} [GeV] < 0.2  ; z; A_{UT}^{sin(#Phi_{h}+#Phi_{s})}");
    Collins_z_frame->Draw();
    graph_Collins_vs_z_1->Draw("P SAME");
    guideLineZ_col->Draw();
    gPad->Update();
    c_Collins_zBin_1->Update();
    c_Collins_zBin_1->Write();
    TCanvas* c_Collins_zBin_2 = new TCanvas("Aut_Collins_vs_z_Pt_2", "Collins Asymmetry vs z", 800, 600);
    Collins_z_frame->SetTitle("Collins Asymmetry vs z | pion+ | 0.2 < P_{hT} [GeV] < 0.4  ; z; A_{UT}^{sin(#Phi_{h}+#Phi_{s})}");
    Collins_z_frame->Draw();
    graph_Collins_vs_z_2->Draw("P SAME");
    guideLineZ_col->Draw();
    gPad->Update();
    c_Collins_zBin_2->Update();
    c_Collins_zBin_2->Write();
    TCanvas* c_Collins_zBin_3 = new TCanvas("Aut_Collins_vs_z_Pt_3", "Collins Asymmetry vs z", 800, 600);
    Collins_z_frame->SetTitle("Collins Asymmetry vs z | pion+ | 0.4 < P_{hT} [GeV] < 0.6  ; z; A_{UT}^{sin(#Phi_{h}+#Phi_{s})}");
    Collins_z_frame->Draw();
    graph_Collins_vs_z_3->Draw("P SAME");
    guideLineZ_col->Draw();
    gPad->Update();
    c_Collins_zBin_3->Update();
    c_Collins_zBin_3->Write();
    TCanvas* c_Collins_zBin_4 = new TCanvas("Aut_Collins_vs_z_Pt_4", "Collins Asymmetry vs z", 800, 600);
    Collins_z_frame->SetTitle("Collins Asymmetry vs z | pion+ | 0.6 < P_{hT} [GeV] < 0.8  ; z; A_{UT}^{sin(#Phi_{h}+#Phi_{s})}");
    Collins_z_frame->Draw();
    graph_Collins_vs_z_4->Draw("P SAME");
    guideLineZ_col->Draw();
    gPad->Update();
    c_Collins_zBin_4->Update();
    c_Collins_zBin_4->Write();
    TCanvas* c_Collins_zBin_5 = new TCanvas("Aut_Collins_vs_z_Pt_5", "Collins Asymmetry vs z", 800, 600);
    Collins_z_frame->SetTitle("Collins Asymmetry vs z | pion+ | 0.8 < P_{hT} [GeV] < 1.2  ; z; A_{UT}^{sin(#Phi_{h}+#Phi_{s})}");
    Collins_z_frame->Draw();
    graph_Collins_vs_z_5->Draw("P SAME");
    guideLineZ_col->Draw();
    gPad->Update();
    c_Collins_zBin_5->Update();
    c_Collins_zBin_5->Write();

    // Pt color scale vs z 4D
    for (int x = 0; x < xB_nBins; x++) {
        TGraphErrors* graph_Collins_vs_z = new TGraphErrors();
        TGraphErrors* graph_Collins_vs_z_1 = new TGraphErrors();
        TGraphErrors* graph_Collins_vs_z_2 = new TGraphErrors();
        TGraphErrors* graph_Collins_vs_z_3 = new TGraphErrors();
        TGraphErrors* graph_Collins_vs_z_4 = new TGraphErrors();
        TGraphErrors* graph_Collins_vs_z_5 = new TGraphErrors();
        int p_idx = 0, p_idx_1 = 0, p_idx_2 = 0, p_idx_3 = 0, p_idx_4 = 0, p_idx_5 = 0;
        for(int z = 0; z < z_nBins; z++){
            //if(z == 0 || z == 7) continue; AAA
            if (vec_pion_xB_4d[x][z].empty()) continue; // Skip empty bins
            if (vec_pion_z_4d[x][z].empty()) continue;
            // Compute mean xB for this bin
            double sum_z = 0.0, sum_eps = 0;
            for (double val : vec_pion_z_4d[x][z]) sum_z += val;
            double mean_z = sum_z / vec_pion_z_4d[x][z].size();
            for (double valeps : vec_pion_epsilon_4d[x][z]) sum_eps += valeps;
            double mean_eps = sum_eps / vec_pion_epsilon_4d[x][z].size();
            // Get A_UT and error from vectors
            double A_UT = A_UT_collins_4d[x][z];
            double A_UT_err = A_UT_collins_err_4d[x][z];
            if(A_UT > 1 || A_UT < -1) continue;
            // Fill the graph
            graph_Collins_vs_z->SetPoint(p_idx, mean_z, A_UT);
            graph_Collins_vs_z->SetPointError(p_idx, 0.0, A_UT_err); // No x error
            p_idx++;
            if(z < 7){
                graph_Collins_vs_z_1->SetPoint(p_idx_1, mean_z, A_UT);
                graph_Collins_vs_z_1->SetPointError(p_idx_1, 0.0, A_UT_err); // No x error
                p_idx_1++;
            } else if (z < 14){
                graph_Collins_vs_z_2->SetPoint(p_idx_2, mean_z, A_UT);
                graph_Collins_vs_z_2->SetPointError(p_idx_2, 0.0, A_UT_err);
                p_idx_2++;
            } else if (z < 20){
                graph_Collins_vs_z_3->SetPoint(p_idx_3, mean_z, A_UT);
                graph_Collins_vs_z_3->SetPointError(p_idx_3, 0.0, A_UT_err);
                p_idx_3++;
            } else if (z < 26){
                graph_Collins_vs_z_4->SetPoint(p_idx_4, mean_z, A_UT);
                graph_Collins_vs_z_4->SetPointError(p_idx_4, 0.0, A_UT_err);
                p_idx_4++;
            } else if (z < 30){
                graph_Collins_vs_z_5->SetPoint(p_idx_5, mean_z, A_UT);
                graph_Collins_vs_z_5->SetPointError(p_idx_5, 0.0, A_UT_err);
                p_idx_5++;
            }
        }

        graph_Collins_vs_z->SetTitle(Form("Collins Asymmetry vs z (z-P_{hT} bin) for bin (%d, x_{B}-Q^{2}}) ; z; A_{UT}^{sin(#Phi_{h}+#Phi_{s})}", x+1));
        graph_Collins_vs_z_1->SetMarkerStyle(20), graph_Collins_vs_z_2->SetMarkerStyle(20), graph_Collins_vs_z_3->SetMarkerStyle(20), graph_Collins_vs_z_4->SetMarkerStyle(20), graph_Collins_vs_z_5->SetMarkerStyle(20);
        graph_Collins_vs_z_1->SetLineColor(kTeal+5);
        graph_Collins_vs_z_1->SetMarkerColor(kTeal+5);
        graph_Collins_vs_z_2->SetLineColor(kAzure+5);
        graph_Collins_vs_z_2->SetMarkerColor(kAzure+5);
        graph_Collins_vs_z_3->SetLineColor(kViolet+5);
        graph_Collins_vs_z_3->SetMarkerColor(kViolet+5);
        graph_Collins_vs_z_4->SetLineColor(kPink+5);
        graph_Collins_vs_z_4->SetMarkerColor(kPink+5);
        graph_Collins_vs_z_5->SetLineColor(kOrange+5);
        graph_Collins_vs_z_5->SetMarkerColor(kOrange+5);

        TCanvas* c_Collins_z_2 = new TCanvas(Form("Aut_Collins_vs_z_4d_xQ2_bin%d", x+1), "Collins Asymmetry vs z for bin (n, x_{B}-Q^{2})", 800, 600);
        graph_Collins_vs_z->Draw("A");
        graph_Collins_vs_z_1->Draw("P SAME");
        graph_Collins_vs_z_2->Draw("P SAME");
        graph_Collins_vs_z_3->Draw("P SAME");
        graph_Collins_vs_z_4->Draw("P SAME");
        graph_Collins_vs_z_5->Draw("P SAME");
        TLine* guideLine2 = new TLine(graph_Collins_vs_z->GetXaxis()->GetXmin(), 0, graph_Collins_vs_z->GetXaxis()->GetXmax(), 0);
        guideLine2->SetLineStyle(2);  
        guideLine2->SetLineColor(kGray+1);
        guideLine2->Draw();

        TLegend* legend = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
        legend->AddEntry(graph_Collins_vs_z_1, "0 < P_{hT} [GeV] < 0.2", "ep");
        legend->AddEntry(graph_Collins_vs_z_2, "0.2 < P_{hT} [GeV] < 0.4", "ep");
        legend->AddEntry(graph_Collins_vs_z_3, "0.4 < P_{hT} [GeV] < 0.6", "ep");
        legend->AddEntry(graph_Collins_vs_z_4, "0.6 < P_{hT} [GeV] < 0.8", "ep");
        legend->AddEntry(graph_Collins_vs_z_5, "0.8 < P_{hT} [GeV] < 1.2", "ep");
        legend->SetFillStyle(0);  // Transparent background
        legend->Draw();
        c_Collins_z_2->Update();
        c_Collins_z_2->Write();
    }

    // ________________________________________________________________________________________________________________________________________________________











    dir_Sivers->cd();
    double zMin = 0.00;
    double zMax = 0.04;
    //vector<vector<float>> A_UT_stat_err(xB_nBins, vector<float> (z_nBins, 0.0));
    for(int x = 0; x < xB_nBins; x++){
        for(int z = 0; z < z_nBins; z++){
            //if(z == 0 || z == 7) continue;
            float entries = hist_dist_xB_4d[x][z]->GetEntries();
            //A_UT_stat_err[x][z] = (1/0.85) * sqrt((1 - pow(0.85*0.1,2))/(1.25*entries));
            if (A_UT_sivers_err_4d[x][z] > 0.04) continue;
            if(z < 7) z_vs_Pt_xQ_bin[x]->SetBinContent(z+1, 1, A_UT_sivers_err_4d[x][z]); // +1 since bin number start from 1 and not 0 as our z counter
            else if(z < 14) z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-7, 2, A_UT_sivers_err_4d[x][z]);
            else if(z < 20){
                z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-14+1, 3, A_UT_sivers_err_4d[x][z]);
                if(z == 14) z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-14, 3, A_UT_sivers_err_4d[x][z]);
            } else if(z < 26){
                z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-20+1, 4, A_UT_sivers_err_4d[x][z]);
                if(z == 20) z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-20, 4, A_UT_sivers_err_4d[x][z]);
            }else if(z < 30){ 
                z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-25+1, 5, A_UT_sivers_err_4d[x][z]);
                if(z == 26){
                    z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-25-1, 5, A_UT_sivers_err_4d[x][z]);
                    z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-25, 5, A_UT_sivers_err_4d[x][z]);
                } else if(z == 29) z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-25+2, 5, A_UT_sivers_err_4d[x][z]);
            }
        }
        z_vs_Pt_xQ_bin[x]->SetMinimum(zMin);
        z_vs_Pt_xQ_bin[x]->SetMaximum(zMax);
        //z_vs_Pt_xQ_bin[x]->Write();

        TCanvas* c = new TCanvas(Form("Aut_Sivers_err_zPt_xQ_bin%d", x+1), Form("Canvas %d", x+1), 600, 600);
        //c->SetLogz();
        z_vs_Pt_xQ_bin[x]->GetZaxis()->SetNoExponent(true);       // Disattiva notazione scientifica
        z_vs_Pt_xQ_bin[x]->Draw("COLZ");
        std::vector<TPolyLine*> gridLines;
        for (int i = 0; i < 7; ++i) {
            for (int j = 0; j < 5; ++j) {
                double z_temp[5];
                // Se j <= 2, utilizziamo Clas_zBins_12, altrimenti Clas_zBins_34
                if (j <= 1) {
                    z_temp[0] = Clas_zBins_12[i][0];
                    z_temp[1] = Clas_zBins_12[i][1];
                    z_temp[2] = Clas_zBins_12[i][1];
                    z_temp[3] = Clas_zBins_12[i][0];
                    z_temp[4] = Clas_zBins_12[i][0];
                } else if(j <= 3) {
                    if(i >= 6) continue;
                    z_temp[0] = Clas_zBins_34[i][0];
                    z_temp[1] = Clas_zBins_34[i][1];
                    z_temp[2] = Clas_zBins_34[i][1];
                    z_temp[3] = Clas_zBins_34[i][0];
                    z_temp[4] = Clas_zBins_34[i][0];
                } else {
                    if(i >= 4) continue;
                    z_temp[0] = Clas_zBins_5[i][0];
                    z_temp[1] = Clas_zBins_5[i][1];
                    z_temp[2] = Clas_zBins_5[i][1];
                    z_temp[3] = Clas_zBins_5[i][0];
                    z_temp[4] = Clas_zBins_5[i][0];
                }
                double p_temp[] = {Clas_PtBins[j][0], Clas_PtBins[j][0], Clas_PtBins[j][1], Clas_PtBins[j][1], Clas_PtBins[j][0]};
                TPolyLine *rect = new TPolyLine(5, z_temp, p_temp);
                //rect->SetLineWidth(2);
                rect->Draw("same");
                gridLines.push_back(rect);
            }
        }
        c->Update();
        c->Write();
    }


    // FULL PLOT OF Z VS PT ERROR
    std::vector<std::string> Q2_labels;
    for (int q = 3; q >= 0; q--) {
        Q2_labels.push_back(Form("%.2f < Q^{2} [GeV^{2}] < %.2f", Clas_Q2Bins[q][0], Clas_Q2Bins[q][1]));
    }
    TCanvas* mainCanvas_siv = new TCanvas("FULL_Sivers_err_zPt_xQ", "Canvas con griglia 5x5", 1500, 1500);
    mainCanvas_siv->Divide(10, 4, 0.0, 0.0); 
    std::vector<int> padPositions = {31, 32, 33, 34, 35, 36, 37, 25, 26, 27, 28, 29, 18, 19, 20, 9, 10}; // Posizioni dei pad
    auto AddQ2Labels_siv = [&](TCanvas* canvas, int nCols, int nRows) {
        for (int row = 0; row < 4; ++row) {
            int padIndex_Pt = row * nCols + 1;  // First pad in the row
            canvas->cd(padIndex_Pt);
            TLatex latex;
            latex.SetNDC();
            //latex.SetTextAlign(21); // Left align vertically, bottom on rotated text
            latex.SetTextFont(2);
            latex.SetTextSize(0.06);
            latex.SetTextAngle(90); // Rotated upwards, positioned at right
            int binIndex = row;
            latex.DrawLatex(0.07, 0.3, Q2_labels[row].c_str());
        }
    };
    // Aggiungere i tuoi 17 grafici nelle rispettive posizioni
    for (int i = 0; i < z_vs_Pt_xQ_bin.size(); ++i) {
        if(i == 0){
            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(2);
            latex.SetTextSize(0.025);
            latex.DrawLatex(0.375, 0.975, "4D bin with a statistical error < 4% | pion+");
        }
        // Seleziona il pad corrispondente (cd() sceglie quale pad usare)
        mainCanvas_siv->cd(padPositions[i]); // i+1 perchè i pad sono indicizzati da 1 a 25
        //gPad->SetLogz();
        // Disegna il grafico nel pad corrente
        z_vs_Pt_xQ_bin[i]->SetStats(false);
        z_vs_Pt_xQ_bin[i]->SetTitle("");
        //if(padPositions[i] == 37 || padPositions[i] == 29 || padPositions[i] == 20 || padPositions[i] == 10) z_vs_Pt_xQ_bin[i]->Draw("COLZ");
        //else z_vs_Pt_xQ_bin[i]->Draw("hist");
        z_vs_Pt_xQ_bin[i]->Draw("hist");
        // Griglia personalizzata (come nel tuo codice)
        std::vector<TPolyLine*> gridLines;
        for (int i = 0; i < 7; ++i) {
            for (int j = 0; j < 5; ++j) {
                double z_temp[5];
                // Se j <= 2, utilizziamo Clas_zBins_12, altrimenti Clas_zBins_34
                if (j <= 1) {
                    z_temp[0] = Clas_zBins_12[i][0];
                    z_temp[1] = Clas_zBins_12[i][1];
                    z_temp[2] = Clas_zBins_12[i][1];
                    z_temp[3] = Clas_zBins_12[i][0];
                    z_temp[4] = Clas_zBins_12[i][0];
                } else if(j <= 3) {
                    if(i >= 6) continue;
                    z_temp[0] = Clas_zBins_34[i][0];
                    z_temp[1] = Clas_zBins_34[i][1];
                    z_temp[2] = Clas_zBins_34[i][1];
                    z_temp[3] = Clas_zBins_34[i][0];
                    z_temp[4] = Clas_zBins_34[i][0];
                } else {
                    if(i >= 4) continue;
                    z_temp[0] = Clas_zBins_5[i][0];
                    z_temp[1] = Clas_zBins_5[i][1];
                    z_temp[2] = Clas_zBins_5[i][1];
                    z_temp[3] = Clas_zBins_5[i][0];
                    z_temp[4] = Clas_zBins_5[i][0];
                }
                double p_temp[] = {Clas_PtBins[j][0], Clas_PtBins[j][0], Clas_PtBins[j][1], Clas_PtBins[j][1], Clas_PtBins[j][0]};
                TPolyLine *rect = new TPolyLine(5, z_temp, p_temp);
                //rect->SetLineWidth(2); // Se vuoi cambiare la larghezza della linea
                rect->Draw("same");
                gridLines.push_back(rect);
            }
        }
        mainCanvas_siv->cd(i + 1); // Ritorno al pad i+1
    }
    AddQ2Labels_siv(mainCanvas_siv, 10, 4);
    mainCanvas_siv->Update();
    mainCanvas_siv->Write();


    // ________________________________________________________________________________________________________________________________________________________








    dir_Collins->cd();
    //vector<vector<float>> A_UT_stat_err(xB_nBins, vector<float> (z_nBins, 0.0));
    for(int x = 0; x < xB_nBins; x++){
        for(int z = 0; z < z_nBins; z++){
            //if(z == 0 || z == 7) continue;
            float entries = hist_dist_xB_4d[x][z]->GetEntries();
            //A_UT_stat_err[x][z] = (1/0.85) * sqrt((1 - pow(0.85*0.1,2))/(1.25*entries));
            if (A_UT_collins_err_4d[x][z] > 0.04) continue;
            if(z < 7) z_vs_Pt_xQ_bin[x]->SetBinContent(z+1, 1, A_UT_collins_err_4d[x][z]); // +1 since bin number start from 1 and not 0 as our z counter
            else if(z < 14) z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-7, 2, A_UT_collins_err_4d[x][z]);
            else if(z < 20){
                z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-14+1, 3, A_UT_collins_err_4d[x][z]);
                if(z == 14) z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-14, 3, A_UT_collins_err_4d[x][z]);
            } else if(z < 26){
                z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-20+1, 4, A_UT_collins_err_4d[x][z]);
                if(z == 20) z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-20, 4, A_UT_collins_err_4d[x][z]);
            }else if(z < 30){ 
                z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-25+1, 5, A_UT_collins_err_4d[x][z]);
                if(z == 26){
                    z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-25-1, 5, A_UT_collins_err_4d[x][z]);
                    z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-25, 5, A_UT_collins_err_4d[x][z]);
                } else if(z == 29) z_vs_Pt_xQ_bin[x]->SetBinContent(z+1-25+2, 5, A_UT_collins_err_4d[x][z]);
            }
        }
        z_vs_Pt_xQ_bin[x]->SetMinimum(zMin);
        z_vs_Pt_xQ_bin[x]->SetMaximum(zMax);
        //z_vs_Pt_xQ_bin[x]->Write();

        TCanvas* c = new TCanvas(Form("Aut_Collins_err_zPt_xQ_bin%d", x+1), Form("Canvas %d", x+1), 600, 600);
        //c->SetLogz();
        z_vs_Pt_xQ_bin[x]->GetZaxis()->SetNoExponent(true);       // Disattiva notazione scientifica
        z_vs_Pt_xQ_bin[x]->Draw("COLZ");
        std::vector<TPolyLine*> gridLines;
        for (int i = 0; i < 7; ++i) {
            for (int j = 0; j < 5; ++j) {
                double z_temp[5];
                // Se j <= 2, utilizziamo Clas_zBins_12, altrimenti Clas_zBins_34
                if (j <= 1) {
                    z_temp[0] = Clas_zBins_12[i][0];
                    z_temp[1] = Clas_zBins_12[i][1];
                    z_temp[2] = Clas_zBins_12[i][1];
                    z_temp[3] = Clas_zBins_12[i][0];
                    z_temp[4] = Clas_zBins_12[i][0];
                } else if(j <= 3) {
                    if(i >= 6) continue;
                    z_temp[0] = Clas_zBins_34[i][0];
                    z_temp[1] = Clas_zBins_34[i][1];
                    z_temp[2] = Clas_zBins_34[i][1];
                    z_temp[3] = Clas_zBins_34[i][0];
                    z_temp[4] = Clas_zBins_34[i][0];
                } else {
                    if(i >= 4) continue;
                    z_temp[0] = Clas_zBins_5[i][0];
                    z_temp[1] = Clas_zBins_5[i][1];
                    z_temp[2] = Clas_zBins_5[i][1];
                    z_temp[3] = Clas_zBins_5[i][0];
                    z_temp[4] = Clas_zBins_5[i][0];
                }
                double p_temp[] = {Clas_PtBins[j][0], Clas_PtBins[j][0], Clas_PtBins[j][1], Clas_PtBins[j][1], Clas_PtBins[j][0]};
                TPolyLine *rect = new TPolyLine(5, z_temp, p_temp);
                //rect->SetLineWidth(2);
                rect->Draw("same");
                gridLines.push_back(rect);
            }
        }
        c->Update();
        c->Write();
    }


    // FULL PLOT OF Z VS PT ERROR
    TCanvas* mainCanvas_col = new TCanvas("FULL_Collins_err_zPt_xQ", "Canvas con griglia 5x5", 1500, 1500);
    mainCanvas_col->Divide(10, 4, 0.0, 0.0); 
    std::vector<int> padPositions_col = {31, 32, 33, 34, 35, 36, 37, 25, 26, 27, 28, 29, 18, 19, 20, 9, 10}; // Posizioni dei pad
    auto AddQ2Labels_col = [&](TCanvas* canvas, int nCols, int nRows) {
        for (int row = 0; row < 4; ++row) {
            int padIndex_Pt = row * nCols + 1;  // First pad in the row
            canvas->cd(padIndex_Pt);
            TLatex latex;
            latex.SetNDC();
            //latex.SetTextAlign(21); // Left align vertically, bottom on rotated text
            latex.SetTextFont(2);
            latex.SetTextSize(0.06);
            latex.SetTextAngle(90); // Rotated upwards, positioned at right
            int binIndex = row;
            latex.DrawLatex(0.07, 0.3, Q2_labels[row].c_str());
        }
    };
    // Aggiungere i tuoi 17 grafici nelle rispettive posizioni
    for (int i = 0; i < z_vs_Pt_xQ_bin.size(); ++i) {
        if(i == 0){
            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(2);
            latex.SetTextSize(0.025);
            latex.DrawLatex(0.375, 0.975, "4D bin with a statistical error < 4% | pion+");
        }
        // Seleziona il pad corrispondente (cd() sceglie quale pad usare)
        mainCanvas_col->cd(padPositions_col[i]); // i+1 perchè i pad sono indicizzati da 1 a 25
        //gPad->SetLogz();
        // Disegna il grafico nel pad corrente
        z_vs_Pt_xQ_bin[i]->SetStats(false);
        z_vs_Pt_xQ_bin[i]->SetTitle("");
        //if(padPositions_col[i] == 37 || padPositions_col[i] == 29 || padPositions_col[i] == 20 || padPositions_col[i] == 10) z_vs_Pt_xQ_bin[i]->Draw("COLZ");
        //else z_vs_Pt_xQ_bin[i]->Draw("hist");
        z_vs_Pt_xQ_bin[i]->Draw("hist");
        // Griglia personalizzata (come nel tuo codice)
        std::vector<TPolyLine*> gridLines;
        for (int i = 0; i < 7; ++i) {
            for (int j = 0; j < 5; ++j) {
                double z_temp[5];
                // Se j <= 2, utilizziamo Clas_zBins_12, altrimenti Clas_zBins_34
                if (j <= 1) {
                    z_temp[0] = Clas_zBins_12[i][0];
                    z_temp[1] = Clas_zBins_12[i][1];
                    z_temp[2] = Clas_zBins_12[i][1];
                    z_temp[3] = Clas_zBins_12[i][0];
                    z_temp[4] = Clas_zBins_12[i][0];
                } else if(j <= 3) {
                    if(i >= 6) continue;
                    z_temp[0] = Clas_zBins_34[i][0];
                    z_temp[1] = Clas_zBins_34[i][1];
                    z_temp[2] = Clas_zBins_34[i][1];
                    z_temp[3] = Clas_zBins_34[i][0];
                    z_temp[4] = Clas_zBins_34[i][0];
                } else {
                    if(i >= 4) continue;
                    z_temp[0] = Clas_zBins_5[i][0];
                    z_temp[1] = Clas_zBins_5[i][1];
                    z_temp[2] = Clas_zBins_5[i][1];
                    z_temp[3] = Clas_zBins_5[i][0];
                    z_temp[4] = Clas_zBins_5[i][0];
                }
                double p_temp[] = {Clas_PtBins[j][0], Clas_PtBins[j][0], Clas_PtBins[j][1], Clas_PtBins[j][1], Clas_PtBins[j][0]};
                TPolyLine *rect = new TPolyLine(5, z_temp, p_temp);
                //rect->SetLineWidth(2); // Se vuoi cambiare la larghezza della linea
                rect->Draw("same");
                gridLines.push_back(rect);
            }
        }
        mainCanvas_col->cd(i + 1); // Ritorno al pad i+1
    }
    AddQ2Labels_col(mainCanvas_col, 10, 4);
    mainCanvas_col->Update();
    mainCanvas_col->Write();



    // ________________________________________________________________________________________________________________________________________________________



    dir_Pretz_ext_xB->cd();
    TGraphErrors* graph_Pretz_vs_xB = new TGraphErrors();
    int p_idp = 0;
    for (int x = 0; x < xB_nBins; x++) {
        if (vec_pion_xB[x].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_xB = 0.0;
        for (double val : vec_pion_xB[x]) sum_xB += val;
        double mean_xB = sum_xB / vec_pion_xB[x].size();
        // Get A_UT and error from vectors
        double A_UT = A_UT_pretz[x];
        double A_UT_err = A_UT_pretz_err[x];
        if(A_UT > 1 || A_UT < -1) continue;
        // Fill the graph
        graph_Pretz_vs_xB->SetPoint(p_idp, mean_xB, A_UT);
        graph_Pretz_vs_xB->SetPointError(p_idp, 0.0, A_UT_err); // No x error
        p_idp++;
    }
    graph_Pretz_vs_xB->SetTitle("Pretzelosity Asymmetry vs x_{B} (x_{B}-Q^{2} binning); x_{B}; A_{UT}^{sin(3#Phi_{h}-#Phi_{s})}");
    graph_Pretz_vs_xB->SetMarkerStyle(20);
    graph_Pretz_vs_xB->SetLineColor(kAzure-7);
    graph_Pretz_vs_xB->SetMarkerColor(kAzure-7);

    TCanvas* c_Pretz_xB = new TCanvas("Aut_Pretz_vs_xB", "Pretzelosity Asymmetry vs x_{B}", 800, 600);
    TLine* guideLine_pretz = new TLine(graph_Pretz_vs_xB->GetXaxis()->GetXmin(), 0, graph_Pretz_vs_xB->GetXaxis()->GetXmax(), 0);
    guideLine_pretz->SetLineStyle(2);  
    guideLine_pretz->SetLineColor(kGray+1);
    graph_Pretz_vs_xB->Draw("AP");
    guideLine_pretz->Draw();
    c_Pretz_xB->Update();
    c_Pretz_xB->Write();


    // Q2 color scale
    TGraphErrors* graph_Pretz_vs_xB_1 = new TGraphErrors();
    TGraphErrors* graph_Pretz_vs_xB_2 = new TGraphErrors();
    TGraphErrors* graph_Pretz_vs_xB_3 = new TGraphErrors();
    TGraphErrors* graph_Pretz_vs_xB_4 = new TGraphErrors();
    int p_idp_1 = 0, p_idp_2 = 0, p_idp_3 = 0, p_idp_4 = 0;
    for (int x = 0; x < xB_nBins; x++) {
        if (vec_pion_xB[x].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_xB = 0.0;
        for (double val : vec_pion_xB[x]) sum_xB += val;
        double mean_xB = sum_xB / vec_pion_xB[x].size();
        // Get A_UT and error from vectors
        double A_UT = A_UT_pretz[x];
        double A_UT_err = A_UT_pretz_err[x];
        if(A_UT > 1 || A_UT < -1) continue;
        // Fill the graph
        if(x < 7){
            graph_Pretz_vs_xB_1->SetPoint(p_idp_1, mean_xB, A_UT);
            graph_Pretz_vs_xB_1->SetPointError(p_idp_1, 0.0, A_UT_err); // No x error
            p_idp_1++;
        } else if (x < 12){
            graph_Pretz_vs_xB_2->SetPoint(p_idp_2, mean_xB, A_UT);
            graph_Pretz_vs_xB_2->SetPointError(p_idp_2, 0.0, A_UT_err);
            p_idp_2++;
        } else if (x < 15){
            graph_Pretz_vs_xB_3->SetPoint(p_idp_3, mean_xB, A_UT);
            graph_Pretz_vs_xB_3->SetPointError(p_idp_3, 0.0, A_UT_err);
            p_idp_3++;
        } else if (x < 18){
            graph_Pretz_vs_xB_4->SetPoint(p_idp_4, mean_xB, A_UT);
            graph_Pretz_vs_xB_4->SetPointError(p_idp_4, 0.0, A_UT_err);
            p_idp_4++;
        }
    }
    graph_Pretz_vs_xB_1->SetTitle("Pretzelosity Asymmetry vs x_{B} (x_{B}-Q^{2} binning); x_{B}; A_{UT}^{sin(3#Phi_{h}-#Phi_{s})}");
    graph_Pretz_vs_xB_1->SetMarkerStyle(20), graph_Pretz_vs_xB_2->SetMarkerStyle(20), graph_Pretz_vs_xB_3->SetMarkerStyle(20), graph_Pretz_vs_xB_4->SetMarkerStyle(20);;
    graph_Pretz_vs_xB_1->SetLineColor(kAzure+5);
    graph_Pretz_vs_xB_1->SetMarkerColor(kAzure+5);
    graph_Pretz_vs_xB_2->SetLineColor(kViolet+5);
    graph_Pretz_vs_xB_2->SetMarkerColor(kViolet+5);
    graph_Pretz_vs_xB_3->SetLineColor(kPink+5);
    graph_Pretz_vs_xB_3->SetMarkerColor(kPink+5);
    graph_Pretz_vs_xB_4->SetLineColor(kOrange+5);
    graph_Pretz_vs_xB_4->SetMarkerColor(kOrange+5);

    TCanvas* c_Pretz_xB_2 = new TCanvas("Aut_Pretz_vs_xB_sepQ2", "Pretzelosity Asymmetry vs x_{B}", 800, 600);
    graph_Pretz_vs_xB->Draw("A");
    graph_Pretz_vs_xB_1->Draw("P SAME");
    graph_Pretz_vs_xB_2->Draw("P SAME");
    graph_Pretz_vs_xB_3->Draw("P SAME");
    graph_Pretz_vs_xB_4->Draw("P SAME");
    guideLine_pretz->Draw();

    TLegend* legend_pretz = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
    legend_pretz->AddEntry(graph_Pretz_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
    legend_pretz->AddEntry(graph_Pretz_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
    legend_pretz->AddEntry(graph_Pretz_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
    legend_pretz->AddEntry(graph_Pretz_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 10", "ep");
    legend_pretz->SetFillStyle(0);  // Transparent background
    legend_pretz->Draw();
    c_Pretz_xB_2->Update();
    c_Pretz_xB_2->Write();


    // Draw the plot of AUT

    // Q2 color scale vs xB 4D
    for (int z = 0; z < z_nBins; z++) {
        //if(z == 0 || z == 7) continue;
        TGraphErrors* graph_Pretz_vs_xB = new TGraphErrors();
        TGraphErrors* graph_Pretz_vs_xB_1 = new TGraphErrors();
        TGraphErrors* graph_Pretz_vs_xB_2 = new TGraphErrors();
        TGraphErrors* graph_Pretz_vs_xB_3 = new TGraphErrors();
        TGraphErrors* graph_Pretz_vs_xB_4 = new TGraphErrors();
        int p_idp = 0, p_idp_1 = 0, p_idp_2 = 0, p_idp_3 = 0, p_idp_4 = 0;
        for(int x = 0; x < xB_nBins; x++){
            if (vec_pion_xB_4d[x][z].empty()) continue; // Skip empty bins
            if (vec_pion_z_4d[x][z].empty()) continue;
            // Compute mean xB for this bin
            double sum_xB = 0.0;
            for (double val : vec_pion_xB_4d[x][z]) sum_xB += val;
            double mean_xB = sum_xB / vec_pion_xB_4d[x][z].size();
            // Get A_UT and error from vectors
            double A_UT = A_UT_pretz_4d[x][z];
            double A_UT_err = A_UT_pretz_err_4d[x][z];
            if(A_UT > 1 || A_UT < -1) continue;
            // Fill the graph
            graph_Pretz_vs_xB->SetPoint(p_idp, mean_xB, A_UT);
            graph_Pretz_vs_xB->SetPointError(p_idp, 0.0, A_UT_err); // No x error
            p_idp++;
            if(x < 7){
                graph_Pretz_vs_xB_1->SetPoint(p_idp_1, mean_xB, A_UT);
                graph_Pretz_vs_xB_1->SetPointError(p_idp_1, 0.0, A_UT_err); // No x error
                p_idp_1++;
            } else if (x < 12){
                graph_Pretz_vs_xB_2->SetPoint(p_idp_2, mean_xB, A_UT);
                graph_Pretz_vs_xB_2->SetPointError(p_idp_2, 0.0, A_UT_err);
                p_idp_2++;
            } else if (x < 15){
                graph_Pretz_vs_xB_3->SetPoint(p_idp_3, mean_xB, A_UT);
                graph_Pretz_vs_xB_3->SetPointError(p_idp_3, 0.0, A_UT_err);
                p_idp_3++;
            } else if (x < 18){
                graph_Pretz_vs_xB_4->SetPoint(p_idp_4, mean_xB, A_UT);
                graph_Pretz_vs_xB_4->SetPointError(p_idp_4, 0.0, A_UT_err);
                p_idp_4++;
            }
        }

        graph_Pretz_vs_xB->SetTitle(Form("Pretzelosity Asymmetry vs x_{B} (x_{B}-Q^{2} bin) for bin (%d, z-P_{hT}) ; x_{B}; A_{UT}^{sin(3#Phi_{h}-#Phi_{s})}", z+1));
        graph_Pretz_vs_xB_1->SetMarkerStyle(20), graph_Pretz_vs_xB_2->SetMarkerStyle(20), graph_Pretz_vs_xB_3->SetMarkerStyle(20), graph_Pretz_vs_xB_4->SetMarkerStyle(20);;
        graph_Pretz_vs_xB_1->SetLineColor(kAzure+5);
        graph_Pretz_vs_xB_1->SetMarkerColor(kAzure+5);
        graph_Pretz_vs_xB_2->SetLineColor(kViolet+5);
        graph_Pretz_vs_xB_2->SetMarkerColor(kViolet+5);
        graph_Pretz_vs_xB_3->SetLineColor(kPink+5);
        graph_Pretz_vs_xB_3->SetMarkerColor(kPink+5);
        graph_Pretz_vs_xB_4->SetLineColor(kOrange+5);
        graph_Pretz_vs_xB_4->SetMarkerColor(kOrange+5);

        TCanvas* c_Pretz_xB_2 = new TCanvas(Form("Aut_Pretz_vs_xB_4d_zPt_bin%d", z+1), "Pretzelosity Asymmetry vs x_{B} for bin (n, z-P_{hT})", 800, 600);
        graph_Pretz_vs_xB->Draw("A");
        graph_Pretz_vs_xB_1->Draw("P SAME");
        graph_Pretz_vs_xB_2->Draw("P SAME");
        graph_Pretz_vs_xB_3->Draw("P SAME");
        graph_Pretz_vs_xB_4->Draw("P SAME");
        TLine* guideLine_pretz2 = new TLine(graph_Pretz_vs_xB->GetXaxis()->GetXmin(), 0, graph_Pretz_vs_xB->GetXaxis()->GetXmax(), 0);
        guideLine_pretz2->SetLineStyle(2);  
        guideLine_pretz2->SetLineColor(kGray+1);
        guideLine_pretz2->Draw();

        TLegend* legend_pretz2 = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
        legend_pretz2->AddEntry(graph_Pretz_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
        legend_pretz2->AddEntry(graph_Pretz_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
        legend_pretz2->AddEntry(graph_Pretz_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
        legend_pretz2->AddEntry(graph_Pretz_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 10", "ep");
        legend_pretz2->SetFillStyle(0);  // Transparent background
        legend_pretz2->Draw();
        c_Pretz_xB_2->Update();
        c_Pretz_xB_2->Write();
    }


 // ______________________________________________________________________________________________________________________________________



 dir_Pretz_ext_z->cd();
    TGraphErrors* graph_Pretz_vs_z = new TGraphErrors();
    int p_idpz = 0;
    for (int z = 0; z < z_nBins; z++) {
        if (vec_pion_z[z].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_z = 0.0;
        for (double val : vec_pion_z[z]) sum_z += val;
        double mean_z = sum_z / vec_pion_z[z].size();
        // Get A_UT and error from vectors
        double A_UT = A_UT_pretz_z[z];
        double A_UT_err = A_UT_pretz_err_z[z];
        if(A_UT > 1 || A_UT < -1) continue;
        // Fill the graph
        graph_Pretz_vs_z->SetPoint(p_idpz, mean_z, A_UT);
        graph_Pretz_vs_z->SetPointError(p_idpz, 0.0, A_UT_err); 
        p_idpz++;
    }
    graph_Pretz_vs_z->SetTitle("Pretzelosity Asymmetry vs z (z-P_{hT} binning); z; A_{UT}^{sin(3#Phi_{h}-#Phi_{s})}");
    graph_Pretz_vs_z->SetMarkerStyle(20);
    graph_Pretz_vs_z->SetLineColor(kAzure-7);
    graph_Pretz_vs_z->SetMarkerColor(kAzure-7);

    TCanvas* c_Pretz_z = new TCanvas("Aut_Pretz_vs_z", "Pretzelosity Asymmetry vs z", 800, 600);
    TLine* guideLineZ_pretz = new TLine(graph_Pretz_vs_z->GetXaxis()->GetXmin(), 0, graph_Pretz_vs_z->GetXaxis()->GetXmax(), 0);
    guideLineZ_pretz->SetLineStyle(2);  
    guideLineZ_pretz->SetLineColor(kGray+1);
    graph_Pretz_vs_z->Draw("AP");
    guideLineZ_pretz->Draw();
    c_Pretz_z->Update();
    c_Pretz_z->Write();

    TGraphErrors* graph_Pretz_vs_z_1 = new TGraphErrors();
    TGraphErrors* graph_Pretz_vs_z_2 = new TGraphErrors();
    TGraphErrors* graph_Pretz_vs_z_3 = new TGraphErrors();
    TGraphErrors* graph_Pretz_vs_z_4 = new TGraphErrors();
    TGraphErrors* graph_Pretz_vs_z_5 = new TGraphErrors();
    int p_idpz_1 = 0, p_idpz_2 = 0, p_idpz_3 = 0, p_idpz_4 = 0, p_idpz_5 = 0;
    for (int z = 0; z < z_nBins; z++) {
        if (vec_pion_z[z].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_z = 0.0;
        for (double val : vec_pion_z[z]) sum_z += val;
        double mean_z = sum_z / vec_pion_z[z].size();
        // Get A_UT and error from vectors
        double A_UT = A_UT_pretz_z[z];
        double A_UT_err = A_UT_pretz_err_z[z];
        if(A_UT > 1 || A_UT < -1) continue;
        // Fill the graph
        if(z < 7){
            graph_Pretz_vs_z_1->SetPoint(p_idpz_1, mean_z, A_UT);
            graph_Pretz_vs_z_1->SetPointError(p_idpz_1, 0.0, A_UT_err); // No x error
            p_idpz_1++;
        } else if (z < 14){
            graph_Pretz_vs_z_2->SetPoint(p_idpz_2, mean_z, A_UT);
            graph_Pretz_vs_z_2->SetPointError(p_idpz_2, 0.0, A_UT_err);
            p_idpz_2++;
        } else if (z < 20){
            graph_Pretz_vs_z_3->SetPoint(p_idpz_3, mean_z, A_UT);
            graph_Pretz_vs_z_3->SetPointError(p_idpz_3, 0.0, A_UT_err);
            p_idpz_3++;
        } else if (z < 26){
            graph_Pretz_vs_z_4->SetPoint(p_idpz_4, mean_z, A_UT);
            graph_Pretz_vs_z_4->SetPointError(p_idpz_4, 0.0, A_UT_err);
            p_idpz_4++;
        } else if (z < 30){
            graph_Pretz_vs_z_5->SetPoint(p_idpz_5, mean_z, A_UT);
            graph_Pretz_vs_z_5->SetPointError(p_idpz_5, 0.0, A_UT_err);
            p_idpz_5++;
        }
    }
    graph_Pretz_vs_z->SetTitle("Pretzelosity Asymmetry vs z (z-P_{hT} bin) ; z; A_{UT}^{sin(3#Phi_{h}-#Phi_{s})}");
    graph_Pretz_vs_z_1->SetMarkerStyle(20), graph_Pretz_vs_z_2->SetMarkerStyle(20), graph_Pretz_vs_z_3->SetMarkerStyle(20), graph_Pretz_vs_z_4->SetMarkerStyle(20), graph_Pretz_vs_z_5->SetMarkerStyle(20);
    graph_Pretz_vs_z_1->SetLineColor(kTeal+5);
    graph_Pretz_vs_z_1->SetMarkerColor(kTeal+5);
    graph_Pretz_vs_z_2->SetLineColor(kAzure+5);
    graph_Pretz_vs_z_2->SetMarkerColor(kAzure+5);
    graph_Pretz_vs_z_3->SetLineColor(kViolet+5);
    graph_Pretz_vs_z_3->SetMarkerColor(kViolet+5);
    graph_Pretz_vs_z_4->SetLineColor(kPink+5);
    graph_Pretz_vs_z_4->SetMarkerColor(kPink+5);
    graph_Pretz_vs_z_5->SetLineColor(kOrange+5);
    graph_Pretz_vs_z_5->SetMarkerColor(kOrange+5);

    TCanvas* c_Pretz_z_2 = new TCanvas("Aut_Pretz_vs_z_sepPt", "Pretzelosity Asymmetry vs z for bin (n, x_{B}-Q^{2})", 800, 600);
    graph_Pretz_vs_z->Draw("A");
    graph_Pretz_vs_z_1->Draw("P SAME");
    graph_Pretz_vs_z_2->Draw("P SAME");
    graph_Pretz_vs_z_3->Draw("P SAME");
    graph_Pretz_vs_z_4->Draw("P SAME");
    graph_Pretz_vs_z_5->Draw("P SAME");
    guideLineZ_pretz->Draw();

    TLegend* legend_z_pretz = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
    legend_z_pretz->AddEntry(graph_Pretz_vs_z_1, "0.0 < P_{hT} [GeV] < 0.2", "ep");
    legend_z_pretz->AddEntry(graph_Pretz_vs_z_2, "0.2 < P_{hT} [GeV] < 0.4", "ep");
    legend_z_pretz->AddEntry(graph_Pretz_vs_z_3, "0.4 < P_{hT} [GeV] < 0.6", "ep");
    legend_z_pretz->AddEntry(graph_Pretz_vs_z_4, "0.6 < P_{hT} [GeV] < 0.8", "ep");
    legend_z_pretz->AddEntry(graph_Pretz_vs_z_5, "0.8 < P_{hT} [GeV] < 1.2", "ep");
    legend_z_pretz->SetFillStyle(0);  // Transparent background
    legend_z_pretz->Draw();
    gPad->Update();
    c_Pretz_z_2->Update();
    c_Pretz_z_2->Write();

    gStyle->SetOptStat(0);
    TH1F* Pretz_z_frame = new TH1F("", "", 1, graph_Pretz_vs_z->GetXaxis()->GetXmin(), graph_Pretz_vs_z->GetXaxis()->GetXmax());
    Pretz_z_frame->SetMinimum(graph_Pretz_vs_z->GetYaxis()->GetXmin());
    Pretz_z_frame->SetMaximum(graph_Pretz_vs_z->GetYaxis()->GetXmax());
    TCanvas* c_Pretz_zBin_1 = new TCanvas("Aut_Pretz_vs_z_Pt_1", "Pretzelosity Asymmetry vs z", 800, 600);
    Pretz_z_frame->SetTitle("Pretzelosity Asymmetry vs z | pion+ | 0.0 < P_{hT} [GeV] < 0.2  ; z; A_{UT}^{sin(3#Phi_{h}-#Phi_{s})}");
    Pretz_z_frame->Draw();
    graph_Pretz_vs_z_1->Draw("P SAME");
    guideLineZ_pretz->Draw();
    gPad->Update();
    c_Pretz_zBin_1->Update();
    c_Pretz_zBin_1->Write();
    TCanvas* c_Pretz_zBin_2 = new TCanvas("Aut_Pretz_vs_z_Pt_2", "Pretzelosity Asymmetry vs z", 800, 600);
    Pretz_z_frame->SetTitle("Pretzelosity Asymmetry vs z | pion+ | 0.2 < P_{hT} [GeV] < 0.4  ; z; A_{UT}^{sin(3#Phi_{h}-#Phi_{s})}");
    Pretz_z_frame->Draw();
    graph_Pretz_vs_z_2->Draw("P SAME");
    guideLineZ_pretz->Draw();
    gPad->Update();
    c_Pretz_zBin_2->Update();
    c_Pretz_zBin_2->Write();
    TCanvas* c_Pretz_zBin_3 = new TCanvas("Aut_Pretz_vs_z_Pt_3", "Pretzelosity Asymmetry vs z", 800, 600);
    Pretz_z_frame->SetTitle("Pretzelosity Asymmetry vs z | pion+ | 0.4 < P_{hT} [GeV] < 0.6  ; z; A_{UT}^{sin(3#Phi_{h}-#Phi_{s})}");
    Pretz_z_frame->Draw();
    graph_Pretz_vs_z_3->Draw("P SAME");
    guideLineZ_pretz->Draw();
    gPad->Update();
    c_Pretz_zBin_3->Update();
    c_Pretz_zBin_3->Write();
    TCanvas* c_Pretz_zBin_4 = new TCanvas("Aut_Pretz_vs_z_Pt_4", "Pretzelosity Asymmetry vs z", 800, 600);
    Pretz_z_frame->SetTitle("Pretzelosity Asymmetry vs z | pion+ | 0.6 < P_{hT} [GeV] < 0.8  ; z; A_{UT}^{sin(3#Phi_{h}-#Phi_{s})}");
    Pretz_z_frame->Draw();
    graph_Pretz_vs_z_4->Draw("P SAME");
    guideLineZ_pretz->Draw();
    gPad->Update();
    c_Pretz_zBin_4->Update();
    c_Pretz_zBin_4->Write();
    TCanvas* c_Pretz_zBin_5 = new TCanvas("Aut_Pretz_vs_z_Pt_5", "Pretzelosity Asymmetry vs z", 800, 600);
    Pretz_z_frame->SetTitle("Pretzelosity Asymmetry vs z | pion+ | 0.8 < P_{hT} [GeV] < 1.2  ; z; A_{UT}^{sin(3#Phi_{h}-#Phi_{s})}");
    Pretz_z_frame->Draw();
    graph_Pretz_vs_z_5->Draw("P SAME");
    guideLineZ_pretz->Draw();
    gPad->Update();
    c_Pretz_zBin_5->Update();
    c_Pretz_zBin_5->Write();

    // Pt color scale vs z 4D
    for (int x = 0; x < xB_nBins; x++) {
        TGraphErrors* graph_Pretz_vs_z = new TGraphErrors();
        TGraphErrors* graph_Pretz_vs_z_1 = new TGraphErrors();
        TGraphErrors* graph_Pretz_vs_z_2 = new TGraphErrors();
        TGraphErrors* graph_Pretz_vs_z_3 = new TGraphErrors();
        TGraphErrors* graph_Pretz_vs_z_4 = new TGraphErrors();
        TGraphErrors* graph_Pretz_vs_z_5 = new TGraphErrors();
        int p_idx = 0, p_idx_1 = 0, p_idx_2 = 0, p_idx_3 = 0, p_idx_4 = 0, p_idx_5 = 0;
        for(int z = 0; z < z_nBins; z++){
            //if(z == 0 || z == 7) continue; AAA
            if (vec_pion_xB_4d[x][z].empty()) continue; // Skip empty bins
            if (vec_pion_z_4d[x][z].empty()) continue;
            // Compute mean xB for this bin
            double sum_z = 0.0;
            for (double val : vec_pion_z_4d[x][z]) sum_z += val;
            double mean_z = sum_z / vec_pion_z_4d[x][z].size();
            // Get A_UT and error from vectors
            double A_UT = A_UT_pretz_4d[x][z];
            double A_UT_err = A_UT_pretz_err_4d[x][z];
            if(A_UT > 1 || A_UT < -1) continue;
            // Fill the graph
            graph_Pretz_vs_z->SetPoint(p_idx, mean_z, A_UT);
            graph_Pretz_vs_z->SetPointError(p_idx, 0.0, A_UT_err); // No x error
            p_idx++;
            if(z < 7){
                graph_Pretz_vs_z_1->SetPoint(p_idx_1, mean_z, A_UT);
                graph_Pretz_vs_z_1->SetPointError(p_idx_1, 0.0, A_UT_err); // No x error
                p_idx_1++;
            } else if (z < 14){
                graph_Pretz_vs_z_2->SetPoint(p_idx_2, mean_z, A_UT);
                graph_Pretz_vs_z_2->SetPointError(p_idx_2, 0.0, A_UT_err);
                p_idx_2++;
            } else if (z < 20){
                graph_Pretz_vs_z_3->SetPoint(p_idx_3, mean_z, A_UT);
                graph_Pretz_vs_z_3->SetPointError(p_idx_3, 0.0, A_UT_err);
                p_idx_3++;
            } else if (z < 26){
                graph_Pretz_vs_z_4->SetPoint(p_idx_4, mean_z, A_UT);
                graph_Pretz_vs_z_4->SetPointError(p_idx_4, 0.0, A_UT_err);
                p_idx_4++;
            } else if (z < 30){
                graph_Pretz_vs_z_5->SetPoint(p_idx_5, mean_z, A_UT);
                graph_Pretz_vs_z_5->SetPointError(p_idx_5, 0.0, A_UT_err);
                p_idx_5++;
            }
        }

        graph_Pretz_vs_z->SetTitle(Form("Pretzelosity Asymmetry vs z (z-P_{hT} bin) for bin (%d, x_{B}-Q^{2}}) ; z; A_{UT}^{sin(3#Phi_{h}-#Phi_{s})}", x+1));
        graph_Pretz_vs_z_1->SetMarkerStyle(20), graph_Pretz_vs_z_2->SetMarkerStyle(20), graph_Pretz_vs_z_3->SetMarkerStyle(20), graph_Pretz_vs_z_4->SetMarkerStyle(20), graph_Pretz_vs_z_5->SetMarkerStyle(20);
        graph_Pretz_vs_z_1->SetLineColor(kTeal+5);
        graph_Pretz_vs_z_1->SetMarkerColor(kTeal+5);
        graph_Pretz_vs_z_2->SetLineColor(kAzure+5);
        graph_Pretz_vs_z_2->SetMarkerColor(kAzure+5);
        graph_Pretz_vs_z_3->SetLineColor(kViolet+5);
        graph_Pretz_vs_z_3->SetMarkerColor(kViolet+5);
        graph_Pretz_vs_z_4->SetLineColor(kPink+5);
        graph_Pretz_vs_z_4->SetMarkerColor(kPink+5);
        graph_Pretz_vs_z_5->SetLineColor(kOrange+5);
        graph_Pretz_vs_z_5->SetMarkerColor(kOrange+5);

        TCanvas* c_Pretz_z_2 = new TCanvas(Form("Aut_Pretz_vs_z_4d_xQ2_bin%d", x+1), "Pretzelosity Asymmetry vs z for bin (n, x_{B}-Q^{2})", 800, 600);
        graph_Pretz_vs_z->Draw("A");
        graph_Pretz_vs_z_1->Draw("P SAME");
        graph_Pretz_vs_z_2->Draw("P SAME");
        graph_Pretz_vs_z_3->Draw("P SAME");
        graph_Pretz_vs_z_4->Draw("P SAME");
        graph_Pretz_vs_z_5->Draw("P SAME");
        TLine* guideLine2 = new TLine(graph_Pretz_vs_z->GetXaxis()->GetXmin(), 0, graph_Pretz_vs_z->GetXaxis()->GetXmax(), 0);
        guideLine2->SetLineStyle(2);  
        guideLine2->SetLineColor(kGray+1);
        guideLine2->Draw();

        TLegend* legend = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
        legend->AddEntry(graph_Pretz_vs_z_1, "0 < P_{hT} [GeV] < 0.2", "ep");
        legend->AddEntry(graph_Pretz_vs_z_2, "0.2 < P_{hT} [GeV] < 0.4", "ep");
        legend->AddEntry(graph_Pretz_vs_z_3, "0.4 < P_{hT} [GeV] < 0.6", "ep");
        legend->AddEntry(graph_Pretz_vs_z_4, "0.6 < P_{hT} [GeV] < 0.8", "ep");
        legend->AddEntry(graph_Pretz_vs_z_5, "0.8 < P_{hT} [GeV] < 1.2", "ep");
        legend->SetFillStyle(0);  // Transparent background
        legend->Draw();
        c_Pretz_z_2->Update();
        c_Pretz_z_2->Write();
    }


    
    // ________________________________________________________________________________________________________________________________________________________

    dir_Sivers_funct->cd();
    TH1F* x_vaaals = new TH1F("x_vals", "x_{B} mean | z>0.2, M_{x}>1.6 GeV; <x_{B}>; counts", nbin, 0, 1);
    TH1F* Pt_vaaals = new TH1F("Pt_vals", "P_{hT} mean | z>0.2, M_{x}>1.6 GeV; <P_{hT}> [GeV]; counts", nbin, 0, 1.2);
    for (int x = 0; x < xB_nBins; x++) {
        if (vec_pion_Pt[x].empty()) continue; // Skip empty bins
        // Compute mean xB for this bin
        double sum_Pt = 0.0;
        for (double val : vec_pion_Pt[x]) sum_Pt += val;
        double mean_Pt = sum_Pt / vec_pion_Pt[x].size();
        vec_pion_Pt_mean.push_back(mean_Pt);
        for (int z = 0; z < z_nBins; z++){
            if (vec_pion_Pt_4d[x][z].empty()) continue;
            double sum_Pt_4d = 0, sum_xB_4d = 0, sum_z_4d = 0, sum_y_4d = 0;
            for (double valPt : vec_pion_Pt_4d[x][z]) sum_Pt_4d += valPt;
            double m_Pt = sum_Pt_4d / vec_pion_Pt_4d[x][z].size();
            Pt_vaaals->Fill(m_Pt);
            for (double valxB : vec_pion_xB_4d[x][z]) sum_xB_4d += valxB;
            double m_xB = sum_xB_4d / vec_pion_xB_4d[x][z].size();
            x_vaaals->Fill(m_xB);
            for (double valy : vec_pion_y_4d[x][z]) sum_y_4d += valy;
            double m_y = sum_y_4d / vec_pion_y_4d[x][z].size();
            for (double valz : vec_pion_z_4d[x][z]) sum_z_4d += valz;
            double m_z = sum_z_4d / vec_pion_z_4d[x][z].size();
            if (A_UT_sivers_4d[x][z] > 0.3 || A_UT_sivers_4d[x][z] < -0.3) continue;
            if (A_UT_sivers_err_4d[x][z] > 0.06) continue;
            vec_pion_xB_4d_mean.push_back(m_xB);
            vec_pion_Pt_4d_mean.push_back(m_Pt);
            vec_pion_z_4d_mean.push_back(m_z);
            vec_pion_y_4d_mean.push_back(m_y);
            A_UT_sivers_single_4d.push_back(A_UT_sivers_4d[x][z]);
            A_UT_sivers_single_err_4d.push_back(A_UT_sivers_err_4d[x][z]);
            if (m_xB <= 0 || m_xB >= 1) cout << "indice x:" << x << " e z: "<< z << "ha xB fuori range" << endl;
            if (m_z <= 0 || m_z >= 1) cout << "indice x:" << x << " e z: "<< z << "ha z fuori range" << endl;
            if (m_Pt <= 0 || m_xB >= 1.2) cout << "indice x:" << x << " e z: "<< z << "ha Pt fuori range" << endl;
        }
    }
    x_vaaals->Write(), Pt_vaaals->Write();
    cout << "Bin used for the Sivers fit (err < 0.06): " << A_UT_sivers_single_4d.size() << endl;
    double Sivers_fit_params_xQ2[8];
    double cov_matrix[8][8];
    FitSiversFunction(vec_pion_xB_4d_mean, vec_pion_Pt_4d_mean, A_UT_sivers_single_4d, A_UT_sivers_single_err_4d, Sivers_fit_params_xQ2, cov_matrix);

    double Pt_fixed = 0.3736;
    TGraph* sivers_model_xQ2 = new TGraph();
    TGraph* sivers_model_fperp = new TGraph();
    TGraph* sivers_model_xQ22 = new TGraph();
    TGraph* sivers_model_xQ23 = new TGraph();
    TGraphErrors* error_band = BuildSiversErrorBand(Pt_fixed, Sivers_fit_params_xQ2, cov_matrix, 1);
    TGraphErrors* error_band_f = BuildSiversErrorBand(Pt_fixed, Sivers_fit_params_xQ2, cov_matrix, -1);
    error_band->SetLineColor(kRed-10), error_band_f->SetLineColor(kRed-10);;
    //double Pt_fixed = 0.3;
    int n_points = 100;
    for (int i = 0; i < n_points; ++i) {
        double x = 0.01 + i * (0.7 - 0.01) / (n_points - 1); // da 0.01 a 0.7
        double x_f = 0.01 + i * (0.999 - 0.01) / (n_points - 1); // da 0.01 a 0.7
        double pt = 0.3736, pt2 = 0.8, pt3 = 0.1;
        double y = SiversModel(x, Pt_fixed, Sivers_fit_params_xQ2);
        double yf = SiversModel(x_f, Pt_fixed, Sivers_fit_params_xQ2);
        double y2 = SiversModel(x, pt2, Sivers_fit_params_xQ2);
        double y3 = SiversModel(x, pt3, Sivers_fit_params_xQ2);
        sivers_model_xQ2->SetPoint(i, x, y);
        sivers_model_fperp->SetPoint(i, x_f, -yf);
        sivers_model_xQ22->SetPoint(i, x, y2);
        sivers_model_xQ23->SetPoint(i, x, y3);
    }
    sivers_model_xQ2->SetLineColor(kRed), sivers_model_fperp->SetLineColor(kRed);
    sivers_model_xQ22->SetLineColor(kBlue);
    sivers_model_xQ23->SetLineColor(kOrange);
    sivers_model_xQ2->SetLineWidth(2), sivers_model_fperp->SetLineWidth(2),sivers_model_xQ23->SetLineWidth(2), sivers_model_xQ22->SetLineWidth(2);

    TCanvas* c_test = new TCanvas("c_SiversModel", "", 800, 600);
    error_band_f->SetTitle("Sivers Model | #Delta^{N} f_{q/p^{#uparrow}}(x) = -f_{1T}^{#perpq}(x) | k_{#perp} fixed; x_{B}; Sivers model");
    error_band_f->GetXaxis()->SetLimits(0, 1), error_band_f->GetXaxis()->SetRangeUser(0,1);
    sivers_model_fperp->GetXaxis()->SetLimits(0, 1), sivers_model_fperp->GetXaxis()->SetRangeUser(0,1);
    error_band_f->GetYaxis()->SetNoExponent(kTRUE); 
    //error_band_f->GetYaxis()->SetMoreLogLabels(kTRUE);
    error_band_f->Draw("A");
    sivers_model_fperp->Draw("L SAME"); 
    TLine* gLine = new TLine(0, 0, 1, 0);
    gLine->SetLineStyle(2), gLine->SetLineColor(kGray+1);
    gLine->Draw();
    c_test->Update();
    c_test->Write();

    TCanvas* c_sivers = new TCanvas("Sivers_binxQ2", "Sivers Fit", 800, 600);
    graph_AUT_vs_xB->GetYaxis()->SetLimits(-0.05, 0.06);
    graph_AUT_vs_xB->SetTitle("Sivers Asymmetry vs x_{B} (x_{B}-Q^{2} binning) | k_{#perp} fixed; x_{B}; A_{UT}^{sin(#Phi_{h}-#Phi_{s})}");
    graph_AUT_vs_xB->Draw("A");  
    error_band->Draw("SAME"); 
    graph_AUT_vs_xB->Draw("P SAME");         
    sivers_model_xQ2->Draw("L SAME");     
    //sivers_model_xQ22->Draw("L SAME"); 
    //sivers_model_xQ23->Draw("L SAME");
    guideLine->Draw();

    TLegend* leg = new TLegend(0.15, 0.72, 0.4, 0.88);
    leg->AddEntry(graph_AUT_vs_xB, "Measured A_{UT}^{Sivers}", "ep");
    //leg->AddEntry(sivers_model_xQ23, "Sivers fit model P_{hT} = 0.1 GeV", "l");
    leg->AddEntry(sivers_model_xQ2, "Sivers fit model P_{hT} = 0.37 GeV", "l");
    leg->AddEntry(sivers_model_xQ2, "PDF: NNPDF31_lo_as_0130 | Fixed Q^{2}");
    leg->AddEntry(error_band, "Sivers model 0.95 CL", "l");
    //leg->AddEntry(sivers_model_xQ22, "Sivers fit model P_{hT} = 1 GeV", "l");
    leg->Draw();

    c_sivers->Update();
    c_sivers->Write();

    cout << "Sivers parameters" << endl;
    cout << "Nq_up       =  " << Sivers_fit_params_xQ2[0] << endl;
    cout << "Nq_down     =  " << Sivers_fit_params_xQ2[1] << endl;
    cout << "Nq_sea      =  " << Sivers_fit_params_xQ2[2] << endl;
    cout << "alpha_u     =  " << Sivers_fit_params_xQ2[3] << endl;
    cout << "alpha_down  =  " << Sivers_fit_params_xQ2[4] << endl;
    cout << "alpha_sea   =  " << Sivers_fit_params_xQ2[5] << endl;
    cout << "beta        =  " << Sivers_fit_params_xQ2[6] << endl;
    cout << "M^2_1       =  " << Sivers_fit_params_xQ2[7]*Sivers_fit_params_xQ2[7] << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << "" << endl;



    TGraph* g_u = new TGraph();
    TGraph* g_d = new TGraph();
    TGraph* g_sea = new TGraph();
    TGraph* band_u = BuildSiversErrorBand_flavour(0.3736, Sivers_fit_params_xQ2, cov_matrix, 0, 1.96);
    TGraph* band_d = BuildSiversErrorBand_flavour(0.3736, Sivers_fit_params_xQ2, cov_matrix, 1, 1.96);
    TGraph* band_sea = BuildSiversErrorBand_flavour(0.3736, Sivers_fit_params_xQ2, cov_matrix, 2, 1.96);
    TGraph* band_u_67 = BuildSiversErrorBand_flavour(0.3736, Sivers_fit_params_xQ2, cov_matrix, 0, 1);
    TGraph* band_d_67 = BuildSiversErrorBand_flavour(0.3736, Sivers_fit_params_xQ2, cov_matrix, 1, 1);
    TGraph* band_sea_67 = BuildSiversErrorBand_flavour(0.3736, Sivers_fit_params_xQ2, cov_matrix, 2, 1);
    band_u->SetLineColor(kBlue-10), band_d->SetLineColor(kGreen-10), band_sea->SetLineColor(kMagenta-10);
    band_u_67->SetLineColor(kBlue-9), band_d_67->SetLineColor(kGreen-9), band_sea_67->SetLineColor(kMagenta-9);

    int n_p = 100;

    for (int i = 0; i < n_p; ++i) {
        double x = 0.01 + i * (0.85 - 0.01) / (n_p - 1);
        double yu = SiversModel_u(x, Pt_fixed, Sivers_fit_params_xQ2);
        double yd = SiversModel_d(x, Pt_fixed, Sivers_fit_params_xQ2);
        double ysea = SiversModel_sea(x, Pt_fixed, Sivers_fit_params_xQ2);
        g_u->SetPoint(i, x, -yu*x);
        g_d->SetPoint(i, x, -yd*x);
        g_sea->SetPoint(i, x, -ysea*x);
    }

    g_u->SetLineColor(kBlue), g_u->SetLineWidth(2), g_u->SetTitle("u quark");
    g_d->SetLineColor(kGreen + 2), g_d->SetLineWidth(2), g_d->SetTitle("d quark");
    g_sea->SetLineColor(kMagenta + 2), g_sea->SetLineWidth(2), g_sea->SetTitle("sea quarks");
    //g_u->GetYaxis()->SetRangeUser(-0.02,0.02);
    //error_band_f->SetTitle("Sivers Model | #Delta^{N} f_{q/p^{#uparrow}}(x) = -f_{1T}^{#perpq}(x); x_{B}; f_{1T}^{#perpq}(x)");
    band_u->GetYaxis()->SetNoExponent(kTRUE); 
    band_u->SetTitle("Sivers f_{1T}^{#perp(1)u} contributions by flavour | k_{#perp} fixed;x_{B};xf_{1T}^{#perp(1)u}(x)");
    band_d->SetTitle("Sivers f_{1T}^{#perp(1)d} contributions by flavour | k_{#perp} fixed;x_{B};xf_{1T}^{#perp(1)d}(x)");
    band_sea->SetTitle("Sivers f_{1T}^{#perp(1)sea} contributions by flavour | k_{#perp} fixed;x_{B};xf_{1T}^{#perp(1)sea}(x)");
    TCanvas* c_sivers_up= new TCanvas("c_Sivers_up", "Sivers by flavour", 800, 600);
    TLine* zero = new TLine(0, 0, 0.934, 0);
    zero->SetLineStyle(2);
    zero->SetLineColor(kGray + 2);
    band_u->Draw("A");
    band_u_67->Draw("SAME");
    g_u->Draw("L SAME");
    TLegend* l_up = new TLegend(0.72, 0.78, 0.88, 0.88);
    l_up->AddEntry(band_u, "CL 0.95", "l");
    l_up->AddEntry(band_u_67, "CL 0.67", "l");
    l_up->Draw();
    zero->Draw();
    c_sivers_up->Update();
    c_sivers_up->Write();
    TCanvas* c_sivers_down= new TCanvas("c_Sivers_down", "Sivers by flavour", 800, 600);
    band_d->Draw("A");
    band_d_67->Draw("SAME");
    g_d->Draw("L SAME");
    TLegend* l_down = new TLegend(0.72, 0.25, 0.88, 0.15);
    l_down->AddEntry(band_d, "CL 0.95", "l");
    l_down->AddEntry(band_d_67, "CL 0.67", "l");
    l_down->Draw();
    zero->Draw();
    c_sivers_down->Update();
    c_sivers_down->Write();
    TCanvas* c_sivers_sea= new TCanvas("c_Sivers_sea", "Sivers by flavour", 800, 600);
    band_sea->Draw("A");
    band_sea_67->Draw("SAME");
    g_sea->Draw("L SAME");
    TLegend* l_sea = new TLegend(0.72, 0.25, 0.88, 0.15);
    l_sea->AddEntry(band_sea, "CL 0.95", "l");
    l_sea->AddEntry(band_sea_67, "CL 0.67", "l");
    l_sea->Draw();
    zero->Draw();
    c_sivers_sea->Update();
    c_sivers_sea->Write();

    TCanvas* c_contributions = new TCanvas("c_Sivers_flavours_xB", "Sivers by flavour", 800, 600);
    band_u->GetYaxis()->SetRangeUser(-0.16,0.19);
    band_u->SetTitle("Sivers f_{1T}^{#perp(1)q} contributions by flavour;x_{B};xf_{1T}^{#perp(1)q}(x)");
    band_u->Draw("A");
    band_d->Draw("SAME");
    band_sea->Draw("SAME");
    g_u->Draw("L SAME");
    g_d->Draw("L SAME");
    g_sea->Draw("L SAME");

    TLegend* lego = new TLegend(0.70, 0.75, 0.88, 0.88);
    lego->AddEntry(g_u, "u quark | 95 CL", "l");
    lego->AddEntry(g_d, "d quark | 95 CL", "l");
    lego->AddEntry(g_sea, "sea quarks | 95 CL", "l");
    lego->Draw();
    zero->Draw();

    c_contributions->Update();
    c_contributions->Write();


    TGraph* g_siv_Pt_total = new TGraph();
    TGraph* g_siv_Pt_up = new TGraph();
    TGraph* g_siv_Pt_down = new TGraph();
    TGraph* g_siv_Pt_sea = new TGraph();
    double x_fixed = 0.2406; // it's our <xB>
    TGraph* band_Pt_u = BuildSiversErrorBand_flavour_Pt(x_fixed, Sivers_fit_params_xQ2, cov_matrix, 0, 1.96);
    TGraph* band_Pt_d = BuildSiversErrorBand_flavour_Pt(x_fixed, Sivers_fit_params_xQ2, cov_matrix, 1, 1.96);
    TGraph* band_Pt_sea = BuildSiversErrorBand_flavour_Pt(x_fixed, Sivers_fit_params_xQ2, cov_matrix, 2, 1.96);
    band_Pt_u->SetLineColor(kBlue-10), band_Pt_d->SetLineColor(kGreen-10), band_Pt_sea->SetLineColor(kMagenta-10);
    int point = 100;
    for (int i = 0; i < point; ++i) {
        double Pt = 0.01 + i * (0.95 - 0.01) / (point -1);
        double f_tot = SiversModel(x_fixed, Pt, Sivers_fit_params_xQ2);  
        double f_up = SiversModel_u(x_fixed, Pt, Sivers_fit_params_xQ2);  
        double f_down = SiversModel_d(x_fixed, Pt, Sivers_fit_params_xQ2); 
        double f_sea = SiversModel_sea(x_fixed, Pt, Sivers_fit_params_xQ2); 
        g_siv_Pt_total->SetPoint(i, Pt, -x_fixed*f_tot);
        g_siv_Pt_up->SetPoint(i, Pt, -x_fixed*f_up);
        g_siv_Pt_down->SetPoint(i, Pt, -x_fixed*f_down);
        g_siv_Pt_sea->SetPoint(i, Pt, -x_fixed*f_sea);
    }
    band_Pt_u->SetTitle("Sivers x#Delta^{N} f^{(1)}_{q/p^{#uparrow}}(x_{mean}, P_{hT}) contributions; P_{hT} [GeV];-x#Delta^{N} f^{(1)}_{q/p^{#uparrow}}(P_{hT})");
    g_siv_Pt_total->SetLineWidth(2), g_siv_Pt_total->SetLineColor(kRed);
    g_siv_Pt_up->SetLineWidth(2), g_siv_Pt_up->SetLineColor(kBlue);
    g_siv_Pt_down->SetLineWidth(2), g_siv_Pt_down->SetLineColor(kGreen +2);
    g_siv_Pt_sea->SetLineWidth(2), g_siv_Pt_sea->SetLineColor(kMagenta+2);
    band_Pt_u->GetYaxis()->SetMaxDigits(3), band_Pt_u->GetYaxis()->SetRangeUser(-0.2, 0.25);
    TCanvas* c_siv_Pt = new TCanvas("c_Sivers_flavours_Pt", "", 800, 600);
    band_Pt_u->Draw();
    band_Pt_d->Draw("SAME");
    band_Pt_sea->Draw("SAME");
    //g_siv_Pt_total->Draw("L SAME");
    g_siv_Pt_up->Draw("L SAME");
    g_siv_Pt_down->Draw("L SAME");
    g_siv_Pt_sea->Draw("L SAME");
    TLegend* leg_Pt = new TLegend(0.70, 0.75, 0.88, 0.88);
    //leg_Pt->AddEntry(g_siv_Pt_total, "all contribution", "l");
    leg_Pt->AddEntry(g_siv_Pt_up, "u quark | 95 CL", "l");
    leg_Pt->AddEntry(g_siv_Pt_down, "d quark | 95 CL", "l");
    leg_Pt->AddEntry(g_siv_Pt_sea, "sea quarks | 95 CL", "l");
    leg_Pt->Draw();
    TLine* zero_pt = new TLine(0, 0, 1.044, 0);
    zero_pt->SetLineStyle(2);
    zero_pt->SetLineColor(kGray + 2);
    zero_pt->Draw();
    c_siv_Pt->Update();
    c_siv_Pt->Write();








    // _____________________________________________________ COLLINS _________________________________________________________________________________________

    dir_Collins_funct->cd();
    for (int x = 0; x < xB_nBins; x++) {
        for (int z = 0; z < z_nBins; z++){
            if (vec_pion_Pt_4d[x][z].empty()) continue;
            if (vec_pion_xB_4d[x][z].empty()) continue;
            double sum_Pt_4d = 0, sum_xB_4d = 0, sum_z_4d = 0, sum_y_4d = 0, sum_eps = 0;
            for (double valPt : vec_pion_Pt_4d[x][z]) sum_Pt_4d += valPt;
            double m_Pt = sum_Pt_4d / vec_pion_Pt_4d[x][z].size();
            for (double valxB : vec_pion_xB_4d[x][z]) sum_xB_4d += valxB;
            double m_xB = sum_xB_4d / vec_pion_xB_4d[x][z].size();
            for (double valy : vec_pion_y_4d[x][z]) sum_y_4d += valy;
            double m_y = sum_y_4d / vec_pion_y_4d[x][z].size();
            for (double valz : vec_pion_z_4d[x][z]) sum_z_4d += valz;
            double m_z = sum_z_4d / vec_pion_z_4d[x][z].size();
            for (double valeps : vec_pion_epsilon_4d[x][z]) sum_eps += valeps;
            double mean_eps = sum_eps / vec_pion_epsilon_4d[x][z].size();
            if (A_UT_collins_4d[x][z] > 0.3 || A_UT_collins_4d[x][z] < -0.3) continue;
            if (A_UT_collins_err_4d[x][z] > 0.042) continue;
            vec_pion_xB_4d_mean_coll.push_back(m_xB);
            vec_pion_Pt_4d_mean_coll.push_back(m_Pt);
            vec_pion_z_4d_mean_coll.push_back(m_z);
            vec_pion_y_4d_mean_coll.push_back(m_y);
            A_UT_collins_single_4d.push_back(A_UT_collins_4d[x][z]);
            A_UT_collins_single_err_4d.push_back(A_UT_collins_err_4d[x][z]);
            if (m_xB <= 0 || m_xB >= 1) cout << "indice x:" << x << " e z: "<< z << "ha xB fuori range" << endl;
            if (m_z <= 0 || m_z >= 1) cout << "indice x:" << x << " e z: "<< z << "ha z fuori range" << endl;
            if (m_Pt <= 0 || m_xB >= 1.2) cout << "indice x:" << x << " e z: "<< z << "ha Pt fuori range" << endl;
        }
    }

    cout << "Bin used for the Collins fit (err < 0.04): " << vec_pion_y_4d_mean_coll.size() << endl;

    double Collins_fit_params_xQ2[10];
    double cov_Coll_matrix[10][10];
    FitCollins(vec_pion_xB_4d_mean_coll, vec_pion_z_4d_mean_coll, vec_pion_Pt_4d_mean_coll, A_UT_collins_single_4d, A_UT_collins_single_err_4d, Collins_fit_params_xQ2, cov_Coll_matrix, vec_pion_y_4d_mean_coll);

    double z_fixed = 0.3805;
    double y_fixed = 0.585;
    TGraph* collins_model_xQ2 = new TGraph();
    TGraph* collins_model_fperp = new TGraph();
    TGraphErrors* error_band_collins = BuildCollinsErrorBand(z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2, cov_Coll_matrix, 1);
    TGraphErrors* error_band_collins_f = BuildCollinsErrorBand(z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2, cov_Coll_matrix, -1);
    error_band_collins->SetLineColor(kRed-10), error_band_collins_f->SetLineColor(kRed-10);;
    //double Pt_fixed = 0.3;
    double n_points_coll = 100;
    double n0 = 1 / n_points_coll;
    for (int i = 0; i < n_points_coll; ++i) {
        double x = n0 + i * (0.7 - n0) / (n_points_coll - 1); // da 0.01 a 0.7
        double x_f = n0 + i * (0.99 - n0) / (n_points_coll - 1); // da 0.01 a 0.99
        double y = CollinsAsymmetry(x, z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2);
        double yf = CollinsAsymmetry(x_f, z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2);
        collins_model_xQ2->SetPoint(i, x, y);
        collins_model_fperp->SetPoint(i, x_f, -yf);
    }
    collins_model_xQ2->SetLineColor(kRed), collins_model_fperp->SetLineColor(kRed);
    collins_model_xQ2->SetLineWidth(2), collins_model_fperp->SetLineWidth(2);

    
    TCanvas* c_test_coll = new TCanvas("c_CollinsModel", "", 800, 600);
    error_band_collins_f->SetTitle("Collins Model | #Delta^{N} D_{q/p^{#uparrow}}(x) = -D_{1T}^{#perpq}(x) | k_{#perp} fixed; x_{B}; Collins model");
    error_band_collins_f->GetYaxis()->SetNoExponent(kTRUE); 
    //error_band_f->GetYaxis()->SetMoreLogLabels(kTRUE);
    error_band_collins_f->Draw("A");
    collins_model_fperp->Draw("L SAME"); 
    gLine->Draw();
    c_test_coll->Update();
    c_test_coll->Write();
    

    TCanvas* c_collins = new TCanvas("Collins_binxQ2", "Collins Fit", 800, 600);
    graph_Collins_vs_xB->GetYaxis()->SetLimits(-0.05, 0.06);
    graph_Collins_vs_xB->SetTitle("Collins Asymmetry vs x_{B} (x_{B}-Q^{2} binning) | k_{#perp} fixed; x_{B}; A_{UT}^{sin(#Phi_{h}+#Phi_{s})}");
    graph_Collins_vs_xB ->Draw("A");  
    error_band_collins->Draw("SAME"); 
    graph_Collins_vs_xB->Draw("P SAME");         
    collins_model_xQ2->Draw("L SAME");     
    guideLine->Draw();

    TLegend* leg_coll = new TLegend(0.15, 0.72, 0.4, 0.88);
    leg_coll->AddEntry(graph_Collins_vs_xB , "Measured A_{UT}^{Collins}", "ep");
    leg_coll->AddEntry(collins_model_xQ2, "Collins fit model <P_{hT},z,y>", "l");
    leg_coll->AddEntry(collins_model_xQ2, "PDF: JAM20-SIDIS_PDF_proton_nlo", "l");
    leg_coll->AddEntry(collins_model_xQ2, "FF: JAM20-SIDIS_FF_pion_nlo", "l");
    leg_coll->AddEntry(error_band_collins, "Collins model 0.95 CL", "l");
    leg_coll->Draw();

    c_collins->Update();
    c_collins->Write();


    cout << "Collins parameters" << endl;
    cout << "N^up_T      =  " << Collins_fit_params_xQ2[0] << endl;
    cout << "N^down_T    =  " << Collins_fit_params_xQ2[1] << endl;
    cout << "N^sea_T     =  " << Collins_fit_params_xQ2[2] << endl;
    cout << "alpha       =  " << Collins_fit_params_xQ2[3] << endl;
    cout << "beta        =  " << Collins_fit_params_xQ2[4] << endl;
    cout << "N_C_fav     =  " << Collins_fit_params_xQ2[5] << endl;
    cout << "N_C_dis     =  " << Collins_fit_params_xQ2[6] << endl;
    cout << "gamma_c     =  " << Collins_fit_params_xQ2[7] << endl;
    cout << "delta_c     =  " << Collins_fit_params_xQ2[8] << endl;
    cout << "M^2_H       =  " << Collins_fit_params_xQ2[9]*Collins_fit_params_xQ2[9] << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << "" << endl;


    TGraph* g_coll_u = new TGraph();
    TGraph* g_coll_d = new TGraph();
    TGraph* g_coll_sea = new TGraph();
    TGraph* band_coll_u = BuildCollinsErrorBand_flavour(z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2, cov_Coll_matrix, 0, 1.96);
    TGraph* band_coll_d = BuildCollinsErrorBand_flavour(z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2, cov_Coll_matrix, 1, 1.96);
    TGraph* band_coll_sea = BuildCollinsErrorBand_flavour(z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2, cov_Coll_matrix, 2, 1.96);
    TGraph* band_coll_u_67 = BuildCollinsErrorBand_flavour(z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2, cov_Coll_matrix, 0, 1);
    TGraph* band_coll_d_67 = BuildCollinsErrorBand_flavour(z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2, cov_Coll_matrix, 1, 1);
    TGraph* band_coll_sea_67 = BuildCollinsErrorBand_flavour(z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2, cov_Coll_matrix, 2, 1);
    band_coll_u->SetLineColor(kBlue-10), band_coll_d->SetLineColor(kGreen-10), band_coll_sea->SetLineColor(kMagenta-10);
    band_coll_u_67->SetLineColor(kBlue-9), band_coll_d_67->SetLineColor(kGreen-9), band_coll_sea_67->SetLineColor(kMagenta-9);

    int n_pc = 100;

    for (int i = 0; i < n_pc; ++i) {
        double x = 0.01 + i * (0.9 - 0.01) / (n_pc - 1);
        double yu = DeltaT_q(x, z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2,"u");
        double yd = DeltaT_q(x, z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2, "d");
        double ysea = DeltaT_q(x, z_fixed, Pt_fixed, y_fixed, Collins_fit_params_xQ2, "sea");
        g_coll_u->SetPoint(i, x, yu*x);
        g_coll_d->SetPoint(i, x, yd*x);
        g_coll_sea->SetPoint(i, x, ysea*x);
    }

    g_coll_u->SetLineColor(kBlue), g_coll_u->SetLineWidth(2), g_coll_u->SetTitle("u quark");
    g_coll_d->SetLineColor(kGreen + 2), g_coll_d->SetLineWidth(2), g_coll_d->SetTitle("d quark");
    g_coll_sea->SetLineColor(kMagenta + 2), g_coll_sea->SetLineWidth(2), g_coll_sea->SetTitle("sea quarks");
    //g_coll_u->GetYaxis()->SetRangeUser(-0.02,0.02);
    //error_band_coll_f->SetTitle("Collins Model | #Delta^{N} f_{q/p^{#uparrow}}(x) = -f_{1T}^{#perpq}(x); x_{B}; f_{1T}^{#perpq}(x)");
    band_coll_u->GetYaxis()->SetNoExponent(kTRUE); 
    band_coll_u->SetTitle("#Delta_{T}q(x_{B},Q_{0}^{2}) = #it{N_{u}}^{T}#frac{1}{2}[#it{f}_{q/p}(x_{B},Q^{2}_{0}) + #Deltaq(x_{B},Q^{2}_{0})] ;x_{B}; x#Delta_{T}q^{u}");
    band_coll_d->SetTitle("#Delta_{T}q(x_{B},Q_{0}^{2}) = #it{N_{d}}^{T}#frac{1}{2}[#it{f}_{q/p}(x_{B},Q^{2}_{0}) + #Deltaq(x_{B},Q^{2}_{0})] ;x_{B}; x#Delta_{T}q^{d}");
    band_coll_sea->SetTitle("#Delta_{T}q(x_{B},Q_{0}^{2}) = #it{N_{sea}}^{T}#frac{1}{2}[#it{f}_{q/p}(x_{B},Q^{2}_{0}) + #Deltaq(x_{B},Q^{2}_{0})] ;x_{B}; x#Delta_{T}q^{sea}");
    TCanvas* c_collins_up= new TCanvas("c_Collins_up", "Collins by flavour", 800, 600);
    TLine* zero_coll = new TLine(0, 0, 0.9898, 0);
    zero_coll->SetLineStyle(2), zero_coll->SetLineColor(kGray + 2);
    band_coll_u->Draw("A");
    band_coll_u_67->Draw("SAME");
    g_coll_u->Draw("L SAME");
    TLegend* l_coll_up = new TLegend(0.72, 0.78, 0.88, 0.88);
    l_coll_up->AddEntry(band_coll_u, "CL 0.95", "l");
    l_coll_up->AddEntry(band_coll_u_67, "CL 0.67", "l");
    l_coll_up->Draw();
    zero_coll->Draw();
    c_collins_up->Update();
    c_collins_up->Write();
    TCanvas* c_collins_down= new TCanvas("c_Collins_down", "Collins by flavour", 800, 600);
    band_coll_d->Draw("A");
    band_coll_d_67->Draw("SAME");
    g_coll_d->Draw("L SAME");
    TLegend* l_coll_down = new TLegend(0.12, 0.25, 0.25, 0.15);
    l_coll_down->AddEntry(band_coll_d, "CL 0.95", "l");
    l_coll_down->AddEntry(band_coll_d_67, "CL 0.67", "l");
    l_coll_down->Draw();
    zero_coll->Draw();
    c_collins_down->Update();
    c_collins_down->Write();
    TCanvas* c_collins_sea= new TCanvas("c_Collins_sea", "Collins by flavour", 800, 600);
    band_coll_sea->Draw("A");
    band_coll_sea_67->Draw("SAME");
    g_coll_sea->Draw("L SAME");
    TLegend* l_coll_sea = new TLegend(0.12, 0.25, 0.25, 0.15);
    l_coll_sea->AddEntry(band_coll_sea, "CL 0.95", "l");
    l_coll_sea->AddEntry(band_coll_sea_67, "CL 0.67", "l");
    l_coll_sea->Draw();
    zero_coll->Draw();
    c_collins_sea->Update();
    c_collins_sea->Write();

    TCanvas* c_contributions_coll = new TCanvas("c_Collins_flavours_xB", "Collins by flavour", 800, 600);
    band_coll_u->GetYaxis()->SetRangeUser(-0.4,0.4);
    band_coll_u->SetTitle("#Delta_{T}q(x_{B},Q_{0}^{2}) = #it{N_{q}}^{T}#frac{1}{2}[#it{f}_{q/p}(x_{B},Q^{2}_{0}) + #Deltaq(x_{B},Q^{2}_{0})] ;x_{B}; x#Delta_{T}q^{q}");
    band_coll_u->Draw("A");
    band_coll_u_67->Draw("SAME");
    g_coll_u->Draw("L SAME");
    band_coll_d_67->Draw("SAME");
    band_coll_d->Draw("SAME");
    g_coll_d->Draw("L SAME");
    band_coll_sea->Draw("SAME");
    band_coll_sea_67->Draw("SAME");
    g_coll_sea->Draw("L SAME");

    TLegend* lego_coll = new TLegend(0.70, 0.75, 0.88, 0.88);
    lego_coll->AddEntry(g_coll_u, "u quark | 95 CL", "l");
    lego_coll->AddEntry(g_coll_d, "d quark | 95 CL", "l");
    lego_coll->AddEntry(g_coll_sea, "sea quarks | 95 CL", "l");
    lego_coll->Draw();
    zero_coll->Draw();

    c_contributions_coll->Update();
    c_contributions_coll->Write();

    




    // ________________________________________________________________________________________________________________________________________________________






    
    //outFile.Write();
    outFile.Close(); 

    cout << "Analysis completed, result save in: " << outputFile << endl;
    cout << "\n";

}