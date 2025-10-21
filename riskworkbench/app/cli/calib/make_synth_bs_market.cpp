#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>

static inline double N(double x){ return 0.5*std::erfc(-x/std::sqrt(2.0)); }
static inline double df(double r,double T){ return std::exp(-r*T); }
static double bs_call(double S,double K,double r,double q,double T,double s){
  if (T<=0.0 || s<=0.0) return std::max(0.0, S*df(q,T)-K*df(r,T));
  double v=s*std::sqrt(T), m=std::log(S/K)+(r-q)*T;
  double d1=(m+0.5*s*s*T)/v, d2=d1-v;
  return S*df(q,T)*N(d1)-K*df(r,T)*N(d2);
}
static double bs_put(double S,double K,double r,double q,double T,double s){
  if (T<=0.0 || s<=0.0) return std::max(0.0, K*df(r,T)-S*df(q,T));
  double v=s*std::sqrt(T), m=std::log(S/K)+(r-q)*T;
  double d1=(m+0.5*s*s*T)/v, d2=d1-v;
  return K*df(r,T)*N(-d2)-S*df(q,T)*N(-d1);
}

static void usage(const char* a){
  std::cerr << "Usage: " << a << " out.csv [S0 r q sigma T1,T2,... Kmin Kmax NK]\n";
}

int main(int argc, char** argv){
  if (argc < 9) { usage(argv[0]); return 1; }
  std::string out = argv[1];
  double S0 = std::atof(argv[2]);
  double r  = std::atof(argv[3]);
  double q  = std::atof(argv[4]);
  double sig= std::atof(argv[5]);
  std::string tlist=argv[6];
  double Kmin= std::atof(argv[7]);
  double Kmax= std::atof(argv[8]);
  int NK = std::atoi(argv[9]);

  // parse maturities
  std::vector<double> T;
  size_t p=0;
  while (p<tlist.size()){
    size_t qpos = tlist.find(',', p);
    std::string tok = tlist.substr(p, (qpos==std::string::npos? tlist.size(): qpos)-p);
    if (!tok.empty()) T.push_back(std::atof(tok.c_str()));
    if (qpos==std::string::npos) break; else p = qpos+1;
  }

  std::ofstream f(out);
  if (!f) { std::cerr << "cannot open " << out << "\n"; return 2; }

  // meta
  f << "underlying,S0,r,q\n";
  f << "SYNTH," << S0 << "," << r << "," << q << "\n";
  f << "K,T,type,price_bid,price_ask,mid,iv_bid,iv_ask,iv_mid\n";

  for (double Tm : T) {
    for (int i=0;i<NK;++i) {
      double u = (NK==1 ? 0.0 : double(i)/(NK-1));
      double K = (1.0-u)*Kmin + u*Kmax;

      // alterne call/put pour l’exemple
      bool is_call = (i%2==0);
      double mid = is_call ? bs_call(S0,K,r,q,Tm,sig) : bs_put(S0,K,r,q,Tm,sig);

      // spread bid/ask fictif ±1%
      double bid = mid*0.995, ask = mid*1.005;
      f << K << "," << Tm << "," << (is_call?"C":"P") << ","
        << bid << "," << ask << "," << mid << ",,,\n";
    }
  }
  std::cerr << "written: " << out << "  (sigma*= " << sig << ")\n";
  return 0;
}
