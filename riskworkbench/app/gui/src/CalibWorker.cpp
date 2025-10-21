#include "CalibWorker.hpp"

#include <QFileInfo>
#include <QDateTime>
#include <QDir>
#include <QDebug>

#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>

#include <rw/io/market_csv.hpp>
#include <rw/pricing/analytic_bs.hpp>
#include <rw/pricing/implied_vol.hpp>

using gui::CalibWorker;

namespace {
inline bool fin(double x){ return std::isfinite(x); }
inline double N(double x){ return 0.5 * std::erfc(-x * M_SQRT1_2); }

// petit BS fermé, robuste aux bords
static double bs_price_cp(double S0,double K,double r,double q,double T,double sigma,bool is_call){
  if (!(fin(S0)&&fin(K)&&fin(r)&&fin(q)&&fin(T)&&fin(sigma)) || S0<=0||K<=0||T<0)
    return std::numeric_limits<double>::quiet_NaN();
  if (T==0.0 || sigma<=0.0){
    const double df_r=std::exp(-r*T), df_q=std::exp(-q*T);
    const double fwd = S0*df_q - K*df_r;
    return is_call? std::max(0.0, fwd) : std::max(0.0, -fwd);
  }
  const double df_r=std::exp(-r*T), df_q=std::exp(-q*T);
  const double vs = sigma*std::sqrt(T);
  const double m  = std::log(S0/K) + (r-q+0.5*sigma*sigma)*T;
  const double d1 = m / vs;
  const double d2 = d1 - vs;
  return is_call ? S0*df_q*N(d1)-K*df_r*N(d2) : K*df_r*N(-d2)-S0*df_q*N(-d1);
}
}

CalibWorker::CalibWorker(QObject* parent): QObject(parent) {}

QString CalibWorker::sanitize(const QString& s){
  QString o; o.reserve(s.size());
  for (auto ch: s){
    if (ch.isLetterOrNumber() || ch=='_' || ch=='-') o.push_back(ch);
    else if (ch==' ' || ch=='/' || ch=='\\') o.push_back('_');
  }
  return o.isEmpty()? "UNKNOWN" : o;
}

void CalibWorker::loadCsv(const QString& path){
  try {
    emit progress("Lecture CSV…", 0, 0);
    ms_ = rw::market::read_market_surface(path.toStdString());
    auto meta = rw::io::read_market_meta_csv(path.toStdString());
    underlying_ = meta && !meta->underlying.empty() ? QString::fromStdString(meta->underlying) : "UNKNOWN";
    emit loadedMarket(static_cast<quint64>(ms_.rows.size()), ms_.S0, ms_.r, ms_.q, underlying_);
    emit message("CSV chargé.");
  } catch (const std::exception& e){
    emit failed(QString("loadCsv: %1").arg(e.what()));
  }
}

void CalibWorker::invertIv() {
  if (ms_.rows.empty()) { emit failed("invertIv: pas de données."); return; }

  qDebug() << "[Calib] invertIv: rows=" << ms_.rows.size();
  emit progress("Inversion IV…", 0, (int)ms_.rows.size());

  std::size_t priced=0, ok=0, bad=0;

  for (std::size_t i=0; i<ms_.rows.size(); ++i) {
    auto& r = ms_.rows[i];

    if (std::isfinite(r.iv_mid) && r.iv_mid > 0.0) { ++priced; ++ok; continue; }
    if (!std::isfinite(r.mid_price)) { ++bad; continue; }

    auto res = rw::bs::implied_vol_cp(ms_.S0, r.K, ms_.r, ms_.q, r.T, r.mid_price, r.is_call);
    ++priced;
    if (res.converged && res.sigma > 0.0) { r.iv_mid = res.sigma; ++ok; }
    else { ++bad; }

    if ((i % 32) == 0) emit progress("Inversion IV…", (int)i, (int)ms_.rows.size());
  }

  emit progress("Inversion IV…", (int)ms_.rows.size(), (int)ms_.rows.size());
  qDebug() << "[Calib] invertIv: priced=" << priced << " ok=" << ok << " bad=" << bad;

  emit inverted((quint64)priced, (quint64)ok, (quint64)bad);
  emit message(QString("Inversion terminée: priced=%1 ok=%2 bad=%3")
               .arg(priced).arg(ok).arg(bad));
}

void CalibWorker::buildSmile() {
  if (ms_.rows.empty()) { emit failed("buildSmile: pas de marché."); return; }
  try {
    qDebug() << "[Calib] buildSmile: start (rows=" << ms_.rows.size() << ")";
    emit progress("Construction du smile…", 0, 0);

    smile_ = rw::calib::build_smile_surface(ms_);

    QVector<double> mats;
    mats.reserve((int)smile_.slices.size());
    for (const auto& s : smile_.slices) mats.push_back(s.T);

    qDebug() << "[Calib] buildSmile: slices=" << mats.size() << " mats=" << mats;
    emit smileBuilt(mats);
    emit message("Smile construit.");
  } catch (const std::exception& e) {
    emit failed(QString("buildSmile: %1").arg(e.what()));
  }
}

void CalibWorker::fitGlobal(double sigma0){
  if (ms_.rows.empty()) { emit failed("fitGlobal: pas de marché."); return; }
  try{
    double sigma = sigma0;
    auto rep = rw::calib::fit_bs_global_sigma(ms_, sigma);
    emit globalFitted(rep, sigma);
    emit message(QString("σ global= %1  RMSE_price=%2  n=%3")
                 .arg(sigma,0,'g',10).arg(rep.rmse_price,0,'g',10).arg((qulonglong)rep.n));
  } catch (const std::exception& e){
    emit failed(QString("fitGlobal: %1").arg(e.what()));
  }
}

void CalibWorker::runQc(double tick, double onlyConsistentTol, bool priceFromIv){
  Q_UNUSED(tick);
  if (ms_.rows.empty() || smile_.slices.empty()) { emit failed("runQc: marché/smile manquant."); return; }
  residuals_.clear();
  emit progress("QC résidus…", 0, 0);
  auto rows = rw::qc::compare_fit(ms_, smile_);

  // éventuel “price-from-iv” et “only-consistent”
  std::vector<rw::qc::ErrorRow> used; used.reserve(rows.size());
  for (auto r : rows){
    double price_iv = std::numeric_limits<double>::quiet_NaN();
    if (fin(r.iv_mkt) && r.iv_mkt>0.0)
      price_iv = bs_price_cp(ms_.S0, r.K, ms_.r, ms_.q, r.T, r.iv_mkt, r.is_call);
    if (priceFromIv && fin(price_iv)) {
      r.price_mkt = price_iv;
      r.err_price = r.price_fit - r.price_mkt;
    }
    bool keep=true;
    if (onlyConsistentTol>0.0 && fin(price_iv) && fin(r.price_mkt))
      keep = (std::fabs(r.price_mkt - price_iv) <= onlyConsistentTol);
    if (keep) used.push_back(r);
  }
  if (used.empty()) used = rows;

  auto rmse = [&](const auto& v, auto getter){
    double s=0; size_t n=0; for (const auto& e: v){ double z=getter(e); if (fin(z)){ s+=z*z; ++n; } }
    return n? std::sqrt(s/n) : std::numeric_limits<double>::quiet_NaN();
  };
  double rmse_p = rmse(used, [](const auto& e){ return e.err_price; });
  double rmse_iv= rmse(used, [](const auto& e){ return e.err_iv;    });

  residuals_ = std::move(used);
  qDebug() << "QC done: rows=" << rows.size();
  emit qcDone(rmse_p, rmse_iv, static_cast<quint64>(residuals_.size()), residuals_);
  emit message("QC terminé.");
}

void CalibWorker::exportSurfaceCsv(const QString& outDirBase){
  if (smile_.slices.empty()) { emit failed("exportSurfaceCsv: smile vide."); return; }

  const QString und = sanitize(underlying_);
  QDir dir(outDirBase);
  dir.mkpath(und);
  QDir out = QDir(dir.filePath(und));

  const double S0 = ms_.S0;
  const double kmin = 0.7*S0, kmax = 1.3*S0;
  const int NK = 41;

  const auto now = QDateTime::currentDateTime();
  const QString fname = QString("ivgrid_%1_%2.csv").arg(und).arg(now.toString("yyyyMMdd-HHmmss"));
  const QString path = out.filePath(fname);

  std::ofstream f(path.toStdString());
  f << "T,K,iv,bs_price\n";
  for (const auto& sl : smile_.slices){
    for (int i=0;i<NK;++i){
      const double u = (double)i/(double)(NK-1);
      const double K = std::exp(std::log(kmin) + u*(std::log(kmax)-std::log(kmin)));
      const double iv = smile_.iv_at(sl.T, K);
      const double p  = bs_price_cp(S0, K, ms_.r, ms_.q, sl.T, iv, /*call?*/ K>=S0);
      f << sl.T << "," << K << "," << iv << "," << p << "\n";
    }
  }
  f.close();
  emit exportedCsv(path);
}

void CalibWorker::saveProjectJson(){
  if (smile_.slices.empty()) { emit failed("saveProjectJson: smile vide."); return; }

  const QString und = sanitize(underlying_);
  QDir dir("data/projects");
  dir.mkpath(".");
  const auto now = QDateTime::currentDateTime();
  const QString fname = QString("calib_%1_%2.json").arg(und).arg(now.toString("yyyyMMdd-HHmmss"));
  const QString path = dir.filePath(fname);

  std::ofstream f(path.toStdString());
  f << "{\n";
  f << "  \"underlying\": \"" << und.toStdString() << "\",\n";
  f << "  \"timestamp\": \"" << now.toString(Qt::ISODate).toStdString() << "\",\n";
  f << "  \"S0\": " << ms_.S0 << ", \"r\": " << ms_.r << ", \"q\": " << ms_.q << ",\n";
  f << "  \"slices\": [\n";
  for (size_t i=0;i<smile_.slices.size();++i){
    const auto& s = smile_.slices[i];
    f << "    { \"T\": " << s.T << ", \"x\": [";
    for (size_t j=0;j<s.x.size();++j){ f << s.x[j]; if (j+1<s.x.size()) f << ","; }
    f << "], \"y\": [";
    for (size_t j=0;j<s.y.size();++j){ f << s.y[j]; if (j+1<s.y.size()) f << ","; }
    f << "] }";
    if (i+1<smile_.slices.size()) f << ",";     
    f << "\n";
  }
  f << "  ]\n";
  f << "}\n";
  f.close();

  emit projectSaved(path);
}
